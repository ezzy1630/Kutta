#include "core/ViewFormat.h"

#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QtEndian>

#include <cmath>
#include <cstring>
#include <limits>

namespace kapp {

namespace {

constexpr quint32 kHeaderBytes = 32;
constexpr char kTriMagic[4] = {'K', 'U', 'T', 'V'};
constexpr char kLineMagic[4] = {'K', 'U', 'T', 'L'};

std::optional<double> optionalNumber(const QJsonValue& v)
{
    if (v.isDouble()) return v.toDouble();
    return std::nullopt;
}

QVector3D vec3(const QJsonValue& v)
{
    const QJsonArray a = v.toArray();
    if (a.size() != 3) return {};
    return QVector3D(float(a[0].toDouble()), float(a[1].toDouble()), float(a[2].toDouble()));
}

QJsonArray vec3Json(const QVector3D& v) { return {v.x(), v.y(), v.z()}; }

Result<FieldDesc> parseField(const QJsonObject& o, bool needOffset)
{
    FieldDesc f;
    f.name = o.value("name").toString();
    f.unit = o.value("unit").toString();
    f.components = o.value("components").toInt(1);
    if (f.name.isEmpty()) return Result<FieldDesc>::failure("field without a name");
    if (f.components != 1 && f.components != 3)
        return Result<FieldDesc>::failure(QStringLiteral("field %1 has %2 components; expected 1 or 3")
                                              .arg(f.name).arg(f.components));
    if (needOffset) {
        const QJsonValue off = o.value("offset_bytes");
        if (!off.isDouble() || off.toDouble() < 0)
            return Result<FieldDesc>::failure(QStringLiteral("field %1 has no offset").arg(f.name));
        f.offsetBytes = qint64(off.toDouble());
    }
    f.min = optionalNumber(o.value("min"));
    f.max = optionalNumber(o.value("max"));
    return f;
}

Result<ViewPart> parsePart(const QJsonObject& o, PartKind kind, int index)
{
    ViewPart p;
    p.kind = kind;
    p.id = o.value("id").toInt(-1);
    p.type = o.value("type").toString();
    p.file = o.value("file").toString();
    p.vertices = quint32(o.value("vertices").toDouble());
    p.indices = quint32(o.value(kind == PartKind::Streamlines ? "polylines" : "indices").toDouble());
    if (p.file.isEmpty() || p.file.contains('/') || p.file.contains('\\') || p.file.startsWith('.'))
        return Result<ViewPart>::failure(QStringLiteral("unsafe or missing file name '%1'").arg(p.file));
    if (kind == PartKind::Slice) {
        p.origin = vec3(o.value("origin"));
        p.normal = vec3(o.value("normal"));
        p.name = o.value("name").toString();
        if (p.name.isEmpty()) p.name = QStringLiteral("slice-%1").arg(index);
    } else {
        p.name = o.value("name").toString();
    }
    for (const QJsonValue& fv : o.value("fields").toArray()) {
        auto f = parseField(fv.toObject(), true);
        if (!f) return Result<ViewPart>::failure(p.file + ": " + f.error());
        p.fields.push_back(std::move(f).value());
    }
    return p;
}

template <class T>
void appendLE(QByteArray& out, T v)
{
    T le = qToLittleEndian(v);
    out.append(reinterpret_cast<const char*>(&le), sizeof(T));
}

} // namespace

void Bounds::expand(const QVector3D& p)
{
    if (!valid) { min = max = p; valid = true; return; }
    min = QVector3D(std::min(min.x(), p.x()), std::min(min.y(), p.y()), std::min(min.z(), p.z()));
    max = QVector3D(std::max(max.x(), p.x()), std::max(max.y(), p.y()), std::max(max.z(), p.z()));
}

void Bounds::expand(const Bounds& b)
{
    if (!b.valid) return;
    expand(b.min);
    expand(b.max);
}

const FieldDesc* ViewPart::field(const QString& n) const
{
    for (const auto& f : fields)
        if (f.name == n) return &f;
    return nullptr;
}

const FieldDesc* ViewManifest::field(const QString& n) const
{
    for (const auto& f : fields)
        if (f.name == n) return &f;
    return nullptr;
}

const MeshData::Field* MeshData::field(const QString& n) const
{
    for (const auto& f : fields)
        if (f.desc.name == n) return &f;
    return nullptr;
}

Bounds MeshData::bounds() const
{
    Bounds b;
    for (std::size_t i = 0; i + 2 < positions.size(); i += 3)
        b.expand(QVector3D(positions[i], positions[i + 1], positions[i + 2]));
    return b;
}

Result<ViewManifest> parseViewManifest(const QJsonObject& json, const QString& dir)
{
    if (json.value("v").toInt() != 1)
        return Result<ViewManifest>::failure(QStringLiteral("unsupported view manifest version %1")
                                                 .arg(json.value("v").toInt()));
    const QJsonObject format = json.value("format").toObject();
    if (format.value("version").toInt(1) != 1 || format.value("header_bytes").toInt(32) != int(kHeaderBytes))
        return Result<ViewManifest>::failure("unsupported view binary format");

    ViewManifest m;
    m.dir = dir;
    m.reference = json.value("reference").toBool();
    m.provenance = json.value("provenance").toString();
    m.synthetic = json.value("synthetic").toBool();
    const QJsonArray b = json.value("bounds").toArray();
    if (b.size() == 2) {
        m.bounds.min = vec3(b[0]);
        m.bounds.max = vec3(b[1]);
        m.bounds.valid = true;
    }
    int i = 0;
    for (const QJsonValue& z : json.value("zones").toArray()) {
        auto p = parsePart(z.toObject(), PartKind::Zone, i++);
        if (!p) return Result<ViewManifest>::failure(p.error());
        m.zones.push_back(std::move(p).value());
    }
    i = 0;
    for (const QJsonValue& s : json.value("slices").toArray()) {
        auto p = parsePart(s.toObject(), PartKind::Slice, i++);
        if (!p) return Result<ViewManifest>::failure(p.error());
        m.slices.push_back(std::move(p).value());
    }
    for (const QJsonValue& fv : json.value("fields").toArray()) {
        auto f = parseField(fv.toObject(), false);
        if (!f) return Result<ViewManifest>::failure(f.error());
        m.fields.push_back(std::move(f).value());
    }
    return m;
}

Result<ViewManifest> loadViewManifest(const QString& manifestPath)
{
    QFile f(manifestPath);
    if (!f.open(QIODevice::ReadOnly))
        return Result<ViewManifest>::failure(QStringLiteral("cannot open %1").arg(manifestPath));
    QJsonParseError err;
    const QJsonDocument doc = QJsonDocument::fromJson(f.readAll(), &err);
    if (!doc.isObject())
        return Result<ViewManifest>::failure(QStringLiteral("%1: %2").arg(manifestPath, err.errorString()));
    return parseViewManifest(doc.object(), QFileInfo(manifestPath).absolutePath());
}

Result<MeshData> readMeshBytes(const QByteArray& bytes, const ViewPart& part)
{
    using R = Result<MeshData>;
    const auto size = quint64(bytes.size());
    if (size < kHeaderBytes) return R::failure(part.file + ": shorter than its header");
    const char* d = bytes.constData();
    const bool lines = std::memcmp(d, kLineMagic, 4) == 0;
    if (!lines && std::memcmp(d, kTriMagic, 4) != 0) return R::failure(part.file + ": not a Kutta view file");
    auto u32 = [d](int off) { return qFromLittleEndian<quint32>(d + off); };
    if (u32(4) != 1) return R::failure(part.file + ": unsupported version");
    if (u32(8) != kHeaderBytes) return R::failure(part.file + ": unexpected header size");
    const quint64 V = u32(12), I = u32(16), F = u32(20);
    if (V != part.vertices || I != part.indices)
        return R::failure(QStringLiteral("%1: counts (%2, %3) disagree with the manifest (%4, %5)")
                              .arg(part.file).arg(V).arg(I).arg(part.vertices).arg(part.indices));
    if (F != part.fields.size()) return R::failure(part.file + ": field count disagrees with the manifest");
    const quint64 indexCount = lines ? I + 1 : I;
    if (!lines && I % 3 != 0) return R::failure(part.file + ": index count is not a multiple of 3");
    const quint64 posEnd = kHeaderBytes + 12 * V;
    const quint64 idxEnd = posEnd + 4 * indexCount;
    if (idxEnd > size) return R::failure(part.file + ": truncated geometry");

    MeshData m;
    m.polylines = lines;
    m.positions.resize(3 * V);
    m.indices.resize(indexCount);
    // Little-endian hosts (every supported target) copy directly.
    static_assert(Q_BYTE_ORDER == Q_LITTLE_ENDIAN, "view files are little endian");
    std::memcpy(m.positions.data(), d + kHeaderBytes, 12 * V);
    std::memcpy(m.indices.data(), d + posEnd, 4 * indexCount);
    if (lines) {
        if (m.indices.front() != 0 || m.indices.back() != V)
            return R::failure(part.file + ": polyline offsets do not span the vertices");
        for (quint64 k = 1; k < indexCount; ++k)
            if (m.indices[k] < m.indices[k - 1]) return R::failure(part.file + ": polyline offsets decrease");
    } else {
        for (auto idx : m.indices)
            if (idx >= V) return R::failure(part.file + ": triangle index out of range");
    }
    quint64 expectedOffset = idxEnd;
    for (const FieldDesc& fd : part.fields) {
        const quint64 n = V * quint64(fd.components);
        if (quint64(fd.offsetBytes) != expectedOffset || expectedOffset + 4 * n > size)
            return R::failure(QStringLiteral("%1: field %2 is outside the file").arg(part.file, fd.name));
        MeshData::Field f{fd, std::vector<float>(n)};
        std::memcpy(f.values.data(), d + expectedOffset, 4 * n);
        expectedOffset += 4 * n;
        m.fields.push_back(std::move(f));
    }
    if (expectedOffset != size) return R::failure(part.file + ": trailing bytes after the last field");
    return m;
}

Result<MeshData> readMeshFile(const QString& path, const ViewPart& part)
{
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) return Result<MeshData>::failure(QStringLiteral("cannot open %1").arg(path));
    return readMeshBytes(f.readAll(), part);
}

QByteArray writeMeshBytes(const MeshData& mesh, ViewPart* part)
{
    const quint32 V = quint32(mesh.vertexCount());
    const quint32 I = mesh.polylines ? quint32(mesh.indices.size() - 1) : quint32(mesh.indices.size());
    QByteArray out;
    out.append(mesh.polylines ? kLineMagic : kTriMagic, 4);
    appendLE<quint32>(out, 1);
    appendLE<quint32>(out, kHeaderBytes);
    appendLE<quint32>(out, V);
    appendLE<quint32>(out, I);
    appendLE<quint32>(out, quint32(mesh.fields.size()));
    appendLE<quint32>(out, 0);
    appendLE<quint32>(out, 0);
    out.append(reinterpret_cast<const char*>(mesh.positions.data()), qsizetype(mesh.positions.size() * 4));
    out.append(reinterpret_cast<const char*>(mesh.indices.data()), qsizetype(mesh.indices.size() * 4));
    part->vertices = V;
    part->indices = I;
    part->fields.clear();
    for (const auto& f : mesh.fields) {
        FieldDesc d = f.desc;
        d.offsetBytes = out.size();
        if (auto r = finiteRange(f)) { d.min = r->first; d.max = r->second; }
        else { d.min.reset(); d.max.reset(); }
        part->fields.push_back(d);
        out.append(reinterpret_cast<const char*>(f.values.data()), qsizetype(f.values.size() * 4));
    }
    return out;
}

QJsonObject toJson(const FieldDesc& f, bool withOffset)
{
    QJsonObject o{{"name", f.name}, {"unit", f.unit}, {"components", f.components}, {"location", "vertex"}};
    if (withOffset) o["offset_bytes"] = double(f.offsetBytes);
    o["min"] = f.min ? QJsonValue(*f.min) : QJsonValue();
    o["max"] = f.max ? QJsonValue(*f.max) : QJsonValue();
    return o;
}

QJsonObject toJson(const ViewPart& p)
{
    QJsonObject o{{"file", p.file}, {"vertices", double(p.vertices)}};
    o[p.kind == PartKind::Streamlines ? "polylines" : "indices"] = double(p.indices);
    if (p.kind == PartKind::Zone) {
        o["id"] = p.id;
        o["name"] = p.name;
        o["type"] = p.type;
    } else {
        o["name"] = p.name;
    }
    if (p.kind == PartKind::Slice) {
        o["origin"] = vec3Json(p.origin);
        o["normal"] = vec3Json(p.normal);
    }
    QJsonArray fields;
    for (const auto& f : p.fields) fields.append(toJson(f, true));
    o["fields"] = fields;
    return o;
}

std::optional<std::pair<double, double>> finiteRange(const MeshData::Field& f)
{
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    const int c = f.desc.components;
    for (std::size_t i = 0; i + c <= f.values.size(); i += c) {
        double v = f.values[i];
        if (c == 3) v = std::sqrt(double(f.values[i]) * f.values[i] + double(f.values[i + 1]) * f.values[i + 1] +
                                  double(f.values[i + 2]) * f.values[i + 2]);
        if (!std::isfinite(v)) continue;
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
    if (lo > hi) return std::nullopt;
    return std::pair{lo, hi};
}

} // namespace kapp
