#include "core/Stl.h"

#include <QFile>
#include <QFileInfo>
#include <QSaveFile>
#include <QTextStream>

#include <cmath>
#include <cstring>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

namespace kapp {

double unitScale(LengthUnit u) { return u == LengthUnit::Millimetre ? 1e-3 : 1.0; }
QString unitName(LengthUnit u) { return u == LengthUnit::Millimetre ? QStringLiteral("mm") : QStringLiteral("m"); }
LengthUnit unitFromName(const QString& s) { return s == QLatin1String("mm") ? LengthUnit::Millimetre : LengthUnit::Metre; }

LengthUnit guessUnit(const Bounds& b)
{
    if (!b.valid) return LengthUnit::Metre;
    const QVector3D e = b.max - b.min;
    return std::max({e.x(), e.y(), e.z()}) > 20.0f ? LengthUnit::Millimetre : LengthUnit::Metre;
}

namespace {

// Welds exactly coincident vertices (CAD exports repeat each vertex per facet).
class Welder {
public:
    explicit Welder(MeshData& m) : m_mesh(m) {}
    std::uint32_t add(float x, float y, float z)
    {
        Key k{bits(x), bits(y), bits(z)};
        auto [it, inserted] = m_index.try_emplace(k, std::uint32_t(m_mesh.positions.size() / 3));
        if (inserted) m_mesh.positions.insert(m_mesh.positions.end(), {x, y, z});
        return it->second;
    }

private:
    static std::uint32_t bits(float f)
    {
        if (f == 0.0f) f = 0.0f;  // -0 and +0 weld
        std::uint32_t u;
        std::memcpy(&u, &f, 4);
        return u;
    }
    struct Key {
        std::uint32_t x, y, z;
        bool operator==(const Key& o) const { return x == o.x && y == o.y && z == o.z; }
    };
    struct Hash {
        std::size_t operator()(const Key& k) const { return (std::size_t(k.x) * 73856093u) ^ (std::size_t(k.y) * 19349663u) ^ (std::size_t(k.z) * 83492791u); }
    };
    MeshData& m_mesh;
    std::unordered_map<Key, std::uint32_t, Hash> m_index;
};

// Gives each smooth patch around a vertex its own copy where faces meet at
// more than 40° (wing endplates, trailing edges, box corners), so shading
// shows the edge instead of averaging the normals across it.
void splitCreases(MeshData& m)
{
    const auto& I = m.indices;
    const std::size_t nv = m.positions.size() / 3, nt = I.size() / 3;
    std::vector<QVector3D> faceNormal(nt);
    for (std::size_t t = 0; t < nt; ++t) {
        auto p = [&](std::uint32_t i) { return QVector3D(m.positions[3 * i], m.positions[3 * i + 1], m.positions[3 * i + 2]); };
        faceNormal[t] = QVector3D::crossProduct(p(I[3 * t + 1]) - p(I[3 * t]), p(I[3 * t + 2]) - p(I[3 * t])).normalized();
    }
    std::vector<std::uint32_t> start(nv + 1, 0), corners(I.size());
    for (std::uint32_t i : I) ++start[i + 1];
    for (std::size_t v = 0; v < nv; ++v) start[v + 1] += start[v];
    std::vector<std::uint32_t> fill(start.begin(), start.end() - 1);
    for (std::size_t c = 0; c < I.size(); ++c) corners[fill[I[c]]++] = std::uint32_t(c);
    const float cosCrease = std::cos(40.0f * 3.14159265f / 180.0f);
    struct Group { QVector3D sum; std::uint32_t index; };
    std::vector<Group> groups;
    for (std::size_t v = 0; v < nv; ++v) {
        groups.clear();
        for (std::uint32_t k = start[v]; k < start[v + 1]; ++k) {
            const std::uint32_t c = corners[k];
            const QVector3D n = faceNormal[c / 3];
            Group* g = nullptr;
            for (auto& x : groups)
                if (QVector3D::dotProduct(x.sum.normalized(), n) >= cosCrease) { g = &x; break; }
            if (!g) {
                std::uint32_t index = std::uint32_t(v);
                if (!groups.empty()) {
                    index = std::uint32_t(m.positions.size() / 3);
                    m.positions.insert(m.positions.end(), {m.positions[3 * v], m.positions[3 * v + 1], m.positions[3 * v + 2]});
                }
                groups.push_back({QVector3D(), index});
                g = &groups.back();
            }
            g->sum += n;
            m.indices[c] = g->index;
        }
    }
}

void finish(StlSolid& s)
{
    splitCreases(s.mesh);
    s.triangles = s.mesh.indices.size() / 3;
    const auto& P = s.mesh.positions;
    double area = 0;
    for (std::size_t t = 0; t + 2 < s.mesh.indices.size(); t += 3) {
        const std::uint32_t a = s.mesh.indices[t], b = s.mesh.indices[t + 1], c = s.mesh.indices[t + 2];
        const QVector3D pa(P[3 * a], P[3 * a + 1], P[3 * a + 2]), pb(P[3 * b], P[3 * b + 1], P[3 * b + 2]), pc(P[3 * c], P[3 * c + 1], P[3 * c + 2]);
        area += 0.5 * double(QVector3D::crossProduct(pb - pa, pc - pa).length());
    }
    s.area = area;
}

bool degenerate(const float* v)
{
    auto same = [](const float* a, const float* b) { return a[0] == b[0] && a[1] == b[1] && a[2] == b[2]; };
    return same(v, v + 3) || same(v + 3, v + 6) || same(v, v + 6);
}

Result<StlFile> readBinary(const QByteArray& bytes, const QString& path, float scale)
{
    using R = Result<StlFile>;
    quint32 count = 0;
    std::memcpy(&count, bytes.constData() + 80, 4);
    if (qint64(count) * 50 + 84 != bytes.size())
        return R::failure(QStringLiteral("%1 is not a complete binary STL (%2 triangles declared, %3 bytes)")
                              .arg(QFileInfo(path).fileName()).arg(count).arg(bytes.size()));
    StlFile f;
    f.path = path;
    f.binary = true;
    StlSolid s;
    s.name = QFileInfo(path).completeBaseName();
    s.mesh.indices.reserve(std::size_t(count) * 3);
    Welder w(s.mesh);
    const char* p = bytes.constData() + 84;
    for (quint32 i = 0; i < count; ++i, p += 50) {
        float v[9];
        std::memcpy(v, p + 12, 36);
        for (float& x : v) {
            if (!std::isfinite(x)) return R::failure(QStringLiteral("%1: triangle %2 has a non-finite coordinate").arg(QFileInfo(path).fileName()).arg(i + 1));
            x *= scale;
        }
        if (degenerate(v)) continue;
        for (int k = 0; k < 3; ++k) s.mesh.indices.push_back(w.add(v[3 * k], v[3 * k + 1], v[3 * k + 2]));
    }
    finish(s);
    f.solids.push_back(std::move(s));
    return f;
}

Result<StlFile> readAscii(const QByteArray& bytes, const QString& path, float scale)
{
    using R = Result<StlFile>;
    StlFile f;
    f.path = path;
    const QString file = QFileInfo(path).fileName();
    std::optional<StlSolid> cur;
    std::unique_ptr<Welder> welder;
    float tri[9];
    int nv = 0;
    int lineNo = 0;
    qsizetype pos = 0;
    while (pos < bytes.size()) {
        qsizetype nl = bytes.indexOf('\n', pos);
        if (nl < 0) nl = bytes.size();
        const QByteArray line = bytes.mid(pos, nl - pos).trimmed();
        pos = nl + 1;
        ++lineNo;
        if (line.isEmpty()) continue;
        if (line.startsWith("solid")) {
            if (cur) return R::failure(QStringLiteral("%1 line %2: solid inside solid %3").arg(file).arg(lineNo).arg(cur->name));
            cur.emplace();
            cur->name = QString::fromUtf8(line.mid(5).trimmed());
            if (cur->name.isEmpty()) cur->name = QFileInfo(path).completeBaseName() + (f.solids.empty() ? QString() : QString::number(f.solids.size() + 1));
            welder = std::make_unique<Welder>(cur->mesh);
        } else if (line.startsWith("vertex")) {
            if (!cur) return R::failure(QStringLiteral("%1 line %2: vertex outside a solid").arg(file).arg(lineNo));
            const QList<QByteArray> parts = line.simplified().split(' ');
            if (parts.size() != 4 || nv >= 3) return R::failure(QStringLiteral("%1 line %2: expected \"vertex x y z\"").arg(file).arg(lineNo));
            for (int k = 0; k < 3; ++k) {
                bool ok = false;
                const float x = parts[k + 1].toFloat(&ok);
                if (!ok || !std::isfinite(x)) return R::failure(QStringLiteral("%1 line %2: unreadable coordinate").arg(file).arg(lineNo));
                tri[3 * nv + k] = x * scale;
            }
            ++nv;
        } else if (line.startsWith("endloop")) {
            if (nv != 3) return R::failure(QStringLiteral("%1 line %2: a facet needs exactly three vertices").arg(file).arg(lineNo));
            if (!degenerate(tri))
                for (int k = 0; k < 3; ++k) cur->mesh.indices.push_back(welder->add(tri[3 * k], tri[3 * k + 1], tri[3 * k + 2]));
            nv = 0;
        } else if (line.startsWith("endsolid")) {
            if (!cur) return R::failure(QStringLiteral("%1 line %2: endsolid without solid").arg(file).arg(lineNo));
            finish(*cur);
            f.solids.push_back(std::move(*cur));
            cur.reset();
            welder.reset();
        } else if (!line.startsWith("facet") && !line.startsWith("outer") && !line.startsWith("endfacet")) {
            return R::failure(QStringLiteral("%1 line %2: unexpected \"%3\"").arg(file).arg(lineNo).arg(QString::fromUtf8(line.left(40))));
        }
    }
    if (cur) {  // tolerate a missing final endsolid
        finish(*cur);
        f.solids.push_back(std::move(*cur));
    }
    if (f.solids.empty()) return R::failure(file + ": no solids");
    return f;
}

} // namespace

Result<StlFile> readStl(const QString& path, LengthUnit unit)
{
    using R = Result<StlFile>;
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) return R::failure(QStringLiteral("cannot read %1: %2").arg(path, file.errorString()));
    const QByteArray bytes = file.readAll();
    if (bytes.size() < 15) return R::failure(QFileInfo(path).fileName() + " is too short to be an STL file");
    const float scale = float(unitScale(unit));
    // Binary files may begin with "solid" too; the declared size decides.
    bool binary = !bytes.startsWith("solid");
    if (!binary && bytes.size() >= 84) {
        quint32 count = 0;
        std::memcpy(&count, bytes.constData() + 80, 4);
        if (qint64(count) * 50 + 84 == bytes.size() && !bytes.left(qMin<qsizetype>(bytes.size(), 512)).contains("facet")) binary = true;
    }
    auto r = binary ? readBinary(bytes, path, scale) : readAscii(bytes, path, scale);
    if (!r) return r;
    std::size_t total = 0;
    for (const auto& s : r.value().solids) total += s.triangles;
    if (total == 0) return R::failure(QFileInfo(path).fileName() + " has no triangles");
    return r;
}

void rescale(StlFile& f, double factor)
{
    for (auto& s : f.solids) {
        for (float& x : s.mesh.positions) x = float(x * factor);
        s.area *= factor * factor;
    }
}

Bounds stlBounds(const StlFile& f)
{
    Bounds b;
    for (const auto& s : f.solids) b.expand(s.mesh.bounds());
    return b;
}

bool writeAsciiStl(const QString& path, const std::vector<std::pair<QString, const MeshData*>>& solids, LengthUnit unit, QString* error)
{
    QSaveFile f(path);
    if (!f.open(QIODevice::WriteOnly)) {
        if (error) *error = f.errorString();
        return false;
    }
    QTextStream out(&f);
    out.setRealNumberPrecision(9);
    const double inv = 1.0 / unitScale(unit);
    for (const auto& [name, mesh] : solids) {
        out << "solid " << name << '\n';
        const auto& P = mesh->positions;
        for (std::size_t t = 0; t + 2 < mesh->indices.size(); t += 3) {
            QVector3D v[3];
            for (int k = 0; k < 3; ++k) {
                const std::uint32_t i = mesh->indices[t + std::size_t(k)];
                v[k] = QVector3D(P[3 * i], P[3 * i + 1], P[3 * i + 2]) * float(inv);
            }
            const QVector3D n = QVector3D::crossProduct(v[1] - v[0], v[2] - v[0]).normalized();
            out << " facet normal " << n.x() << ' ' << n.y() << ' ' << n.z() << "\n  outer loop\n";
            for (const auto& p : v) out << "   vertex " << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
            out << "  endloop\n endfacet\n";
        }
        out << "endsolid " << name << '\n';
    }
    out.flush();
    if (!f.commit()) {
        if (error) *error = f.errorString();
        return false;
    }
    return true;
}

} // namespace kapp
