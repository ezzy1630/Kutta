#include "core/SceneData.h"

#include <QDir>
#include <QFile>
#include <QJsonArray>
#include <QLocale>

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

namespace kapp {

PreparedMesh::PreparedMesh(MeshData data, bool flatNormal, QVector3D planeNormal) : m_data(std::move(data))
{
    const std::size_t V = m_data.vertexCount();
    m_pn.assign(V * 6, 0.0f);
    for (std::size_t i = 0; i < V; ++i) {
        m_pn[6 * i] = m_data.positions[3 * i];
        m_pn[6 * i + 1] = m_data.positions[3 * i + 1];
        m_pn[6 * i + 2] = m_data.positions[3 * i + 2];
    }
    m_bounds = m_data.bounds();
    if (m_data.polylines) return;
    if (flatNormal) {
        const QVector3D n = planeNormal.normalized();
        for (std::size_t i = 0; i < V; ++i) {
            m_pn[6 * i + 3] = n.x();
            m_pn[6 * i + 4] = n.y();
            m_pn[6 * i + 5] = n.z();
        }
        return;
    }
    // Area-weighted vertex normals.
    const auto& P = m_data.positions;
    for (std::size_t t = 0; t + 2 < m_data.indices.size(); t += 3) {
        const std::uint32_t a = m_data.indices[t], b = m_data.indices[t + 1], c = m_data.indices[t + 2];
        const QVector3D pa(P[3 * a], P[3 * a + 1], P[3 * a + 2]), pb(P[3 * b], P[3 * b + 1], P[3 * b + 2]),
            pc(P[3 * c], P[3 * c + 1], P[3 * c + 2]);
        const QVector3D n = QVector3D::crossProduct(pb - pa, pc - pa);
        for (std::uint32_t v : {a, b, c}) {
            m_pn[6 * v + 3] += n.x();
            m_pn[6 * v + 4] += n.y();
            m_pn[6 * v + 5] += n.z();
        }
    }
    for (std::size_t i = 0; i < V; ++i) {
        QVector3D n(m_pn[6 * i + 3], m_pn[6 * i + 4], m_pn[6 * i + 5]);
        n = n.isNull() ? QVector3D(0, 0, 1) : n.normalized();
        m_pn[6 * i + 3] = n.x();
        m_pn[6 * i + 4] = n.y();
        m_pn[6 * i + 5] = n.z();
    }
}

bool PreparedMesh::hasField(const FieldRef& f) const
{
    const auto* fd = m_data.field(f.manifestName());
    if (!fd) return false;
    return f.component < 0 || fd->desc.components == 3;
}

QStringList PreparedMesh::fieldNames() const
{
    QStringList out;
    for (const auto& f : m_data.fields) out << f.desc.name;
    return out;
}

std::shared_ptr<const std::vector<float>> PreparedMesh::scalar(const FieldRef& f) const
{
    const std::string key = f.key().toStdString();
    std::lock_guard lock(m_mutex);
    if (auto it = m_scalars.find(key); it != m_scalars.end()) return it->second;
    const auto* fd = m_data.field(f.manifestName());
    if (!fd) return nullptr;
    const std::size_t V = m_data.vertexCount();
    auto out = std::make_shared<std::vector<float>>(V);
    const int c = fd->desc.components;
    for (std::size_t i = 0; i < V; ++i) {
        float v;
        if (c == 1) v = fd->values[i];
        else if (f.component >= 0) v = fd->values[3 * i + std::size_t(f.component)];
        else {
            const float x = fd->values[3 * i], y = fd->values[3 * i + 1], z = fd->values[3 * i + 2];
            v = std::sqrt(x * x + y * y + z * z);
        }
        (*out)[i] = std::isfinite(v) ? v : kMissingSentinel;
    }
    m_scalars.emplace(key, out);
    return out;
}

std::optional<std::pair<double, double>> PreparedMesh::scalarRange(const FieldRef& f) const
{
    auto s = scalar(f);
    if (!s) return std::nullopt;
    double lo = INFINITY, hi = -INFINITY;
    for (float v : *s) {
        if (v >= 1e37f) continue;
        lo = std::min(lo, double(v));
        hi = std::max(hi, double(v));
    }
    if (lo > hi) return std::nullopt;
    return std::pair{lo, hi};
}

void PreparedMesh::setFaceEdges(std::vector<std::uint32_t> lines)
{
    std::lock_guard lock(m_mutex);
    m_edges = std::move(lines);
    m_edgesBuilt = true;
    m_faceEdges = true;
}

const std::vector<std::uint32_t>& PreparedMesh::edges() const
{
    std::lock_guard lock(m_mutex);
    if (m_edgesBuilt) return m_edges;
    m_edgesBuilt = true;
    if (m_data.polylines) return m_edges;
    std::unordered_set<std::uint64_t> seen;
    seen.reserve(m_data.indices.size());
    for (std::size_t t = 0; t + 2 < m_data.indices.size(); t += 3) {
        for (int e = 0; e < 3; ++e) {
            std::uint32_t a = m_data.indices[t + std::size_t(e)], b = m_data.indices[t + std::size_t((e + 1) % 3)];
            if (a > b) std::swap(a, b);
            if (seen.insert((std::uint64_t(a) << 32) | b).second) {
                m_edges.push_back(a);
                m_edges.push_back(b);
            }
        }
    }
    return m_edges;
}

std::optional<double> PreparedMesh::sampleNear(const QVector3D& p, const FieldRef& f) const
{
    auto s = scalar(f);
    if (!s || m_data.polylines) return std::nullopt;
    const auto& P = m_data.positions;
    double best = INFINITY;
    std::optional<double> value;
    for (std::size_t t = 0; t + 2 < m_data.indices.size(); t += 3) {
        const std::uint32_t ia = m_data.indices[t], ib = m_data.indices[t + 1], ic = m_data.indices[t + 2];
        const QVector3D a(P[3 * ia], P[3 * ia + 1], P[3 * ia + 2]), b(P[3 * ib], P[3 * ib + 1], P[3 * ib + 2]),
            c(P[3 * ic], P[3 * ic + 1], P[3 * ic + 2]);
        const QVector3D centroid = (a + b + c) / 3.0f;
        const double d = (centroid - p).lengthSquared();
        if (d >= best) continue;
        // Barycentric weights of the projection onto the triangle plane.
        const QVector3D v0 = b - a, v1 = c - a, v2 = p - a;
        const float d00 = QVector3D::dotProduct(v0, v0), d01 = QVector3D::dotProduct(v0, v1),
                    d11 = QVector3D::dotProduct(v1, v1), d20 = QVector3D::dotProduct(v2, v0),
                    d21 = QVector3D::dotProduct(v2, v1);
        const float den = d00 * d11 - d01 * d01;
        if (std::abs(den) < 1e-30f) continue;
        float w1 = (d11 * d20 - d01 * d21) / den, w2 = (d00 * d21 - d01 * d20) / den;
        w1 = std::clamp(w1, 0.0f, 1.0f);
        w2 = std::clamp(w2, 0.0f, 1.0f - w1);
        const float w0 = 1 - w1 - w2;
        const float va = (*s)[ia], vb = (*s)[ib], vc = (*s)[ic];
        if (va >= 1e37f || vb >= 1e37f || vc >= 1e37f) continue;
        best = d;
        value = double(w0 * va + w1 * vb + w2 * vc);
    }
    return value;
}

const SceneItem* SceneData::item(const QString& id) const
{
    for (const auto& i : items)
        if (i.id == id) return &i;
    return nullptr;
}

std::vector<FieldRef> SceneData::availableFields() const
{
    std::vector<FieldRef> out;
    auto add = [&out](const FieldRef& f) {
        for (const auto& e : out)
            if (e == f) return;
        out.push_back(f);
    };
    for (const auto& it : items) {
        if (!it.mesh) continue;
        for (const auto& f : it.mesh->data().fields) {
            FieldRef r = FieldRef::parse(f.desc.name);
            add(r);
            if (f.desc.components == 3)
                for (int c = 0; c < 3; ++c) { FieldRef rc = r; rc.component = c; add(rc); }
        }
    }
    std::stable_sort(out.begin(), out.end(), [](const FieldRef& a, const FieldRef& b) {
        if (a.role != b.role) return int(a.role) < int(b.role);
        if (a.quantity != b.quantity) return quantityOrder(a.quantity) < quantityOrder(b.quantity);
        return a.component < b.component;
    });
    return out;
}

std::optional<std::pair<double, double>> SceneData::range(const FieldRef& f, const QStringList& visibleIds) const
{
    double lo = INFINITY, hi = -INFINITY;
    for (const auto& it : items) {
        if (!it.mesh || (!visibleIds.isEmpty() && !visibleIds.contains(it.id))) continue;
        if (auto r = it.mesh->scalarRange(f)) {
            lo = std::min(lo, r->first);
            hi = std::max(hi, r->second);
        }
    }
    if (lo > hi) return std::nullopt;
    return std::pair{lo, hi};
}

QString SceneData::unitOf(const QString& quantity) const
{
    for (const auto& it : items) {
        if (!it.mesh) continue;
        for (const auto& f : it.mesh->data().fields)
            if (FieldRef::parse(f.desc.name).quantity == quantity) return f.desc.unit;
    }
    return {};
}

QString itemLabelForSlice(const ViewPart& s)
{
    const QLocale c = QLocale::c();
    const QVector3D n = s.normal.normalized();
    auto axisPlane = [&](int axis, const char* name) -> QString {
        const float comps[3] = {n.x(), n.y(), n.z()};
        const float o[3] = {s.origin.x(), s.origin.y(), s.origin.z()};
        if (std::abs(std::abs(comps[axis]) - 1) < 1e-4f)
            return QStringLiteral("%1 = %2 m").arg(name, c.toString(o[axis], 'g', 4));
        return {};
    };
    for (int a = 0; a < 3; ++a)
        if (QString l = axisPlane(a, a == 0 ? "x" : a == 1 ? "y" : "z"); !l.isEmpty()) return l;
    return QStringLiteral("plane through (%1, %2, %3)")
        .arg(c.toString(s.origin.x(), 'g', 3), c.toString(s.origin.y(), 'g', 3), c.toString(s.origin.z(), 'g', 3));
}

Result<std::shared_ptr<SceneData>> loadScene(const ViewManifest& m, const QString& label)
{
    using R = Result<std::shared_ptr<SceneData>>;
    auto scene = std::make_shared<SceneData>();
    scene->label = label;
    scene->provenance = m.provenance;
    scene->synthetic = m.synthetic;
    scene->reference = m.reference;
    scene->geometryOnly = m.geometryOnly;
    const QDir dir(m.dir);
    int zoneIndex = 0;
    for (const auto& z : m.zones) {
        auto mesh = readMeshFile(dir.filePath(z.file), z);
        if (!mesh) return R::failure(mesh.error());
        SceneItem it;
        it.id = "zone:" + z.name;
        it.label = z.name;
        it.kind = ItemKind::Zone;
        it.zoneType = z.type;
        it.color = zoneColor(zoneIndex++);
        // Walls first: inlets, outlets and symmetry planes usually hide the body.
        it.defaultVisible = z.type == "wall";
        auto prepared = std::make_shared<PreparedMesh>(std::move(mesh).value());
        // Polygon edges are optional: without them (or if unreadable) the
        // viewer falls back to triangle edges rather than failing the view.
        if (z.edges) {
            QFile f(dir.filePath(z.edges->file));
            if (f.open(QIODevice::ReadOnly))
                if (auto lines = readEdgeLines(f.readAll(), *z.edges, z.vertices)) prepared->setFaceEdges(std::move(lines).value());
        }
        it.mesh = std::move(prepared);
        scene->bounds.expand(it.mesh->bounds());
        scene->items.push_back(std::move(it));
    }
    int n = 0;
    for (const auto& s : m.slices) {
        auto mesh = readMeshFile(dir.filePath(s.file), s);
        if (!mesh) return R::failure(mesh.error());
        SceneItem it;
        it.id = QStringLiteral("slice:%1").arg(n++);
        it.label = itemLabelForSlice(s);
        it.kind = ItemKind::Slice;
        it.mesh = std::make_shared<PreparedMesh>(std::move(mesh).value(), true, s.normal);
        scene->items.push_back(std::move(it));
    }
    auto addPost = [&](const std::vector<ViewPart>& parts, ItemKind kind) -> QString {
        for (const auto& p : parts) {
            auto mesh = readMeshFile(dir.filePath(p.file), p);
            if (!mesh) return mesh.error();
            SceneItem it;
            it.id = "post:" + p.name;
            it.label = p.name;
            it.kind = kind;
            it.mesh = std::make_shared<PreparedMesh>(std::move(mesh).value());
            scene->items.push_back(std::move(it));
        }
        return {};
    };
    for (const QString& err : {addPost(m.isosurfaces, ItemKind::IsoSurface), addPost(m.streamlines, ItemKind::Streamlines),
                               addPost(m.probes, ItemKind::Probe)})
        if (!err.isEmpty()) return R::failure(err);
    if (!scene->bounds.valid) scene->bounds = m.bounds;
    return scene;
}

QColor zoneColor(int index)
{
    // Okabe–Ito hues, softened for shaded 3D surfaces, then lighter variants.
    static const QColor base[] = {{86, 158, 214}, {230, 159, 0}, {0, 158, 115}, {204, 121, 167},
                                  {213, 94, 0}, {240, 228, 66}, {120, 190, 233}, {170, 140, 210}};
    constexpr int n = int(std::size(base));
    const QColor c = base[index % n];
    return index / n % 2 ? c.lighter(135) : c;
}

bool enclosesDomain(const SceneItem& it, const Bounds& domain)
{
    if (it.kind != ItemKind::Zone || !it.mesh || !domain.valid) return false;
    const Bounds b = it.mesh->bounds();
    const QVector3D e = domain.max - domain.min, z = b.max - b.min;
    // Flat within 0.1% of the domain: float noise or a slight slope in a
    // ground mesh must not make it a body.
    int span = 0, flat = 0;
    for (int a = 0; a < 3; ++a) {
        span += e[a] > 0 && z[a] >= 0.98f * e[a];
        flat += z[a] <= 1e-3f * e.length();
    }
    return span >= 2 && flat >= 1;
}

bool isEnclosingWall(const SceneItem& it, const SceneData& s)
{
    return it.kind == ItemKind::Zone && it.zoneType == "wall" && (enclosesDomain(it, s.bounds) || it.id.contains("ground"));
}

bool isExternalFlow(const SceneData& s)
{
    bool enclosure = false, body = false;
    for (const auto& it : s.items) {
        if (it.kind != ItemKind::Zone || it.zoneType != "wall") continue;
        (isEnclosingWall(it, s) ? enclosure : body) = true;
    }
    return enclosure && body;
}

QString planarFace(const SceneData& s)
{
    if (!s.planar || !s.bounds.valid) return {};
    QString best;
    float top = -std::numeric_limits<float>::infinity();
    const float depth = s.bounds.max.z() - s.bounds.min.z();
    for (const auto& it : s.items) {
        if (it.kind != ItemKind::Zone || it.zoneType != "symmetry" || !it.mesh) continue;
        const Bounds b = it.mesh->bounds();
        // A cap spans the plane and is flat through the depth.
        if (b.max.z() - b.min.z() > 1e-3f * std::max(depth, 1e-30f) || !enclosesDomain(it, s.bounds)) continue;
        if (b.max.z() > top) {
            top = b.max.z();
            best = it.id;
        }
    }
    return best;
}

void makePlanar(SceneData& s)
{
    s.planar = true;
    const QString face = planarFace(s);
    if (face.isEmpty()) return;
    const float depth = s.bounds.max.z() - s.bounds.min.z();
    std::erase_if(s.items, [&](const SceneItem& it) {
        if (it.id == face || it.kind != ItemKind::Zone || it.zoneType != "symmetry" || !it.mesh) return false;
        const Bounds b = it.mesh->bounds();
        return b.max.z() - b.min.z() <= 1e-3f * std::max(depth, 1e-30f) && enclosesDomain(it, s.bounds);
    });
    // The plane is the whole domain: the mesh reads best on the neutral
    // surface colour, as in Fluent, not flooded with a zone colour.
    for (auto& it : s.items)
        if (it.id == face) {
            it.label = QStringLiteral("2D domain");
            it.color = QColor();
        }
}

QStringList defaultMeshItems(const SceneData& s)
{
    if (const QString face = planarFace(s); !face.isEmpty()) return {face};
    QStringList body, all;
    for (const auto& it : s.items) {
        if (it.kind != ItemKind::Zone || !it.mesh) continue;
        all << it.id;
        if (it.zoneType == "symmetry" || it.zoneType == "interior") continue;
        if (enclosesDomain(it, s.bounds)) continue;
        body << it.id;
    }
    return body.isEmpty() ? all : body;
}

std::optional<int> commonPlaneAxis(const SceneData& s, const QStringList& ids)
{
    if (ids.isEmpty() || !s.bounds.valid) return std::nullopt;
    const float tol = 1e-3f * (s.bounds.max - s.bounds.min).length();
    std::optional<int> axis;
    for (const QString& id : ids) {
        const SceneItem* it = s.item(id);
        if (!it || !it->mesh) return std::nullopt;
        const Bounds b = it->mesh->bounds();
        const QVector3D e = b.max - b.min;
        std::optional<int> flat;
        for (int a = 0; a < 3; ++a)
            if (e[a] <= tol) flat = a;
        if (!flat || (axis && *axis != *flat)) return std::nullopt;
        axis = flat;
    }
    return axis;
}

std::optional<Bounds> planarFocus(const SceneData& s)
{
    if (!s.planar || !s.bounds.valid || !isExternalFlow(s)) return std::nullopt;
    Bounds body;
    for (const auto& it : s.items)
        if (it.kind == ItemKind::Zone && it.zoneType == "wall" && it.mesh && !isEnclosingWall(it, s)) body.expand(it.mesh->bounds());
    if (!body.valid) return std::nullopt;
    const QVector3D e = body.max - body.min;
    const float grow = std::max(e.x(), e.y());
    Bounds out;
    out.min = QVector3D(std::max(body.min.x() - grow, s.bounds.min.x()), std::max(body.min.y() - grow, s.bounds.min.y()), s.bounds.min.z());
    out.max = QVector3D(std::min(body.max.x() + grow, s.bounds.max.x()), std::min(body.max.y() + grow, s.bounds.max.y()), s.bounds.max.z());
    out.valid = true;
    return out;
}

QStringList defaultVisibleItems(const SceneData& s)
{
    // External aero: the body's walls, without the ground or tunnel. Internal
    // flow: the planes through it (a slice, else the symmetry plane), since
    // walls enclose the flow. A 2D case: its plane.
    if (const QString face = planarFace(s); !face.isEmpty()) return {face};
    QStringList walls, slices, symmetry;
    for (const auto& it : s.items) {
        if (it.kind == ItemKind::Slice) slices << it.id;
        if (it.kind != ItemKind::Zone) continue;
        if (it.zoneType == "wall" && !isEnclosingWall(it, s)) walls << it.id;
        else if (it.zoneType == "symmetry") symmetry << it.id;
    }
    if (isExternalFlow(s)) return walls;
    if (slices.isEmpty() && symmetry.isEmpty()) {
        QStringList all;
        for (const auto& it : s.items)
            if (it.kind == ItemKind::Zone && it.zoneType == "wall") all << it.id;
        return all;
    }
    return slices.isEmpty() ? symmetry : slices;
}

FieldRef defaultField(const SceneData& s)
{
    const auto fields = s.availableFields();
    const bool external = isExternalFlow(s);
    const QStringList prefs = external ? QStringList{"pressure_coefficient", "pressure", "velocity_magnitude"}
                                       : QStringList{"velocity_magnitude", "velocity", "pressure"};
    const FieldRole role = s.runDelta ? FieldRole::Run : FieldRole::Kutta;
    for (const QString& q : prefs)
        for (const auto& f : fields)
            if (f.quantity == q && f.role == role && f.component < 0) return f;
    return fields.empty() ? FieldRef{} : fields.front();
}

std::shared_ptr<SceneData> runDifference(const SceneData& a, const SceneData& b, QStringList* skipped)
{
    auto out = std::make_shared<SceneData>();
    out->label = a.label + " − " + b.label;
    out->provenance = a.provenance.isEmpty() ? b.provenance : a.provenance;
    out->synthetic = a.synthetic || b.synthetic;
    out->runDelta = true;
    out->planar = a.planar;
    out->bounds = a.bounds;
    for (const auto& ia : a.items) {
        const SceneItem* ib = b.item(ia.id);
        const bool lines = ia.mesh && ia.mesh->data().polylines;
        if (!ib || !ia.mesh || !ib->mesh || lines) {
            if (skipped && !lines) *skipped << ia.label;
            continue;
        }
        const MeshData& ma = ia.mesh->data();
        const MeshData& mb = ib->mesh->data();
        if (ma.positions != mb.positions || ma.indices != mb.indices) {
            if (skipped) *skipped << ia.label;
            continue;
        }
        MeshData d;
        d.positions = ma.positions;
        d.indices = ma.indices;
        for (const auto& fa : ma.fields) {
            const FieldRef r = FieldRef::parse(fa.desc.name);
            if (r.role != FieldRole::Kutta) continue;
            const auto* fb = mb.field(fa.desc.name);
            if (!fb || fb->desc.components != fa.desc.components) continue;
            MeshData::Field keep = fa;   // run A's own values stay available
            MeshData::Field diff{fa.desc, fa.values};
            diff.desc.name = r.quantity + "@run";
            for (std::size_t k = 0; k < diff.values.size(); ++k) diff.values[k] = fa.values[k] - fb->values[k];
            d.fields.push_back(std::move(keep));
            d.fields.push_back(std::move(diff));
        }
        SceneItem it = ia;
        it.mesh = std::make_shared<PreparedMesh>(std::move(d), ia.kind == ItemKind::Slice,
                                                 ia.kind == ItemKind::Slice
                                                     ? QVector3D(ia.mesh->positionsNormals()[3], ia.mesh->positionsNormals()[4],
                                                                 ia.mesh->positionsNormals()[5])
                                                     : QVector3D());
        out->items.push_back(std::move(it));
    }
    return out;
}

} // namespace kapp
