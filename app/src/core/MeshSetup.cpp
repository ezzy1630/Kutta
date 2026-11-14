#include "core/MeshSetup.h"
#include "core/MeshCli.h"

#include <QJsonArray>
#include <QObject>
#include <QSet>

#include <algorithm>
#include <cmath>

namespace kapp {

namespace {

QJsonArray vec(const QVector3D& v, double k) { return {double(v.x()) * k, double(v.y()) * k, double(v.z()) * k}; }
QVector3D vec(const QJsonValue& v, double k)
{
    const QJsonArray a = v.toArray();
    return a.size() == 3 ? QVector3D(float(a[0].toDouble() * k), float(a[1].toDouble() * k), float(a[2].toDouble() * k)) : QVector3D();
}

// Round to a tidy number for defaults the user will read: 1, 2, 2.5 or 5 × 10^n.
double tidy(double x)
{
    if (!(x > 0)) return x;
    const double p = std::pow(10.0, std::floor(std::log10(x)));
    const double m = x / p;
    return (m < 1.5 ? 1 : m < 2.25 ? 2 : m < 3.5 ? 2.5 : m < 7.5 ? 5 : 10) * p;
}

double extent(const Bounds& b, int axis) { return b.valid ? double(b.max[axis] - b.min[axis]) : 0.0; }

} // namespace

std::pair<int, int> parseAxis(const QString& axis)
{
    const int sign = axis.startsWith('-') ? -1 : 1;
    const QChar c = axis.isEmpty() ? QChar('x') : axis.back();
    return {c == 'y' ? 1 : c == 'z' ? 2 : 0, sign};
}

QJsonObject MeshSetup::toJson() const
{
    const double k = 1.0;  // mesh.json lengths are metres
    QJsonArray boxesJson;
    for (const auto& b : boxes)
        boxesJson.append(QJsonObject{{"name", b.name}, {"min", vec(b.min, k)}, {"max", vec(b.max, k)}, {"size", b.size * k}});
    return QJsonObject{
        {"v", 1},
        {"units", unitName(units)},
        {"flow", flow},
        {"up", up},
        {"tunnel", QJsonObject{{"upstream", tunnel.upstream}, {"downstream", tunnel.downstream}, {"side", tunnel.side}, {"top", tunnel.top},
                               {"ground", tunnel.ground}, {"symmetry", tunnel.symmetry.isEmpty() ? QJsonValue() : QJsonValue(tunnel.symmetry)}}},
        {"base_size", baseSize * k},
        {"surface", [&] {
             QJsonObject sf{{"min_size", minSize * k}, {"max_size", maxSize * k}, {"curvature_degrees", curvatureDegrees}, {"proximity_cells", proximityCells}};
             if (proximityMinSize > 0) sf["proximity_min_size"] = proximityMinSize * k;
             return sf;
         }()},
        {"refinement_boxes", boxesJson},
        {"layers", QJsonObject{{"count", layers}, {"first_height", firstHeight * k}, {"growth", growth},
                               {"zones", QJsonArray::fromStringList(layerZones)}}},
        {"max_cells", double(maxCells)},
    };
}

MeshSetup MeshSetup::fromJson(const QJsonObject& o, bool lengthsInUnits)
{
    MeshSetup s;
    s.units = unitFromName(o.value("units").toString("m"));
    const double k = lengthsInUnits ? unitScale(s.units) : 1.0;
    s.flow = o.value("flow").toString(s.flow);
    s.up = o.value("up").toString(s.up);
    const QJsonObject t = o.value("tunnel").toObject();
    s.tunnel.upstream = t.value("upstream").toDouble(s.tunnel.upstream);
    s.tunnel.downstream = t.value("downstream").toDouble(s.tunnel.downstream);
    s.tunnel.side = t.value("side").toDouble(s.tunnel.side);
    s.tunnel.top = t.value("top").toDouble(s.tunnel.top);
    s.tunnel.ground = t.value("ground").toBool(s.tunnel.ground);
    s.tunnel.symmetry = t.value("symmetry").toString();
    s.baseSize = o.value("base_size").toDouble(s.baseSize / k) * k;
    const QJsonObject sf = o.value("surface").toObject();
    s.minSize = sf.value("min_size").toDouble(s.minSize / k) * k;
    s.maxSize = sf.value("max_size").toDouble(s.maxSize / k) * k;
    s.curvatureDegrees = sf.value("curvature_degrees").toDouble(s.curvatureDegrees);
    s.proximityCells = sf.value("proximity_cells").toInt(s.proximityCells);
    s.proximityMinSize = sf.value("proximity_min_size").toDouble(0) * k;
    for (const QJsonValue& v : o.value("refinement_boxes").toArray()) {
        const QJsonObject b = v.toObject();
        s.boxes.push_back({b.value("name").toString(), vec(b.value("min"), k), vec(b.value("max"), k), b.value("size").toDouble() * k});
    }
    const QJsonObject l = o.value("layers").toObject();
    s.layers = l.value("count").toInt(s.layers);
    s.firstHeight = l.value("first_height").toDouble(s.firstHeight / k) * k;
    s.growth = l.value("growth").toDouble(s.growth);
    if (l.contains("zones")) {
        s.layerZones.clear();
        for (const QJsonValue& z : l.value("zones").toArray()) s.layerZones << z.toString();
    }
    s.maxCells = qint64(o.value("max_cells").toDouble(double(s.maxCells)));
    return s;
}

// The ranges kutta mesh v1 accepts (src/kutta/mesh/mesher/settings.cpp), so
// Generate is refused here, with the field's name, rather than by the mesher.
QStringList MeshSetup::problems() const
{
    QStringList p;
    auto positive = [&](double v, const QString& what) { if (!(v > 0) || !std::isfinite(v)) p << QObject::tr("%1 must be greater than zero.").arg(what); };
    positive(baseSize, QObject::tr("Base size"));
    positive(minSize, QObject::tr("Minimum surface size"));
    positive(maxSize, QObject::tr("Maximum surface size"));
    if (minSize > maxSize) p << QObject::tr("Minimum surface size is larger than the maximum.");
    if (maxSize > baseSize) p << QObject::tr("Maximum surface size is larger than the base size.");
    if (proximityMinSize < 0 || proximityMinSize > minSize) p << QObject::tr("Smallest gap cell must be between zero (automatic) and the minimum surface size.");
    if (!(curvatureDegrees >= 0.1 && curvatureDegrees <= 180)) p << QObject::tr("Curvature normal angle must be between 0.1 and 180°.");
    if (proximityCells < 1 || proximityCells > 20) p << QObject::tr("Cells per gap must be between 1 and 20.");
    if (layers < 0 || layers > 100) p << QObject::tr("Boundary layers must be between 0 and 100.");
    if (layers > 0) {
        positive(firstHeight, QObject::tr("First layer height"));
        if (!(growth > 1 && growth <= 3)) p << QObject::tr("Growth rate must be above 1 and at most 3.");
    }
    if (maxCells < 1000 || maxCells > 100'000'000) p << QObject::tr("Maximum cells must be between 1,000 and 100 million.");
    if (flow != "+x" || up != "+z") p << QObject::tr("kutta mesh meshes with the flow along +x and +z up; rotate the geometry in CAD to match.");
    for (double d : {tunnel.upstream, tunnel.downstream, tunnel.side, tunnel.top})
        if (!(d >= 0.01 && d <= 100)) {
            p << QObject::tr("Tunnel distances must be between 0.01 and 100 car lengths.");
            break;
        }
    QSet<QString> names;
    for (const auto& b : boxes) {
        if (b.name.trimmed().isEmpty()) p << QObject::tr("A refinement region has no name.");
        else if (names.contains(b.name)) p << QObject::tr("Two refinement regions are called %1; names must be unique.").arg(b.name);
        names.insert(b.name);
        if (!(b.min.x() < b.max.x() && b.min.y() < b.max.y() && b.min.z() < b.max.z()))
            p << QObject::tr("Refinement region %1: each minimum must be below its maximum.").arg(b.name);
        if (!(b.size > 0)) p << QObject::tr("Refinement region %1: size must be greater than zero.").arg(b.name);
    }
    return p;
}

MeshSetup defaultMeshSetup(const Bounds& g, LengthUnit units)
{
    MeshSetup s;
    s.units = units;
    if (!g.valid) return s;
    const double L = std::max({extent(g, 0), 1e-3});
    const double W = std::max(extent(g, 1), 1e-3), H = std::max(extent(g, 2), 1e-3);
    s.baseSize = tidy(L / 6);
    s.maxSize = tidy(L / 80);
    s.minSize = tidy(L / 1500);
    // y+ ~ 30 at the first cell's centre (half its height) at 20 m/s over the
    // body length: wall functions on the first cell.
    s.firstHeight = tidy(2 * firstHeightForYPlus(30, 20, L));
    const QVector3D lo = g.min, hi = g.max;
    // Around the body, then the near wake: the regions FSAE teams refine first.
    s.boxes.push_back({QObject::tr("body"), QVector3D(float(lo.x() - 0.1 * L), float(lo.y() - 0.1 * W), lo.z()),
                       QVector3D(float(hi.x() + 0.15 * L), float(hi.y() + 0.1 * W), float(hi.z() + 0.25 * H)), tidy(L / 120)});
    s.boxes.push_back({QObject::tr("wake"), QVector3D(hi.x(), float(lo.y() - 0.05 * W), lo.z()),
                       QVector3D(float(hi.x() + 1.0 * L), float(hi.y() + 0.05 * W), float(hi.z() + 0.2 * H)), tidy(L / 60)});
    return s;
}

Bounds tunnelBounds(const MeshSetup& s, const Bounds& g)
{
    Bounds t;
    if (!g.valid) return t;
    const float L = std::max(g.max.x() - g.min.x(), 1e-6f);
    t.min = QVector3D(g.min.x() - float(s.tunnel.upstream) * L, s.tunnel.symmetry == "y" ? 0.0f : g.min.y() - float(s.tunnel.side) * L,
                      s.tunnel.ground ? 0.0f : g.min.z() - float(s.tunnel.side) * L);
    t.max = QVector3D(g.max.x() + float(s.tunnel.downstream) * L, g.max.y() + float(s.tunnel.side) * L, g.max.z() + float(s.tunnel.top) * L);
    if (s.baseSize > 0)
        for (int a = 0; a < 3; ++a) {
            const double cells = std::ceil(double(t.max[a] - t.min[a]) / s.baseSize - 1e-4);  // float bounds
            if (cells > 0 && cells < 1e5) t.max[a] = float(t.min[a] + cells * s.baseSize);
        }
    t.valid = true;
    return t;
}

std::vector<TunnelFace> tunnelFaces(const MeshSetup& s)
{
    const bool half = s.tunnel.symmetry == "y";
    return {{"inlet", "velocity-inlet", 0, false},
            {"outlet", "pressure-outlet", 0, true},
            {half ? "symmetry" : "side-low", "symmetry", 1, false},
            {"side-high", "symmetry", 1, true},
            {s.tunnel.ground ? "ground" : "bottom", s.tunnel.ground ? "moving-ground" : "symmetry", 2, false},
            {"top", "symmetry", 2, true}};
}

const QStringList& reservedZoneNames()
{
    static const QStringList names{"inlet", "outlet", "ground", "top", "bottom", "symmetry", "side-low", "side-high"};
    return names;
}

double CellEstimate::solverBytes() const { return total() * kSolverBytesPerCell; }

// Calibrated on kutta mesh's coarse Ahmed (144,877 cells) and FSAE
// (1,312,313) meshes, each within 5 % given the layered wall area; the app
// passes the whole STL area, which reads 141k and 1.41M (+8 %). On the
// default FSAE setup (4.42M) it gives 3.5M. Other setups are untested: a
// planning figure.
CellEstimate estimateCells(const MeshSetup& s, const Bounds& g, double area)
{
    CellEstimate e;
    if (!g.valid || !(s.minSize > 0) || !(s.maxSize > 0) || !(s.baseSize > 0)) return e;
    // Octree sizes are the base size halved: a requested size rounds down to one.
    auto level = [&](double h) { return s.baseSize / std::pow(2.0, std::ceil(std::log2(s.baseSize / std::max(h, 1e-6)) - 1e-9)); };
    // Next to layered walls, cells are at most four last-layer thicknesses.
    double wall = s.maxSize;
    if (s.layers > 0) wall = std::min(wall, 4 * s.firstHeight * std::pow(s.growth, s.layers - 1));
    const double hs = level(wall);
    // Curvature refines a share of the surface part-way to the minimum size;
    // the sharpest edges and narrowest gaps (a few percent) reach the minimum.
    const double curved = std::clamp(18.0 / std::max(1.0, s.curvatureDegrees) * 0.25, 0.05, 0.8);
    const double sharp = s.proximityCells > 0 ? 0.03 : 0.015;
    const double hc = std::min(hs, std::max(s.minSize, s.maxSize / 8)), hm = std::min(hs, s.minSize);
    // The castellated wall has about twice the faces of a smooth one (steps, hanging faces).
    const double perArea = 2.0 * (sharp / (hm * hm) + curved / (hc * hc) + (1 - curved - sharp) / (hs * hs));
    e.surfaceFaces = area * perArea;
    // A band about 2.6 cells thick at each level from the wall out, sizes
    // doubling: 2.6 · 4/3 cells per wall face.
    e.octreeCells = 2.6 * 4.0 / 3.0 * e.surfaceFaces;
    e.layerCells = s.layers > 0 ? e.surfaceFaces * s.layers : 0;
    const Bounds t = tunnelBounds(s, g);
    double boxVolume = 0;
    for (const auto& b : s.boxes) {
        const QVector3D d = b.max - b.min;
        const double v = std::max(0.0f, d.x()) * double(std::max(0.0f, d.y())) * std::max(0.0f, d.z());
        boxVolume += v;
        if (b.size > 0) e.boxCells += v / std::pow(level(b.size), 3);
    }
    const QVector3D d = t.max - t.min;
    const double domain = double(d.x()) * d.y() * d.z();
    e.farFieldCells = std::max(0.0, domain - boxVolume) / (s.baseSize * s.baseSize * s.baseSize);
    e.overLimit = e.total() > double(s.maxCells);
    return e;
}

double firstHeightForYPlus(double yPlus, double speed, double length, double density, double viscosity)
{
    const double re = density * speed * length / viscosity;
    const double cf = std::pow(2 * std::log10(re) - 0.65, -2.3);
    const double uTau = std::sqrt(0.5 * cf) * speed;
    return yPlus * viscosity / (density * uTau);
}

double yPlusForFirstHeight(double height, double speed, double length, double density, double viscosity)
{
    return height / firstHeightForYPlus(1.0, speed, length, density, viscosity);
}

} // namespace kapp
