#include "core/MeshSetup.h"

#include <QJsonArray>
#include <QObject>

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
    const double k = 1.0 / unitScale(units);
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
        {"surface", QJsonObject{{"min_size", minSize * k}, {"max_size", maxSize * k}, {"curvature_degrees", curvatureDegrees},
                                {"proximity_cells", proximityCells}}},
        {"refinement_boxes", boxesJson},
        {"layers", QJsonObject{{"count", layers}, {"first_height", firstHeight * k}, {"growth", growth},
                               {"zones", QJsonArray::fromStringList(layerZones)}}},
        {"max_cells", double(maxCells)},
    };
}

MeshSetup MeshSetup::fromJson(const QJsonObject& o)
{
    MeshSetup s;
    s.units = unitFromName(o.value("units").toString("m"));
    const double k = unitScale(s.units);
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

QStringList MeshSetup::problems() const
{
    QStringList p;
    auto positive = [&](double v, const QString& what) { if (!(v > 0) || !std::isfinite(v)) p << QObject::tr("%1 must be greater than zero.").arg(what); };
    positive(baseSize, QObject::tr("Base size"));
    positive(minSize, QObject::tr("Minimum surface size"));
    positive(maxSize, QObject::tr("Maximum surface size"));
    if (minSize > maxSize) p << QObject::tr("Minimum surface size is larger than the maximum.");
    if (maxSize > baseSize) p << QObject::tr("Maximum surface size is larger than the base size.");
    if (!(curvatureDegrees > 0 && curvatureDegrees < 90)) p << QObject::tr("Curvature normal angle must be between 0 and 90°.");
    if (proximityCells < 0) p << QObject::tr("Cells per gap cannot be negative.");
    if (layers < 0 || layers > 60) p << QObject::tr("Boundary layers must be between 0 and 60.");
    if (layers > 0) {
        positive(firstHeight, QObject::tr("First layer height"));
        if (!(growth >= 1 && growth <= 3)) p << QObject::tr("Growth rate must be between 1 and 3.");
    }
    if (maxCells < 1000) p << QObject::tr("Maximum cells must be at least 1,000.");
    const auto [fa, fs] = parseAxis(flow);
    const auto [ua, us] = parseAxis(up);
    if (fa == ua) p << QObject::tr("Flow and up directions must be different axes.");
    for (double d : {tunnel.upstream, tunnel.downstream, tunnel.side, tunnel.top})
        if (!(d > 0)) {
            p << QObject::tr("Tunnel distances must be greater than zero.");
            break;
        }
    for (const auto& b : boxes) {
        if (b.name.trimmed().isEmpty()) p << QObject::tr("A refinement region has no name.");
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
    // y+ ~ 30 at 20 m/s over the body length: wall functions on the first cell.
    s.firstHeight = tidy(firstHeightForYPlus(30, 20, L));
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
    const auto [fa, fs] = parseAxis(s.flow);
    const auto [ua, us] = parseAxis(s.up);
    const int sa = 3 - fa - ua;
    const float L = std::max(g.max[fa] - g.min[fa], 1e-3f);
    t.min = g.min;
    t.max = g.max;
    // Upstream is against the flow.
    const float up = float(s.tunnel.upstream) * L, down = float(s.tunnel.downstream) * L;
    if (fs > 0) { t.min[fa] -= up; t.max[fa] += down; }
    else { t.max[fa] += up; t.min[fa] -= down; }
    const float top = float(s.tunnel.top) * L, side = float(s.tunnel.side) * L;
    if (us > 0) {
        t.max[ua] += top;
        if (!s.tunnel.ground) t.min[ua] -= top;
    } else {
        t.min[ua] -= top;
        if (!s.tunnel.ground) t.max[ua] += top;
    }
    t.min[sa] -= side;
    t.max[sa] += side;
    if (s.tunnel.symmetry == "y" && sa == 1) t.min[1] = std::max(t.min[1], 0.0f);
    t.valid = true;
    return t;
}

std::vector<TunnelFace> tunnelFaces(const MeshSetup& s)
{
    const auto [fa, fs] = parseAxis(s.flow);
    const auto [ua, us] = parseAxis(s.up);
    const int sa = 3 - fa - ua;
    std::vector<TunnelFace> f;
    f.push_back({"inlet", "velocity-inlet", fa, fs < 0});
    f.push_back({"outlet", "pressure-outlet", fa, fs > 0});
    f.push_back({s.tunnel.ground ? "ground" : "bottom", s.tunnel.ground ? "moving-ground" : "symmetry", ua, us < 0});
    f.push_back({"top", "symmetry", ua, us > 0});
    const bool half = s.tunnel.symmetry == "y" && sa == 1;
    f.push_back({half ? "symmetry" : "sides", "symmetry", sa, false});
    f.push_back({"sides", "symmetry", sa, true});
    return f;
}

double CellEstimate::solverBytes() const { return total() * kPlanningBytesPerCell; }

CellEstimate estimateCells(const MeshSetup& s, const Bounds& g, double area)
{
    CellEstimate e;
    if (!g.valid || !(s.minSize > 0) || !(s.maxSize > 0) || !(s.baseSize > 0)) return e;
    // Curvature refines a share of the surface part-way to the minimum size;
    // the sharpest edges and narrowest gaps (a few percent) reach the minimum;
    // the rest sits near the maximum. Faces per area add up over the bands.
    const double curved = std::clamp(18.0 / std::max(1.0, s.curvatureDegrees) * 0.25, 0.05, 0.8);
    const double sharp = s.proximityCells > 0 ? 0.03 : 0.015;
    const double hc = std::max(s.minSize, s.maxSize / 8);
    const double perArea = sharp / (s.minSize * s.minSize) + curved / (hc * hc) + (1 - curved - sharp) / (s.maxSize * s.maxSize);
    e.surfaceFaces = area * perArea * 2.0;  // trimmed cells cut into about two faces each
    // A band about four cells thick at each octree level, sizes doubling to
    // the base size: 4A/h² summed over levels ≈ 5.3 A/h².
    e.octreeCells = 5.3 * area * perArea;
    e.layerCells = s.layers > 0 ? e.surfaceFaces / 2.0 * s.layers : 0;
    const Bounds t = tunnelBounds(s, g);
    double boxVolume = 0;
    for (const auto& b : s.boxes) {
        const QVector3D d = b.max - b.min;
        const double v = std::max(0.0f, d.x()) * double(std::max(0.0f, d.y())) * std::max(0.0f, d.z());
        boxVolume += v;
        if (b.size > 0) e.boxCells += v / (b.size * b.size * b.size);
    }
    const QVector3D d = t.max - t.min;
    const double domain = double(d.x()) * d.y() * d.z();
    // The far field grades from the finest region to the base size: about
    // twice the cells of a uniform base-size grid over what is left.
    e.farFieldCells = 2.0 * std::max(0.0, domain - boxVolume) / (s.baseSize * s.baseSize * s.baseSize);
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
