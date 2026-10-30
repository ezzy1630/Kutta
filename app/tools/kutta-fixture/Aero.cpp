#include "Aero.h"

#include <QSet>

#include <cmath>

namespace fixture {

namespace {

struct V { double x = 0, y = 0, z = 0; };
V operator+(V a, V b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
V operator-(V a, V b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
V operator*(V a, double s) { return {a.x * s, a.y * s, a.z * s}; }
double dot(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
V cross(V a, V b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
double norm(V a) { return std::sqrt(dot(a, a)); }

std::optional<V> vec(const QJsonValue& v)
{
    const QJsonArray a = v.toArray();
    if (a.size() != 3 || !a[0].isDouble() || !a[1].isDouble() || !a[2].isDouble()) return std::nullopt;
    return V{a[0].toDouble(), a[1].toDouble(), a[2].toDouble()};
}
QJsonArray json(V v) { return {v.x, v.y, v.z}; }
QJsonValue opt(std::optional<double> v) { return v ? QJsonValue(*v) : QJsonValue(); }

struct Loads { V pressure, viscous, pressureMoment, viscousMoment; };

} // namespace

QJsonObject aeroFromLoads(const QJsonObject& loads, const QJsonObject& spec, QString* error)
{
    auto fail = [&](const QString& why) {
        if (error) *error = why;
        return QJsonObject{};
    };
    static const QSet<QString> keys{"v", "axes", "moment_centre", "car_zones", "groups", "vehicle"};
    for (const QString& k : spec.keys())
        if (!keys.contains(k)) return fail("aero spec: unknown key " + k);
    if (spec.contains("v") && spec.value("v") != QJsonValue(1)) return fail("aero spec: v must be 1");
    const QJsonObject axes = spec.value("axes").toObject();
    V drag = vec(axes.value("drag")).value_or(V{1, 0, 0}), lift = vec(axes.value("lift")).value_or(V{0, 0, 1});
    if (norm(drag) == 0 || norm(lift) == 0) return fail("aero spec: axes must be nonzero");
    drag = drag * (1 / norm(drag));
    lift = lift * (1 / norm(lift));
    if (std::abs(dot(drag, lift)) > 1e-6) return fail("aero spec: drag and lift axes must be orthogonal");
    const V side = cross(lift, drag);
    const V centre = vec(spec.value("moment_centre")).value_or(V{});

    std::map<QString, Loads> zones;
    QStringList order;
    for (const QJsonValue& zv : loads.value("zones").toArray()) {
        const QJsonObject z = zv.toObject();
        const QString name = z.value("name").toString();
        Loads l{*vec(z.value("pressure")), *vec(z.value("viscous")), *vec(z.value("pressure_moment")), *vec(z.value("viscous_moment"))};
        zones[name] = l;
        order << name;
    }
    QStringList car;
    for (const QJsonValue& z : spec.value("car_zones").toArray()) car << z.toString();
    if (car.isEmpty())
        for (const QJsonValue& z : loads.value("report_zones").toArray()) car << z.toString();
    if (car.isEmpty()) car = order;
    for (const QString& z : car)
        if (!zones.count(z)) return fail("aero spec: unknown wall zone " + z);

    struct Group { QString name; QStringList zones; };
    std::vector<Group> groups;
    QSet<QString> grouped, names;
    for (const QJsonValue& gv : spec.value("groups").toArray()) {
        const QJsonObject g = gv.toObject();
        Group x{g.value("name").toString(), {}};
        if (x.name.isEmpty() || names.contains(x.name)) return fail("aero spec: group names must be unique and nonempty");
        names.insert(x.name);
        for (const QJsonValue& z : g.value("zones").toArray()) {
            const QString n = z.toString();
            if (!car.contains(n)) return fail("aero spec: group zone is not a car zone: " + n);
            if (grouped.contains(n)) return fail("aero spec: zone in two groups: " + n);
            grouped.insert(n);
            x.zones << n;
        }
        groups.push_back(x);
    }

    const QJsonObject ref = loads.value("reference").toObject();
    auto positive = [&](const char* k) { const double v = ref.value(k).toDouble(); return ref.value(k).isDouble() && v > 0 ? std::optional<double>(v) : std::nullopt; };
    const auto area = positive("area"), density = positive("density"), velocity = positive("velocity"), length = positive("length");
    const std::optional<double> q = density && velocity ? std::optional<double>(0.5 * *density * *velocity * *velocity) : std::nullopt;

    auto sum = [&](const QStringList& names) {
        Loads t;
        for (const QString& n : names) {
            const Loads& l = zones.at(n);
            t.pressure = t.pressure + l.pressure;
            t.viscous = t.viscous + l.viscous;
            t.pressureMoment = t.pressureMoment + l.pressureMoment;
            t.viscousMoment = t.viscousMoment + l.viscousMoment;
        }
        return t;
    };
    const Loads total = sum(car);
    const double totalDownforce = -dot(total.pressure + total.viscous, lift);
    auto row = [&](const QString& name, const QString& kind, const QStringList& zs) {
        const Loads l = sum(zs);
        const V f = l.pressure + l.viscous;
        // Moments about the moment centre: M_C = M_O + (O - C) x F, O the origin.
        const V m = l.pressureMoment + l.viscousMoment + cross(V{} - centre, f);
        auto comp = [&](V x) { return QJsonObject{{"drag", dot(x, drag)}, {"lift", dot(x, lift)}, {"side", dot(x, side)}}; };
        // Every key is present; null when the references needed are missing.
        const std::optional<double> qa = q && area ? std::optional<double>(*q * *area) : std::nullopt;
        const std::optional<double> qal = qa && length ? std::optional<double>(*qa * *length) : std::nullopt;
        auto ratio = [](double v, std::optional<double> d) { return d ? QJsonValue(v / *d) : QJsonValue(); };
        const QJsonObject coeffs{{"cd", ratio(dot(f, drag), qa)}, {"cl", ratio(dot(f, lift), qa)}, {"cs", ratio(dot(f, side), qa)},
                                 {"cd_pressure", ratio(dot(l.pressure, drag), qa)}, {"cd_viscous", ratio(dot(l.viscous, drag), qa)},
                                 {"cl_pressure", ratio(dot(l.pressure, lift), qa)}, {"cl_viscous", ratio(dot(l.viscous, lift), qa)},
                                 {"cla", ratio(dot(f, lift), q)}, {"cda", ratio(dot(f, drag), q)},
                                 {"cm_pitch", ratio(dot(m, side), qal)}, {"cm_roll", ratio(dot(m, drag), qal)}, {"cm_yaw", ratio(dot(m, lift), qal)}};
        const double down = -dot(f, lift);
        return QJsonObject{{"name", name}, {"kind", kind}, {"zones", QJsonArray::fromStringList(zs)},
                           {"force", QJsonObject{{"pressure", comp(l.pressure)}, {"viscous", comp(l.viscous)}, {"total", comp(f)}}},
                           {"downforce", down},
                           {"moment", QJsonObject{{"pitch", dot(m, side)}, {"roll", dot(m, drag)}, {"yaw", dot(m, lift)}}},
                           {"coefficients", coeffs},
                           {"downforce_share", totalDownforce != 0 ? QJsonValue(down / totalDownforce) : QJsonValue()}};
    };
    QJsonArray rows;
    for (const QString& z : order)
        if (car.contains(z)) rows.append(row(z, "zone", {z}));
    for (const Group& g : groups) rows.append(row(g.name, "group", g.zones));
    rows.append(row("Total", "total", car));
    QStringList ungrouped;
    for (const QString& z : car)
        if (!grouped.contains(z)) ungrouped << z;

    const V f = total.pressure + total.viscous;
    const double totalDrag = dot(f, drag);
    QJsonObject balance{{"status", "needs_vehicle"}, {"front_axle", QJsonValue()}, {"wheelbase", QJsonValue()}, {"front_load", QJsonValue()},
                        {"rear_load", QJsonValue()}, {"front_percent", QJsonValue()}, {"cop_from_front_axle", QJsonValue()}};
    const QJsonObject vehicle = spec.value("vehicle").toObject();
    if (!vehicle.isEmpty()) {
        const auto axle = vec(vehicle.value("front_axle"));
        const double wb = vehicle.value("wheelbase").toDouble();
        if (!axle || !(wb > 0)) return fail("aero spec: vehicle needs front_axle and a positive wheelbase");
        const V mFront = total.pressureMoment + total.viscousMoment + cross(V{} - *axle, f);
        const double rear = dot(mFront, side) / wb, front = totalDownforce - rear;
        const bool down = totalDownforce > 0;
        balance = {{"status", down ? "defined" : "no_downforce"}, {"front_axle", json(*axle)}, {"wheelbase", wb},
                   {"front_load", front}, {"rear_load", rear},
                   {"front_percent", down ? QJsonValue(100 * front / totalDownforce) : QJsonValue()},
                   {"cop_from_front_axle", down ? QJsonValue(wb * rear / totalDownforce) : QJsonValue()}};
    }
    if (error) error->clear();
    return {{"v", 1}, {"provenance", loads.value("provenance")}, {"approximate", loads.value("approximate")},
            {"modified", loads.value("modified")}, {"iterations", loads.value("iterations")},
            {"reference", QJsonObject{{"area", opt(area)}, {"density", opt(density)}, {"velocity", opt(velocity)}, {"length", opt(length)},
                                      {"dynamic_pressure", opt(q)}}},
            {"axes", QJsonObject{{"drag", json(drag)}, {"lift", json(lift)}, {"side", json(side)}}},
            {"moment_centre", json(centre)},
            {"rows", rows},
            {"ungrouped", QJsonArray::fromStringList(ungrouped)},
            {"lift_to_drag", totalDrag != 0 ? QJsonValue(totalDownforce / totalDrag) : QJsonValue()},
            {"balance", balance}};
}

QJsonObject ConvergenceConfig::toJson() const
{
    QJsonObject res;
    for (auto it = residuals.begin(); it != residuals.end(); ++it) res[it.key()] = it.value();
    QJsonArray reps;
    for (const ReportCondition& r : reports)
        reps.append(QJsonObject{{"report", r.report}, {"stop_criterion", r.criterion}, {"previous_values", r.previous}, {"ignore_before", r.ignoreBefore}});
    return {{"enabled", enabled}, {"condition", any ? "any" : "all"}, {"residuals", res}, {"reports", reps}};
}

ConvergenceConfig convergenceFrom(const QJsonArray& settings, const QJsonObject& overrides)
{
    QMap<QString, QJsonValue> v;
    for (const QJsonValue& s : settings) {
        const QString p = s.toObject().value("path").toString();
        if (p.startsWith("monitors/")) v[p] = s.toObject().value("value");
    }
    for (auto it = overrides.begin(); it != overrides.end(); ++it)
        if (it.key().startsWith("monitors/")) v[it.key()] = it.value();
    ConvergenceConfig c;
    c.enabled = v.value("monitors/check-convergence").toBool();
    c.any = v.value("monitors/condition").toString() == "any";
    for (auto it = v.begin(); it != v.end(); ++it) {
        const QStringList parts = it.key().split('/');
        if (parts.size() == 4 && parts[1] == "residual" && parts[3] == "check" && it.value().toBool()) {
            const QJsonValue crit = v.value("monitors/residual/" + parts[2] + "/absolute-criterion");
            if (crit.isDouble()) c.residuals[parts[2]] = crit.toDouble();
        }
        if (parts.size() == 4 && parts[1] == "report" && parts[3] == "active" && it.value().toBool()) {
            const QString base = "monitors/report/" + parts[2] + "/";
            c.reports.push_back({parts[2], v.value(base + "stop-criterion").toDouble(1e-3), v.value(base + "previous-values").toInt(5),
                                 v.value(base + "ignore-before").toInt(100)});
        }
    }
    return c;
}

QJsonObject ConvergenceCheck::sample(int iter, const QMap<QString, double>& residuals, const QMap<QString, double>& reportForces)
{
    QJsonArray conditions;
    QStringList met;
    int count = 0;
    if (!m_config.residuals.isEmpty()) {
        ++count;
        bool all = true;
        QJsonObject eqs;
        for (auto it = m_config.residuals.begin(); it != m_config.residuals.end(); ++it) {
            const double r = residuals.value(it.key(), INFINITY);
            const bool ok = r <= it.value();
            all &= ok;
            eqs[it.key()] = QJsonObject{{"value", r}, {"criterion", it.value()}, {"met", ok}};
        }
        conditions.append(QJsonObject{{"kind", "residuals"}, {"met", all}, {"equations", eqs}});
        if (all) met << "residuals";
    }
    for (const ReportCondition& r : m_config.reports) {
        ++count;
        auto& h = m_history[r.report];
        if (iter > r.ignoreBefore && reportForces.contains(r.report)) h.push_back(reportForces.value(r.report));
        std::optional<double> value;
        if (int(h.size()) >= r.previous + 1 && h.back() != 0) {
            double worst = 0;
            const std::size_t n = h.size() - 1;
            for (int k = 1; k <= r.previous; ++k)
                worst = std::max(worst, std::abs(h[n - std::size_t(k) + 1] - h[n - std::size_t(k)]) / std::abs(h.back()));
            value = worst;
        }
        const bool ok = value && *value < r.criterion;
        conditions.append(QJsonObject{{"kind", "report"}, {"report", r.report}, {"value", opt(value)}, {"criterion", r.criterion}, {"met", ok}});
        if (ok) met << "report:" + r.report;
    }
    m_met = count > 0 && (m_config.any ? !met.isEmpty() : met.size() == count);
    m_metBy = m_met ? met : QStringList{};
    return {{"met", m_met}, {"conditions", conditions}};
}

} // namespace fixture
