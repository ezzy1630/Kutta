// kutta-fixture: a stand-in for the `kutta` CLI that speaks docs/app-protocol.md
// (plus the app's proposed `post` and `--live-view` additions) from recorded and
// synthetic data. The app talks to it through the same process boundary as the
// real CLI; switching is a path change.
//
//   mixing_elbow.cas.h5   replays a REAL recorded Kutta run (Ansys example, MIT)
//   fsae_car*.cas.h5      a SYNTHETIC car; every output says "synthetic": true
//
// Environment:
//   KUTTA_FIXTURE_DIR         fixture data directory
//   KUTTA_FIXTURE_SPEED       replay speed multiplier (default 6)
//   KUTTA_FIXTURE_HOLD_AT     pause after this iteration until stdin closes
//   KUTTA_FIXTURE_DEVICES     "macos" | "windows" device inventory
//   KUTTA_FIXTURE_CAR_DETAIL  car tessellation multiplier (1 ~ 0.4M triangles)

#include "SyntheticCar.h"

#include "core/CliClient.h"
#include "core/Geometry.h"
#include "core/Protocol.h"
#include "core/SettingsModel.h"
#include "core/ViewFormat.h"

#include <QCoreApplication>
#include <QDateTime>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QThread>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <thread>

namespace {

using namespace fixture;

QTextStream& out()
{
    static QTextStream s(stdout);
    return s;
}

void emitJson(const QJsonObject& o)
{
    out() << QJsonDocument(o).toJson(QJsonDocument::Compact) << '\n';
    out().flush();
}

void eprint(const QString& s)
{
    std::fputs((s + "\n").toUtf8().constData(), stderr);
    std::fflush(stderr);
}

[[noreturn]] void usage(bool jsonl)
{
    const QString msg = "Usage: kutta-fixture devices --json | inspect CASE [--data DAT] --json | solve CASE ... --out DIR "
                        "--progress jsonl | export-view DIR [--slice ox,oy,oz,nx,ny,nz]... | post DIR --request FILE";
    if (jsonl) emitJson({{"v", 1}, {"event", "error"}, {"message", msg}, {"issues", QJsonArray{}}});
    eprint(msg);
    std::exit(2);
}

[[noreturn]] void fail(const QString& message, const QJsonArray& issues = {}, bool jsonl = true)
{
    if (jsonl) emitJson({{"v", 1}, {"event", "error"}, {"message", message}, {"issues", issues}});
    eprint("kutta-fixture: " + message);
    std::exit(1);
}

QString fixtureDir()
{
    const QString env = qEnvironmentVariable("KUTTA_FIXTURE_DIR");
    if (!env.isEmpty()) return env;
    const QDir exe(QCoreApplication::applicationDirPath());
    for (const QString& rel : {QStringLiteral("fixtures"), QStringLiteral("../Resources/fixtures")})
        if (QFileInfo(exe.filePath(rel)).isDir()) return QDir(exe.filePath(rel)).absolutePath();
#ifdef KUTTA_FIXTURE_DIR
    return QStringLiteral(KUTTA_FIXTURE_DIR);
#else
    return exe.filePath("fixtures");
#endif
}

// Plain file, or the qCompress'd ".qz" variant.
QByteArray readFixture(const QString& rel)
{
    const QString base = QDir(fixtureDir()).filePath(rel);
    QFile f(base);
    if (f.open(QIODevice::ReadOnly)) return f.readAll();
    QFile z(base + ".qz");
    if (z.open(QIODevice::ReadOnly)) return qUncompress(z.readAll());
    return {};
}

QJsonObject readJson(const QString& path)
{
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) return {};
    return QJsonDocument::fromJson(f.readAll()).object();
}

void writeFile(const QString& path, const QByteArray& data)
{
    QDir().mkpath(QFileInfo(path).absolutePath());
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly)) fail("cannot write " + path, {}, false);
    f.write(data);
}

void writeJson(const QString& path, const QJsonObject& o) { writeFile(path, QJsonDocument(o).toJson(QJsonDocument::Indented)); }

struct Args {
    QStringList positional;
    QMap<QString, QStringList> opts;
    bool has(const QString& k) const { return opts.contains(k); }
    QString value(const QString& k, const QString& d = {}) const { return opts.value(k).value(0, d); }
};

Args parseArgs(const QStringList& args)
{
    Args a;
    for (int i = 0; i < args.size(); ++i) {
        const QString& s = args[i];
        if (s == "--json") a.opts["json"] << "1";
        else if (s.startsWith("--")) {
            if (i + 1 >= args.size()) usage(false);
            a.opts[s.mid(2)] << args[++i];
        } else a.positional << s;
    }
    return a;
}

enum class CaseKind { Elbow, Car };

CaseKind caseKind(const QString& casePath, bool jsonl)
{
    const QString name = QFileInfo(casePath).fileName();
    if (name.startsWith("mixing_elbow")) return CaseKind::Elbow;
    if (name.startsWith("fsae_car")) return CaseKind::Car;
    fail(QStringLiteral("The fixture CLI has no data for %1. It serves mixing_elbow.cas.h5 (recorded) and "
                        "fsae_car.cas.h5 (synthetic); point the app at the real kutta CLI for other cases.").arg(name),
         {}, jsonl);
}

QJsonArray rawElbowSettings()
{
    const QJsonObject elbow = QJsonDocument::fromJson(readFixture("cases/mixing_elbow/inspect.json")).object();
    QJsonArray raw;
    for (const QJsonValue& v : elbow.value("settings").toArray())
        if (v.toObject().value("path").toString().startsWith("case/")) raw.append(v);
    return raw;
}

QJsonObject inspectJson(const QString& casePath, const QString& dataPath, const CarParams& car = {})
{
    const bool data = !dataPath.isEmpty();
    QJsonObject o;
    if (caseKind(casePath, false) == CaseKind::Elbow) {
        o = QJsonDocument::fromJson(readFixture(data ? "cases/mixing_elbow/inspect-data.json" : "cases/mixing_elbow/inspect.json")).object();
        if (o.isEmpty()) fail("fixture data missing: cases/mixing_elbow/inspect.json", {}, false);
    } else {
        o = carInspect(car, data, rawElbowSettings());
    }
    QJsonObject file = o.value("file").toObject();
    file["path"] = QFileInfo(casePath).absoluteFilePath();
    o["file"] = file;
    if (data) {
        QJsonObject d = o.value("data").toObject();
        d["file"] = QFileInfo(dataPath).absoluteFilePath();
        o["data"] = d;
    }
    return o;
}

QJsonObject devicesJson()
{
    QString which = qEnvironmentVariable("KUTTA_FIXTURE_DEVICES");
    if (which.isEmpty()) {
#ifdef Q_OS_WIN
        which = "windows";
#else
        which = "macos";
#endif
    }
    return QJsonDocument::fromJson(readFixture(QStringLiteral("devices-%1.json").arg(which))).object();
}

// ---- stop handling: stdin EOF ------------------------------------------------------------

std::atomic<bool> g_stop{false};

void watchStdin()
{
    std::thread([] {
        char buf[256];
        while (std::fread(buf, 1, sizeof buf, stdin) > 0) {}
        g_stop = true;
    }).detach();
}

void sleepMs(double ms)
{
    if (ms > 0) QThread::usleep(static_cast<unsigned long>(ms * 1000));
}

// ---- solve -----------------------------------------------------------------------------------

struct RunSpec {
    CaseKind kind;
    QString casePath, dataPath, outDir;
    QStringList approximations;
    int iterations = 1000, monitor = 10, liveView = 0;
    QJsonObject overrides;
    CarParams car;
    QJsonObject device;
};

QString approximationNotice(const RunSpec& r)
{
    if (r.approximations.isEmpty()) return {};
    QString n = "APPROXIMATE RESULT: not a faithful Fluent reproduction; wfgc: recon/cell-lsf? Fluent #t -> Kutta run without "
                "warped-face gradient correction (WFGC off, #f); V6 CPU hybrid diagnostic versus Fluent: velocity L2 3.601459%, "
                "pressure L2 2.961922%, k L2 15.478919%, omega L2 8.474918%, pressure-drop error 1.185077%; combined solver "
                "discrepancy, not isolated WFGC error or a general bound; evidence: src/kutta/fv/VERIFICATION.md#v6-initialization-comparison";
    return n;
}

QString provenance(const RunSpec& r)
{
    QStringList parts;
    if (QString a = approximationNotice(r); !a.isEmpty()) parts << a;
    if (!r.overrides.isEmpty()) {
        QStringList keys = r.overrides.keys();
        parts << "MODIFIED RESULT: case settings overridden (" + keys.join(", ") + ")";
    }
    if (r.kind == CaseKind::Car) parts << "SYNTHETIC FIXTURE: invented geometry and fields for app development; not a CFD result";
    else if (!r.overrides.isEmpty()) parts << "FIXTURE REPLAY: recorded mixing-elbow fields; overrides did not change these numbers";
    return parts.join(" | ");
}

void applyCarOverrides(RunSpec& r)
{
    const QJsonObject& o = r.overrides;
    auto num = [&](const char* k, double& dst) { if (o.contains(k)) dst = o.value(k).toDouble(); };
    num("boundary/inlet/velocity-magnitude", r.car.speed);
    num("materials/air/density", r.car.density);
    num("materials/air/viscosity", r.car.viscosity);
    num("reference/area", r.car.refArea);
    num("reference/density", r.car.refDensity);
    num("reference/velocity", r.car.refVelocity);
    num("reference/pressure", r.car.refPressure);
}

QJsonObject runDescriptor(const RunSpec& r)
{
    return {{"case", r.kind == CaseKind::Elbow ? "mixing_elbow" : "fsae_car"},
            {"variant", r.dataPath.isEmpty() ? "cold" : "seeded"},
            {"reference", !r.dataPath.isEmpty()},
            {"provenance", provenance(r)},
            {"synthetic", r.kind == CaseKind::Car},
            {"params", QJsonObject{{"speed", r.car.speed}, {"density", r.car.density}, {"viscosity", r.car.viscosity},
                                   {"ref_area", r.car.refArea}, {"ref_density", r.car.refDensity},
                                   {"ref_velocity", r.car.refVelocity}, {"ref_pressure", r.car.refPressure}}}};
}

CarParams paramsFrom(const QJsonObject& d)
{
    const QJsonObject p = d.value("params").toObject();
    CarParams c;
    c.speed = p.value("speed").toDouble(c.speed);
    c.density = p.value("density").toDouble(c.density);
    c.viscosity = p.value("viscosity").toDouble(c.viscosity);
    c.refArea = p.value("ref_area").toDouble(c.refArea);
    c.refDensity = p.value("ref_density").toDouble(c.refDensity);
    c.refVelocity = p.value("ref_velocity").toDouble(c.refVelocity);
    c.refPressure = p.value("ref_pressure").toDouble(c.refPressure);
    c.detail = std::max(1, qEnvironmentVariableIntValue("KUTTA_FIXTURE_CAR_DETAIL"));
    return c;
}

// Write a car view (walls + flow boundaries + slices) into dir.
void writeCarView(const QString& dir, const CarParams& prm, bool reference, const QString& prov,
                  const std::vector<kapp::Plane>& slices, bool wallsOnly)
{
    QDir().mkpath(dir);
    QJsonArray zones, slicesJson;
    std::map<QString, std::pair<std::optional<double>, std::optional<double>>> ranges;
    std::map<QString, std::pair<QString, int>> meta;
    kapp::Bounds bounds;
    auto collect = [&](const kapp::ViewPart& part) {
        for (const auto& f : part.fields) {
            auto& r = ranges[f.name];
            if (f.min) r.first = r.first ? std::min(*r.first, *f.min) : *f.min;
            if (f.max) r.second = r.second ? std::max(*r.second, *f.max) : *f.max;
            meta[f.name] = {f.unit, f.components};
        }
    };
    int ordinal = 0;
    for (const CarZone& z : carZones()) {
        if (wallsOnly && (z.type != "wall" || z.name == "ground")) continue;
        kapp::MeshData m = carZoneMesh(z.name, prm, reference);
        bounds.expand(m.bounds());
        kapp::ViewPart part;
        part.kind = kapp::PartKind::Zone;
        part.id = z.id;
        part.name = z.name;
        part.type = z.type;
        part.file = QStringLiteral("zone-%1.surface.bin").arg(ordinal++);
        writeFile(QDir(dir).filePath(part.file), kapp::writeMeshBytes(m, &part));
        QJsonObject zj = kapp::toJson(part);
        zj["faces"] = double(m.indices.size() / 3);
        zones.append(zj);
        collect(part);
    }
    int n = 0;
    for (const kapp::Plane& s : slices) {
        kapp::MeshData m = carSlice(s.origin, s.normal, prm, reference);
        kapp::ViewPart part;
        part.kind = kapp::PartKind::Slice;
        part.origin = s.origin;
        part.normal = s.normal.normalized();
        part.file = QStringLiteral("slice-%1.bin").arg(n++);
        writeFile(QDir(dir).filePath(part.file), kapp::writeMeshBytes(m, &part));
        QJsonObject sj = kapp::toJson(part);
        sj.remove("name");
        slicesJson.append(sj);
        collect(part);
    }
    QJsonArray fields;
    for (const auto& [name, r] : ranges)
        fields.append(QJsonObject{{"name", name}, {"unit", meta[name].first}, {"components", meta[name].second}, {"location", "vertex"},
                                  {"min", r.first ? QJsonValue(*r.first) : QJsonValue()}, {"max", r.second ? QJsonValue(*r.second) : QJsonValue()}});
    writeJson(QDir(dir).filePath("manifest.json"),
              {{"v", 1}, {"format", QJsonObject{{"name", "kutta-view"}, {"version", 1}, {"endian", "little"}, {"header_bytes", 32}}},
               {"reference", reference}, {"provenance", prov}, {"synthetic", true},
               {"bounds", QJsonArray{QJsonArray{bounds.min.x(), bounds.min.y(), bounds.min.z()}, QJsonArray{bounds.max.x(), bounds.max.y(), bounds.max.z()}}},
               {"zones", zones}, {"slices", slicesJson}, {"fields", fields}});
}

void copyFixtureDir(const QString& rel, const QString& dest)
{
    QDir().mkpath(dest);
    const QDir src(QDir(fixtureDir()).filePath(rel));
    for (const QFileInfo& fi : src.entryInfoList(QDir::Files)) {
        QString name = fi.fileName();
        if (name.endsWith(".qz")) name.chop(3);
        writeFile(QDir(dest).filePath(name), readFixture(rel + "/" + name));
    }
}

QString csvNotice(const RunSpec& r)
{
    const QString p = provenance(r);
    if (p.isEmpty()) return {};
    QString lines;
    for (const QString& part : p.split(" | ")) lines += "# " + part + "\n";
    return lines;
}

int solve(const Args& a)
{
    const bool jsonl = a.value("progress") == "jsonl";
    if (a.positional.size() != 1 || !a.has("out") || !jsonl) usage(jsonl);
    RunSpec r;
    r.casePath = a.positional[0];
    r.kind = caseKind(r.casePath, true);
    r.dataPath = a.value("data");
    r.outDir = a.value("out");
    r.approximations = a.value("approximate").split(',', Qt::SkipEmptyParts);
    for (const QString& n : r.approximations)
        if (n != "wfgc") fail(QStringLiteral("unknown approximation '%1'; the only allowed name is wfgc").arg(n));
    if (r.approximations.removeDuplicates() > 0) fail("duplicate approximation names");
    r.iterations = a.has("iterations") ? a.value("iterations").toInt() : 1000;
    r.liveView = a.value("live-view", "0").toInt();
    if (r.iterations < 1) fail("--iterations must be at least 1");

    auto summary = kapp::parseCaseSummary(inspectJson(r.casePath, r.dataPath));
    if (!summary) fail(summary.error());

    if (a.has("setup")) {
        QFile f(a.value("setup"));
        if (!f.open(QIODevice::ReadOnly) || f.size() > (1 << 20)) fail("cannot read setup file (at most 1 MiB)");
        QJsonParseError pe;
        const QJsonDocument doc = QJsonDocument::fromJson(f.readAll(), &pe);
        const QJsonObject setup = doc.object();
        if (!doc.isObject() || setup.value("v").toInt() != 1 || !setup.value("overrides").isObject())
            fail("setup file must be {\"v\":1,\"overrides\":{...}}");
        for (const QString& k : setup.keys())
            if (k != "v" && k != "overrides") fail("unknown setup property " + k);
        kapp::SettingsModel model;
        model.reset(summary->settings);
        QJsonArray issues;
        const QJsonObject ov = setup.value("overrides").toObject();
        for (auto it = ov.begin(); it != ov.end(); ++it)
            if (QString e = model.validate(it.key(), it.value()); !e.isEmpty())
                issues.append(QJsonObject{{"path", it.key()}, {"setting", it.key()}, {"label", it.key()}, {"class", "unsupported"}, {"explanation", e}});
        if (!issues.isEmpty()) fail("Invalid setup overrides", issues);
        r.overrides = ov;
        if (ov.contains("run/iterations")) r.iterations = ov.value("run/iterations").toInt();
        if (ov.contains("run/monitor-interval")) r.monitor = ov.value("run/monitor-interval").toInt();
    }

    QJsonArray blocking;
    for (const auto* i : summary->blocking(r.approximations)) {
        blocking.append(QJsonObject{{"path", i->path}, {"setting", i->setting}, {"label", i->label},
                                    {"class", i->cls == kapp::IssueClass::Unsupported ? "unsupported" : "approximable"},
                                    {"approximation", i->approximation}, {"explanation", i->explanation}, {"evidence", i->evidence}});
    }
    if (!blocking.isEmpty()) {
        QString msg = "Fluent setup rejected:\n";
        for (const QJsonValue& b : blocking) {
            const QJsonObject o = b.toObject();
            msg += "  " + o.value("setting").toString() + ": " + o.value("label").toString();
            if (!o.value("approximation").toString().isEmpty())
                msg += "; opt in with --approximate " + o.value("approximation").toString();
            msg += "\n";
        }
        fail(msg, blocking);
    }

    const QJsonObject devs = devicesJson();
    const QString backend = a.value("backend", "auto");
    for (const QJsonValue& d : devs.value("devices").toArray()) {
        const QJsonObject o = d.toObject();
        if ((backend == "auto" && o.value("default").toBool()) || o.value("backend").toString() == backend) { r.device = o; break; }
    }
    if (r.device.isEmpty()) fail(QStringLiteral("backend %1 is not available on this machine").arg(backend));

    if (r.kind == CaseKind::Car) applyCarOverrides(r);
    watchStdin();
    QDir().mkpath(r.outDir);
    const double speed = std::max(0.1, qEnvironmentVariable("KUTTA_FIXTURE_SPEED", "6").toDouble());
    const int holdAt = qEnvironmentVariableIntValue("KUTTA_FIXTURE_HOLD_AT");
    const QString prov = provenance(r);

    emitJson({{"v", 1}, {"event", "setup"}, {"cells", double(r.kind == CaseKind::Elbow ? 17822 : 18412903)},
              {"backend", r.device.value("backend")}, {"device", r.device.value("name")},
              {"approximations", QJsonArray::fromStringList(r.approximations)},
              {"init", r.dataPath.isEmpty() ? "hybrid" : "fluent-data"}, {"overrides", r.overrides}, {"modified", !r.overrides.isEmpty()}});
    if (QString n = approximationNotice(r); !n.isEmpty()) emitJson({{"v", 1}, {"event", "warning"}, {"message", n}});
    if (r.kind == CaseKind::Car)
        emitJson({{"v", 1}, {"event", "warning"}, {"message", "SYNTHETIC FIXTURE: this run is generated by kutta-fixture for app development; no solver ran."}});
    else
        emitJson({{"v", 1}, {"event", "warning"}, {"message", "FIXTURE REPLAY: progress and fields were recorded from a real Kutta run of this case on Apple M4 Pro (Metal); no solver is running now."}});
    if (r.kind == CaseKind::Elbow && !r.overrides.isEmpty())
        emitJson({{"v", 1}, {"event", "warning"}, {"message", "FIXTURE REPLAY: overrides were validated, but the recorded numbers do not change."}});
    emitJson({{"v", 1}, {"event", "warning"}, {"message", "Case convergence criteria are not exposed by the adapter; stopping at the requested iteration limit."}});

    QElapsedTimer clock;
    clock.start();
    QStringList residualRows;
    QStringList residualNames;
    int lastIter = 0;
    bool stopped = false;
    auto snapshot = [&](int iter, double convergence) {
        const QString dir = QDir(r.outDir).filePath(QStringLiteral("live/iter-%1").arg(iter, 4, 10, QChar('0')));
        if (r.kind == CaseKind::Elbow) {
            // Recorded snapshots exist for 40, 120 and 300 iterations.
            const QString rel = QStringLiteral("cases/mixing_elbow/live/iter-%1").arg(iter, 4, 10, QChar('0'));
            if (!QFileInfo(QDir(fixtureDir()).filePath(rel)).isDir()) return;
            copyFixtureDir(rel, dir);
        } else {
            CarParams prm = r.car;
            prm.convergence = convergence;
            writeCarView(dir, prm, false, prov, {}, true);
        }
        emitJson({{"v", 1}, {"event", "snapshot"}, {"iter", iter}, {"view", QDir(dir).absoluteFilePath("manifest.json")}});
    };
    auto onIteration = [&](QJsonObject ev) {
        lastIter = ev.value("iter").toInt();
        if (residualNames.isEmpty()) residualNames = kapp::residualOrder(ev.value("residuals").toObject().keys());
        QStringList row{QString::number(lastIter)};
        for (const QString& n : residualNames) row << QString::number(ev.value("residuals").toObject().value(n).toDouble(), 'g', 17);
        row << QString::number(ev.value("sec_per_iter").toDouble(), 'g', 17);
        residualRows << row.join(',');
        emitJson(ev);
        if (holdAt > 0 && lastIter >= holdAt)
            while (!g_stop) sleepMs(20);
    };

    std::vector<std::pair<QString, QVector3D>> finalForce;  // car reports
    if (r.kind == CaseKind::Elbow) {
        const QString variant = r.dataPath.isEmpty() ? "cold" : "seeded";
        const QList<QByteArray> lines = readFixture("cases/mixing_elbow/progress-" + variant + ".jsonl").split('\n');
        double lastElapsed = 0;
        const int snaps[] = {40, 120, 300};
        std::size_t nextSnap = 0;
        for (const QByteArray& line : lines) {
            const QJsonObject ev = QJsonDocument::fromJson(line).object();
            if (ev.value("event").toString() != "iteration") continue;
            const int iter = ev.value("iter").toInt();
            if (iter > r.iterations) break;
            const double el = ev.value("elapsed").toDouble();
            sleepMs((el - lastElapsed) * 1000 / speed);
            lastElapsed = el;
            onIteration(ev);
            while (r.liveView > 0 && nextSnap < 3 && iter >= snaps[nextSnap]) snapshot(snaps[nextSnap++], 1);
            if (g_stop) { stopped = true; break; }
        }
        const QString run = "cases/mixing_elbow/run-" + variant;
        for (const char* f : {"forces.csv", "fluxes.csv"}) {
            QByteArray csv = readFixture(run + "/" + f);
            // Recorded outputs keep their own headers; add the replay notice.
            csv.prepend(QStringLiteral("# FIXTURE REPLAY: recorded final values from a 600-iteration Kutta run\n%1")
                            .arg(r.overrides.isEmpty() ? QString() : "# MODIFIED RESULT: case settings overridden; recorded values did not change\n").toUtf8());
            writeFile(QDir(r.outDir).filePath(f), csv);
        }
    } else {
        const auto forces = carForces(r.car);
        QVector3D total, fw, rw;
        for (const ZoneForce& f : forces) {
            total += f.pressure + f.viscous;
            if (f.zone == "front-wing") fw = f.pressure + f.viscous;
            if (f.zone == "rear-wing") rw = f.pressure + f.viscous;
        }
        const double qa = 0.5 * r.car.refDensity * r.car.refVelocity * r.car.refVelocity * r.car.refArea;
        finalForce = {{"downforce", total}, {"drag", total}, {"front-wing-downforce", fw}, {"rear-wing-downforce", rw}};
        const double stepMs = 1000.0 / (150.0 * speed / 6.0);
        QElapsedTimer pace;
        pace.start();
        for (int it = r.monitor; it <= r.iterations; it += r.monitor) {
            const double t = it;
            const double osc = 1 + 0.22 * std::sin(t / 13.0) * std::exp(-t / 600);
            auto res = [&](double a0, double tau, double floor) { return (a0 * std::exp(-t / tau) + floor * (1 + 0.15 * std::sin(t / 7.0))) * osc; };
            QJsonObject resid{{"continuity", res(0.9, 70, 6e-5)}, {"x-velocity", res(0.08, 90, 4e-6)}, {"y-velocity", res(0.05, 95, 3e-6)},
                              {"z-velocity", res(0.06, 85, 3.5e-6)}, {"k", res(0.3, 110, 2e-5)}, {"omega", res(0.2, 120, 1.5e-5)}};
            const double settle = 1 - std::exp(-t / 110) * std::cos(t / 45);
            QJsonObject reports;
            for (const auto& [name, f] : finalForce) {
                const double proj = name.contains("downforce") ? -double(f.z()) : double(f.x());
                const double val = proj * settle * (1 + 0.004 * std::sin(t / 9.0) * std::exp(-t / 400));
                reports[name] = QJsonObject{{"force", val}, {"coefficient", val / qa}};
            }
            sleepMs(stepMs * r.monitor - pace.nsecsElapsed() / 1e6);
            pace.restart();
            const double el = clock.nsecsElapsed() / 1e9;
            onIteration({{"v", 1}, {"event", "iteration"}, {"iter", it}, {"residuals", resid}, {"reports", reports},
                         {"sec_per_iter", el / it}, {"elapsed", el}});
            if (r.liveView > 0 && it % r.liveView == 0) snapshot(it, 1 - std::exp(-t / 150));
            if (g_stop) { stopped = true; break; }
        }
        // forces.csv in the CLI's layout.
        QString csv = csvNotice(r);
        csv += "kind,name,zone_id,pressure_x,pressure_y,pressure_z,viscous_x,viscous_y,viscous_z,total_x,total_y,total_z,direction_x,direction_y,direction_z,pressure_projection,viscous_projection,total_projection,coefficient\n";
        const double settle = 1 - std::exp(-double(lastIter) / 110) * std::cos(lastIter / 45.0);
        for (const ZoneForce& f : forces) {
            const CarZone* z = nullptr;
            for (const auto& cz : carZones()) if (cz.name == f.zone) z = &cz;
            const QVector3D p = f.pressure * float(settle), v = f.viscous * float(settle), t = p + v;
            csv += QStringLiteral("wall,\"%1\",%2,%3,%4,%5,%6,%7,%8,%9,%10,%11,,,,,,,\n").arg(f.zone).arg(z ? z->id : -1)
                       .arg(p.x()).arg(p.y()).arg(p.z()).arg(v.x()).arg(v.y()).arg(v.z()).arg(t.x()).arg(t.y()).arg(t.z());
        }
        for (const auto& [name, f0] : finalForce) {
            const QVector3D f = f0 * float(settle);
            const QVector3D dir = name.contains("downforce") ? QVector3D(0, 0, -1) : QVector3D(1, 0, 0);
            const double proj = QVector3D::dotProduct(f, dir);
            csv += QStringLiteral("report,\"%1\",,%2,%3,%4,0,0,0,%2,%3,%4,%5,%6,%7,%8,0,%8,%9\n").arg(name)
                       .arg(f.x()).arg(f.y()).arg(f.z()).arg(dir.x()).arg(dir.y()).arg(dir.z()).arg(proj).arg(proj / qa);
        }
        writeFile(QDir(r.outDir).filePath("forces.csv"), csv.toUtf8());
        writeFile(QDir(r.outDir).filePath("fluxes.csv"),
                  (csvNotice(r) + "name,zone_id,mass_flow_outward\n\"inlet\",2," + QString::number(-r.car.density * r.car.speed * 32, 'g', 10) +
                   "\n\"outlet\",3," + QString::number(r.car.density * r.car.speed * 32, 'g', 10) + "\n").toUtf8());
    }
    QString res = csvNotice(r) + "iteration," + residualNames.join(',') + ",seconds_per_iteration\n" + residualRows.join('\n') + "\n";
    writeFile(QDir(r.outDir).filePath("residuals.csv"), res.toUtf8());
    writeJson(QDir(r.outDir).filePath("fixture-run.json"), runDescriptor(r));

    const double elapsed = clock.nsecsElapsed() / 1e9;
    const QString stopReason = stopped ? "stopped" : "iteration_limit";
    QJsonObject sum;
    if (r.kind == CaseKind::Elbow)
        sum = QJsonDocument::fromJson(readFixture(QStringLiteral("cases/mixing_elbow/run-%1/summary.json").arg(r.dataPath.isEmpty() ? "cold" : "seeded"))).object();
    sum["v"] = 1;
    sum["fixture"] = true;
    sum["synthetic"] = r.kind == CaseKind::Car;
    sum["modified"] = !r.overrides.isEmpty();
    sum["overrides"] = r.overrides;
    sum["approximate"] = !r.approximations.isEmpty();
    if (r.kind == CaseKind::Car && !r.approximations.isEmpty()) {
        sum["approximation_notice"] = approximationNotice(r);
        sum["approximations"] = QJsonArray{QJsonObject{{"name", "wfgc"}, {"setting", "recon/cell-lsf?"}, {"fluent_value", "#t"},
                                                       {"kutta_action", "run without warped-face gradient correction (WFGC off, #f)"},
                                                       {"effect_documentation", "src/kutta/fv/VERIFICATION.md#v6-initialization-comparison"},
                                                       {"measured_effect", "not measured on this synthetic case"}}};
    }
    sum["case_file"] = QFileInfo(r.casePath).absoluteFilePath();
    sum["data_file"] = r.dataPath.isEmpty() ? QJsonValue() : QJsonValue(QFileInfo(r.dataPath).absoluteFilePath());
    sum["backend"] = r.device.value("backend");
    sum["device"] = r.device.value("name");
    sum["iterations"] = lastIter;
    sum["stop_reason"] = stopReason;
    sum["converged"] = false;
    sum["replay_seconds"] = elapsed;
    if (r.kind == CaseKind::Car) {
        sum["cells"] = 18412903;
        sum["wall_time_seconds"] = elapsed;
        sum["kutta_version"] = "fixture";
    }
    writeJson(QDir(r.outDir).filePath("summary.json"), sum);
    emitJson({{"v", 1}, {"event", "done"}, {"converged", false}, {"iterations", lastIter}, {"elapsed", elapsed},
              {"out", r.outDir}, {"stop_reason", stopReason}});
    return 0;
}

std::vector<kapp::Plane> parseSlices(const QStringList& args)
{
    std::vector<kapp::Plane> out;
    for (const QString& s : args) {
        const QStringList v = s.split(',');
        if (v.size() != 6) fail("--slice takes ox,oy,oz,nx,ny,nz", {}, false);
        double x[6];
        for (int i = 0; i < 6; ++i) {
            bool ok = false;
            x[i] = v[i].toDouble(&ok);
            if (!ok || !std::isfinite(x[i])) fail("--slice values must be finite numbers", {}, false);
        }
        const QVector3D n{float(x[3]), float(x[4]), float(x[5])};
        if (n.isNull()) fail("--slice normal must be nonzero", {}, false);
        out.push_back({QVector3D(float(x[0]), float(x[1]), float(x[2])), n.normalized()});
    }
    return out;
}

int exportView(const Args& a)
{
    if (a.positional.size() != 1) usage(false);
    const QString dir = a.positional[0];
    const QJsonObject desc = readJson(QDir(dir).filePath("fixture-run.json"));
    if (desc.isEmpty()) fail(dir + " is not a finished fixture run (no summary)", {}, false);
    const auto slices = parseSlices(a.opts.value("slice"));
    const QString view = QDir(dir).filePath("view");
    QDir(view).removeRecursively();
    if (desc.value("case").toString() == "mixing_elbow") {
        copyFixtureDir("cases/mixing_elbow/run-" + desc.value("variant").toString() + "/view", view);
        if (!slices.empty())
            eprint("kutta-fixture: recorded elbow view has one slice (z = 0.025 m); requested slices need the real CLI");
        QJsonObject m = readJson(QDir(view).filePath("manifest.json"));
        if (!desc.value("provenance").toString().isEmpty()) m["provenance"] = desc.value("provenance").toString();
        writeJson(QDir(view).filePath("manifest.json"), m);
    } else {
        writeCarView(view, paramsFrom(desc), desc.value("reference").toBool(), desc.value("provenance").toString(), slices, false);
    }
    eprint("Wrote " + view);
    return 0;
}

QVector3D vec(const QJsonValue& v)
{
    const QJsonArray a = v.toArray();
    return a.size() == 3 ? QVector3D(float(a[0].toDouble()), float(a[1].toDouble()), float(a[2].toDouble())) : QVector3D();
}

int post(const Args& a)
{
    if (a.positional.size() != 1 || !a.has("request")) usage(false);
    const QString dir = a.positional[0];
    const QJsonObject desc = readJson(QDir(dir).filePath("fixture-run.json"));
    if (desc.isEmpty()) fail(dir + " is not a finished fixture run", {}, false);
    const QJsonObject req = readJson(a.value("request"));
    if (req.value("v").toInt() != 1) fail("post request must have v: 1", {}, false);
    const bool car = desc.value("case").toString() == "fsae_car";
    const CarParams prm = paramsFrom(desc);
    const QString prov = desc.value("provenance").toString();
    // Elbow sections/integrals use the exported view surfaces.
    auto viewManifest = kapp::loadViewManifest(QDir(dir).filePath("view/manifest.json"));
    for (const QJsonValue& iv : req.value("items").toArray()) {
        const QJsonObject item = iv.toObject();
        const QString id = item.value("id").toString(), kind = item.value("kind").toString();
        if (id.isEmpty() || id.contains('/') || id.contains('\\') || id.startsWith('.')) fail("post item needs a plain id", {}, false);
        const QString outDir = QDir(dir).filePath("post/" + id);
        QDir(outDir).removeRecursively();
        QDir().mkpath(outDir);
        QJsonObject manifest{{"v", 1}, {"id", id}, {"kind", kind}, {"provenance", prov}, {"synthetic", car}};
        auto needsCar = [&] {
            if (!car) fail(kind + " on the recorded elbow needs the real kutta CLI (the fixture has no volume fields)", {}, false);
        };
        if (kind == "streamlines") {
            needsCar();
            const QJsonObject seed = item.value("seed").toObject();
            kapp::MeshData m = carStreamlines(vec(seed.value("from")), vec(seed.value("to")), seed.value("count").toInt(30),
                                              item.value("max_length").toDouble(8), prm);
            kapp::ViewPart part;
            part.kind = kapp::PartKind::Streamlines;
            part.name = item.value("name").toString("Pathlines");
            part.file = "streamlines.bin";
            writeFile(QDir(outDir).filePath(part.file), kapp::writeMeshBytes(m, &part));
            QJsonObject pj = kapp::toJson(part);
            pj["kind"] = "streamlines";
            manifest["parts"] = QJsonArray{pj};
        } else if (kind == "iso-surface") {
            needsCar();
            const QString field = item.value("field").toString("total_pressure_coefficient");
            kapp::MeshData m = carIsoSurface(field, item.value("value").toDouble(0), prm);
            kapp::ViewPart part;
            part.kind = kapp::PartKind::IsoSurface;
            part.name = item.value("name").toString(QStringLiteral("%1 = %2").arg(field).arg(item.value("value").toDouble()));
            part.file = "iso.bin";
            writeFile(QDir(outDir).filePath(part.file), kapp::writeMeshBytes(m, &part));
            QJsonObject pj = kapp::toJson(part);
            pj["kind"] = "iso-surface";
            manifest["parts"] = QJsonArray{pj};
        } else if (kind == "line-probe") {
            needsCar();
            QStringList fields;
            for (const QJsonValue& f : item.value("fields").toArray()) fields << f.toString();
            if (fields.isEmpty()) fields << "velocity_magnitude";
            const ProbeResult pr = carLineProbe(vec(item.value("from")), vec(item.value("to")), item.value("samples").toInt(200), fields, prm);
            QJsonArray series, x, pts;
            for (double d : pr.distance) x.append(d);
            for (const auto& p : pr.points) pts.append(QJsonArray{p.x(), p.y(), p.z()});
            for (const auto& s : pr.series) {
                QJsonArray y;
                for (double v : s.y) y.append(std::isfinite(v) ? QJsonValue(v) : QJsonValue());
                series.append(QJsonObject{{"name", s.name}, {"unit", s.unit}, {"x", x}, {"y", y}});
            }
            manifest["x_label"] = "Distance along line";
            manifest["x_unit"] = "m";
            manifest["series"] = series;
            manifest["points"] = pts;
        } else if (kind == "section") {
            if (!viewManifest) fail("export the view before requesting sections: " + viewManifest.error(), {}, false);
            const QString field = item.value("field").toString("pressure_coefficient");
            const QVector3D origin = vec(item.value("origin")), normal = vec(item.value("normal"));
            const int axis = item.value("x_axis").toInt(0);
            QJsonArray series;
            for (const QJsonValue& zv : item.value("zones").toArray()) {
                const QString zone = zv.toString();
                for (const auto& part : viewManifest->zones) {
                    if (part.name != zone) continue;
                    auto mesh = kapp::readMeshFile(QDir(viewManifest->dir).filePath(part.file), part);
                    if (!mesh) fail(mesh.error(), {}, false);
                    const auto pts = kapp::planeSection(mesh.value(), field, -1, origin, normal);
                    QJsonArray x, y;
                    for (const auto& p : pts) {
                        x.append(p.position[axis]);
                        y.append(std::isfinite(p.value) ? QJsonValue(p.value) : QJsonValue());
                    }
                    const auto* fd = part.field(field);
                    series.append(QJsonObject{{"name", zone}, {"unit", fd ? fd->unit : ""}, {"field", field}, {"x", x}, {"y", y}, {"segments", true}});
                }
            }
            manifest["x_label"] = QString("Position (%1)").arg(QChar('x' + axis));
            manifest["x_unit"] = "m";
            manifest["series"] = series;
        } else if (kind == "surface-integral") {
            if (!viewManifest) fail("export the view before requesting integrals: " + viewManifest.error(), {}, false);
            const QString field = item.value("field").toString("pressure"), report = item.value("report").toString("area-weighted-average");
            const kapp::IntegralKind ik = report == "area" ? kapp::IntegralKind::Area : report == "integral" ? kapp::IntegralKind::Integral
                                        : report == "minimum" ? kapp::IntegralKind::Minimum : report == "maximum" ? kapp::IntegralKind::Maximum
                                                              : kapp::IntegralKind::AreaWeightedAverage;
            QJsonArray rows;
            QString unit;
            double totalArea = 0, totalInt = 0;
            for (const QJsonValue& zv : item.value("zones").toArray()) {
                for (const auto& part : viewManifest->zones) {
                    if (part.name != zv.toString()) continue;
                    auto mesh = kapp::readMeshFile(QDir(viewManifest->dir).filePath(part.file), part);
                    if (!mesh) fail(mesh.error(), {}, false);
                    if (const auto* fd = part.field(field)) unit = fd->unit;
                    const auto r = kapp::integrate(mesh.value(), field, ik);
                    const auto area = kapp::integrate(mesh.value(), {}, kapp::IntegralKind::Area);
                    rows.append(QJsonObject{{"zone", part.name}, {"value", r ? QJsonValue(r->value) : QJsonValue()}, {"area", area ? area->value : 0}});
                    if (r && area) { totalArea += area->value; totalInt += ik == kapp::IntegralKind::AreaWeightedAverage ? r->value * area->value : r->value; }
                }
            }
            manifest["field"] = field;
            manifest["report"] = report;
            manifest["unit"] = ik == kapp::IntegralKind::Area ? "m^2" : ik == kapp::IntegralKind::Integral ? unit + " m^2" : unit;
            manifest["rows"] = rows;
            manifest["total"] = ik == kapp::IntegralKind::AreaWeightedAverage ? (totalArea > 0 ? totalInt / totalArea : 0) : totalInt;
        } else {
            fail("unknown post-processing kind " + kind, {}, false);
        }
        writeJson(QDir(outDir).filePath("manifest.json"), manifest);
    }
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    QCoreApplication app(argc, argv);
    const QStringList args = app.arguments().mid(1);
    if (args.isEmpty()) usage(false);
    const QString cmd = args[0];
    const Args a = parseArgs(args.mid(1));
    if (cmd == "devices") {
        if (!a.has("json")) usage(false);
        out() << QJsonDocument(devicesJson()).toJson(QJsonDocument::Compact) << '\n';
        return 0;
    }
    if (cmd == "inspect") {
        if (a.positional.size() != 1 || !a.has("json")) usage(false);
        caseKind(a.positional[0], false);
        out() << QJsonDocument(inspectJson(a.positional[0], a.value("data"))).toJson(QJsonDocument::Compact) << '\n';
        return 0;
    }
    if (cmd == "solve") return solve(a);
    if (cmd == "export-view") return exportView(a);
    if (cmd == "post") return post(a);
    usage(false);
}
