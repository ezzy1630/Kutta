// kutta-fixture: a stand-in for the `kutta` CLI that speaks docs/app-protocol.md
// (including export-view --post and solve --live-view) from recorded and
// synthetic data. The app talks to it through the same process boundary as the
// real CLI; switching is a path change.
//
//   <name>.cas[.h5]       replays a REAL recorded Kutta run of that case when
//                         fixtures/cases/<name>/ exists (app/fixtures/README.md)
//   fsae_car*.cas.h5      a SYNTHETIC car; every output says "synthetic": true
//
// Environment:
//   KUTTA_FIXTURE_DIR         fixture data directory
//   KUTTA_FIXTURE_SPEED       replay speed multiplier (default 6)
//   KUTTA_FIXTURE_HOLD_AT     pause after this iteration until stdin closes
//   KUTTA_FIXTURE_FAIL_AT     fail after this iteration: an error event with
//                             KUTTA_FIXTURE_FAIL_MESSAGE and exit 1, or, with
//                             no message, abort as a crash would
//   KUTTA_FIXTURE_DEVICES     "macos" | "windows" device inventory
//   KUTTA_FIXTURE_CAR_DETAIL  car tessellation multiplier (1 ~ 0.4M triangles)

#include "Aero.h"
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

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <map>
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
                        "--progress jsonl [--live-view DIR [--live-every N]] | export-view DIR [--post FILE] [--slice ox,oy,oz,nx,ny,nz]... "
                        "| export-geometry CASE --out DIR [--json] | aero DIR [--aero SPEC] [--json]";
    if (jsonl) emitJson({{"v", 1}, {"event", "error"}, {"message", msg}, {"issues", QJsonArray{}}});
    eprint(msg);
    std::exit(2);
}

// Invalid request syntax: exit 2, like the real CLI.
[[noreturn]] void badRequest(const QString& message)
{
    eprint("kutta-fixture: " + message);
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

// The pressure-velocity coupling a run uses, as the CLI's setup event reports
// it: the override or case value, except that Coupled without the
// pseudo-transient formulation runs SIMPLEC.
QString couplingRun(const QJsonObject& inspected, const QJsonObject& overrides)
{
    QString scheme;
    bool pseudo = true;
    for (const QJsonValue& v : inspected.value("settings").toArray()) {
        const QJsonObject o = v.toObject();
        if (o.value("path") == "methods/pressure-velocity-coupling") scheme = o.value("value").toString();
        if (o.value("path") == "methods/pseudo-transient") pseudo = o.value("value").toBool(true);
    }
    scheme = overrides.value("methods/pressure-velocity-coupling").toString(scheme);
    return scheme == "coupled" && !pseudo ? QStringLiteral("simplec") : scheme;
}

// Run controls (iterations, monitor interval, convergence) never make a result MODIFIED.
QJsonObject caseOverrides(const QJsonObject& overrides)
{
    QJsonObject o;
    for (auto it = overrides.begin(); it != overrides.end(); ++it)
        if (!kapp::SettingsModel::isRunControl(it.key())) o.insert(it.key(), it.value());
    return o;
}
bool caseModified(const QJsonObject& overrides) { return !caseOverrides(overrides).isEmpty(); }

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

enum class CaseKind { Recorded, Car };

// "elbow_param.cas.h5" -> "elbow_param": the recording's directory name.
QString caseStem(const QString& casePath)
{
    QString name = QFileInfo(casePath).fileName();
    for (const char* ext : {".h5", ".gz", ".cas", ".dat"})
        if (name.endsWith(QLatin1String(ext))) name.chop(int(qstrlen(ext)));
    return name;
}

QStringList recordedCases()
{
    return QDir(QDir(fixtureDir()).filePath("cases")).entryList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
}

CaseKind caseKind(const QString& casePath, bool jsonl)
{
    const QString name = QFileInfo(casePath).fileName();
    if (recordedCases().contains(caseStem(casePath))) return CaseKind::Recorded;
    if (name.startsWith("fsae_car")) return CaseKind::Car;
    fail(QStringLiteral("The fixture CLI has no data for %1. It replays recorded runs of %2 and generates "
                        "fsae_car.cas.h5 (synthetic); point the app at the real kutta CLI for other cases.")
             .arg(name, recordedCases().join(", ")),
         {}, jsonl);
}

// A recorded case's files: cases/<stem>/<rel>.
QString recorded(const QString& stem, const QString& rel) { return QStringLiteral("cases/%1/%2").arg(stem, rel); }

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
    if (caseKind(casePath, false) == CaseKind::Recorded) {
        const QString rel = recorded(caseStem(casePath), data ? "inspect-data.json" : "inspect.json");
        o = QJsonDocument::fromJson(readFixture(rel)).object();
        if (o.isEmpty())
            fail(data ? QStringLiteral("The fixture has no recording of %1 with a data file; open it without one.").arg(QFileInfo(casePath).fileName())
                      : "fixture data missing: " + rel, {}, false);
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
    QString casePath, dataPath, outDir, stem;
    QStringList approximations;
    int iterations = 1000, monitor = 10, liveView = 0;  // liveView: --live-every, 0 when off
    QString liveDir;
    QJsonObject overrides;
    CarParams car;
    QJsonObject device;
};

// The recorded run's summary: what the real CLI wrote for this variant.
QJsonObject recordedSummary(const RunSpec& r)
{
    return QJsonDocument::fromJson(readFixture(recorded(r.stem, QStringLiteral("run-%1/summary.json").arg(r.dataPath.isEmpty() ? "cold" : "seeded")))).object();
}

QString approximationNotice(const RunSpec& r)
{
    if (r.approximations.isEmpty()) return {};
    // A recording ran with exactly the approximations its case needs.
    if (r.kind == CaseKind::Recorded) return recordedSummary(r).value("approximation_notice").toString();
    QString n = "APPROXIMATE RESULT: not a faithful Fluent reproduction; wfgc: recon/cell-lsf? Fluent #t -> Kutta run without "
                "warped-face gradient correction (WFGC off, #f); V6 CPU hybrid diagnostic versus Fluent: velocity L2 3.366748%, "
                "pressure L2 2.665514%, k L2 14.722368%, omega L2 2.805466%, pressure-drop error 1.973008%; combined solver "
                "discrepancy, not isolated WFGC error or a general bound; evidence: src/kutta/fv/VERIFICATION.md#v6-initialization-comparison";
    return n;
}

QString provenance(const RunSpec& r)
{
    QStringList parts;
    if (QString a = approximationNotice(r); !a.isEmpty()) parts << a;
    if (caseModified(r.overrides)) {
        QStringList keys = caseOverrides(r.overrides).keys();
        parts << "MODIFIED RESULT: case settings overridden (" + keys.join(", ") + ")";
    }
    if (r.kind == CaseKind::Car) parts << "SYNTHETIC FIXTURE: invented geometry and fields for app development; not a CFD result";
    else if (caseModified(r.overrides)) parts << "FIXTURE REPLAY: recorded " + r.stem + " fields; overrides did not change these numbers";
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
    num("reference/length", r.car.refLength);
}

QJsonObject runDescriptor(const RunSpec& r)
{
    return {{"case", r.kind == CaseKind::Recorded ? r.stem : QStringLiteral("fsae_car")},
            {"recorded", r.kind == CaseKind::Recorded},
            {"variant", r.dataPath.isEmpty() ? "cold" : "seeded"},
            {"reference", !r.dataPath.isEmpty()},
            {"provenance", provenance(r)},
            {"synthetic", r.kind == CaseKind::Car},
            {"params", QJsonObject{{"speed", r.car.speed}, {"density", r.car.density}, {"viscosity", r.car.viscosity},
                                   {"ref_area", r.car.refArea}, {"ref_density", r.car.refDensity},
                                   {"ref_velocity", r.car.refVelocity}, {"ref_pressure", r.car.refPressure},
                                   {"ref_length", r.car.refLength}}}};
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
    c.refLength = p.value("ref_length").toDouble(c.refLength);
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
    r.stem = caseStem(r.casePath);
    r.dataPath = a.value("data");
    r.outDir = a.value("out");
    r.approximations = a.value("approximate").split(',', Qt::SkipEmptyParts);
    static const QStringList allowed{"wfgc", "compressibility", "ke-options", "reverse-flow"};
    for (const QString& n : r.approximations)
        if (!allowed.contains(n)) fail(QStringLiteral("unknown --approximate setting '%1'; allowed: %2").arg(n, allowed.join(',')));
    if (r.approximations.removeDuplicates() > 0) fail("duplicate approximation names");
    r.iterations = a.has("iterations") ? a.value("iterations").toInt() : 1000;
    r.liveDir = a.value("live-view");
    if (a.has("live-every") && r.liveDir.isEmpty()) usage(true);
    r.liveView = r.liveDir.isEmpty() ? 0 : a.value("live-every", "20").toInt();
    if (!r.liveDir.isEmpty() && r.liveView < 1) fail("--live-every must be a positive integer");
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
    const QJsonObject inspected = inspectJson(r.casePath, r.dataPath, r.car);
    const ConvergenceConfig conv = convergenceFrom(inspected.value("settings").toArray(), r.overrides);
    ConvergenceCheck convergence(conv);
    watchStdin();
    QDir().mkpath(r.outDir);
    const double speed = std::max(0.1, qEnvironmentVariable("KUTTA_FIXTURE_SPEED", "6").toDouble());
    const int holdAt = qEnvironmentVariableIntValue("KUTTA_FIXTURE_HOLD_AT");
    const int failAt = qEnvironmentVariableIntValue("KUTTA_FIXTURE_FAIL_AT");
    const QString prov = provenance(r);

    const QString liveRunId = QStringLiteral("run-%1-0").arg(QDateTime::currentMSecsSinceEpoch());
    const QString liveRoot = r.liveDir.isEmpty() ? QString() : QDir(r.liveDir).absoluteFilePath("live");
    QJsonObject setup{{"v", 1}, {"event", "setup"}, {"cells", double(r.kind == CaseKind::Recorded ? summary->cells : 18412903)},
                      {"backend", r.device.value("backend")}, {"device", r.device.value("name")},
                      {"approximations", QJsonArray::fromStringList(r.approximations)},
                      {"init", r.dataPath.isEmpty() ? "hybrid" : "fluent-data"}, {"overrides", r.overrides}, {"modified", caseModified(r.overrides)},
                      {"convergence", conv.toJson()}, {"pressure_velocity_coupling", couplingRun(inspected, r.overrides)}};
    if (r.liveView > 0)
        setup["live_view"] = QJsonObject{{"directory", liveRoot}, {"run_id", liveRunId}, {"every", r.liveView},
                                         {"unit", "iterations"}, {"budget_fraction", 0.05}};
    emitJson(setup);
    if (QString n = approximationNotice(r); !n.isEmpty()) emitJson({{"v", 1}, {"event", "warning"}, {"message", n}});
    if (r.kind == CaseKind::Car)
        emitJson({{"v", 1}, {"event", "warning"}, {"message", "SYNTHETIC FIXTURE: this run is generated by kutta-fixture for app development; no solver ran."}});
    else
        emitJson({{"v", 1}, {"event", "warning"}, {"message", "FIXTURE REPLAY: progress and fields were recorded from a real Kutta run of this case on Apple M4 Pro (Metal); no solver is running now."}});
    if (r.kind == CaseKind::Recorded && !r.overrides.isEmpty())
        emitJson({{"v", 1}, {"event", "warning"}, {"message", "FIXTURE REPLAY: overrides were validated, but the recorded numbers do not change."}});
    if (!conv.enabled)
        emitJson({{"v", 1}, {"event", "warning"}, {"message", "Convergence checking is off (monitors/check-convergence); stopping at the requested iteration limit."}});
    else if (conv.residuals.isEmpty() && conv.reports.empty())
        emitJson({{"v", 1}, {"event", "warning"}, {"message", "Convergence checking is on but no residual is checked and no report condition is active; stopping at the requested iteration limit."}});

    QElapsedTimer clock;
    clock.start();
    QStringList residualRows;
    QStringList residualNames;
    int lastIter = 0;
    bool stopped = false, converged = false;
    int liveSequence = 0;
    // Publish like the CLI: complete the snapshot directory, then swap the pointer.
    auto snapshot = [&](int iter, double convergence) -> void {
        QElapsedTimer cost;
        cost.start();
        const QString rel = QStringLiteral("%1/snapshot-%2").arg(liveRunId).arg(liveSequence + 1);
        const QString dir = QDir(liveRoot).filePath(rel), staging = dir + ".staging";
        QDir(staging).removeRecursively();
        if (r.kind == CaseKind::Recorded) {
            const QString src = recorded(r.stem, QStringLiteral("live/iter-%1").arg(iter, 4, 10, QChar('0')));
            if (!QFileInfo(QDir(fixtureDir()).filePath(src)).isDir()) return;
            copyFixtureDir(src, staging);
        } else {
            CarParams prm = r.car;
            prm.convergence = convergence;
            writeCarView(staging, prm, false, prov, {}, true);
        }
        if (!QDir().rename(staging, dir)) fail("cannot publish live snapshot " + dir);
        ++liveSequence;
        const double elapsed = clock.nsecsElapsed() / 1e9;
        const QString pointer = QDir(liveRoot).filePath("manifest.json");
        writeJson(pointer + ".tmp", {{"v", 1}, {"run_id", liveRunId}, {"sequence", liveSequence}, {"iter", iter}, {"elapsed", elapsed},
                                     {"manifest", rel + "/manifest.json"}});
        QFile::remove(pointer);
        QFile::rename(pointer + ".tmp", pointer);
        emitJson({{"v", 1}, {"event", "snapshot"}, {"iter", iter}, {"manifest", QDir(dir).absoluteFilePath("manifest.json")},
                  {"run_id", liveRunId}, {"sequence", liveSequence}, {"snapshot_seconds", cost.nsecsElapsed() / 1e9}, {"elapsed", elapsed}});
    };
    // sample: a monitor sample. Iterations 1-5 between samples carry residuals
    // only and are not convergence samples, as in the CLI.
    auto onIteration = [&](QJsonObject ev, bool sample) {
        lastIter = ev.value("iter").toInt();
        if (conv.enabled && sample) {
            QMap<QString, double> res, reps;
            const QJsonObject ro = ev.value("residuals").toObject(), po = ev.value("reports").toObject();
            for (auto it = ro.begin(); it != ro.end(); ++it) res[it.key()] = it.value().toDouble();
            for (auto it = po.begin(); it != po.end(); ++it) reps[it.key()] = it.value().toObject().value("force").toDouble();
            ev["convergence"] = convergence.sample(lastIter, res, reps);
            converged = convergence.met();
        }
        if (residualNames.isEmpty()) residualNames = kapp::residualOrder(ev.value("residuals").toObject().keys());
        QStringList row{QString::number(lastIter)};
        for (const QString& n : residualNames) row << QString::number(ev.value("residuals").toObject().value(n).toDouble(), 'g', 17);
        row << QString::number(ev.value("sec_per_iter").toDouble(), 'g', 17);
        residualRows << row.join(',');
        emitJson(ev);
        if (holdAt > 0 && lastIter >= holdAt)
            while (!g_stop) sleepMs(20);
        if (failAt > 0 && lastIter >= failAt) {
            const QString msg = qEnvironmentVariable("KUTTA_FIXTURE_FAIL_MESSAGE");
            if (msg.isEmpty()) std::abort();
            emitJson({{"v", 1}, {"event", "error"}, {"message", msg}, {"issues", QJsonArray{}}});
            std::exit(1);
        }
    };

    std::vector<std::pair<QString, QVector3D>> finalForce;  // car reports
    if (r.kind == CaseKind::Recorded) {
        const QString variant = r.dataPath.isEmpty() ? "cold" : "seeded";
        const QList<QByteArray> lines = readFixture(recorded(r.stem, "progress-" + variant + ".jsonl")).split('\n');
        const int recordedIterations = recordedSummary(r).value("iterations").toInt();
        double lastElapsed = 0;
        // Recorded live snapshots, in iteration order.
        std::vector<int> snaps;
        for (const QString& d : QDir(QDir(fixtureDir()).filePath(recorded(r.stem, "live"))).entryList({"iter-*"}, QDir::Dirs, QDir::Name))
            snaps.push_back(d.mid(5).toInt());
        std::size_t nextSnap = 0;
        int interval = 0;  // the recording's monitor interval: its first sample after iteration 5
        for (const QByteArray& line : lines) {
            const QJsonObject ev = QJsonDocument::fromJson(line).object();
            if (ev.value("event").toString() == "iteration" && ev.value("iter").toInt() > 5) { interval = ev.value("iter").toInt(); break; }
        }
        for (const QByteArray& line : lines) {
            const QJsonObject ev = QJsonDocument::fromJson(line).object();
            if (ev.value("event").toString() != "iteration") continue;
            const int iter = ev.value("iter").toInt();
            if (iter > r.iterations) break;
            const double el = ev.value("elapsed").toDouble();
            sleepMs((el - lastElapsed) * 1000 / speed);
            lastElapsed = el;
            onIteration(ev, iter > 5 || (interval > 0 && iter % interval == 0) || iter == recordedIterations);
            while (r.liveView > 0 && nextSnap < snaps.size() && iter >= snaps[nextSnap]) snapshot(snaps[nextSnap++], 1);
            if (converged) break;
            if (g_stop) { stopped = true; break; }
        }
        const QString run = recorded(r.stem, "run-" + variant);
        for (const char* f : {"forces.csv", "fluxes.csv"}) {
            QByteArray csv = readFixture(run + "/" + f);
            // Recorded outputs keep their own headers; add the replay notice.
            csv.prepend(QStringLiteral("# FIXTURE REPLAY: recorded final values from a %1-iteration Kutta run\n%2").arg(recordedIterations)
                            .arg(r.overrides.isEmpty() ? QString() : "# MODIFIED RESULT: case settings overridden; recorded values did not change\n").toUtf8());
            writeFile(QDir(r.outDir).filePath(f), csv);
        }
        // The recorded wall loads, so kutta aero works on a replayed run.
        QJsonObject loads = QJsonDocument::fromJson(readFixture(run + "/loads.json")).object();
        if (!loads.isEmpty()) {
            loads["provenance"] = (QStringLiteral("FIXTURE REPLAY: recorded wall loads of a %1-iteration Kutta run; ").arg(recordedIterations) + loads.value("provenance").toString()).trimmed();
            loads["modified"] = caseModified(r.overrides);
            // The run is the replay: it stopped at lastIter, as summary.json and
            // the done event say. The provenance names the recording's length.
            loads["iterations"] = lastIter;
            writeJson(QDir(r.outDir).filePath("loads.json"), loads);
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
        // Iterations 1-5 (residuals only, as the CLI reports them), then each monitor sample.
        std::vector<int> iters;
        for (int it = 1; it <= std::min(5, r.iterations); ++it)
            if (it % r.monitor != 0) iters.push_back(it);
        for (int it = r.monitor; it <= r.iterations; it += r.monitor) iters.push_back(it);
        for (const int it : iters) {
            const bool sample = it % r.monitor == 0;
            const double t = it;
            const double osc = 1 + 0.22 * std::sin(t / 13.0) * std::exp(-t / 600);
            auto res = [&](double a0, double tau, double floor) { return (a0 * std::exp(-t / tau) + floor * (1 + 0.15 * std::sin(t / 7.0))) * osc; };
            QJsonObject resid{{"continuity", res(0.9, 70, 6e-5)}, {"x-velocity", res(0.08, 90, 4e-6)}, {"y-velocity", res(0.05, 95, 3e-6)},
                              {"z-velocity", res(0.06, 85, 3.5e-6)}, {"k", res(0.3, 110, 2e-5)}, {"omega", res(0.2, 120, 1.5e-5)}};
            const double settle = 1 - std::exp(-t / 110) * std::cos(t / 45);
            QJsonObject reports;
            if (sample) for (const auto& [name, f] : finalForce) {
                const double proj = name.contains("downforce") ? -double(f.z()) : double(f.x());
                const double val = proj * settle * (1 + 0.004 * std::sin(t / 9.0) * std::exp(-t / 400));
                reports[name] = QJsonObject{{"force", val}, {"coefficient", val / qa}};
            }
            sleepMs(stepMs * (sample ? r.monitor : 1) - pace.nsecsElapsed() / 1e6);
            pace.restart();
            const double el = clock.nsecsElapsed() / 1e9;
            onIteration({{"v", 1}, {"event", "iteration"}, {"iter", it}, {"residuals", resid}, {"reports", reports},
                         {"sec_per_iter", el / it}, {"elapsed", el}}, sample);
            if (r.liveView > 0 && it % r.liveView == 0) snapshot(it, 1 - std::exp(-t / 150));
            if (converged) break;
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
        // loads.json: the same wall forces with their moments about the origin.
        QJsonArray zl;
        auto arr = [](const QVector3D& v) { return QJsonArray{double(v.x()), double(v.y()), double(v.z())}; };
        for (const ZoneForce& f : forces) {
            int id = -1;
            for (const auto& cz : carZones()) if (cz.name == f.zone) id = cz.id;
            const float k = float(settle);
            zl.append(QJsonObject{{"name", f.zone}, {"zone_id", id}, {"area", f.area}, {"pressure", arr(f.pressure * k)}, {"viscous", arr(f.viscous * k)},
                                  {"pressure_moment", arr(f.pressureMoment * k)}, {"viscous_moment", arr(f.viscousMoment * k)}});
        }
        writeJson(QDir(r.outDir).filePath("loads.json"),
                  {{"v", 1}, {"synthetic", true}, {"origin", QJsonArray{0, 0, 0}},
                   {"units", QJsonObject{{"force", "N"}, {"moment", "N m"}, {"area", "m^2"}}},
                   {"reference", QJsonObject{{"area", r.car.refArea}, {"density", r.car.refDensity}, {"velocity", r.car.refVelocity},
                                             {"length", r.car.refLength}, {"pressure", r.car.refPressure}}},
                   {"report_zones", QJsonArray{"monocoque", "sidepods", "front-wing", "rear-wing", "tyre-front", "tyre-rear"}},
                   {"zones", zl}, {"provenance", prov}, {"approximate", !r.approximations.isEmpty()},
                   {"modified", caseModified(r.overrides)}, {"iterations", lastIter}});
        writeFile(QDir(r.outDir).filePath("fluxes.csv"),
                  (csvNotice(r) + "name,zone_id,mass_flow_outward\n\"inlet\",2," + QString::number(-r.car.density * r.car.speed * 32, 'g', 10) +
                   "\n\"outlet\",3," + QString::number(r.car.density * r.car.speed * 32, 'g', 10) + "\n").toUtf8());
    }
    QString res = csvNotice(r) + "iteration," + residualNames.join(',') + ",seconds_per_iteration\n" + residualRows.join('\n') + "\n";
    writeFile(QDir(r.outDir).filePath("residuals.csv"), res.toUtf8());
    writeJson(QDir(r.outDir).filePath("fixture-run.json"), runDescriptor(r));

    const double elapsed = clock.nsecsElapsed() / 1e9;
    const QString stopReason = converged ? "converged" : stopped ? "stopped" : "iteration_limit";
    QJsonObject convergenceOut = conv.toJson();
    convergenceOut["met_by"] = QJsonArray::fromStringList(converged ? convergence.metBy() : QStringList{});
    convergenceOut["iteration"] = converged ? QJsonValue(lastIter) : QJsonValue();
    QJsonObject sum;
    if (r.kind == CaseKind::Recorded) sum = recordedSummary(r);
    sum["v"] = 1;
    sum["fixture"] = true;
    sum["synthetic"] = r.kind == CaseKind::Car;
    sum["modified"] = caseModified(r.overrides);
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
    sum["converged"] = converged;
    sum["convergence"] = convergenceOut;
    sum["convergence_criteria_available"] = true;
    sum["pressure_velocity_coupling"] = couplingRun(inspected, r.overrides);
    sum["replay_seconds"] = elapsed;
    if (r.kind == CaseKind::Car) {
        sum["cells"] = 18412903;
        sum["wall_time_seconds"] = elapsed;
        sum["kutta_version"] = "fixture";
    }
    writeJson(QDir(r.outDir).filePath("summary.json"), sum);
    emitJson({{"v", 1}, {"event", "done"}, {"converged", converged}, {"iterations", lastIter}, {"elapsed", elapsed},
              {"out", r.outDir}, {"stop_reason", stopReason}, {"convergence", convergenceOut}});
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

QVector3D vec(const QJsonValue& v)
{
    const QJsonArray a = v.toArray();
    if (a.size() != 3) badRequest("expected an [x, y, z] array");
    return QVector3D(float(a[0].toDouble()), float(a[1].toDouble()), float(a[2].toDouble()));
}

// The --post batch, checked like the CLI: v:1, known keys, unique names.
QJsonObject readBatch(const QString& path)
{
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly) || f.size() > (1 << 20)) badRequest("cannot read post specification (at most 1 MiB)");
    QJsonParseError pe;
    const QJsonDocument doc = QJsonDocument::fromJson(f.readAll(), &pe);
    if (!doc.isObject()) badRequest("post specification: " + pe.errorString());
    const QJsonObject b = doc.object();
    if (b.value("v").toInt() != 1) badRequest("post specification must have v: 1");
    static const QStringList keys{"v", "fields", "fields_policy", "streamlines", "isosurfaces", "probes", "integrals", "reference"};
    QStringList names;
    for (const QString& k : b.keys()) {
        if (!keys.contains(k)) badRequest("unknown post specification key " + k);
        if (k == "fields_policy" && b.value(k).toString() != "strict" && b.value(k).toString() != "available")
            badRequest("fields_policy must be \"strict\" or \"available\"");
        if (k == "v" || k == "fields" || k == "fields_policy" || k == "reference") continue;
        for (const QJsonValue& op : b.value(k).toArray()) {
            const QString n = op.toObject().value("name").toString();
            if (n.isEmpty() || names.contains(n)) badRequest("operation names must be present and unique: '" + n + "'");
            names << n;
        }
    }
    return b;
}

QStringList strings(const QJsonValue& v)
{
    QStringList out;
    for (const QJsonValue& x : v.toArray()) out << x.toString();
    return out;
}

// A polyline part beside the manifest; descriptor extras go in `meta`.
QJsonObject writePolylines(const QString& dir, const QString& file, const QString& name, kapp::PartKind kind,
                           const kapp::MeshData& m, QJsonObject meta)
{
    kapp::ViewPart part;
    part.kind = kind;
    part.name = name;
    part.file = file;
    part.meta = std::move(meta);
    writeFile(QDir(dir).filePath(file), kapp::writeMeshBytes(m, &part));
    return kapp::toJson(part);
}

// Probe CSV in the CLI's column layout.
void writeProbeCsv(const QString& path, const kapp::MeshData& m)
{
    QString csv = "sample,x,y,z";
    for (const auto& f : m.fields) csv += ",\"" + f.desc.name + "\"";
    csv += "\n";
    for (std::size_t v = 0; v < m.vertexCount(); ++v) {
        csv += QStringLiteral("%1,%2,%3,%4").arg(v).arg(m.positions[3 * v]).arg(m.positions[3 * v + 1]).arg(m.positions[3 * v + 2]);
        for (const auto& f : m.fields) {
            const float x = f.values[v * std::size_t(f.desc.components)];
            csv += "," + (std::isfinite(x) ? QString::number(x, 'g', 9) : QString());
        }
        csv += "\n";
    }
    writeFile(path, csv.toUtf8());
}

kapp::MeshData lineProbe(const QJsonObject& op, const CarParams& prm)
{
    QStringList fields = strings(op.value("fields"));
    if (fields.isEmpty()) fields << "velocity_magnitude";
    const ProbeResult pr = carLineProbe(vec(op.value("start")), vec(op.value("end")), op.value("count").toInt(100), fields, prm);
    kapp::MeshData m;
    m.polylines = true;
    for (const QVector3D& p : pr.points) m.positions.insert(m.positions.end(), {p.x(), p.y(), p.z()});
    m.indices = {0, quint32(pr.points.size())};
    std::vector<float> valid(pr.points.size(), 1.0f);
    for (const auto& s : pr.series) {
        kapp::MeshData::Field f{{s.name, s.unit, 1, -1, {}, {}}, {}};
        for (std::size_t k = 0; k < s.y.size(); ++k) {
            f.values.push_back(float(s.y[k]));
            if (!std::isfinite(s.y[k])) valid[k] = 0;
        }
        m.fields.push_back(std::move(f));
    }
    m.fields.push_back({{"arc_length", "m", 1, -1, {}, {}}, std::vector<float>(pr.distance.begin(), pr.distance.end())});
    m.fields.push_back({{"valid", "1", 1, -1, {}, {}}, valid});
    return m;
}

// Section through exported surface zones: one curve per connected piece.
kapp::MeshData sectionProbe(const QJsonObject& op, const kapp::ViewManifest& view, QJsonArray* lines)
{
    const QVector3D origin = vec(op.value("origin")), normal = vec(op.value("normal"));
    if (normal.isNull()) badRequest("section normal must be nonzero");
    const QStringList fields = strings(op.value("fields"));
    kapp::MeshData out;
    out.polylines = true;
    out.indices = {0};
    std::vector<std::vector<float>> values(fields.size());
    std::vector<float> arc;
    const float tol = std::max(1e-9f, (view.bounds.max - view.bounds.min).length() * 1e-6f);
    for (const QString& zone : strings(op.value("zones"))) {
        const auto zp = std::find_if(view.zones.begin(), view.zones.end(), [&](const kapp::ViewPart& z) { return z.name == zone; });
        if (zp == view.zones.end()) fail("unknown section zone " + zone, {}, false);
        auto mesh = kapp::readMeshFile(QDir(view.dir).filePath(zp->file), *zp);
        if (!mesh) fail(mesh.error(), {}, false);
        // Curves are traced on the first field; every field samples the same segments.
        std::vector<std::vector<std::vector<kapp::SectionPoint>>> perField;
        for (const QString& f : fields) {
            if (!zp->field(f)) fail("unknown or unavailable post-processing field: " + f, {}, false);
            perField.push_back(kapp::sectionCurves(kapp::planeSection(mesh.value(), f, -1, origin, normal.normalized()), tol));
        }
        if (perField.empty()) continue;
        for (std::size_t c = 0; c < perField[0].size(); ++c) {
            const auto& curve = perField[0][c];
            double s = 0;
            for (std::size_t k = 0; k < curve.size(); ++k) {
                if (k) s += (curve[k].position - curve[k - 1].position).length();
                const QVector3D p = curve[k].position;
                out.positions.insert(out.positions.end(), {p.x(), p.y(), p.z()});
                arc.push_back(float(s));
                for (std::size_t f = 0; f < fields.size(); ++f) values[f].push_back(float(perField[f][c][k].value));
            }
            out.indices.push_back(quint32(out.vertexCount()));
            const bool closed = curve.size() > 2 && (curve.front().position - curve.back().position).length() <= tol;
            lines->append(QJsonObject{{"zone", zone}, {"closed", closed}, {"samples", int(curve.size())}, {"length", s}});
        }
    }
    for (std::size_t f = 0; f < fields.size(); ++f) {
        const auto* fd = view.field(fields[int(f)]);
        out.fields.push_back({{fields[int(f)], fd ? fd->unit : QString(), 1, -1, {}, {}}, std::move(values[f])});
    }
    out.fields.push_back({{"arc_length", "m", 1, -1, {}, {}}, arc});
    out.fields.push_back({{"valid", "1", 1, -1, {}, {}}, std::vector<float>(arc.size(), 1.0f)});
    return out;
}

// Area-weighted surface report over exported zones. The fixture has no face
// mass flux, so mass flows and mass-weighted averages are null, as the CLI
// reports undefined quantities.
QJsonObject surfaceReport(const QJsonObject& op, const kapp::ViewManifest& view)
{
    QStringList zones = strings(op.value("zones"));
    if (zones.isEmpty())
        for (const auto& z : view.zones) zones << z.name;
    const QStringList fields = strings(op.value("fields"));
    auto entry = [&](const QString& name, int id, double area, const std::map<QString, std::pair<double, double>>& sums,
                     const std::map<QString, QString>& units) {
        QJsonObject averages;
        for (const QString& f : fields) {
            auto it = sums.find(f);
            const bool have = it != sums.end() && it->second.second > 0;
            averages[f] = QJsonObject{{"area_integral", have ? QJsonValue(it->second.first) : QJsonValue()},
                                      {"area_weighted_average", have ? QJsonValue(it->second.first / it->second.second) : QJsonValue()},
                                      {"area_coverage", area > 0 && have ? it->second.second / area : 0},
                                      {"mass_flow_weighted_average", QJsonValue()}, {"mass_coverage", 0},
                                      {"unit", units.count(f) ? units.at(f) : QString()}};
        }
        return QJsonObject{{"zone", name}, {"zone_id", id}, {"area", area}, {"mass_flow", QJsonValue()}, {"absolute_mass_flow", QJsonValue()},
                           {"mass_flow_source", "unavailable"}, {"averages", averages}};
    };
    QJsonArray rows;
    double totalArea = 0;
    std::map<QString, std::pair<double, double>> totals;
    std::map<QString, QString> units;
    for (const QString& zone : zones) {
        const auto zp = std::find_if(view.zones.begin(), view.zones.end(), [&](const kapp::ViewPart& z) { return z.name == zone; });
        if (zp == view.zones.end()) fail("unknown report zone " + zone, {}, false);
        auto mesh = kapp::readMeshFile(QDir(view.dir).filePath(zp->file), *zp);
        if (!mesh) fail(mesh.error(), {}, false);
        const double area = kapp::integrate(mesh.value(), {}, kapp::IntegralKind::Area).value_or(kapp::SurfaceIntegral{}).value;
        std::map<QString, std::pair<double, double>> sums;
        for (const QString& f : fields) {
            const auto* fd = zp->field(f);
            if (!fd) continue;  // boundary field absent on this zone: zero coverage
            units[f] = fd->unit;
            if (auto r = kapp::integrate(mesh.value(), f, kapp::IntegralKind::Integral)) {
                sums[f] = {r->value, r->area};
                totals[f].first += r->value;
                totals[f].second += r->area;
            }
        }
        totalArea += area;
        rows.append(entry(zone, zp->id, area, sums, units));
    }
    return {{"name", op.value("name").toString()}, {"zones", rows}, {"total", entry("selected-zones", -1, totalArea, totals, units)},
            {"units", QJsonObject{{"area", "m^2"}, {"mass_flow", "kg/s"}, {"force", "N"}, {"moment", "N m"}, {"centre_of_pressure", "m"}}}};
}

int exportView(const Args& a)
{
    if (a.positional.size() != 1) usage(false);
    const QString dir = a.positional[0];
    const QJsonObject desc = readJson(QDir(dir).filePath("fixture-run.json"));
    if (desc.isEmpty()) fail(dir + " is not a finished fixture run (no summary)", {}, false);
    const auto slices = parseSlices(a.opts.value("slice"));
    const QJsonObject batch = a.has("post") ? readBatch(a.value("post")) : QJsonObject{};
    const bool car = !desc.value("recorded").toBool(desc.value("case").toString() == "mixing_elbow");
    const QString stem = desc.value("case").toString();
    const CarParams prm = paramsFrom(desc);
    // Build beside the old view and swap at the end: a failed export keeps it.
    const QString view = QDir(dir).filePath("view"), staging = QDir(dir).filePath("view.staging");
    QDir(staging).removeRecursively();
    if (!car) {
        copyFixtureDir(recorded(stem, "run-" + desc.value("variant").toString() + "/view"), staging);
        if (!slices.empty())
            eprint("kutta-fixture: a recorded view keeps its recorded slices; requested slices need the real CLI");
        QJsonObject m = readJson(QDir(staging).filePath("manifest.json"));
        if (!desc.value("provenance").toString().isEmpty()) m["provenance"] = desc.value("provenance").toString();
        writeJson(QDir(staging).filePath("manifest.json"), m);
    } else {
        writeCarView(staging, prm, desc.value("reference").toBool(), desc.value("provenance").toString(), slices, false);
    }
    if (!batch.isEmpty()) {
        auto base = kapp::loadViewManifest(QDir(staging).filePath("manifest.json"));
        if (!base) fail(base.error(), {}, false);
        QJsonObject manifest = readJson(QDir(staging).filePath("manifest.json"));
        auto needsCar = [&](const QString& what) {
            if (!car) fail(what + " on a recorded case needs the real kutta CLI (the fixture has no volume fields)", {}, false);
        };
        QJsonArray traces, isos, probes, reports, available;
        int n = 0;
        for (const QJsonValue& v : batch.value("streamlines").toArray()) {
            const QJsonObject op = v.toObject();
            needsCar("streamlines");
            const QJsonObject seeds = op.value("seeds").toObject();
            if (seeds.value("type").toString() != "rake") fail("the fixture traces rake seeds only", {}, false);
            const kapp::MeshData m = carStreamlines(vec(seeds.value("start")), vec(seeds.value("end")), seeds.value("count").toInt(20),
                                                   op.value("max_length").toDouble(8), prm);
            traces.append(writePolylines(staging, QStringLiteral("streamlines-%1.bin").arg(n++), op.value("name").toString(),
                                         kapp::PartKind::Streamlines, m, {{"mode", op.value("mode").toString("streamline")}, {"velocity_time_model", "steady"}}));
        }
        n = 0;
        for (const QJsonValue& v : batch.value("isosurfaces").toArray()) {
            const QJsonObject op = v.toObject();
            needsCar("iso-surfaces");
            const QString field = op.value("field").toString();
            const kapp::MeshData m = carIsoSurface(field, op.value("value").toDouble(), prm);
            kapp::ViewPart part;
            part.kind = kapp::PartKind::IsoSurface;
            part.name = op.value("name").toString();
            part.file = QStringLiteral("iso-%1.bin").arg(n++);
            part.meta = {{"field", field}, {"value", op.value("value").toDouble()}, {"component", -1}};
            writeFile(QDir(staging).filePath(part.file), kapp::writeMeshBytes(m, &part));
            isos.append(kapp::toJson(part));
        }
        n = 0;
        for (const QJsonValue& v : batch.value("probes").toArray()) {
            const QJsonObject op = v.toObject();
            const QString type = op.value("type").toString();
            const QString file = QStringLiteral("probe-%1.bin").arg(n), csv = QStringLiteral("probe-%1.csv").arg(n);
            ++n;
            QJsonArray lines;
            kapp::MeshData m;
            if (type == "line") {
                needsCar("line probes");
                m = lineProbe(op, prm);
                lines.append(QJsonObject{{"closed", false}, {"samples", int(m.vertexCount())}});
            } else if (type == "section") {
                m = sectionProbe(op, base.value(), &lines);
            } else {
                badRequest("probe type must be line or section");
            }
            writeProbeCsv(QDir(staging).filePath(csv), m);
            probes.append(writePolylines(staging, file, op.value("name").toString(), kapp::PartKind::Probe, m,
                                         {{"type", type}, {"csv", csv}, {"lines", lines}}));
        }
        for (const QJsonValue& v : batch.value("integrals").toArray()) reports.append(surfaceReport(v.toObject(), base.value()));
        // Every exported field is sampleable here; the car has a synthetic volume.
        for (const auto& f : base->fields) {
            const bool wall = f.name.startsWith("wall_shear") || f.name.startsWith("y_plus") || f.name.startsWith("skin_friction");
            available.append(QJsonObject{{"name", f.name}, {"unit", f.unit}, {"components", f.components},
                                         {"location", car && !wall ? "volume" : "boundary"}});
        }
        manifest["streamlines"] = traces;
        manifest["isosurfaces"] = isos;
        manifest["probes"] = probes;
        manifest["reports"] = reports;
        manifest["available_fields"] = available;
        if (!reports.isEmpty()) {
            manifest["reports_file"] = "reports.json";
            writeJson(QDir(staging).filePath("reports.json"), {{"v", 1}, {"reports", reports}});
        }
        writeJson(QDir(staging).filePath("manifest.json"), manifest);
    }
    QDir(view).removeRecursively();
    if (!QDir().rename(staging, view)) fail("cannot publish " + view, {}, false);
    eprint("Wrote " + view);
    return 0;
}

// aero: coefficients, breakdown and balance from a run's loads.json.
int aero(const Args& a)
{
    const bool jsonl = a.has("json");
    if (a.positional.size() != 1) usage(jsonl);
    const QString dir = a.positional[0];
    QJsonObject spec;
    if (a.has("aero")) {
        QFile f(a.value("aero"));
        QJsonParseError pe;
        if (!f.open(QIODevice::ReadOnly) || f.size() > (1 << 20)) {
            if (jsonl) emitJson({{"v", 1}, {"event", "error"}, {"message", "cannot read aero spec (at most 1 MiB)"}, {"issues", QJsonArray{}}});
            badRequest("cannot read aero spec (at most 1 MiB)");
        }
        const QJsonDocument doc = QJsonDocument::fromJson(f.readAll(), &pe);
        if (!doc.isObject()) {
            if (jsonl) emitJson({{"v", 1}, {"event", "error"}, {"message", "aero spec: " + pe.errorString()}, {"issues", QJsonArray{}}});
            badRequest("aero spec: " + pe.errorString());
        }
        spec = doc.object();
    }
    const QJsonObject loads = readJson(QDir(dir).filePath("loads.json"));
    if (loads.isEmpty())
        fail(QStringLiteral("%1 has no loads.json: the run was solved before aero loads were written, or did not finish; solve it again.").arg(dir), {}, jsonl);
    QString error;
    QJsonObject out = aeroFromLoads(loads, spec, &error);
    if (!error.isEmpty()) {
        if (jsonl) emitJson({{"v", 1}, {"event", "error"}, {"message", error}, {"issues", QJsonArray{}}});
        badRequest(error);
    }
    writeJson(QDir(dir).filePath("aero.json.tmp"), out);
    QFile::remove(QDir(dir).filePath("aero.json"));
    QFile::rename(QDir(dir).filePath("aero.json.tmp"), QDir(dir).filePath("aero.json"));
    if (jsonl) {
        out["event"] = "aero";
        emitJson(out);
    } else eprint("Wrote " + QDir(dir).filePath("aero.json"));
    return 0;
}

// export-geometry: the case's boundary zones before any run. A recorded case
// reuses its recorded view's zones; the car is generated. Fields in zone files are
// ignored by readers of a geometry-only manifest.
int exportGeometry(const Args& a)
{
    if (a.positional.size() != 1 || !a.has("out")) usage(false);
    const QString casePath = a.positional[0];
    const bool jsonl = a.has("json");
    const bool car = caseKind(casePath, jsonl) == CaseKind::Car;
    const QString root = QDir(a.value("out")).filePath("geometry"), staging = QDir(a.value("out")).filePath("geometry.staging");
    QDir(staging).removeRecursively();
    if (car) writeCarView(staging, paramsFrom({}), false, {}, {}, false);
    else copyFixtureDir(recorded(caseStem(casePath), "run-cold/view"), staging);
    for (const QString& s : QDir(staging).entryList({"slice-*"})) QFile::remove(QDir(staging).filePath(s));
    QJsonObject m = readJson(QDir(staging).filePath("manifest.json"));
    m["geometry_only"] = true;
    m["reference"] = false;
    m["provenance"] = QString();
    m["slices"] = QJsonArray{};
    m["fields"] = QJsonArray{};
    writeJson(QDir(staging).filePath("manifest.json"), m);
    QDir(root).removeRecursively();
    if (!QDir().rename(staging, root)) fail("cannot publish " + root, {}, jsonl);
    if (jsonl)
        emitJson({{"v", 1}, {"event", "geometry"}, {"manifest", QFileInfo(QDir(root).filePath("manifest.json")).absoluteFilePath()},
                  {"zones", int(m.value("zones").toArray().size())}});
    else eprint("Wrote " + root);
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
    if (cmd == "export-geometry") return exportGeometry(a);
    if (cmd == "aero") return aero(a);
    usage(false);
}
