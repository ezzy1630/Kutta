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
//   KUTTA_FIXTURE_MESH_HOLD   `mesh` pauses at this stage until stdin closes
//
//   mesh GEOM.stl...          a SYNTHETIC stand-in for `kutta mesh`: it meshes
//                             nothing. It accepts STL solids named after the
//                             synthetic car's walls (sample-stl writes them),
//                             writes a placeholder mesh.msh the other commands
//                             read as that car (each solid's size and position
//                             carried over), and a synthetic mesh-report.json.
//   sample-stl DIR            the synthetic car's walls as STL files in mm

#include "Aero.h"
#include "SyntheticCar.h"

#include "core/CliClient.h"
#include "core/MeshSetup.h"
#include "core/Stl.h"
#include "core/Geometry.h"
#include "core/Protocol.h"
#include "core/SettingsModel.h"
#include "core/ViewFormat.h"

#include <QCoreApplication>
#include <QCryptographicHash>
#include <QDateTime>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QRegularExpression>
#include <QThread>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <map>
#include <thread>

namespace {

using namespace fixture;
using kapp::LengthUnit;

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
                        "| export-geometry CASE --out DIR [--json] | aero DIR [--aero SPEC] [--json] "
                        "| new-case MESH [--speed U] [--role ZONE=ROLE]... [--out CASE.kutta.json [--force]] [--json] "
                        "| convert CASE.kutta.json [--setup FILE] --out CASE.kutta.json [--force] [--json]";
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
        else if (s == "--force") a.opts["force"] << "1";
        else if (s == "--case") a.opts["case"] << "1";
        else if (s.startsWith("--")) {
            if (i + 1 >= args.size()) usage(false);
            a.opts[s.mid(2)] << args[++i];
        } else a.positional << s;
    }
    return a;
}

enum class CaseKind { Recorded, Car };

bool isKuttaCase(const QString& path) { return path.endsWith(".json", Qt::CaseInsensitive); }

// A Kutta case's mesh, absolute (relative paths are relative to the case file).
QString kuttaMesh(const QString& casePath)
{
    QFile f(casePath);
    if (!f.open(QIODevice::ReadOnly)) return {};
    const QString mesh = QJsonDocument::fromJson(f.readAll()).object().value("mesh").toString();
    if (mesh.isEmpty()) return {};
    return QDir::cleanPath(QFileInfo(mesh).isAbsolute() ? mesh : QDir(QFileInfo(casePath).absolutePath()).absoluteFilePath(mesh));
}

// The settings the fixture carries in a Kutta case, at the JSON path where the
// real CLI writes them (docs/kutta-case.md); empty for any other setting.
// Saving writes these, and inspecting a Kutta case reads them back over the
// recording's inventory, so Save Case round-trips them without a solver.
QStringList kuttaCasePath(const QString& setting, const QJsonObject& kase)
{
    static const QStringList reference{"area", "length", "density", "velocity", "pressure", "viscosity", "depth"};
    const QStringList p = setting.split('/');
    if (p.size() == 2 && p[0] == "reference" && reference.contains(p[1])) return {"reference_values", p[1]};
    if (setting == "general/operating-pressure") return {"operating_pressure"};
    if (p.size() == 3 && p[0] == "boundary" &&
        kase.value("boundaries").toObject().value(p[1]).toObject().value("type") == "velocity-inlet") {
        if (p[2] == "velocity-magnitude") return {"boundaries", p[1], "speed"};
        if (p[2] == "turb-intensity") return {"boundaries", p[1], "turbulence", "intensity"};
        if (p[2] == "turb-viscosity-ratio") return {"boundaries", p[1], "turbulence", "viscosity_ratio"};
    }
    return {};
}

// Fluent's defaults for what a Kutta case leaves out (docs/kutta-case.md, Defaults).
QJsonValue kuttaCaseDefault(const QStringList& path)
{
    static const QMap<QString, double> d{{"reference_values/area", 1},     {"reference_values/length", 1},
                                         {"reference_values/density", 1.225}, {"reference_values/velocity", 1},
                                         {"reference_values/pressure", 0},  {"reference_values/viscosity", 1.7894e-05},
                                         {"reference_values/depth", 1},     {"operating_pressure", 101325}};
    const auto it = d.find(path.join('/'));
    return it == d.end() ? QJsonValue() : QJsonValue(*it);
}

QJsonValue jsonAt(const QJsonObject& o, const QStringList& path)
{
    QJsonValue v = o;
    for (const QString& k : path) v = v.toObject().value(k);
    return v;
}

void setJsonAt(QJsonObject& o, QStringList path, const QJsonValue& value)
{
    const QString key = path.takeFirst();
    if (path.isEmpty()) {
        o[key] = value;
        return;
    }
    QJsonObject child = o.value(key).toObject();
    setJsonAt(child, path, value);
    o[key] = child;
}

// ---- fixture meshes (`mesh`) ---------------------------------------------------------

const QByteArray kFixtureMeshMarker = "Kutta fixture synthetic mesh";

// The placeholder mesh.msh `mesh` writes: a marker line, then JSON.
QJsonObject fixtureMeshInfo(const QString& meshPath)
{
    QFile f(meshPath);
    if (!f.open(QIODevice::ReadOnly)) return {};
    const QByteArray head = f.readLine();
    if (!head.startsWith(kFixtureMeshMarker)) return {};
    return QJsonDocument::fromJson(f.readAll()).object();
}

// A case or mesh path's fixture mesh, if it is one.
QJsonObject fixtureMeshOf(const QString& path)
{
    if (isKuttaCase(path)) return fixtureMeshInfo(kuttaMesh(path));
    return fixtureMeshInfo(path);
}

QMap<QString, SolidTransform> meshTransforms(const QString& path)
{
    QMap<QString, SolidTransform> out;
    const QJsonObject t = fixtureMeshOf(path).value("transforms").toObject();
    for (auto it = t.begin(); it != t.end(); ++it) {
        const QJsonObject o = it.value().toObject();
        const QJsonArray p = o.value("pivot").toArray(), d = o.value("offset").toArray();
        SolidTransform x;
        x.pivot = QVector3D(float(p[0].toDouble()), float(p[1].toDouble()), float(p[2].toDouble()));
        x.offset = QVector3D(float(d[0].toDouble()), float(d[1].toDouble()), float(d[2].toDouble()));
        x.scale = float(o.value("scale").toDouble(1));
        out[it.key()] = x;
    }
    return out;
}

QJsonObject transformsJson(const QMap<QString, SolidTransform>& t)
{
    QJsonObject o;
    for (auto it = t.begin(); it != t.end(); ++it)
        o[it.key()] = QJsonObject{{"pivot", QJsonArray{it->pivot.x(), it->pivot.y(), it->pivot.z()}},
                                  {"offset", QJsonArray{it->offset.x(), it->offset.y(), it->offset.z()}},
                                  {"scale", it->scale}};
    return o;
}

// "elbow_param.cas.h5" -> "elbow_param": the recording's directory name. A
// Kutta case replays its mesh's recording.
QString caseStem(const QString& casePath)
{
    if (isKuttaCase(casePath))
        if (const QString mesh = kuttaMesh(casePath); !mesh.isEmpty()) return caseStem(mesh);
    QString name = QFileInfo(casePath).fileName();
    for (const char* ext : {".h5", ".gz", ".cas", ".dat", ".msh"})
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
    if (!fixtureMeshOf(casePath).isEmpty()) return CaseKind::Car;
    if (recordedCases().contains(caseStem(casePath))) return CaseKind::Recorded;
    if (name.startsWith("fsae_car") || caseStem(casePath).startsWith("fsae_car")) return CaseKind::Car;
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
        if (const QJsonObject fm = fixtureMeshOf(casePath); fm.contains("cells")) {
            // A fixture mesh states its own (synthetic) size.
            QJsonObject mesh = o.value("mesh").toObject();
            for (const char* k : {"cells", "faces", "nodes", "cell_types"}) mesh[k] = fm.value(k);
            o["mesh"] = mesh;
        }
    }
    QJsonObject file = o.value("file").toObject();
    file["path"] = QFileInfo(casePath).absoluteFilePath();
    if (isKuttaCase(casePath)) {
        // The recording's settings, with the case file's identity and boundary types.
        QFile f(casePath);
        const QByteArray text = f.open(QIODevice::ReadOnly) ? f.readAll() : QByteArray();
        const QJsonObject kutta = QJsonDocument::fromJson(text).object();
        file["format"] = "kutta";
        file["fluent_version"] = QJsonValue::Null;
        file["setup_sha256"] = QString::fromLatin1(QCryptographicHash::hash(text, QCryptographicHash::Sha256).toHex());
        file["sha256"] = file["setup_sha256"];
        file["mesh_file"] = QJsonObject{{"path", kuttaMesh(casePath)}, {"sha256", QString()}};
        const QJsonObject boundaries = kutta.value("boundaries").toObject();
        QJsonArray zones;
        for (const QJsonValue& zv : o.value("zones").toArray()) {
            QJsonObject z = zv.toObject();
            const auto boundary = boundaries.value(z.value("name").toString()).toObject();
            const QString type = boundary.value("type").toString();
            if (!type.isEmpty()) z["type"] = type;
            if (type == "wall") {
                auto condition = z.value("condition").toObject();
                const QString motion = boundary.value("motion").toString("stationary");
                condition["wall_motion"] = motion == "moving" ? "translating" : motion;
                for (const auto* key : {"angular_speed", "axis_origin", "axis_direction"}) condition[key] = boundary.value(key);
                z["condition"] = condition;
            }
            zones.append(z);
        }
        o["zones"] = zones;
        QJsonArray settings;
        for (const QJsonValue& sv : o.value("settings").toArray()) {
            QJsonObject st = sv.toObject();
            const auto path = st.value("path").toString().split('/');
            if (path.size() == 3 && path[0] == "boundary") {
                const auto b = boundaries.value(path[1]).toObject();
                if (path[2] == "wall-motion") st["value"] = b.value("motion").toString("stationary");
                if (path[2] == "rotation-speed") st["value"] = b.value("angular_speed");
            }
            if (const QStringList at = kuttaCasePath(st.value("path").toString(), kutta); !at.isEmpty()) {
                const QJsonValue v = jsonAt(kutta, at);
                if (!v.isUndefined()) st["value"] = v;
                else if (const QJsonValue d = kuttaCaseDefault(at); !d.isUndefined()) st["value"] = d;
            }
            settings.append(st);
        }
        o["settings"] = settings;
    }
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
                                   {"ref_length", r.car.refLength}, {"transforms", transformsJson(r.car.transforms)}}}};
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
    const QJsonObject t = p.value("transforms").toObject();
    for (auto it = t.begin(); it != t.end(); ++it) {
        const QJsonObject o = it.value().toObject();
        const QJsonArray pv = o.value("pivot").toArray(), d = o.value("offset").toArray();
        SolidTransform x;
        x.pivot = QVector3D(float(pv[0].toDouble()), float(pv[1].toDouble()), float(pv[2].toDouble()));
        x.offset = QVector3D(float(d[0].toDouble()), float(d[1].toDouble()), float(d[2].toDouble()));
        x.scale = float(o.value("scale").toDouble(1));
        c.transforms[it.key()] = x;
    }
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

    if (r.kind == CaseKind::Car) {
        applyCarOverrides(r);
        r.car.transforms = meshTransforms(r.casePath);
    }
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
    if (car) {
        CarParams prm = paramsFrom({});
        prm.transforms = meshTransforms(casePath);
        writeCarView(staging, prm, false, {}, {}, false);
    }
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

// kutta new-case: the recorded reply for a mesh, with the requested speed and
// roles applied the way the CLI's FSAE template applies them.
QJsonObject starterBoundary(const QString& role, double speed)
{
    const QJsonArray x{1, 0, 0};
    const QJsonObject turbulence{{"intensity", 0.01}, {"viscosity_ratio", 10}};
    if (role == "velocity-inlet") return {{"type", role}, {"speed", speed}, {"direction", x}, {"turbulence", turbulence}};
    if (role == "pressure-outlet") return {{"type", role}, {"gauge_pressure", 0}, {"turbulence", turbulence}};
    if (role == "moving-ground") return {{"type", "wall"}, {"motion", "moving"}, {"speed", speed}, {"direction", x}};
    return {{"type", role}};
}

// The synthetic car's mesh (fsae_car*.msh) as kutta new-case would read it: its
// zones, roles suggested by the CLI's rules (src/kutta/setup/kutta_case.cpp
// suggest_role), groups by the app's aero grouping, and the FSAE template. It
// gives the fixture a 3D starter with wheels, and so the CLI's wheel note.
QJsonObject carNewCaseReply()
{
    const QJsonObject inspect = carInspect(paramsFrom({}), false, {});
    QStringList names;
    for (const QJsonValue& zv : inspect.value("zones").toArray()) names << zv.toObject().value("name").toString();
    QMap<QString, QString> groupOf;
    for (const auto& g : kapp::suggestAeroGroups(names))
        for (const QString& z : g.zones) groupOf[z] = g.name;
    QJsonArray zones;
    qint64 faces = 0;
    for (const QJsonValue& zv : inspect.value("zones").toArray()) {
        const QJsonObject z = zv.toObject();
        const QString name = z.value("name").toString(), type = z.value("type").toString();
        const QStringList words = name.toLower().split(QRegularExpression("[^a-z0-9]+"), Qt::SkipEmptyParts);
        QString role = "wall", reason = "a wall of the car";
        if (words.contains("inlet")) role = "velocity-inlet", reason = "name contains \"inlet\"";
        else if (words.contains("outlet")) role = "pressure-outlet", reason = "name contains \"outlet\"";
        else if (words.contains("ground")) role = "moving-ground", reason = "name contains \"ground\"; the road moves with the freestream";
        else if (type == "symmetry") role = "symmetry", reason = "mesh zone type symmetry";
        faces += z.value("faces").toInteger();
        zones.append(QJsonObject{{"id", z.value("id")}, {"name", name}, {"mesh_type", type}, {"kind", "boundary"},
                                 {"faces", z.value("faces")}, {"role", role}, {"suggested_role", role}, {"reason", reason},
                                 {"group", role == "wall" ? QJsonValue(groupOf.value(name, "Body")) : QJsonValue()}});
    }
    const QJsonObject m = inspect.value("mesh").toObject();
    const QJsonArray x{1, 0, 0};
    const QJsonObject kase{
        {"kutta_case", 1},
        {"models", QJsonObject{{"viscous", "k-omega-sst"}}},
        {"material", QJsonObject{{"name", "air"}, {"density", 1.225}, {"viscosity", 1.7894e-05}}},
        {"methods", QJsonObject{{"turbulent_kinetic_energy", "second-order-upwind"}, {"specific_dissipation_rate", "second-order-upwind"}}},
        {"reference_values", QJsonObject{{"area", 1}, {"length", 1}, {"density", 1.225}}},
        {"reports", QJsonObject{{"drag", QJsonObject{{"type", "drag"}, {"direction", x}, {"coefficient", true}}},
                                {"downforce", QJsonObject{{"type", "lift"}, {"direction", QJsonArray{0, 0, -1}}, {"coefficient", true}}}}}};
    return {{"v", 1}, {"template", "fsae-external"}, {"speed", 15},
            {"mesh", QJsonObject{{"dimension", 3}, {"cells", m.value("cells")}, {"faces", m.value("faces").toInteger(faces)},
                                 {"nodes", m.value("nodes")}, {"polyhedral", false}}},
            {"zones", zones},
            {"roles", QJsonArray{"velocity-inlet", "pressure-outlet", "moving-ground", "wall", "symmetry", "interior"}},
            {"notes", QJsonArray{"Reference area is Fluent's default, 1 m^2; set your car's frontal area (reference_values/area) for meaningful coefficients."}},
            {"case", kase}};
}

int newCase(const Args& a)
{
    const bool json = a.has("json");
    if (a.positional.size() != 1) usage(false);
    const QString mesh = a.positional[0];
    QJsonObject o = caseStem(mesh).startsWith("fsae_car") || !fixtureMeshInfo(mesh).isEmpty() ? carNewCaseReply()
                                                          : QJsonDocument::fromJson(readFixture(recorded(caseStem(mesh), "new-case.json"))).object();
    if (o.isEmpty())
        fail(QStringLiteral("The fixture CLI reads naca4412_ground.msh (recorded) and fsae_car.msh (synthetic) only; point the app at the real "
                            "kutta CLI to set up %1.")
                 .arg(QFileInfo(mesh).fileName()), {}, json);
    if (!QFileInfo(mesh).isFile()) fail(mesh + ": no such file", {}, json);
    double speed = o.value("speed").toDouble();
    if (a.has("speed")) {
        bool ok = false;
        speed = a.value("speed").toDouble(&ok);
        if (!ok || !(speed > 0)) fail("speed: --speed must be a positive number of m/s", {}, json);
    }
    QMap<QString, QString> roles;
    const QJsonArray known = o.value("roles").toArray();
    for (const QString& r : a.opts.value("role")) {
        const qsizetype eq = r.lastIndexOf('=');
        if (eq < 0 || !known.contains(r.mid(eq + 1))) fail("role: --role " + r + ": expected ZONE=ROLE", {}, json);
        roles[r.left(eq)] = r.mid(eq + 1);
    }
    QJsonObject kase = o.value("case").toObject(), boundaries;
    QMap<QString, QJsonObject> rotations;
    for (const auto& spec : a.opts.value("wheel")) {
        const auto eq = spec.lastIndexOf('=');
        const auto parts = spec.mid(eq + 1).split(',');
        if (eq <= 0 || parts.size() != 7) fail("wheel: expected ZONE=RADIUS,X,Y,Z,AX,AY,AZ", {}, json);
        std::array<double, 7> values;
        for (int i = 0; i < 7; ++i) {
            bool ok = false;
            values[i] = parts[i].toDouble(&ok);
            if (!ok || !std::isfinite(values[i])) fail("wheel: geometry must be finite", {}, json);
        }
        if (!(values[0] > 0) || (values[4] == 0 && values[5] == 0 && values[6] == 0)) fail("wheel: positive radius and nonzero axis required", {}, json);
        const QString name = spec.left(eq);
        if (rotations.contains(name)) fail("wheel: duplicate wheel", {}, json);
        rotations[name] = {{"type", "wall"}, {"motion", "rotating"}, {"angular_speed", speed / values[0]}, {"relative", false},
                           {"axis_origin", QJsonArray{values[1], values[2], values[3]}}, {"axis_direction", QJsonArray{values[4], values[5], values[6]}}};
    }
    QJsonArray zones, walls;
    bool inlet = false, outlet = false, wheels = false;
    for (const QJsonValue& zv : o.value("zones").toArray()) {
        QJsonObject z = zv.toObject();
        const QString name = z.value("name").toString();
        z["role"] = roles.value(name, z.value("suggested_role").toString());
        const QString role = z.value("role").toString();
        boundaries[name] = starterBoundary(role, speed);
        if (rotations.contains(name)) {
            if (role != "wall") fail("wheel: must name a wall boundary zone", {}, json);
            boundaries[name] = rotations.value(name);
            rotations.remove(name);
        }
        if (role == "wall") walls.append(name);
        inlet = inlet || role == "velocity-inlet";
        outlet = outlet || role == "pressure-outlet";
        wheels = wheels || (role == "wall" && z.value("group") == "Wheels" && boundaries.value(name).toObject().value("motion") != "rotating");
        zones.append(z);
    }
    if (!rotations.isEmpty()) fail("wheel: unknown wall boundary zone", {}, json);
    kase["boundaries"] = boundaries;
    // The CLI's notes that follow the roles, in its words (case_commands.cpp).
    QJsonArray notes;
    for (const QJsonValue& n : o.value("notes").toArray())
        if (!n.toString().startsWith("No zone is a") && !n.toString().startsWith("Wheels are stationary")) notes.append(n);
    if (!inlet) notes.append("No zone is a velocity inlet; assign the inflow boundary its role.");
    if (!outlet) notes.append("No zone is a pressure outlet; assign the outflow boundary its role.");
    if (wheels) notes.append("Wheels are stationary walls in this starter; make them rotating walls (axis, origin, angular speed) if the car's wheels turn.");
    o["notes"] = notes;
    QJsonObject reports = kase.value("reports").toObject();
    for (const QString& name : reports.keys()) {
        QJsonObject r = reports.value(name).toObject();
        r["zones"] = walls;
        reports[name] = r;
    }
    kase["reports"] = reports;
    QJsonObject reference = kase.value("reference_values").toObject();
    reference["velocity"] = speed;
    kase["reference_values"] = reference;
    const QString target = a.value("out");
    kase["mesh"] = target.isEmpty() ? mesh : QDir(QFileInfo(target).absolutePath()).relativeFilePath(QFileInfo(mesh).absoluteFilePath());
    if (!target.isEmpty()) {
        if (QFileInfo::exists(target) && !a.has("force")) fail(target + " exists; pass --force to replace it", {}, json);
        writeJson(target, kase);
    }
    QJsonObject meshInfo = o.value("mesh").toObject();
    meshInfo["path"] = QFileInfo(mesh).absoluteFilePath();
    if (const QJsonObject fm = fixtureMeshInfo(mesh); fm.contains("cells"))
        for (const char* k : {"cells", "faces", "nodes"}) meshInfo[k] = fm.value(k);
    o["mesh"] = meshInfo;
    o["zones"] = zones;
    o["speed"] = speed;
    o["case"] = kase;
    o["out"] = target.isEmpty() ? QJsonValue(QJsonValue::Null) : QJsonValue(QFileInfo(target).absoluteFilePath());
    if (json) out() << QJsonDocument(o).toJson(QJsonDocument::Compact) << '\n';
    else out() << QJsonDocument(kase).toJson(QJsonDocument::Indented);
    return 0;
}

// kutta convert: a Kutta case copied to a new place (its mesh path follows),
// with the setting changes kuttaCasePath knows written in; it refuses other
// changes rather than drop them.
int convert(const Args& a)
{
    const bool json = a.has("json");
    if (a.positional.size() != 1 || !a.has("out")) usage(false);
    const QString input = a.positional[0], target = a.value("out");
    if (!isKuttaCase(input))
        fail("The fixture CLI saves Kutta cases only; point the app at the real kutta CLI to save a Fluent case as a Kutta case.", {}, json);
    if (QFileInfo::exists(target) && !a.has("force")) fail(target + " exists; pass --force to replace it", {}, json);
    QJsonObject kase = readJson(input);
    if (kase.isEmpty()) fail(input + ": not a Kutta case", {}, json);
    const QJsonObject changes = a.has("setup") ? caseOverrides(readJson(a.value("setup")).value("overrides").toObject()) : QJsonObject();
    QStringList unwritable;
    for (auto it = changes.begin(); it != changes.end(); ++it) {
        if (const QStringList at = kuttaCasePath(it.key(), kase); !at.isEmpty()) setJsonAt(kase, at, it.value());
        else unwritable << it.key();
    }
    if (!unwritable.isEmpty())
        fail(QStringLiteral("The fixture CLI saves reference values, the operating pressure and inlet speed and turbulence only; point the app "
                            "at the real kutta CLI to save %1.").arg(unwritable.join(", ")), {}, json);
    const QString mesh = kuttaMesh(input);
    const bool same = QFileInfo(input).canonicalFilePath() == QFileInfo(target).canonicalFilePath();
    kase["mesh"] = QDir(QFileInfo(target).absolutePath()).relativeFilePath(mesh);
    writeJson(target, kase);
    if (json)
        out() << QJsonDocument(QJsonObject{{"v", 1}, {"out", QFileInfo(target).absoluteFilePath()}, {"replaced_input", same},
                                           {"omitted", QJsonArray{}}, {"run_controls_not_saved", QJsonArray{}}})
                     .toJson(QJsonDocument::Compact)
              << '\n';
    else eprint("Wrote " + target);
    return 0;
}


// ---- mesh: a SYNTHETIC stand-in for `kutta mesh` ------------------------------------

int mesh(const Args& a)
{
    const bool jsonl = a.value("progress") == "jsonl";
    if (a.positional.isEmpty() || !a.has("out")) usage(jsonl);
    kapp::MeshSetup setup;
    bool haveSetup = false;
    if (a.has("setup")) {
        const QJsonObject o = readJson(a.value("setup"));
        if (o.isEmpty()) fail("cannot read the mesh setup " + a.value("setup"), {}, jsonl);
        if (o.value("v").toInt() != 1) fail("mesh setup: \"v\" must be 1", {}, jsonl);
        setup = kapp::MeshSetup::fromJson(o);
        haveSetup = true;
    }
    if (const QStringList p = setup.problems(); !p.isEmpty()) fail("mesh setup: " + p.join(" "), {}, jsonl);
    const QString outDir = a.value("out");
    if (QFileInfo::exists(QDir(outDir).filePath("mesh.msh")) && !a.has("force")) fail(outDir + " already holds a mesh; pass --force to replace it", {}, jsonl);
    // Cancelled by an interrupt, as the app escalates after closing stdin
    // (stdin is not watched: from a shell it is often closed already).
    std::signal(SIGINT, [](int) { g_stop = true; });
    std::signal(SIGTERM, [](int) { g_stop = true; });
    const double speed = std::max(0.1, qEnvironmentVariable("KUTTA_FIXTURE_SPEED", "6").toDouble());
    const QString hold = qEnvironmentVariable("KUTTA_FIXTURE_MESH_HOLD");
    QElapsedTimer clock;
    clock.start();
    auto stage = [&](const QString& name, double from, double to, const QString& message, double ms) {
        const int steps = 4;
        for (int i = 0; i <= steps; ++i) {
            if (jsonl)
                emitJson({{"v", 1}, {"event", "mesh-progress"}, {"stage", name}, {"fraction", from + (to - from) * i / steps}, {"message", message}});
            if (name == hold && i == 1)
                while (!g_stop) sleepMs(20);
            if (g_stop) {
                if (jsonl) emitJson({{"v", 1}, {"event", "error"}, {"message", "mesh generation cancelled"}, {"advice", "Nothing was written."}});
                std::exit(130);
            }
            sleepMs(ms / steps / speed);
        }
    };

    // read: the solids must be the synthetic car's walls.
    LengthUnit units = setup.units;
    std::vector<kapp::StlSolid> solids;
    stage("read", 0.0, 0.08, QStringLiteral("Reading %1 STL file(s)").arg(a.positional.size()), 600);
    kapp::Bounds raw;
    for (const QString& path : a.positional) {
        auto r = kapp::readStl(path, haveSetup ? units : LengthUnit::Metre);
        if (!r) fail(r.error(), {}, jsonl);
        for (auto& sld : r.value().solids) {
            raw.expand(sld.mesh.bounds());
            solids.push_back(std::move(sld));
        }
    }
    if (!haveSetup) {
        units = kapp::guessUnit(raw);
        const double k = kapp::unitScale(units);
        for (auto& sld : solids) {
            for (float& x : sld.mesh.positions) x = float(x * k);
            sld.area *= k * k;
        }
    }
    const QStringList walls = carWallZones();
    QStringList names;
    for (const auto& sld : solids) names << sld.name;
    QStringList missing, unknown;
    for (const QString& w : walls) if (!names.contains(w)) missing << w;
    for (const QString& n : names) if (!walls.contains(n)) unknown << n;
    if (!missing.isEmpty() || !unknown.isEmpty() || names.removeDuplicates() > 0)
        emitJson({{"v", 1}, {"event", "error"},
                  {"message", QStringLiteral("The fixture CLI meshes only its synthetic car: one solid each named %1 (got %2).")
                                  .arg(walls.join(", "), names.join(", "))},
                  {"advice", "Point the app at the real kutta CLI to mesh other geometry, or write the car with `kutta-fixture sample-stl DIR`."}}),
        std::exit(1);
    kapp::Bounds g;
    double area = 0;
    QMap<QString, SolidTransform> transforms;
    QJsonArray solidsJson;
    for (const auto& sld : solids) {
        const kapp::Bounds b = sld.mesh.bounds();
        g.expand(b);
        area += sld.area;
        const SolidTransform t = fitTransform(carZoneBounds(sld.name), b);
        // A baseline part, re-tessellated: no change.
        if (std::abs(t.scale - 1) > 0.02f || t.offset.length() > 0.01f * (b.max - b.min).length()) transforms[sld.name] = t;
        solidsJson.append(QJsonObject{{"name", sld.name}, {"triangles", double(sld.triangles)}, {"area", sld.area}});
    }
    if (!haveSetup) setup = kapp::defaultMeshSetup(g, units);
    stage("repair", 0.08, 0.15, "Checking surfaces for holes and overlaps", 400);
    const kapp::CellEstimate est = kapp::estimateCells(setup, g, area);
    stage("octree", 0.15, 0.45, QStringLiteral("Refining towards %1 mm on the surface").arg(setup.minSize * 1000, 0, 'g', 3), 1000);
    stage("trim", 0.45, 0.62, "Trimming octree cells at the surfaces", 600);
    stage("layers", 0.62, 0.85, QStringLiteral("%1 prism layers, first height %2 mm").arg(setup.layers).arg(setup.firstHeight * 1000, 0, 'g', 3), 900);
    stage("quality", 0.85, 0.94, "Orthogonal quality, skewness, volume ratio", 300);
    stage("write", 0.94, 1.0, "Writing mesh.msh", 300);

    const qint64 cells = std::max<qint64>(1000, qint64(est.total()));
    const qint64 prisms = qint64(est.layerCells);
    const qint64 rest = cells - prisms;
    const QJsonObject types{{"hex", double(qint64(rest * 0.86))}, {"poly", double(qint64(rest * 0.12))}, {"prism", double(prisms)},
                            {"tet", double(qint64(rest * 0.015))}, {"pyramid", double(rest - qint64(rest * 0.86) - qint64(rest * 0.12) - qint64(rest * 0.015))}};
    QJsonArray zones;
    const QMap<QString, double> share{{"inlet", 0.004}, {"outlet", 0.004}, {"ground", 0.03}, {"tunnel", 0.012}};
    for (const CarZone& z : carZones()) {
        const qint64 faces = walls.contains(z.name) ? qint64(est.surfaceFaces * 0.16) : qint64(cells * share.value(z.name, 0.01));
        zones.append(QJsonObject{{"name", z.name}, {"type", z.type == "wall" && z.name == "ground" ? "wall" : z.type}, {"faces", double(faces)}});
    }
    auto histogram = [](const std::vector<double>& edges, const std::vector<double>& weights, qint64 total) {
        QJsonArray e, c;
        double sum = 0;
        for (double w : weights) sum += w;
        for (double x : edges) e.append(x);
        for (double w : weights) c.append(double(qint64(total * w / sum)));
        return QJsonObject{{"edges", e}, {"counts", c}};
    };
    const std::vector<double> unit{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    const std::vector<double> ratio{1, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20};
    QJsonObject quality{
        {"orthogonal_quality", QJsonObject{{"min", 0.118}, {"mean", 0.884}, {"max", 1.0},
                                           {"histogram", histogram(unit, {0.00002, 0.0001, 0.0006, 0.002, 0.006, 0.015, 0.04, 0.09, 0.22, 0.63}, cells)}}},
        {"skewness", QJsonObject{{"min", 0.0}, {"mean", 0.121}, {"max", 0.874},
                                 {"histogram", histogram(unit, {0.62, 0.22, 0.09, 0.04, 0.016, 0.007, 0.003, 0.0012, 0.0004, 0.00003}, cells)}}},
        {"volume_ratio", QJsonObject{{"min", 1.0}, {"mean", 1.43}, {"max", 11.8},
                                     {"histogram", histogram(ratio, {0.71, 0.2, 0.06, 0.02, 0.006, 0.002, 0.001, 0.0003, 0.0001, 0.00001}, cells)}}}};
    auto worstAt = [&](const QString& zone, const QString& metric, double value) {
        kapp::Bounds b;
        for (const auto& sld : solids) if (sld.name == zone) b = sld.mesh.bounds();
        const QVector3D p(b.max.x(), b.center().y(), b.center().z());
        return QJsonObject{{"metric", metric}, {"value", value}, {"centroid", QJsonArray{p.x(), p.y(), p.z()}}, {"zone", zone}};
    };
    const QJsonArray worst{worstAt("rear-wing", "orthogonal_quality", 0.118), worstAt("front-wing", "orthogonal_quality", 0.142),
                           worstAt("tyre-front", "skewness", 0.874), worstAt("sidepods", "volume_ratio", 11.8)};
    const double L = g.valid ? g.max.x() - g.min.x() : 1.0;
    const QJsonObject report{
        {"v", 1}, {"synthetic", true},
        {"provenance", "SYNTHETIC FIXTURE: kutta-fixture did not mesh anything; counts follow the app's planning estimate and the quality figures are invented."},
        {"cells", QJsonObject{{"total", double(cells)}}}, {"cell_types", types},
        {"faces", double(qint64(cells * 3.4))}, {"nodes", double(qint64(cells * 1.15))},
        {"zones", zones}, {"quality", quality}, {"worst_cells", worst},
        {"repairs", QJsonArray{QJsonObject{{"solid", "sidepods"}, {"kind", "closed small gaps"}, {"count", 2}},
                               QJsonObject{{"solid", "front-wing"}, {"kind", "removed sliver triangles"}, {"count", 14}}}},
        {"y_plus", QJsonObject{{"estimate", kapp::yPlusForFirstHeight(setup.firstHeight, 16, L)}, {"speed", 16}}},
        {"wall_seconds", clock.elapsed() / 1000.0}, {"peak_memory_bytes", double(cells) * 620}};
    QDir().mkpath(outDir);
    QFile m(QDir(outDir).filePath("mesh.msh"));
    if (!m.open(QIODevice::WriteOnly)) fail("cannot write " + m.fileName(), {}, jsonl);
    m.write(kFixtureMeshMarker + ": not a Fluent mesh. kutta-fixture reads it as its synthetic car.\n");
    m.write(QJsonDocument(QJsonObject{{"fixture_mesh", 1}, {"solids", solidsJson}, {"transforms", transformsJson(transforms)},
                                      {"setup", setup.toJson()}, {"cells", report.value("cells").toObject().value("total")},
                                      {"faces", report.value("faces")}, {"nodes", report.value("nodes")}, {"cell_types", types}})
                .toJson(QJsonDocument::Indented));
    m.close();
    writeJson(QDir(outDir).filePath("mesh-report.json"), report);
    if (a.has("case")) {
        QJsonObject kase = carNewCaseReply().value("case").toObject(), boundaries;
        for (const CarZone& z : carZones())
            boundaries[z.name] = starterBoundary(z.name == "inlet" ? "velocity-inlet" : z.name == "outlet" ? "pressure-outlet" : z.name == "ground" ? "moving-ground"
                                                 : z.name == "tunnel" ? "symmetry" : "wall", 15);
        kase["boundaries"] = boundaries;
        kase["mesh"] = "mesh.msh";
        writeJson(QDir(outDir).filePath("case.kutta.json"), kase);
    }
    const QJsonObject done{{"v", 1}, {"event", "mesh-done"}, {"mesh", QFileInfo(m.fileName()).absoluteFilePath()},
                           {"report", QFileInfo(QDir(outDir).filePath("mesh-report.json")).absoluteFilePath()},
                           {"cells", double(cells)}, {"seconds", clock.elapsed() / 1000.0}, {"synthetic", true}};
    if (jsonl || a.has("json")) emitJson(done);
    else eprint(QStringLiteral("Wrote %1 (%2 cells, SYNTHETIC)").arg(m.fileName()).arg(cells));
    return 0;
}

// sample-stl DIR: the synthetic car's walls as ASCII STL in millimetres, for
// trying New Case from Geometry with the fixture: fsae-car.stl (every wall
// but the rear wing), rear-wing.stl, and rear-wing-v2.stl, a larger rear wing
// set higher and further back (the same solid name, for a swap).
int sampleStl(const Args& a)
{
    if (a.positional.size() != 1) usage(false);
    const QString dir = a.positional[0];
    QDir().mkpath(dir);
    const double resolution = 0.3;
    std::vector<std::pair<QString, kapp::MeshData>> body;
    for (const QString& z : carWallZones())
        if (z != "rear-wing") body.push_back({z, carZoneSurface(z, resolution)});
    auto write = [&](const QString& name, const std::vector<std::pair<QString, kapp::MeshData>>& solids) {
        std::vector<std::pair<QString, const kapp::MeshData*>> refs;
        for (const auto& [n, m] : solids) refs.push_back({n, &m});
        QString err;
        if (!kapp::writeAsciiStl(QDir(dir).filePath(name), refs, LengthUnit::Millimetre, &err)) fail(err, {}, false);
        out() << QDir(dir).filePath(name) << '\n';
    };
    write("fsae-car.stl", body);
    kapp::MeshData wing = carZoneSurface("rear-wing", resolution);
    write("rear-wing.stl", {{"rear-wing", wing}});
    const kapp::Bounds b = wing.bounds();
    SolidTransform t;
    t.pivot = b.center();
    t.scale = 1.12f;
    t.offset = QVector3D(0.03f, 0, 0.04f);
    for (std::size_t i = 0; i + 2 < wing.positions.size(); i += 3) {
        const QVector3D p = t.apply(QVector3D(wing.positions[i], wing.positions[i + 1], wing.positions[i + 2]));
        wing.positions[i] = p.x(), wing.positions[i + 1] = p.y(), wing.positions[i + 2] = p.z();
    }
    write("rear-wing-v2.stl", {{"rear-wing", wing}});
    out().flush();
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
    if (cmd == "new-case") return newCase(a);
    if (cmd == "convert") return convert(a);
    if (cmd == "mesh") return mesh(a);
    if (cmd == "sample-stl") return sampleStl(a);
    usage(false);
}
