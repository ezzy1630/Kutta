#include "core/MeshCli.h"

#include <QDir>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QObject>
#include <QVector3D>

#include <algorithm>
#include <cmath>

#ifdef Q_OS_UNIX
#include <csignal>
#endif

namespace kapp {

// ---- The contract's names (keep them here) ------------------------------------

QStringList meshArguments(const MeshRequest& r)
{
    QStringList a{"mesh"};
    a << r.geometry;
    if (!r.setupFile.isEmpty()) a << "--setup" << r.setupFile;
    a << "--out" << r.outDir;
    if (r.speed > 0) a << "--speed" << QString::number(r.speed, 'g', 10);
    a << "--progress" << "jsonl" << "--force";
    return a;
}

QString meshFileIn(const QString& outDir) { return QDir(outDir).filePath("mesh.msh"); }
QString meshReportFileIn(const QString& outDir) { return QDir(outDir).filePath("mesh-report.json"); }

QStringList meshStages() { return {"read", "repair", "octree", "trim", "layers", "quality", "write"}; }

QString meshStageLabel(const QString& stage)
{
    static const QMap<QString, QString> labels{
        {"read", QObject::tr("Import Geometry")},       {"repair", QObject::tr("Repair Surfaces")},
        {"octree", QObject::tr("Generate the Octree")}, {"trim", QObject::tr("Trim to the Geometry")},
        {"layers", QObject::tr("Add Boundary Layers")}, {"quality", QObject::tr("Check Mesh Quality")},
        {"write", QObject::tr("Write the Mesh")}};
    return labels.value(stage, stage);
}

Result<MeshEvent> parseMeshEvent(const QByteArray& line)
{
    using R = Result<MeshEvent>;
    QJsonParseError pe;
    const QJsonDocument doc = QJsonDocument::fromJson(line, &pe);
    if (!doc.isObject()) return R::failure("not JSON: " + pe.errorString());
    const QJsonObject o = doc.object();
    if (o.value("v").toInt() != 1) return R::failure("unsupported protocol version");
    const QString ev = o.value("event").toString();
    if (ev == "mesh-progress")
        return MeshEvent(MeshProgressEvent{o.value("stage").toString(), std::clamp(o.value("fraction").toDouble(), 0.0, 1.0), o.value("message").toString()});
    if (ev == "mesh-done") return MeshEvent(MeshDoneEvent{o});
    if (ev == "error") return MeshEvent(MeshErrorEvent{o.value("message").toString(), o.value("advice").toString()});
    return R::failure("unknown event " + ev);
}

namespace {

// Strict reads with the offending key's path: a present key of the wrong type
// is an error; absent or null is "not reported".
struct Reader {
    QStringList errors;
    void bad(const QString& path, const char* what) { errors << QStringLiteral("%1: expected %2").arg(path, QString::fromLatin1(what)); }
    static bool absent(const QJsonValue& v) { return v.isUndefined() || v.isNull(); }
    std::optional<double> number(const QJsonObject& o, const QString& key, const QString& path)
    {
        const QJsonValue v = o.value(key);
        if (absent(v)) return std::nullopt;
        if (!v.isDouble()) return bad(path + "/" + key, "a number"), std::nullopt;
        return v.toDouble();
    }
    std::optional<qint64> count(const QJsonObject& o, const QString& key, const QString& path)
    {
        const auto d = number(o, key, path);
        if (!d) return std::nullopt;
        if (*d < 0 || *d != std::floor(*d)) return bad(path + "/" + key, "a whole number"), std::nullopt;
        return qint64(*d);
    }
    std::optional<QJsonObject> object(const QJsonObject& o, const QString& key, const QString& path)
    {
        const QJsonValue v = o.value(key);
        if (absent(v)) return std::nullopt;
        if (!v.isObject()) return bad(path + "/" + key, "an object"), std::nullopt;
        return v.toObject();
    }
    std::optional<QJsonArray> array(const QJsonObject& o, const QString& key, const QString& path)
    {
        const QJsonValue v = o.value(key);
        if (absent(v)) return std::nullopt;
        if (!v.isArray()) return bad(path + "/" + key, "an array"), std::nullopt;
        return v.toArray();
    }
    QString string(const QJsonObject& o, const QString& key, const QString& path)
    {
        const QJsonValue v = o.value(key);
        if (absent(v)) return {};
        if (!v.isString()) return bad(path + "/" + key, "a string"), QString();
        return v.toString();
    }
    std::optional<QVector3D> point(const QJsonObject& o, const QString& key, const QString& path)
    {
        const auto a = array(o, key, path);
        if (!a) return std::nullopt;
        if (a->size() != 3 || !(*a)[0].isDouble() || !(*a)[1].isDouble() || !(*a)[2].isDouble()) return bad(path + "/" + key, "[x, y, z]"), std::nullopt;
        return QVector3D(float((*a)[0].toDouble()), float((*a)[1].toDouble()), float((*a)[2].toDouble()));
    }
    // [{"from":a,"to":b,"cells":n}, ...]: the report counts faces under "cells" too.
    std::vector<MeshBin> bins(const QJsonObject& o, const QString& key, const QString& path)
    {
        std::vector<MeshBin> out;
        const auto a = array(o, key, path);
        if (!a) return out;
        for (qsizetype i = 0; i < a->size(); ++i) {
            const QString p = QStringLiteral("%1/%2/%3").arg(path, key).arg(i);
            if (!(*a)[i].isObject()) {
                bad(p, "an object");
                continue;
            }
            const QJsonObject b = (*a)[i].toObject();
            const auto from = number(b, "from", p), to = number(b, "to", p);
            const auto n = count(b, "cells", p);
            if (!from || !to || !n) {
                bad(p, "from, to and cells");
                continue;
            }
            out.push_back({*from, *to, *n});
        }
        return out;
    }
};

std::optional<MeshDistribution> distribution(Reader& rd, const QJsonObject& q, const QString& key)
{
    const QString path = "quality/" + key;
    const auto o = rd.object(q, key, "quality");
    if (!o) return std::nullopt;
    MeshDistribution d;
    d.max = rd.number(*o, "max", path);
    d.p50 = rd.number(*o, "p50", path);
    d.p99 = rd.number(*o, "p99", path);
    d.p999 = rd.number(*o, "p999", path);
    d.faces = rd.count(*o, "faces", path);
    d.histogram = rd.bins(*o, "histogram", path);
    return d;
}

} // namespace

std::vector<MeshGap> MeshReport::closedGaps() const
{
    std::vector<MeshGap> out;
    if (gaps)
        for (const auto& g : *gaps)
            if (g.closed()) out.push_back(g);
    return out;
}

Result<MeshReport> parseMeshReport(const QJsonObject& o)
{
    using R = Result<MeshReport>;
    if (o.isEmpty()) return R::failure(QObject::tr("the mesh report is empty"));
    if (o.value("v").toInt() != 1 || !o.value("v").isDouble())
        return R::failure(o.contains("v") ? QObject::tr("mesh report version %1 is not supported").arg(o.value("v").toVariant().toString())
                                          : QObject::tr("the mesh report has no version (v)"));
    Reader rd;
    MeshReport r;
    r.raw = o;
    // Size and zones are required.
    const auto cells = rd.count(o, "cells", ""), faces = rd.count(o, "faces", ""), nodes = rd.count(o, "nodes", "");
    if (!cells || !faces || !nodes) rd.errors << QObject::tr("cells, faces and nodes are required");
    r.cells = cells.value_or(0);
    r.faces = faces.value_or(0);
    r.nodes = nodes.value_or(0);
    if (const auto z = rd.array(o, "zones", "")) {
        for (qsizetype i = 0; i < z->size(); ++i) {
            const QString p = QStringLiteral("zones/%1").arg(i);
            const QJsonObject zo = (*z)[i].toObject();
            const QString name = rd.string(zo, "name", p);
            if (name.isEmpty()) rd.bad(p + "/name", "a zone name");
            r.zones.push_back({name, rd.string(zo, "type", p), rd.count(zo, "faces", p).value_or(0)});
        }
    } else {
        rd.errors << QObject::tr("zones is required");
    }
    if (const auto t = rd.object(o, "cell_types", "")) {
        r.cellTypes.hex = rd.count(*t, "hex", "cell_types");
        r.cellTypes.hexWithoutHangingNodes = rd.count(*t, "hex_without_hanging_nodes", "cell_types");
        r.cellTypes.poly = rd.count(*t, "poly", "cell_types");
        r.cellTypes.prism = rd.count(*t, "prism", "cell_types");
        r.cellTypes.hexFractionOfCore = rd.number(*t, "hex_fraction_of_core", "cell_types");
    }
    if (const auto q = rd.object(o, "quality", "")) {
        MeshQualityReport m;
        m.orthogonalMin = rd.number(*q, "orthogonal_min", "quality");
        m.orthogonalMean = rd.number(*q, "orthogonal_mean", "quality");
        m.skewnessMax = rd.number(*q, "skewness_max", "quality");
        m.volumeMin = rd.number(*q, "volume_min", "quality");
        m.orthogonalHistogram = rd.bins(*q, "orthogonal_histogram", "quality");
        m.cellsBelow005 = rd.count(*q, "cells_below_0_05", "quality");
        m.cellsBelow015 = rd.count(*q, "cells_below_0_15", "quality");
        m.negativeCells = rd.count(*q, "negative_cells", "quality");
        m.nonconvexCells = rd.count(*q, "nonconvex_cells", "quality");
        if (const auto w = rd.object(*q, "worst_cell", "quality")) m.worstCell = rd.point(*w, "centroid", "quality/worst_cell");
        m.volumeRatioCore = distribution(rd, *q, "adjacent_volume_ratio_core");
        m.volumeRatioLayers = distribution(rd, *q, "adjacent_volume_ratio_layers");
        r.quality = m;
    }
    if (const auto c = rd.object(o, "castellation", ""))
        if (const auto g = rd.array(*c, "gaps", "castellation")) {
            std::vector<MeshGap> gaps;
            for (qsizetype i = 0; i < g->size(); ++i) {
                const QString p = QStringLiteral("castellation/gaps/%1").arg(i);
                const QJsonObject go = (*g)[i].toObject();
                const auto between = rd.array(go, "between", p);
                const auto across = rd.count(go, "min_cells_across", p);
                const auto width = rd.number(go, "width_m", p);
                if (!between || between->size() != 2 || !across || !width) {
                    rd.bad(p, "between [zone, zone], min_cells_across and width_m");
                    continue;
                }
                const auto contact = go.value("contact");
                if (!Reader::absent(contact) && !contact.isBool()) {
                    rd.bad(p + "/contact", "a boolean");
                    continue;
                }
                gaps.push_back({(*between)[0].toString(), (*between)[1].toString(), *across, *width,
                                contact.toBool(*width == 0)});
            }
            r.gaps = gaps;
        }
    if (const auto l = rd.object(o, "layers", "")) {
        MeshLayers m;
        m.requested = rd.count(*l, "requested", "layers");
        m.generated = rd.count(*l, "generated", "layers");
        m.cells = rd.count(*l, "cells", "layers");
        m.stacks = rd.count(*l, "stacks", "layers");
        m.collapseRounds = rd.count(*l, "collapse_rounds", "layers");
        if (const auto c = rd.object(*l, "coverage", "layers")) {
            m.coverageFull = rd.number(*c, "full", "layers/coverage");
            m.coverageThinned = rd.number(*c, "thinned", "layers/coverage");
            m.coverageNone = rd.number(*c, "none", "layers/coverage");
            m.wallArea = rd.number(*c, "wall_area_m2", "layers/coverage");
        }
        if (const auto f = rd.object(*l, "first_height_achieved", "layers")) {
            m.firstHeightMean = rd.number(*f, "mean", "layers/first_height_achieved");
            m.firstHeightMin = rd.number(*f, "min", "layers/first_height_achieved");
        }
        if (const auto t = rd.object(*l, "terminations", "layers")) {
            m.splitPoints = rd.count(*t, "split_points", "layers/terminations");
            m.liftedPoints = rd.count(*t, "lifted_points", "layers/terminations");
        }
        r.layers = m;
    }
    if (const auto sn = rd.object(o, "snap", "")) {
        MeshSnap m;
        if (const auto d = rd.object(*sn, "point_distance_m", "snap")) {
            m.distanceMax = rd.number(*d, "max", "snap/point_distance_m");
            m.distanceMean = rd.number(*d, "mean", "snap/point_distance_m");
            m.distanceP95 = rd.number(*d, "p95", "snap/point_distance_m");
            m.distanceP99 = rd.number(*d, "p99", "snap/point_distance_m");
        }
        if (const auto d = rd.object(*sn, "point_distance_over_cell_size", "snap")) {
            m.relativeMax = rd.number(*d, "max", "snap/point_distance_over_cell_size");
            m.relativeP99 = rd.number(*d, "p99", "snap/point_distance_over_cell_size");
        }
        m.wallPoints = rd.count(*sn, "wall_points", "snap");
        m.partiallySnapped = rd.count(*sn, "partially_snapped", "snap");
        m.unsnapped = rd.count(*sn, "unsnapped", "snap");
        m.featurePoints = rd.count(*sn, "feature_points", "snap");
        r.snap = m;
    }
    if (const auto y = rd.object(o, "estimated_first_cell_y_plus", "")) {
        const QString p = "estimated_first_cell_y_plus";
        MeshYPlus m;
        m.min = rd.number(*y, "min", p);
        m.p05 = rd.number(*y, "p05", p);
        m.p50 = rd.number(*y, "p50", p);
        m.p95 = rd.number(*y, "p95", p);
        m.max = rd.number(*y, "max", p);
        m.speed = rd.number(*y, "speed_m_s", p);
        m.below5 = rd.number(*y, "area_fraction_below_5", p);
        m.from5to30 = rd.number(*y, "area_fraction_5_to_30", p);
        m.from30to300 = rd.number(*y, "area_fraction_30_to_300", p);
        m.above300 = rd.number(*y, "area_fraction_above_300", p);
        m.method = rd.string(*y, "method", p);
        r.yPlus = m;
    }
    if (const auto rep = rd.object(o, "repairs", ""))
        for (auto it = rep->begin(); it != rep->end(); ++it)
            if (it.value().isDouble() && it.value().toDouble() == std::floor(it.value().toDouble())) r.repairs[it.key()] = qint64(it.value().toDouble());
    if (const auto m = rd.array(o, "missing_surface_zones", "")) {
        QStringList names;
        for (const QJsonValue& v : *m) names << v.toString();
        r.missingSurfaceZones = names;
    }
    if (const auto w = rd.array(o, "warnings", ""))
        for (const QJsonValue& v : *w) r.warnings << v.toString();
    r.wallSeconds = rd.number(o, "wall_seconds", "");
    r.peakMemoryBytes = rd.count(o, "peak_memory_bytes", "");
    r.kuttaCommit = rd.string(o, "kutta_commit", "");
    if (!rd.errors.isEmpty()) {
        for (QString& e : rd.errors) if (e.startsWith('/')) e.remove(0, 1);
        return R::failure(QObject::tr("the mesh report does not match kutta mesh v1: %1").arg(rd.errors.join("; ")));
    }
    return r;
}

Result<MeshReport> loadMeshReport(const QString& path)
{
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) return Result<MeshReport>::failure(QObject::tr("no mesh report at %1").arg(QDir::toNativeSeparators(path)));
    QJsonParseError pe;
    const QJsonDocument doc = QJsonDocument::fromJson(f.readAll(), &pe);
    if (!doc.isObject()) return Result<MeshReport>::failure(QObject::tr("unreadable mesh report: %1").arg(pe.errorString()));
    return parseMeshReport(doc.object());
}

// ---- Process --------------------------------------------------------------------

MeshProcess::MeshProcess(const CliInfo& cli, const MeshRequest& req, QObject* parent) : QObject(parent), m_cli(cli), m_req(req)
{
    m_proc.setProgram(cli.program);
    m_proc.setArguments(arguments());
    m_proc.setProcessChannelMode(QProcess::SeparateChannels);
    connect(&m_proc, &QProcess::readyReadStandardOutput, this, &MeshProcess::readStdout);
    connect(&m_proc, &QProcess::readyReadStandardError, this, [this] { emit stderrText(cliMessage(QString::fromUtf8(m_proc.readAllStandardError()))); });
    connect(&m_proc, &QProcess::finished, this, [this](int code, QProcess::ExitStatus st) {
        m_escalate.stop();
        readStdout();
        if (!m_partial.trimmed().isEmpty()) {
            if (auto ev = parseMeshEvent(m_partial.trimmed())) emit meshEvent(ev.value());
            m_partial.clear();
        }
        emit finished(st == QProcess::CrashExit ? -1 : code, m_cancelled);
    });
    connect(&m_proc, &QProcess::errorOccurred, this, [this](QProcess::ProcessError e) {
        if (e != QProcess::FailedToStart) return;
        emit protocolError(QObject::tr("Could not start %1: %2").arg(m_cli.program, m_proc.errorString()));
        emit finished(-1, m_cancelled);
    });
    m_escalate.setSingleShot(true);
    connect(&m_escalate, &QTimer::timeout, this, [this] {
        if (!running()) return;
        if (m_escalation++ == 0) {
#ifdef Q_OS_UNIX
            ::kill(pid_t(m_proc.processId()), SIGINT);
#else
            m_proc.terminate();
#endif
            m_escalate.start(3000);
        } else {
            m_proc.kill();
        }
    });
}

MeshProcess::~MeshProcess()
{
    if (running()) {
        m_proc.kill();
        m_proc.waitForFinished(2000);
    }
}

void MeshProcess::start()
{
    m_clock.start();
    m_proc.start(QIODevice::ReadWrite);
}

void MeshProcess::cancel()
{
    if (!running() || m_cancelled) return;
    m_cancelled = true;
    m_proc.closeWriteChannel();
    m_escalate.start(1500);
}

void MeshProcess::readStdout()
{
    m_partial += m_proc.readAllStandardOutput();
    qsizetype nl;
    while ((nl = m_partial.indexOf('\n')) >= 0) {
        const QByteArray line = m_partial.left(nl).trimmed();
        m_partial.remove(0, nl + 1);
        if (line.isEmpty()) continue;
        auto ev = parseMeshEvent(line);
        if (ev) emit meshEvent(ev.value());
        else emit protocolError(QObject::tr("Unreadable mesh progress line (%1): %2").arg(ev.error(), QString::fromUtf8(line.left(200))));
    }
}

} // namespace kapp
