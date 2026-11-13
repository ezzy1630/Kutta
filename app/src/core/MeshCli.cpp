#include "core/MeshCli.h"

#include <QDir>
#include <algorithm>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QObject>

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
    a << "--out" << r.outDir << "--progress" << "jsonl" << "--force";
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

std::optional<double> number(const QJsonObject& o, std::initializer_list<const char*> keys)
{
    for (const char* k : keys)
        if (o.value(k).isDouble()) return o.value(k).toDouble();
    return std::nullopt;
}

QString qualityLabel(const QString& key)
{
    static const QMap<QString, QString> l{{"orthogonal_quality", QObject::tr("Orthogonal quality")},
                                          {"skewness", QObject::tr("Skewness")},
                                          {"volume_ratio", QObject::tr("Volume ratio")},
                                          {"aspect_ratio", QObject::tr("Aspect ratio")}};
    return l.value(key, QString(key).replace('_', ' '));
}

} // namespace

const MeshQuality* MeshReport::metric(const QString& key) const
{
    for (const auto& q : quality)
        if (q.key == key) return &q;
    return nullptr;
}

Result<MeshReport> parseMeshReport(const QJsonObject& o)
{
    using R = Result<MeshReport>;
    if (o.isEmpty()) return R::failure(QObject::tr("the mesh report is empty"));
    if (o.contains("v") && o.value("v").toInt() != 1) return R::failure(QObject::tr("mesh report version %1 is not supported").arg(o.value("v").toInt()));
    MeshReport r;
    r.raw = o;
    const QJsonObject mesh = o.value("mesh").isObject() ? o.value("mesh").toObject() : o;
    // Cells: a total, with the types beside it or in their own object.
    const QJsonValue cells = mesh.value("cells");
    QJsonObject types = mesh.value("cell_types").toObject();
    if (cells.isObject()) {
        const QJsonObject c = cells.toObject();
        r.cells = qint64(c.value("total").toDouble());
        for (auto it = c.begin(); it != c.end(); ++it)
            if (it.key() != "total" && it.value().isDouble()) types.insert(it.key(), it.value());
    } else {
        r.cells = qint64(cells.toDouble());
    }
    for (auto it = types.begin(); it != types.end(); ++it) r.cellTypes[it.key()] = qint64(it.value().toDouble());
    if (r.cells == 0)
        for (qint64 n : r.cellTypes) r.cells += n;
    r.faces = qint64(mesh.value("faces").toDouble());
    r.nodes = qint64(mesh.value("nodes").toDouble());
    for (const QJsonValue& v : o.value("zones").toArray()) {
        const QJsonObject z = v.toObject();
        r.zones.push_back({z.value("name").toString(), z.value("type").toString(), qint64(z.value("faces").toDouble())});
    }
    // Quality: {"orthogonal_quality": {"min":..,"mean":..}, ...} or a list of {"name"/"metric", ...}.
    auto addQuality = [&](const QString& key, const QJsonObject& q) {
        MeshQuality m;
        m.key = key;
        m.label = qualityLabel(key);
        m.min = number(q, {"min"});
        m.mean = number(q, {"mean", "average"});
        m.max = number(q, {"max"});
        m.higherIsBetter = key == "orthogonal_quality";
        const QJsonObject h = q.value("histogram").toObject();
        for (const QJsonValue& e : h.value("edges").toArray()) m.edges.push_back(e.toDouble());
        for (const QJsonValue& c : h.value("counts").toArray()) m.counts.push_back(qint64(c.toDouble()));
        if (m.edges.size() != m.counts.size() + 1) m.edges.clear(), m.counts.clear();
        r.quality.push_back(std::move(m));
    };
    const QJsonValue quality = o.value("quality");
    if (quality.isObject()) {
        const QJsonObject q = quality.toObject();
        for (const char* k : {"orthogonal_quality", "skewness", "volume_ratio", "aspect_ratio"})
            if (q.value(k).isObject()) addQuality(k, q.value(k).toObject());
        for (auto it = q.begin(); it != q.end(); ++it)
            if (it.value().isObject() && !r.metric(it.key())) addQuality(it.key(), it.value().toObject());
    } else {
        for (const QJsonValue& v : quality.toArray()) {
            const QJsonObject q = v.toObject();
            addQuality(q.value("name").toString(q.value("metric").toString()), q);
        }
    }
    for (const QJsonValue& v : o.value("worst_cells").toArray()) {
        const QJsonObject w = v.toObject();
        MeshWorstCell c;
        c.metric = w.value("metric").toString();
        c.value = w.value("value").toDouble();
        const QJsonArray p = w.value(w.contains("centroid") ? "centroid" : "position").toArray();
        if (p.size() == 3) c.position = QVector3D(float(p[0].toDouble()), float(p[1].toDouble()), float(p[2].toDouble()));
        c.zone = w.value("zone").toString();
        r.worst.push_back(c);
    }
    for (const QJsonValue& v : o.value("repairs").toArray()) {
        if (v.isString()) r.repairs << v.toString();
        else {
            const QJsonObject x = v.toObject();
            QString text = x.value("message").toString(x.value("description").toString(x.value("kind").toString()));
            if (x.value("count").isDouble()) text = QObject::tr("%1 (%2)").arg(text).arg(qint64(x.value("count").toDouble()));
            if (!x.value("solid").toString().isEmpty()) text = x.value("solid").toString() + ": " + text;
            r.repairs << text;
        }
    }
    const QJsonValue yp = o.value("y_plus");
    if (yp.isDouble()) r.yPlus = yp.toDouble();
    else if (yp.isObject()) {
        r.yPlus = number(yp.toObject(), {"estimate", "mean", "typical"});
        r.yPlusSpeed = number(yp.toObject(), {"speed"});
    } else r.yPlus = number(o, {"first_cell_y_plus", "y_plus_estimate"});
    r.wallSeconds = number(o, {"wall_seconds", "wall_time_s", "seconds"});
    if (auto m = number(o, {"peak_memory_bytes", "peak_memory"})) r.peakMemoryBytes = qint64(*m);
    r.synthetic = o.value("synthetic").toBool();
    r.provenance = o.value("provenance").toString();
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
