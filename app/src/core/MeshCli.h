#pragma once

// The one adapter between the app and `kutta mesh` (STL to volume mesh).
// Every command-line flag, event name, file name and report key of that
// contract is spelled in MeshCli.cpp only, so a rename in the CLI is a change
// to that file. Contract as coded (assumptions in docs/app.md, "New Case from
// Geometry"):
//
//   kutta mesh GEOM.stl [MORE.stl ...] --setup mesh.json --out DIR
//              [--progress jsonl] [--json] [--force]
//   DIR/mesh.msh, DIR/mesh-report.json
//   {"v":1,"event":"mesh-progress","stage":"read|repair|octree|trim|layers|quality|write","fraction":0..1,"message":"..."}
//   {"v":1,"event":"mesh-done",...summary...}
//   {"v":1,"event":"error","message":"...","advice":"..."}

#include "core/CliClient.h"

#include <QElapsedTimer>
#include <QJsonObject>
#include <QMap>
#include <QProcess>
#include <QTimer>

#include <optional>
#include <variant>
#include <vector>

namespace kapp {

struct MeshRequest {
    QStringList geometry;  // STL files
    QString setupFile;     // mesh.json
    QString outDir;
};
QStringList meshArguments(const MeshRequest& r);
QString meshFileIn(const QString& outDir);
QString meshReportFileIn(const QString& outDir);

// Stages in the order the mesher runs them, and Fluent Meshing's words for them.
QStringList meshStages();
QString meshStageLabel(const QString& stage);

struct MeshProgressEvent { QString stage; double fraction = 0; QString message; };
struct MeshDoneEvent { QJsonObject summary; };
struct MeshErrorEvent { QString message, advice; };
using MeshEvent = std::variant<MeshProgressEvent, MeshDoneEvent, MeshErrorEvent>;
Result<MeshEvent> parseMeshEvent(const QByteArray& line);

// DIR/mesh-report.json, read tolerantly: anything the report lacks stays empty
// and the page says so rather than inventing it.
struct MeshQuality {
    QString key;     // "orthogonal_quality", "skewness", "volume_ratio", ...
    QString label;   // Fluent's words
    std::optional<double> min, mean, max;
    bool higherIsBetter = true;
    // Optional histogram (bin edges, counts); empty when the report has none.
    std::vector<double> edges;
    std::vector<qint64> counts;
};
struct MeshWorstCell {
    QString metric;
    double value = 0;
    std::optional<QVector3D> position;
    QString zone;
};
struct MeshReportZone { QString name, type; qint64 faces = 0; };
struct MeshReport {
    qint64 cells = 0, faces = 0, nodes = 0;
    QMap<QString, qint64> cellTypes;   // hex, prism, poly, tet, pyramid, ...
    std::vector<MeshReportZone> zones;
    std::vector<MeshQuality> quality;
    std::vector<MeshWorstCell> worst;
    QStringList repairs;
    std::optional<double> yPlus;       // estimated first-cell y+ (mean or typical)
    std::optional<double> yPlusSpeed;  // the speed it assumes, when stated
    std::optional<double> wallSeconds;
    std::optional<qint64> peakMemoryBytes;
    bool synthetic = false;            // fixture data
    QString provenance;
    QJsonObject raw;
    const MeshQuality* metric(const QString& key) const;
};
Result<MeshReport> parseMeshReport(const QJsonObject& o);
Result<MeshReport> loadMeshReport(const QString& path);

// One `kutta mesh` run. Cancel closes stdin, then interrupts, then kills; the
// mesher is assumed to publish DIR's files only on success.
class MeshProcess : public QObject {
    Q_OBJECT
public:
    MeshProcess(const CliInfo& cli, const MeshRequest& req, QObject* parent = nullptr);
    ~MeshProcess() override;
    void start();
    void cancel();
    bool running() const { return m_proc.state() != QProcess::NotRunning; }
    bool cancelled() const { return m_cancelled; }
    QStringList arguments() const { return meshArguments(m_req); }
    const MeshRequest& request() const { return m_req; }
    double elapsedSeconds() const { return m_clock.isValid() ? m_clock.elapsed() / 1000.0 : 0; }

signals:
    void meshEvent(const kapp::MeshEvent& e);
    void stderrText(const QString& text);
    void protocolError(const QString& message);
    // exitCode -1 on a crash or a failure to start.
    void finished(int exitCode, bool cancelled);

private:
    void readStdout();
    CliInfo m_cli;
    MeshRequest m_req;
    QProcess m_proc;
    QByteArray m_partial;
    QTimer m_escalate;
    QElapsedTimer m_clock;
    bool m_cancelled = false;
    int m_escalation = 0;
};

} // namespace kapp
