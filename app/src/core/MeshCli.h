#pragma once

// The one adapter between the app and `kutta mesh` (STL to volume mesh).
// Every command-line flag, event name, file name and report key of that
// contract is spelled in MeshCli.cpp only, so a rename in the CLI is a change
// to that file. The contract is docs/app-protocol.md, "Automatic STL meshing":
//
//   kutta mesh GEOM.stl [MORE.stl ...] --setup mesh.json --out DIR
//              [--speed M/S] [--progress jsonl] [--json] [--force]
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
    double speed = 0;      // m/s for the y+ estimate; 0: the mesher's default
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

// DIR/mesh-report.json (docs/app-protocol.md, "Automatic STL meshing"),
// parsed strictly: a key of the wrong type fails with its path; an optional
// key the report lacks (or reports as null) stays empty, and the page says
// "not reported" rather than inventing a value. All lengths are metres.
struct MeshBin { double from = 0, to = 0; qint64 count = 0; };
struct MeshDistribution {      // adjacent-cell volume ratio over a set of faces
    std::optional<double> max, p50, p99, p999;
    std::optional<qint64> faces;
    std::vector<MeshBin> histogram;
};
struct MeshQualityReport {
    std::optional<double> orthogonalMin, orthogonalMean, skewnessMax, volumeMin;
    std::vector<MeshBin> orthogonalHistogram;
    std::optional<qint64> cellsBelow005, cellsBelow015, negativeCells, nonconvexCells;
    std::optional<QVector3D> worstCell;  // centroid of the lowest orthogonal quality cell
    std::optional<MeshDistribution> volumeRatioCore, volumeRatioLayers;
};
struct MeshGap {               // castellation/gaps: a slot between two zones (or a zone and the ground)
    QString a, b;
    qint64 cellsAcross = 0;    // fluid cells across at its narrowest; 0: the slot closed
    double width = 0;          // m
    bool contact = false;
    std::optional<std::pair<QVector3D, QVector3D>> sample;
    bool closed() const { return cellsAcross == 0 && !contact; }
};
struct MeshLayers {
    std::optional<qint64> requested, generated, cells, stacks, collapseRounds, splitPoints, liftedPoints;
    std::optional<double> coverageFull, coverageThinned, coverageNone, wallArea;  // area fractions of the layered walls; m²
    std::optional<double> firstHeightMean, firstHeightMin;                         // m
};
struct MeshSnap {
    std::optional<double> distanceMax, distanceMean, distanceP95, distanceP99;  // m, point to STL surface
    std::optional<double> relativeMax, relativeP99;                             // over the local cell size
    std::optional<qint64> wallPoints, partiallySnapped, unsnapped, featurePoints;
};
struct MeshYPlus {             // estimated_first_cell_y_plus: flat-plate estimate, not solved y+
    std::optional<double> min, p05, p50, p95, max, speed;
    std::optional<double> below5, from5to30, from30to300, above300;  // area fractions
    QString method;
};
struct MeshCellTypes {
    std::optional<qint64> hex, hexWithoutHangingNodes, poly, prism;
    std::optional<double> hexFractionOfCore;
};
struct MeshReportZone { QString name, type; qint64 faces = 0; };
struct MeshReport {
    qint64 cells = 0, faces = 0, nodes = 0;
    MeshCellTypes cellTypes;
    std::vector<MeshReportZone> zones;
    std::optional<MeshQualityReport> quality;
    std::optional<std::vector<MeshGap>> gaps;   // castellation/gaps
    std::optional<MeshLayers> layers;
    std::optional<MeshSnap> snap;
    std::optional<MeshYPlus> yPlus;
    QMap<QString, qint64> repairs;              // repairs' counts (open_edges, reoriented_triangles, ...)
    std::optional<QStringList> missingSurfaceZones;
    QStringList warnings;
    std::optional<double> wallSeconds;
    std::optional<qint64> peakMemoryBytes;      // the mesher's own footprint, not the solver's
    QString kuttaCommit;
    QJsonObject raw;
    std::vector<MeshGap> closedGaps() const;
};
Result<MeshReport> parseMeshReport(const QJsonObject& o);
Result<MeshReport> loadMeshReport(const QString& path);

// Solver device memory for a mesh: about 2.6 kB per cell, measured on Metal
// for steady SST on mesher meshes (2,528-2,589 B per cell; mesher REPORT.md).
// An estimate, not a measurement of this mesh.
constexpr double kSolverBytesPerCell = 2600;
// A 12 GB GPU holds about 4.3 million cells at that rate once the CUDA
// context is allowed for (the 4.4M-cell default FSAE mesh used 11.4 GB).
constexpr double kTwelveGbUsableBytes = 4.3e6 * kSolverBytesPerCell;
inline double solverBytesEstimate(qint64 cells) { return double(cells) * kSolverBytesPerCell; }

// One `kutta mesh` run. Cancel closes stdin, then interrupts, then kills; the
// mesher stages its work and publishes DIR's files only on success.
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
