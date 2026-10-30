#pragma once

// Application state with no widgets: the open case, its settings and
// overrides, the run in progress, run history and loaded view datasets.
// Widgets observe it through signals; scenarios and tests drive it directly.

#include "core/CliClient.h"
#include "core/RunHistory.h"
#include "core/SceneData.h"
#include "core/SettingsModel.h"
#include "core/Study.h"

#include <QElapsedTimer>
#include <QHash>
#include <QObject>
#include <QPointer>

#include <memory>
#include <optional>

namespace kapp {

enum class RunState { Idle, Starting, Running, Stopping };
enum class ConsoleLevel { Command, Info, Warning, Error, Table, Success };

// export-view --post operation families, in the batch's key order.
enum class PostKind { Streamlines, IsoSurface, Probe, Integral };

// What one post-processing operation produced in the republished run view.
struct PostOutput {
    QString name;
    ViewPart part;        // geometry operations (not integrals)
    QJsonObject report;   // integrals: the reports.json entry
    QString provenance;
    bool synthetic = false;
    std::shared_ptr<const SceneData> scene;  // the run's reloaded dataset
};

struct LiveRun {
    QString runId;
    int requestedIterations = 0;
    std::optional<SetupEvent> setup;
    RunProgress progress;
    std::optional<IterationEvent> last;
    int lastSnapshotIter = 0;
    int lastSnapshotSequence = 0;
    QElapsedTimer clock;
};

class Session : public QObject {
    Q_OBJECT
public:
    Session(const CliInfo& cli, const QString& dataRoot, QObject* parent = nullptr);
    ~Session() override;

    CliClient* cli() { return m_cli; }
    const CliInfo& cliInfo() const { return m_cli->info(); }

    // Devices
    void refreshDevices();
    const std::vector<Device>& devices() const { return m_devices; }
    QString backend() const { return m_backend; }
    void setBackend(const QString& b);

    // Case
    void openCase(const QString& casePath, const QString& dataPath = {});
    bool hasCase() const { return m_case.has_value(); }
    const CaseSummary& caseSummary() const { return *m_case; }
    QString casePath() const { return m_casePath; }
    QString dataPath() const { return m_dataPath; }
    QString caseName() const;
    SettingsModel* settings() { return &m_settings; }
    QStringList recentCases() const;
    // The last case that failed to open and why; cleared when a case opens.
    struct OpenError { QString casePath, message; };
    const OpenError& openError() const { return m_openError; }
    void clearOpenError() { m_openError = {}; }

    // Approximation opt-in: explicit per approximation, reset for each case.
    QStringList approvedApproximations() const { return m_approved; }
    void setApproximationApproved(const QString& name, bool approved);
    int liveViewInterval() const { return m_liveView; }
    void setLiveViewInterval(int n) { m_liveView = n; }
    int requestedIterations() const;

    // Run
    bool canRun(QString* reason = nullptr) const;
    QString run();  // returns the new run id, or empty when refused (see runRefused)
    void stop();
    RunState state() const { return m_state; }
    const LiveRun& live() const { return m_live; }
    bool running() const { return m_state != RunState::Idle; }

    // Design studies: variants of the case run back to back.
    const std::vector<Study>& studies() const { return m_studies; }
    const Study* study(const QString& id) const;
    // Saves a new study made by makeSweep; returns its id, or empty with why.
    QString addStudy(Study s, QString* error = nullptr);
    void deleteStudy(const QString& id);
    // Runs the study's queued variants one after another (one GPU).
    bool startStudy(const QString& id, QString* why = nullptr);
    // Stops the running variant and cancels the queued ones.
    void stopStudy();
    QString activeStudy() const { return m_activeStudy; }

    // Aero post-processing inputs, kept with the case's runs (aero-spec.json).
    const AeroSpec& aeroSpec() const { return m_aeroSpec; }
    bool aeroSpecSaved() const { return m_aeroSpecSaved; }
    void setAeroSpec(const AeroSpec& spec);
    // Coefficients, breakdown and balance of a finished run (kutta aero, cached
    // until the spec changes).
    void aero(const QString& runId, std::function<void(Result<AeroResult>)> done);
    std::optional<AeroResult> aeroIfComputed(const QString& runId) const;
    // Car walls: the case's force-report zones, else every wall but the ground.
    QStringList carZones() const;

    // History
    RunHistory& history() { return m_history; }
    const RunHistory& history() const { return m_history; }
    void deleteRun(const QString& id);
    void restoreRunSettings(const QString& id);

    // Datasets for viewports: "run:<id>", "live", "diff:<a>:<b>", and the open
    // case's boundary mesh, available as soon as the case is read.
    QString meshSceneKey() const;
    using SceneCallback = std::function<void(const QString& key, Result<std::shared_ptr<SceneData>>)>;
    std::shared_ptr<SceneData> sceneIfLoaded(const QString& key) const { return m_scenes.value(key); }
    void loadScene(const QString& key, SceneCallback done);
    QString latestSceneKey() const;
    std::vector<Plane> slices(const QString& runId) const;
    void setSlices(const QString& runId, const std::vector<Plane>& planes, SceneCallback done);
    std::vector<Plane> defaultSlices() const;
    // Post-processing. Each run keeps one batch (app-post.json); adding an
    // operation replaces any of the same name and republishes the whole view,
    // so traces, iso-surfaces and probes persist with the run.
    void postProcess(const QString& runId, PostKind kind, const QJsonObject& op, std::function<void(Result<PostOutput>)> done);
    // The run's last request of that name (forms reopen with it), or empty.
    QJsonObject postOperation(const QString& runId, const QString& name) const;
    // Fields export-view can sample for a run whose view has been published.
    std::vector<AvailableField> availableFields(const QString& runId) const;

signals:
    void devicesChanged();
    void caseOpening(const QString& path);
    void caseOpened();
    void caseFailed(const QString& message);
    void runStateChanged(kapp::RunState state);
    void runStarted(const QString& id);
    void runFinished(const QString& id);
    void runRefused(const QString& message);
    void iteration(const kapp::IterationEvent& e);
    void liveSceneChanged(int iter);
    void sceneChanged(const QString& key);
    void historyChanged();
    void studiesChanged();
    void aeroSpecChanged();
    void aeroComputed(const QString& runId);
    void approvalsChanged();
    void consoleLine(kapp::ConsoleLevel level, const QString& text);

private:
    QString start(const QJsonObject& caseOverrides, const QJsonObject& runControls, const QStringList& approximations,
                  const QString& studyId);
    void startNextVariant();
    void saveStudy(const Study& s);
    void failOpen(const QString& casePath, const QString& message);
    void onProgress(const ProgressEvent& e);
    void finishRun(int exitCode, bool crashed);
    void printResidualRow(const IterationEvent& e);
    void loadSceneNow(const QString& key, SceneCallback done);
    void exportAndLoad(const QString& runId, const std::vector<Plane>& planes, SceneCallback done);
    void exportNow(const QString& runId, const std::vector<Plane>& planes, bool retried, SceneCallback done);
    QJsonObject postBatch(const QString& runId) const;
    void writePostBatch(const QString& runId, const QJsonObject& batch) const;

    CliClient* m_cli;
    QString m_dataRoot;
    OpenError m_openError;
    std::vector<Device> m_devices;
    QString m_backend = "auto";
    std::optional<CaseSummary> m_case;
    QString m_casePath, m_dataPath;
    SettingsModel m_settings;
    QStringList m_approved;
    int m_liveView = 0;
    RunState m_state = RunState::Idle;
    QPointer<SolveProcess> m_proc;
    std::optional<RunRecord> m_current;
    LiveRun m_live;
    RunHistory m_history;
    QHash<QString, std::shared_ptr<SceneData>> m_scenes;
    QHash<QString, std::vector<SceneCallback>> m_loading;
    QHash<QString, ViewManifest> m_views;  // run id -> last published view
    // One export-view per run directory at a time: it republishes view/.
    QHash<QString, std::vector<std::function<void()>>> m_exportQueue;
    std::vector<Study> m_studies;
    StudyStore m_studyStore;
    QString m_activeStudy;
    AeroSpec m_aeroSpec;
    bool m_aeroSpecSaved = false;
    QHash<QString, AeroResult> m_aero;
    QHash<QString, std::vector<std::function<void(Result<AeroResult>)>>> m_aeroPending;
    int m_tableRows = 0;
    QStringList m_tableNames;
};

} // namespace kapp
