#pragma once

// Application state with no widgets: the open case, its settings and
// overrides, the run in progress, run history and loaded view datasets.
// Widgets observe it through signals; scenarios and tests drive it directly.

#include "core/CliClient.h"
#include "core/RunHistory.h"
#include "core/SceneData.h"
#include "core/SettingsModel.h"

#include <QElapsedTimer>
#include <QHash>
#include <QObject>
#include <QPointer>

#include <memory>
#include <optional>

namespace kapp {

enum class RunState { Idle, Starting, Running, Stopping };
enum class ConsoleLevel { Command, Info, Warning, Error, Table, Success };

struct LiveRun {
    QString runId;
    int requestedIterations = 0;
    std::optional<SetupEvent> setup;
    RunProgress progress;
    std::optional<IterationEvent> last;
    int lastSnapshotIter = 0;
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

    // History
    RunHistory& history() { return m_history; }
    const RunHistory& history() const { return m_history; }
    void deleteRun(const QString& id);
    void restoreRunSettings(const QString& id);

    // Datasets for viewports: "run:<id>", "live", "diff:<a>:<b>".
    using SceneCallback = std::function<void(const QString& key, Result<std::shared_ptr<SceneData>>)>;
    std::shared_ptr<SceneData> sceneIfLoaded(const QString& key) const { return m_scenes.value(key); }
    void loadScene(const QString& key, SceneCallback done);
    QString latestSceneKey() const;
    std::vector<Plane> slices(const QString& runId) const;
    void setSlices(const QString& runId, const std::vector<Plane>& planes, SceneCallback done);
    std::vector<Plane> defaultSlices() const;
    // Post-processing; results come back as a manifest. Surface parts are
    // merged into the run's dataset.
    void postProcess(const QString& runId, QJsonObject item, std::function<void(Result<QJsonObject>)> done);

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
    void approvalsChanged();
    void consoleLine(kapp::ConsoleLevel level, const QString& text);

private:
    void onProgress(const ProgressEvent& e);
    void finishRun(int exitCode, bool crashed);
    void printResidualRow(const IterationEvent& e);
    void loadSceneNow(const QString& key, SceneCallback done);
    void exportAndLoad(const QString& runId, const std::vector<Plane>& planes, SceneCallback done);

    CliClient* m_cli;
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
    int m_tableRows = 0;
    QStringList m_tableNames;
};

} // namespace kapp
