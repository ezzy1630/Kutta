#pragma once

// Runs the `kutta` CLI as a child process (docs/app-protocol.md). The solver is
// never linked into the app: a crash, hang or GPU reset in the solver cannot
// take the workbench down, and the CLI stays the single source of truth.

#include "core/Protocol.h"

#include <QObject>
#include <QPointer>
#include <QProcess>
#include <QTimer>
#include <QVector3D>

#include <functional>

namespace kapp {

struct CliInfo {
    QString program;     // absolute path to kutta / kutta-fixture
    bool fixture = false;
    QString fixtureDir;  // fixtures served by the fixture CLI
    QString describe() const;
};

// Where to find the CLI. A real solver always wins over the fixture: explicit
// path, KUTTA_CLI, the packaged copy, the in-tree build, `kutta` on PATH; only
// then a packaged or development fixture CLI.
CliInfo locateCli(const QString& explicitPath = {}, const QString& devRealCli = {}, const QString& devFixture = {});

struct Plane { QVector3D origin, normal; QString arg() const; };

struct SolveRequest {
    QString casePath, dataPath;
    QStringList approximations;
    QString backend = "auto";
    int iterations = 0;          // 0: leave to the case/setup file
    QString outDir;
    QString setupFile;           // written by the app from the overrides
    int liveViewInterval = 0;    // --live-view OUT --live-every N; 0 disables
};

class SolveProcess : public QObject {
    Q_OBJECT
public:
    SolveProcess(const CliInfo& cli, const SolveRequest& req, QObject* parent = nullptr);
    ~SolveProcess() override;
    void start();
    // Close stdin: the CLI finishes the current iteration, writes outputs and emits done.
    void stop();
    bool running() const { return m_proc.state() != QProcess::NotRunning; }
    const SolveRequest& request() const { return m_req; }
    QStringList arguments() const;

signals:
    void progress(const kapp::ProgressEvent& e);
    void stderrText(const QString& text);
    void protocolError(const QString& message);
    void finished(int exitCode, bool crashed);

private:
    void readStdout();
    CliInfo m_cli;
    SolveRequest m_req;
    QProcess m_proc;
    QByteArray m_partial;
    QTimer m_stopEscalation;
};

class CliClient : public QObject {
    Q_OBJECT
public:
    explicit CliClient(CliInfo info, QObject* parent = nullptr);
    const CliInfo& info() const { return m_info; }

    void devices(std::function<void(Result<std::vector<Device>>)> done);
    void inspect(const QString& casePath, const QString& dataPath, std::function<void(Result<CaseSummary>)> done);
    // Publishes DIR/view atomically. `postFile` (optional) is an export-view
    // --post batch: derived fields, traces, iso-surfaces, probes and reports.
    void exportView(const QString& runDir, const std::vector<Plane>& slices, const QString& postFile,
                    std::function<void(Result<ViewManifest>)> done);

    // The case's boundary zones without a run (export-geometry), published
    // atomically under outDir/geometry/.
    void exportGeometry(const QString& casePath, const QString& outDir, std::function<void(Result<ViewManifest>)> done);

    SolveProcess* solve(const SolveRequest& req, QObject* parent);

signals:
    // Every command line the app runs, for the console.
    void commandStarted(const QString& commandLine);
    void commandOutput(const QString& stderrText);

private:
    using RawDone = std::function<void(int exitCode, const QByteArray& out, const QByteArray& err)>;
    void run(const QStringList& args, RawDone done);
    CliInfo m_info;
};

QString shellQuote(const QStringList& argv);

} // namespace kapp
