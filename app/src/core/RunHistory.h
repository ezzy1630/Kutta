#pragma once

// Every solve is kept: the CLI's output directory plus the app's own record
// (app-run.json) and the progress events it received (progress.jsonl). Runs
// belong to a case, identified by the case file's SHA-256.

#include "core/Protocol.h"

#include <QDateTime>
#include <QJsonObject>
#include <QString>

#include <optional>
#include <vector>

namespace kapp {

enum class RunStatus { Running, Completed, Stopped, Failed, Interrupted };
QString runStatusLabel(RunStatus s);

struct RunRecord {
    QString id;            // "run-0003"
    int number = 0;
    QString dir;           // absolute: the CLI's --out directory
    QString note;          // user label
    QString casePath, caseSha, dataPath;
    QJsonObject overrides;
    QStringList approximations;
    QString backend, device;
    int iterationsRequested = 0;
    int liveViewInterval = 0;
    QDateTime started, finished;
    RunStatus status = RunStatus::Running;
    int exitCode = 0;
    QString error;
    int iterationsDone = 0;
    QString stopReason;
    double elapsed = 0;
    bool fixture = false, synthetic = false;
    QString cli;

    QString title() const;  // "Run 3" or "Run 3 · note"
    bool approximate() const { return !approximations.isEmpty(); }
    bool modified() const { return !overrides.isEmpty(); }
    QJsonObject toJson() const;
    static std::optional<RunRecord> fromJson(const QJsonObject& o, const QString& dir);
};

// What happened to a run that did not complete and what to do next, in
// plain words (disk full, out of memory, crash, stopped); empty when it completed.
QString runOutcomeAdvice(const RunRecord& r);

// Final values of each report in a run's progress log.
struct ReportSeries { QString name; std::vector<int> iter; std::vector<double> force, coefficient; };
struct RunProgress {
    ResidualHistory residuals;
    std::vector<ReportSeries> reports;
    QStringList warnings;
    std::optional<SetupEvent> setup;
    std::optional<DoneEvent> done;
    std::optional<ErrorEvent> error;
};
RunProgress parseProgressLog(const QByteArray& jsonl);
void appendIteration(RunProgress& p, const IterationEvent& e);

// A row of a settings comparison.
struct SettingDiff { QString path, label, unit; QJsonValue a, b; };
// Overrides relative to the case: what differs between two runs.
std::vector<SettingDiff> diffRuns(const RunRecord& a, const RunRecord& b, const std::vector<Setting>& inventory);

class RunHistory {
public:
    // Runs for one case live under root/<case name>-<sha prefix>/.
    void setRoot(const QString& root) { m_root = root; }
    void setCase(const QString& casePath, const QString& sha);
    QString caseDir() const { return m_caseDir; }

    void reload();
    const std::vector<RunRecord>& runs() const { return m_runs; }
    const RunRecord* find(const QString& id) const;
    const RunRecord* latestFinished() const;

    RunRecord create();                   // allocates the next number and directory
    bool save(const RunRecord& r);        // writes app-run.json and updates the list; false if the file could not be written
    bool remove(const QString& id);       // deletes the run directory
    void setNote(const QString& id, const QString& note);

    RunProgress progress(const QString& id) const;
    std::optional<ForcesTable> forces(const QString& id) const;
    std::optional<RunSummary> summary(const QString& id) const;

private:
    QString m_root, m_caseDir;
    std::vector<RunRecord> m_runs;
};

} // namespace kapp
