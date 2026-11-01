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
struct RunRecord;
// The status as the user reads it: a completed run says "Converged" when it
// met its convergence criteria and "Iteration limit" when it ran out of
// iterations first.
QString runStatusText(const RunRecord& r);

struct RunRecord {
    QString id;            // "run-0003"
    int number = 0;
    QString dir;           // absolute: the CLI's --out directory
    QString note;          // user label
    QString casePath, caseSha, dataPath;
    QJsonObject overrides;     // case-changing overrides (MODIFIED)
    QJsonObject runControls;   // iterations, monitor interval, convergence criteria
    QStringList approximations;
    QStringList approximationLabels;  // badge words per approximation ("WFGC off"), parallel
    QStringList convergedBy;   // conditions met when it stopped converged
    QString studyId;
    QString backend, device;
    QString coupling;      // pressure-velocity coupling that ran ("simplec"); empty from older CLIs
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
    // "APPROXIMATE · WFGC off": which approximations, in the CLI's words.
    QString approximateBadge() const;
    bool modified() const { return !overrides.isEmpty(); }
    bool converged() const { return stopReason == QLatin1String("converged"); }
    QJsonObject toJson() const;
    static std::optional<RunRecord> fromJson(const QJsonObject& o, const QString& dir);
};

// Why a finished run stopped, in a sentence: the conditions it converged on,
// a stop by the user, or the iteration limit (with checking on or off).
// `stopReason` is the CLI's; `conv` the run's convergence setup.
QString stopReasonText(const QString& stopReason, const QStringList& metBy, const ConvergenceSetup& conv);

// The run's poor-mesh-numerics cells in a sentence; empty when there were
// none (or the CLI predates the count).
QString poorMeshNumericsText(const RunSummary& s);

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
    // Each condition's state at the last monitor sample (convergence checking on).
    std::optional<ConvergenceState> convergence;
    int convergenceIter = 0;
};
RunProgress parseProgressLog(const QByteArray& jsonl);
void appendIteration(RunProgress& p, const IterationEvent& e);

// What a residual chart marks for one run: its residual criteria (when
// convergence checking was on; equal criteria share one line) and the
// iteration it stopped at, with why. 2D cases omit z-velocity.
struct ConvergenceMark {
    enum class Kind { Criterion, Converged, Stopped, IterationLimit };
    Kind kind = Kind::Criterion;
    double value = 0;  // criterion, or the stopping iteration
    QString label;
};
std::vector<ConvergenceMark> convergenceMarks(const RunProgress& p, bool planar);

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
