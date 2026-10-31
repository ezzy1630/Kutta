#include "core/RunHistory.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QRegularExpression>
#include <QSaveFile>

#include <algorithm>

namespace kapp {

QString runStatusLabel(RunStatus s)
{
    switch (s) {
    case RunStatus::Running: return QStringLiteral("Running");
    case RunStatus::Completed: return QStringLiteral("Completed");
    case RunStatus::Stopped: return QStringLiteral("Stopped");
    case RunStatus::Failed: return QStringLiteral("Failed");
    case RunStatus::Interrupted: return QStringLiteral("Interrupted");
    }
    return {};
}

QString runStatusText(const RunRecord& r)
{
    return r.status == RunStatus::Completed && r.converged() ? QStringLiteral("Converged") : runStatusLabel(r.status);
}

QString stopReasonText(const QString& stopReason, const QStringList& metBy, const ConvergenceSetup& conv)
{
    if (stopReason == "converged") {
        QStringList met;
        for (const QString& id : metBy) met << describeCondition(id, conv);
        return met.isEmpty() ? QStringLiteral("Converged") : QStringLiteral("Converged: %1").arg(met.join("; "));
    }
    if (stopReason == "stopped") return QStringLiteral("Stopped by you; not converged");
    if (stopReason == "iteration_limit")
        return conv.enabled ? QStringLiteral("Iteration limit; convergence criteria not met")
                            : QStringLiteral("Iteration limit; convergence checking was off");
    return {};
}

QString runOutcomeAdvice(const RunRecord& r)
{
    const QString e = r.error.toLower();
    auto has = [&e](std::initializer_list<const char*> words) {
        return std::any_of(words.begin(), words.end(), [&e](const char* w) { return e.contains(QLatin1String(w)); });
    };
    switch (r.status) {
    case RunStatus::Running:
    case RunStatus::Completed: return {};
    case RunStatus::Stopped:
        return QStringLiteral("Stopped by you at iteration %1. The results so far are kept and can be post-processed; run again to continue from the case.")
            .arg(r.iterationsDone);
    case RunStatus::Interrupted:
        return QStringLiteral("Kutta closed while this run was going. Its folder keeps what kutta wrote until then; run again for complete results.");
    case RunStatus::Failed: break;
    }
    // Messages differ by OS and runtime: POSIX strerror, Windows, CUDA.
    if (has({"no space left", "not enough space on the disk", "disk full", "enospc", "disk quota"}))
        return QStringLiteral("The disk is full. Free space on the drive that holds the run folder, then run again; finished runs are kept.");
    if (has({"out of memory", "bad_alloc", "cudaerrormemoryallocation", "outofdevicememory", "insufficient memory"}))
        return QStringLiteral("The device ran out of memory. Close other programs using the GPU, or choose another device, then run again.");
    if (has({"permission denied", "access is denied", "read-only file system"}))
        return QStringLiteral("Kutta may not write to the run folder. Check its permissions, then run again.");
    if (has({"stopped unexpectedly"}))
        return QStringLiteral("kutta stopped before finishing, so this run has no final results. The console shows its last messages; if it happens again, keep the run folder for a bug report.");
    return QStringLiteral("kutta reported the error above; the console has its full messages. Fix the cause, then run again.");
}

namespace {

QString statusKey(RunStatus s)
{
    switch (s) {
    case RunStatus::Running: return "running";
    case RunStatus::Completed: return "completed";
    case RunStatus::Stopped: return "stopped";
    case RunStatus::Failed: return "failed";
    case RunStatus::Interrupted: return "interrupted";
    }
    return "failed";
}

RunStatus statusFrom(const QString& k)
{
    if (k == "running") return RunStatus::Running;
    if (k == "completed") return RunStatus::Completed;
    if (k == "stopped") return RunStatus::Stopped;
    if (k == "interrupted") return RunStatus::Interrupted;
    return RunStatus::Failed;
}

QByteArray readAll(const QString& path)
{
    QFile f(path);
    return f.open(QIODevice::ReadOnly) ? f.readAll() : QByteArray();
}

} // namespace

QString RunRecord::title() const
{
    const QString base = QStringLiteral("Run %1").arg(number);
    return note.isEmpty() ? base : base + " · " + note;
}

QJsonObject RunRecord::toJson() const
{
    return {{"v", 1},
            {"id", id},
            {"number", number},
            {"note", note},
            {"case_path", casePath},
            {"case_sha", caseSha},
            {"data_path", dataPath},
            {"overrides", overrides},
            {"run_controls", runControls},
            {"approximations", QJsonArray::fromStringList(approximations)},
            {"converged_by", QJsonArray::fromStringList(convergedBy)},
            {"study", studyId},
            {"backend", backend},
            {"device", device},
            {"coupling", coupling},
            {"iterations_requested", iterationsRequested},
            {"live_view", liveViewInterval},
            {"started", started.toString(Qt::ISODateWithMs)},
            {"finished", finished.isValid() ? finished.toString(Qt::ISODateWithMs) : QString()},
            {"status", statusKey(status)},
            {"exit_code", exitCode},
            {"error", error},
            {"iterations_done", iterationsDone},
            {"stop_reason", stopReason},
            {"elapsed", elapsed},
            {"fixture", fixture},
            {"synthetic", synthetic},
            {"cli", cli}};
}

std::optional<RunRecord> RunRecord::fromJson(const QJsonObject& o, const QString& dir)
{
    if (o.value("v").toInt() != 1 || o.value("id").toString().isEmpty()) return std::nullopt;
    RunRecord r;
    r.id = o.value("id").toString();
    r.number = o.value("number").toInt();
    r.dir = dir;
    r.note = o.value("note").toString();
    r.casePath = o.value("case_path").toString();
    r.caseSha = o.value("case_sha").toString();
    r.dataPath = o.value("data_path").toString();
    r.overrides = o.value("overrides").toObject();
    r.runControls = o.value("run_controls").toObject();
    for (const QJsonValue& a : o.value("approximations").toArray()) r.approximations << a.toString();
    for (const QJsonValue& a : o.value("converged_by").toArray()) r.convergedBy << a.toString();
    r.studyId = o.value("study").toString();
    r.backend = o.value("backend").toString();
    r.device = o.value("device").toString();
    r.coupling = o.value("coupling").toString();
    r.iterationsRequested = o.value("iterations_requested").toInt();
    r.liveViewInterval = o.value("live_view").toInt();
    r.started = QDateTime::fromString(o.value("started").toString(), Qt::ISODateWithMs);
    r.finished = QDateTime::fromString(o.value("finished").toString(), Qt::ISODateWithMs);
    r.status = statusFrom(o.value("status").toString());
    r.exitCode = o.value("exit_code").toInt();
    r.error = o.value("error").toString();
    r.iterationsDone = o.value("iterations_done").toInt();
    r.stopReason = o.value("stop_reason").toString();
    r.elapsed = o.value("elapsed").toDouble();
    r.fixture = o.value("fixture").toBool();
    r.synthetic = o.value("synthetic").toBool();
    r.cli = o.value("cli").toString();
    return r;
}

void appendIteration(RunProgress& p, const IterationEvent& it)
{
    auto& res = p.residuals;
    for (auto r = it.residuals.begin(); r != it.residuals.end(); ++r)
        if (!res.names.contains(r.key())) {
            res.names << r.key();
            res.values.emplace_back(res.iterations.size(), std::nan(""));
        }
    res.iterations.push_back(it.iter);
    for (int n = 0; n < res.names.size(); ++n)
        res.values[std::size_t(n)].push_back(it.residuals.value(res.names[n], std::nan("")));
    res.secPerIter.push_back(it.secPerIter.value_or(std::nan("")));
    if (it.convergence) {
        p.convergence = it.convergence;
        p.convergenceIter = it.iter;
    }
    for (auto r = it.reports.begin(); r != it.reports.end(); ++r) {
        auto s = std::find_if(p.reports.begin(), p.reports.end(), [&](const ReportSeries& x) { return x.name == r.key(); });
        if (s == p.reports.end()) {
            p.reports.push_back({r.key(), {}, {}, {}});
            s = p.reports.end() - 1;
        }
        s->iter.push_back(it.iter);
        s->force.push_back(r.value().force.value_or(std::nan("")));
        s->coefficient.push_back(r.value().coefficient.value_or(std::nan("")));
    }
}

std::vector<ConvergenceMark> convergenceMarks(const RunProgress& p, bool planar)
{
    std::vector<ConvergenceMark> out;
    if (p.setup && p.setup->convergence.enabled) {
        QMap<double, QStringList> byValue;
        const auto& res = p.setup->convergence.residuals;
        for (auto it = res.begin(); it != res.end(); ++it)
            if (residualShown(it.key(), planar)) byValue[it.value()] << it.key();
        for (auto it = byValue.begin(); it != byValue.end(); ++it) {
            const QString v = QString::number(it.key(), 'g', 3);
            out.push_back({ConvergenceMark::Kind::Criterion, it.key(),
                           byValue.size() == 1 ? QStringLiteral("criterion %1").arg(v)
                                               : QStringLiteral("criterion %1 · %2").arg(v, residualOrder(it.value()).join(", "))});
        }
    }
    if (p.done && !p.residuals.iterations.empty()) {
        const QString why = p.done->stopReason;
        const int at = p.done->iterations;
        if (why == "converged") out.push_back({ConvergenceMark::Kind::Converged, double(at), QStringLiteral("Converged · %1").arg(at)});
        else if (why == "stopped") out.push_back({ConvergenceMark::Kind::Stopped, double(at), QStringLiteral("Stopped by you · %1").arg(at)});
        else if (why == "iteration_limit") out.push_back({ConvergenceMark::Kind::IterationLimit, double(at), QStringLiteral("Iteration limit · %1").arg(at)});
    }
    return out;
}

RunProgress parseProgressLog(const QByteArray& jsonl)
{
    RunProgress p;
    for (const QByteArray& line : jsonl.split('\n')) {
        if (line.trimmed().isEmpty()) continue;
        auto ev = parseProgressLine(line);
        if (!ev) continue;
        std::visit([&](const auto& e) {
            using T = std::decay_t<decltype(e)>;
            if constexpr (std::is_same_v<T, SetupEvent>) p.setup = e;
            else if constexpr (std::is_same_v<T, IterationEvent>) appendIteration(p, e);
            else if constexpr (std::is_same_v<T, WarningEvent>) p.warnings << e.message;
            else if constexpr (std::is_same_v<T, DoneEvent>) p.done = e;
            else if constexpr (std::is_same_v<T, ErrorEvent>) p.error = e;
        }, ev.value());
    }
    return p;
}

std::vector<SettingDiff> diffRuns(const RunRecord& a, const RunRecord& b, const std::vector<Setting>& inventory)
{
    QStringList paths = a.overrides.keys();
    for (const QString& k : b.overrides.keys())
        if (!paths.contains(k)) paths << k;
    std::vector<SettingDiff> out;
    for (const QString& p : paths) {
        const Setting* s = nullptr;
        for (const auto& x : inventory)
            if (x.path == p) { s = &x; break; }
        const QJsonValue caseValue = s ? s->value : QJsonValue();
        const QJsonValue va = a.overrides.contains(p) ? a.overrides.value(p) : caseValue;
        const QJsonValue vb = b.overrides.contains(p) ? b.overrides.value(p) : caseValue;
        if (va == vb) continue;
        out.push_back({p, s ? s->label : p, s ? s->unit : QString(), va, vb});
    }
    auto pseudo = [&](const QString& path, const QString& label, const QJsonValue& va, const QJsonValue& vb) {
        if (va != vb) out.push_back({path, label, {}, va, vb});
    };
    pseudo("run/approximations", "Approximations", a.approximations.join(", "), b.approximations.join(", "));
    pseudo("run/device", "Device", a.device, b.device);
    pseudo("run/data", "Initial data", QFileInfo(a.dataPath).fileName(), QFileInfo(b.dataPath).fileName());
    pseudo("run/iterations-done", "Iterations run", a.iterationsDone, b.iterationsDone);
    return out;
}

void RunHistory::setCase(const QString& casePath, const QString& sha)
{
    QString name = QFileInfo(casePath).fileName();
    name.remove(QRegularExpression("\\.(cas|dat)(\\.h5)?$"));
    name.replace(QRegularExpression("[^A-Za-z0-9._-]"), "_");
    const QString key = sha.isEmpty() ? name : name + "-" + sha.left(12);
    m_caseDir = QDir(m_root).filePath(key);
    QDir().mkpath(m_caseDir);
    reload();
}

void RunHistory::reload()
{
    m_runs.clear();
    if (m_caseDir.isEmpty()) return;
    const QDir d(m_caseDir);
    for (const QFileInfo& fi : d.entryInfoList({"run-*"}, QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name)) {
        const QJsonObject o = QJsonDocument::fromJson(readAll(QDir(fi.absoluteFilePath()).filePath("app-run.json"))).object();
        auto r = RunRecord::fromJson(o, fi.absoluteFilePath());
        if (!r) continue;
        // A record still marked running from an earlier session was interrupted.
        if (r->status == RunStatus::Running) r->status = RunStatus::Interrupted;
        m_runs.push_back(*r);
    }
    std::sort(m_runs.begin(), m_runs.end(), [](const RunRecord& a, const RunRecord& b) { return a.number < b.number; });
}

const RunRecord* RunHistory::find(const QString& id) const
{
    for (const auto& r : m_runs)
        if (r.id == id) return &r;
    return nullptr;
}

const RunRecord* RunHistory::latestFinished() const
{
    for (auto it = m_runs.rbegin(); it != m_runs.rend(); ++it)
        if (it->status == RunStatus::Completed || it->status == RunStatus::Stopped) return &*it;
    return nullptr;
}

RunRecord RunHistory::create()
{
    int next = 1;
    for (const auto& r : m_runs) next = std::max(next, r.number + 1);
    RunRecord r;
    r.number = next;
    r.id = QStringLiteral("run-%1").arg(next, 4, 10, QChar('0'));
    r.dir = QDir(m_caseDir).filePath(r.id);
    QDir().mkpath(r.dir);
    r.started = QDateTime::currentDateTime();
    return r;
}

bool RunHistory::save(const RunRecord& r)
{
    QSaveFile f(QDir(r.dir).filePath("app-run.json"));
    const bool written = f.open(QIODevice::WriteOnly) && f.write(QJsonDocument(r.toJson()).toJson(QJsonDocument::Indented)) >= 0 && f.commit();
    auto it = std::find_if(m_runs.begin(), m_runs.end(), [&](const RunRecord& x) { return x.id == r.id; });
    if (it == m_runs.end()) m_runs.push_back(r);
    else *it = r;
    return written;
}

bool RunHistory::remove(const QString& id)
{
    const RunRecord* r = find(id);
    if (!r || r->status == RunStatus::Running) return false;
    // Only ever delete a directory this history created inside its case folder.
    if (!QFileInfo(r->dir).absolutePath().startsWith(QFileInfo(m_caseDir).absoluteFilePath())) return false;
    QDir(r->dir).removeRecursively();
    std::erase_if(m_runs, [&](const RunRecord& x) { return x.id == id; });
    return true;
}

void RunHistory::setNote(const QString& id, const QString& note)
{
    for (auto& r : m_runs)
        if (r.id == id) {
            r.note = note;
            save(r);
        }
}

RunProgress RunHistory::progress(const QString& id) const
{
    const RunRecord* r = find(id);
    if (!r) return {};
    return parseProgressLog(readAll(QDir(r->dir).filePath("progress.jsonl")));
}

std::optional<ForcesTable> RunHistory::forces(const QString& id) const
{
    const RunRecord* r = find(id);
    if (!r) return std::nullopt;
    auto t = parseForcesCsv(readAll(QDir(r->dir).filePath("forces.csv")));
    if (!t) return std::nullopt;
    return t.value();
}

std::optional<RunSummary> RunHistory::summary(const QString& id) const
{
    const RunRecord* r = find(id);
    if (!r) return std::nullopt;
    const QJsonObject o = QJsonDocument::fromJson(readAll(QDir(r->dir).filePath("summary.json"))).object();
    auto s = parseRunSummary(o);
    if (!s) return std::nullopt;
    return s.value();
}

} // namespace kapp
