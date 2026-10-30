#include "core/Session.h"

#include "core/ChartScale.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLocale>
#include <QSettings>

#include <cmath>

namespace kapp {

Session::Session(const CliInfo& cli, const QString& dataRoot, QObject* parent)
    : QObject(parent), m_cli(new CliClient(cli, this))
{
    m_history.setRoot(QDir(dataRoot).filePath("runs"));
    connect(m_cli, &CliClient::commandStarted, this, [this](const QString& cmd) { emit consoleLine(ConsoleLevel::Command, cmd); });
    connect(m_cli, &CliClient::commandOutput, this, [this](const QString& text) {
        for (const QString& l : text.split('\n', Qt::SkipEmptyParts)) emit consoleLine(ConsoleLevel::Info, l);
    });
}

Session::~Session()
{
    if (m_proc) m_proc->stop();
}

void Session::refreshDevices()
{
    m_cli->devices([this](Result<std::vector<Device>> r) {
        if (!r) {
            emit consoleLine(ConsoleLevel::Error, "Device list unavailable: " + r.error());
            return;
        }
        m_devices = r.value();
        emit devicesChanged();
    });
}

void Session::setBackend(const QString& b)
{
    if (b == m_backend) return;
    m_backend = b;
    emit devicesChanged();
}

QString Session::caseName() const
{
    QString n = QFileInfo(m_casePath).fileName();
    for (const char* ext : {".cas.h5", ".cas", ".msh.h5", ".msh"})
        if (n.endsWith(QLatin1String(ext))) { n.chop(int(qstrlen(ext))); break; }
    return n;
}

QStringList Session::recentCases() const { return QSettings().value("recentCases").toStringList(); }

void Session::openCase(const QString& casePath, const QString& dataPath)
{
    if (running()) {
        emit caseFailed("Stop the running calculation before opening another case.");
        return;
    }
    emit caseOpening(casePath);
    emit consoleLine(ConsoleLevel::Info, QStringLiteral("Reading %1%2").arg(QDir::toNativeSeparators(casePath),
                     dataPath.isEmpty() ? QString() : " with data " + QDir::toNativeSeparators(dataPath)));
    m_cli->inspect(casePath, dataPath, [this, casePath, dataPath](Result<CaseSummary> r) {
        if (!r) {
            emit consoleLine(ConsoleLevel::Error, r.error());
            emit caseFailed(r.error());
            return;
        }
        m_case = std::move(r).value();
        m_casePath = casePath;
        m_dataPath = dataPath;
        m_approved.clear();
        m_scenes.clear();
        m_settings.reset(m_case->settings);
        m_history.setCase(casePath, m_case->sha256);
        QSettings s;
        QStringList recent = s.value("recentCases").toStringList();
        recent.removeAll(casePath);
        recent.prepend(casePath);
        s.setValue("recentCases", recent.mid(0, 10));
        const auto& c = *m_case;
        emit consoleLine(ConsoleLevel::Info,
                         QStringLiteral("Case %1: %2 cells, %3 faces, %4 zones, Fluent %5").arg(caseName(),
                             QLocale().toString(c.cells), QLocale().toString(c.faces), QString::number(c.zones.size()),
                             c.fluentVersion.isEmpty() ? "version unknown" : c.fluentVersion));
        for (const auto& i : c.issues)
            if (i.cls != IssueClass::Informational)
                emit consoleLine(i.cls == IssueClass::Unsupported ? ConsoleLevel::Error : ConsoleLevel::Warning,
                                 QStringLiteral("%1 (%2): %3").arg(i.label, i.setting, i.explanation));
        emit caseOpened();
        emit historyChanged();
        emit approvalsChanged();
    });
}

void Session::setApproximationApproved(const QString& name, bool approved)
{
    if (approved && !m_approved.contains(name)) m_approved << name;
    if (!approved) m_approved.removeAll(name);
    emit approvalsChanged();
}

int Session::requestedIterations() const
{
    const QJsonValue v = m_settings.value("run/iterations");
    return v.isDouble() ? v.toInt() : 1000;
}

bool Session::canRun(QString* reason) const
{
    auto no = [&](const QString& r) { if (reason) *reason = r; return false; };
    if (!m_case) return no("Open a Fluent case first.");
    if (running()) return no("A calculation is already running.");
    const auto blocking = m_case->blocking(m_approved);
    if (!blocking.empty()) {
        QStringList names;
        for (const auto* i : blocking) names << i->label;
        const bool approximable = std::all_of(blocking.begin(), blocking.end(),
                                              [](const SupportIssue* i) { return i->cls == IssueClass::Approximable; });
        return no(approximable ? QStringLiteral("Needs your approval to run approximately: %1.").arg(names.join(", "))
                               : QStringLiteral("Kutta cannot run this case: %1.").arg(names.join(", ")));
    }
    return true;
}

QString Session::run()
{
    QString why;
    if (!canRun(&why)) {
        emit runRefused(why);
        return {};
    }
    RunRecord rec = m_history.create();
    rec.casePath = m_casePath;
    rec.caseSha = m_case->sha256;
    rec.dataPath = m_dataPath;
    rec.overrides = m_settings.overrides();
    rec.approximations = m_approved;
    rec.iterationsRequested = requestedIterations();
    rec.liveViewInterval = m_liveView;
    rec.fixture = cliInfo().fixture;
    rec.synthetic = m_case->synthetic;
    rec.cli = cliInfo().program;
    rec.backend = m_backend;
    m_history.save(rec);

    SolveRequest req;
    req.casePath = m_casePath;
    req.dataPath = m_dataPath;
    req.approximations = m_approved;
    req.backend = m_backend;
    req.outDir = rec.dir;
    req.liveViewInterval = m_liveView;
    if (!rec.overrides.isEmpty()) {
        req.setupFile = QDir(rec.dir).filePath("setup.json");
        QFile f(req.setupFile);
        if (f.open(QIODevice::WriteOnly)) f.write(QJsonDocument(m_settings.setupJson()).toJson(QJsonDocument::Indented));
    } else {
        req.iterations = rec.iterationsRequested;
    }
    m_current = rec;
    m_live = LiveRun{};
    m_live.runId = rec.id;
    m_live.requestedIterations = rec.iterationsRequested;
    m_live.clock.start();
    m_tableRows = 0;
    m_tableNames.clear();
    m_scenes.remove("live");

    auto* proc = m_cli->solve(req, this);
    m_proc = proc;
    connect(proc, &SolveProcess::progress, this, &Session::onProgress);
    connect(proc, &SolveProcess::stderrText, this, [this](const QString& t) {
        for (const QString& l : t.split('\n', Qt::SkipEmptyParts)) emit consoleLine(ConsoleLevel::Info, l);
    });
    connect(proc, &SolveProcess::protocolError, this, [this](const QString& m) { emit consoleLine(ConsoleLevel::Error, m); });
    connect(proc, &SolveProcess::finished, this, &Session::finishRun);
    m_state = RunState::Starting;
    emit runStateChanged(m_state);
    emit runStarted(rec.id);
    emit historyChanged();
    proc->start();
    return rec.id;
}

void Session::stop()
{
    if (!m_proc || m_state == RunState::Stopping || m_state == RunState::Idle) return;
    m_state = RunState::Stopping;
    emit consoleLine(ConsoleLevel::Info, "Stop requested: kutta finishes the current iteration and writes its outputs.");
    emit runStateChanged(m_state);
    m_proc->stop();
}

void Session::printResidualRow(const IterationEvent& e)
{
    // Fluent-style residual table in the console.
    if (m_tableNames.isEmpty()) m_tableNames = residualOrder(e.residuals.keys());
    if (m_tableRows % 20 == 0) {
        QString head = QStringLiteral("  iter");
        for (const QString& n : m_tableNames) head += QStringLiteral("  %1").arg(n.left(11), 11);
        head += QStringLiteral("  %1").arg("time/iter", 9);
        emit consoleLine(ConsoleLevel::Table, head);
    }
    QString row = QStringLiteral("%1").arg(e.iter, 6);
    for (const QString& n : m_tableNames)
        row += QStringLiteral("  %1").arg(e.residuals.contains(n) ? QString::number(e.residuals.value(n), 'e', 4) : QString("—"), 11);
    row += QStringLiteral("  %1").arg(e.secPerIter ? QString::number(*e.secPerIter, 'f', 3) + "s" : QString("—"), 9);
    emit consoleLine(ConsoleLevel::Table, row);
    ++m_tableRows;
}

void Session::onProgress(const ProgressEvent& ev)
{
    if (!m_current) return;
    QFile log(QDir(m_current->dir).filePath("progress.jsonl"));
    // Keep the CLI's own events: history plots read them back.
    auto record = [&](const QJsonObject& o) {
        if (log.open(QIODevice::Append)) log.write(QJsonDocument(o).toJson(QJsonDocument::Compact) + "\n");
    };
    std::visit([&](const auto& e) {
        using T = std::decay_t<decltype(e)>;
        if constexpr (std::is_same_v<T, SetupEvent>) {
            m_live.setup = e;
            m_current->backend = e.backend;
            m_current->device = e.device;
            m_history.save(*m_current);
            m_state = RunState::Running;
            emit runStateChanged(m_state);
            emit consoleLine(ConsoleLevel::Info, QStringLiteral("Solving on %1 (%2), %3 cells, initialization: %4")
                                 .arg(e.device, e.backend, QLocale().toString(e.cells), e.init));
            record({{"v", 1}, {"event", "setup"}, {"cells", double(e.cells)}, {"backend", e.backend}, {"device", e.device},
                    {"approximations", QJsonArray::fromStringList(e.approximations)}, {"init", e.init},
                    {"overrides", e.overrides}, {"modified", e.modified}});
        } else if constexpr (std::is_same_v<T, IterationEvent>) {
            appendIteration(m_live.progress, e);
            m_live.last = e;
            m_current->iterationsDone = e.iter;
            printResidualRow(e);
            QJsonObject res, reps;
            for (auto r = e.residuals.begin(); r != e.residuals.end(); ++r) res[r.key()] = r.value();
            for (auto r = e.reports.begin(); r != e.reports.end(); ++r) {
                QJsonObject v;
                if (r.value().force) v["force"] = *r.value().force;
                if (r.value().coefficient) v["coefficient"] = *r.value().coefficient;
                reps[r.key()] = v;
            }
            QJsonObject o{{"v", 1}, {"event", "iteration"}, {"iter", e.iter}, {"residuals", res}, {"reports", reps}};
            if (e.secPerIter) o["sec_per_iter"] = *e.secPerIter;
            if (e.elapsed) o["elapsed"] = *e.elapsed;
            record(o);
            emit iteration(e);
        } else if constexpr (std::is_same_v<T, WarningEvent>) {
            m_live.progress.warnings << e.message;
            emit consoleLine(ConsoleLevel::Warning, e.message);
            record({{"v", 1}, {"event", "warning"}, {"message", e.message}});
        } else if constexpr (std::is_same_v<T, DoneEvent>) {
            m_live.progress.done = e;
            m_current->stopReason = e.stopReason;
            m_current->elapsed = e.elapsed;
            m_current->iterationsDone = e.iterations;
            emit consoleLine(ConsoleLevel::Success,
                             QStringLiteral("Done: %1 iterations in %2 s (%3). Not converged: this case exposes no convergence criteria.")
                                 .arg(e.iterations).arg(e.elapsed, 0, 'f', 1)
                                 .arg(e.stopReason == "stopped" ? "stopped by you" : "iteration limit"));
            record({{"v", 1}, {"event", "done"}, {"converged", e.converged}, {"iterations", e.iterations},
                    {"elapsed", e.elapsed}, {"out", e.out}, {"stop_reason", e.stopReason}});
        } else if constexpr (std::is_same_v<T, ErrorEvent>) {
            m_live.progress.error = e;
            m_current->error = e.message;
            emit consoleLine(ConsoleLevel::Error, e.message.trimmed());
            for (const auto& i : e.issues)
                if (i.cls != IssueClass::Informational)
                    emit consoleLine(ConsoleLevel::Error, QStringLiteral("  %1: %2").arg(i.setting, i.explanation));
            record({{"v", 1}, {"event", "error"}, {"message", e.message}});
        } else if constexpr (std::is_same_v<T, SnapshotEvent>) {
            auto m = loadViewManifest(e.manifest);
            if (!m) {
                emit consoleLine(ConsoleLevel::Warning, "Live snapshot unreadable: " + m.error());
                return;
            }
            auto scene = kapp::loadScene(m.value(), QStringLiteral("Live · iteration %1").arg(e.iter));
            if (!scene) return;
            m_scenes.insert("live", scene.value());
            m_live.lastSnapshotIter = e.iter;
            emit liveSceneChanged(e.iter);
        }
    }, ev);
}

void Session::finishRun(int exitCode, bool crashed)
{
    if (!m_current) return;
    RunRecord& r = *m_current;
    r.exitCode = exitCode;
    r.finished = QDateTime::currentDateTime();
    if (crashed) {
        r.status = RunStatus::Failed;
        if (r.error.isEmpty()) r.error = "kutta stopped unexpectedly (crash or killed); outputs may be incomplete.";
        emit consoleLine(ConsoleLevel::Error, r.error);
    } else if (exitCode != 0 || m_live.progress.error) {
        r.status = RunStatus::Failed;
        if (r.error.isEmpty()) r.error = QStringLiteral("kutta exited with code %1").arg(exitCode);
    } else {
        r.status = r.stopReason == "stopped" ? RunStatus::Stopped : RunStatus::Completed;
    }
    m_history.save(r);
    const QString id = r.id;
    m_current.reset();
    m_state = RunState::Idle;
    if (m_proc) m_proc->deleteLater();
    emit runStateChanged(m_state);
    emit historyChanged();
    emit runFinished(id);
}

void Session::deleteRun(const QString& id)
{
    if (m_history.remove(id)) {
        for (const QString& k : m_scenes.keys())
            if (k.contains(id)) m_scenes.remove(k);
        emit historyChanged();
    }
}

void Session::restoreRunSettings(const QString& id)
{
    const RunRecord* r = m_history.find(id);
    if (!r) return;
    const QString err = m_settings.applyOverrides(r->overrides, QStringLiteral("Restore settings from %1").arg(r->title()));
    if (!err.isEmpty()) emit consoleLine(ConsoleLevel::Warning, "Some settings could not be restored: " + err);
    m_approved = r->approximations;
    emit approvalsChanged();
}

QString Session::latestSceneKey() const
{
    if (const RunRecord* r = m_history.latestFinished()) return "run:" + r->id;
    return {};
}

std::vector<Plane> Session::defaultSlices() const
{
    if (!m_case || !m_case->bounds.valid) return {};
    // External aero (a ground or moving wall): the car's symmetry plane, which
    // is y = 0 whenever the domain spans it. Otherwise slice the domain through
    // its thinnest direction, where internal flows usually vary least.
    const Bounds dom = m_case->bounds;
    bool external = false;
    for (const auto& z : m_case->zones)
        if (z.name.contains("ground") || z.condition.value("wall_motion").toString() == "moving") external = true;
    if (external) {
        QVector3D o = dom.center();
        if (dom.min.y() <= 0 && dom.max.y() >= 0) o.setY(0);
        return {{o, QVector3D(0, 1, 0)}};
    }
    const QVector3D ext = dom.max - dom.min;
    int thin = 0;
    for (int a = 1; a < 3; ++a)
        if (ext[a] < ext[thin]) thin = a;
    QVector3D n;
    n[thin] = 1;
    return {{dom.center(), n}};
}

std::vector<Plane> Session::slices(const QString& runId) const
{
    const RunRecord* r = m_history.find(runId);
    if (!r) return {};
    QFile f(QDir(r->dir).filePath("app-view.json"));
    if (!f.open(QIODevice::ReadOnly)) return defaultSlices();
    std::vector<Plane> out;
    for (const QJsonValue& v : QJsonDocument::fromJson(f.readAll()).object().value("slices").toArray()) {
        const QJsonArray a = v.toArray();
        if (a.size() == 6)
            out.push_back({QVector3D(float(a[0].toDouble()), float(a[1].toDouble()), float(a[2].toDouble())),
                           QVector3D(float(a[3].toDouble()), float(a[4].toDouble()), float(a[5].toDouble()))});
    }
    return out;
}

void Session::exportAndLoad(const QString& runId, const std::vector<Plane>& planes, SceneCallback done)
{
    const RunRecord* r = m_history.find(runId);
    if (!r) return done("run:" + runId, Result<std::shared_ptr<SceneData>>::failure("unknown run " + runId));
    const QString dir = r->dir, title = r->title();
    QJsonArray arr;
    for (const Plane& p : planes) arr.append(QJsonArray{p.origin.x(), p.origin.y(), p.origin.z(), p.normal.x(), p.normal.y(), p.normal.z()});
    m_cli->exportView(dir, planes, [this, runId, dir, title, arr, done](Result<ViewManifest> m) {
        const QString key = "run:" + runId;
        if (!m) return done(key, Result<std::shared_ptr<SceneData>>::failure(m.error()));
        QFile f(QDir(dir).filePath("app-view.json"));
        if (f.open(QIODevice::WriteOnly)) f.write(QJsonDocument(QJsonObject{{"slices", arr}}).toJson());
        auto scene = kapp::loadScene(m.value(), title);
        if (scene) {
            m_scenes.insert(key, scene.value());
            emit sceneChanged(key);
        }
        done(key, scene);
    });
}

void Session::loadScene(const QString& key, SceneCallback done)
{
    if (auto s = m_scenes.value(key)) return done(key, s);
    // Several viewports often ask for the same dataset at once: export it once.
    if (auto it = m_loading.find(key); it != m_loading.end()) {
        it->push_back(std::move(done));
        return;
    }
    m_loading.insert(key, {std::move(done)});
    loadSceneNow(key, [this](const QString& k, Result<std::shared_ptr<SceneData>> r) {
        const auto waiting = m_loading.take(k);
        for (const auto& cb : waiting) cb(k, r);
    });
}

void Session::loadSceneNow(const QString& key, SceneCallback done)
{
    if (key.startsWith("run:")) {
        const QString id = key.mid(4);
        const RunRecord* r = m_history.find(id);
        if (!r) return done(key, Result<std::shared_ptr<SceneData>>::failure("unknown run " + id));
        const QString manifest = QDir(r->dir).filePath("view/manifest.json");
        if (QFileInfo::exists(manifest)) {
            auto m = loadViewManifest(manifest);
            if (m) {
                auto scene = kapp::loadScene(m.value(), r->title());
                if (scene) m_scenes.insert(key, scene.value());
                return done(key, scene);
            }
        }
        return exportAndLoad(id, slices(id), done);
    }
    if (key.startsWith("diff:")) {
        const QStringList parts = key.split(':');
        if (parts.size() != 3) return done(key, Result<std::shared_ptr<SceneData>>::failure("bad comparison key"));
        loadScene("run:" + parts[1], [this, key, parts, done](const QString&, Result<std::shared_ptr<SceneData>> a) {
            if (!a) return done(key, a);
            loadScene("run:" + parts[2], [this, key, a, done](const QString&, Result<std::shared_ptr<SceneData>> b) {
                if (!b) return done(key, b);
                QStringList skipped;
                auto d = runDifference(*a.value(), *b.value(), &skipped);
                if (!skipped.isEmpty())
                    emit consoleLine(ConsoleLevel::Warning, "Not compared (different geometry or fields): " + skipped.join(", "));
                m_scenes.insert(key, d);
                done(key, d);
            });
        });
        return;
    }
    done(key, Result<std::shared_ptr<SceneData>>::failure("no dataset " + key));
}

void Session::setSlices(const QString& runId, const std::vector<Plane>& planes, SceneCallback done)
{
    m_scenes.remove("run:" + runId);
    for (const QString& k : m_scenes.keys())
        if (k.startsWith("diff:") && k.contains(runId)) m_scenes.remove(k);
    exportAndLoad(runId, planes, done);
}

void Session::postProcess(const QString& runId, QJsonObject item, std::function<void(Result<QJsonObject>)> done)
{
    const RunRecord* r = m_history.find(runId);
    if (!r) return done(Result<QJsonObject>::failure("unknown run"));
    const QString dir = r->dir;
    auto go = [this, runId, dir, item, done] {
        m_cli->post(dir, item, [this, runId, dir, done](Result<QJsonObject> m) {
            if (!m) {
                emit consoleLine(ConsoleLevel::Error, "Post-processing failed: " + m.error());
                return done(m);
            }
            const QJsonObject manifest = m.value();
            if (manifest.contains("parts")) {
                if (auto scene = m_scenes.value("run:" + runId)) {
                    const QString err = addPostParts(*scene, manifest, QDir(dir).filePath("post/" + manifest.value("id").toString()));
                    if (!err.isEmpty()) emit consoleLine(ConsoleLevel::Error, err);
                    emit sceneChanged("run:" + runId);
                }
            }
            done(m);
        });
    };
    // Post-processing reads the exported view; make sure it exists first.
    if (!m_scenes.contains("run:" + runId))
        loadScene("run:" + runId, [go, done](const QString&, Result<std::shared_ptr<SceneData>> s) {
            if (!s) return done(Result<QJsonObject>::failure(s.error()));
            go();
        });
    else
        go();
}

} // namespace kapp
