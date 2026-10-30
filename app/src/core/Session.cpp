#include "core/Session.h"

#include "core/ChartScale.h"

#include <QCryptographicHash>
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
    : QObject(parent), m_cli(new CliClient(cli, this)), m_dataRoot(dataRoot)
{
    m_history.setRoot(QDir(dataRoot).filePath("runs"));
    connect(m_cli, &CliClient::commandStarted, this, [this](const QString& cmd) { emit consoleLine(ConsoleLevel::Command, cmd); });
    connect(m_cli, &CliClient::commandOutput, this, [this](const QString& text) {
        for (const QString& l : text.split('\n', Qt::SkipEmptyParts)) emit consoleLine(ConsoleLevel::Info, l);
    });
}

Session::~Session()
{
    // The process's destructor closes stdin and waits for kutta to write its
    // outputs; nothing may call back into this half-destroyed session meanwhile.
    if (m_proc) {
        m_proc->disconnect(this);
        delete m_proc.data();
    }
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
    rec.overrides = m_settings.caseOverrides();
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
    req.iterations = rec.iterationsRequested;
    if (m_settings.needsSetupFile()) {
        req.setupFile = QDir(rec.dir).filePath("setup.json");
        QFile f(req.setupFile);
        if (f.open(QIODevice::WriteOnly)) f.write(QJsonDocument(m_settings.setupJson()).toJson(QJsonDocument::Indented));
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
            // Snapshots are immutable and sequenced; skip another run's or an older one.
            const QString expected = m_live.setup ? m_live.setup->liveRunId : QString();
            if (!expected.isEmpty() && e.runId != expected) return;
            if (e.sequence > 0 && e.sequence <= m_live.lastSnapshotSequence) return;
            m_live.lastSnapshotSequence = e.sequence;
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

QString Session::meshSceneKey() const
{
    // Keyed by content: a geometry export finishing after another case opened
    // can never be shown as the new case's mesh.
    if (!m_case) return {};
    const QString id = !m_case->sha256.isEmpty() ? m_case->sha256
                                                 : QString::fromLatin1(QCryptographicHash::hash(m_casePath.toUtf8(), QCryptographicHash::Sha256).toHex());
    return "mesh:" + id.left(16);
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

namespace {

// Derived fields added to every exported surface and slice. Each needs
// reference values or SST wall data; with fields_policy "available" kutta
// skips the ones a run cannot provide and lists them in the manifest.
const QStringList kSurfaceFields{"pressure_coefficient", "total_pressure_coefficient", "skin_friction",
                                 "vorticity_magnitude", "velocity_x"};

const char* batchKey(PostKind k)
{
    switch (k) {
    case PostKind::Streamlines: return "streamlines";
    case PostKind::IsoSurface: return "isosurfaces";
    case PostKind::Probe: return "probes";
    case PostKind::Integral: return "integrals";
    }
    return "";
}

} // namespace

QJsonObject Session::postBatch(const QString& runId) const
{
    const RunRecord* r = m_history.find(runId);
    QFile f(r ? QDir(r->dir).filePath("app-post.json") : QString());
    if (r && f.open(QIODevice::ReadOnly)) {
        const QJsonObject o = QJsonDocument::fromJson(f.readAll()).object();
        if (o.value("v").toInt() == 1) return o;
    }
    return {{"v", 1}, {"fields", QJsonArray::fromStringList(kSurfaceFields)}, {"fields_policy", "available"}};
}

void Session::writePostBatch(const QString& runId, const QJsonObject& batch) const
{
    const RunRecord* r = m_history.find(runId);
    if (!r) return;
    QFile f(QDir(r->dir).filePath("app-post.json"));
    if (f.open(QIODevice::WriteOnly)) f.write(QJsonDocument(batch).toJson(QJsonDocument::Indented));
}

QJsonObject Session::postOperation(const QString& runId, const QString& name) const
{
    const QJsonObject batch = postBatch(runId);
    for (const char* k : {"streamlines", "isosurfaces", "probes", "integrals"})
        for (const QJsonValue& v : batch.value(k).toArray())
            if (v.toObject().value("name").toString() == name) return v.toObject();
    return {};
}

std::vector<AvailableField> Session::availableFields(const QString& runId) const
{
    auto it = m_views.find(runId);
    return it == m_views.end() ? std::vector<AvailableField>{} : it->availableFields;
}

void Session::exportAndLoad(const QString& runId, const std::vector<Plane>& planes, SceneCallback done)
{
    auto job = [this, runId, planes, done] {
        exportNow(runId, planes, false, [this, runId, done](const QString& key, Result<std::shared_ptr<SceneData>> s) {
            auto& q = m_exportQueue[runId];
            q.erase(q.begin());
            const auto next = q.empty() ? std::function<void()>{} : q.front();
            if (q.empty()) m_exportQueue.remove(runId);
            done(key, s);
            if (next) next();  // a copy: the queue may change while it runs
        });
    };
    auto& queue = m_exportQueue[runId];
    queue.push_back(job);
    if (queue.size() == 1) job();
}

void Session::exportNow(const QString& runId, const std::vector<Plane>& planes, bool retried, SceneCallback done)
{
    const QString key = "run:" + runId;
    const RunRecord* r = m_history.find(runId);
    if (!r) return done(key, Result<std::shared_ptr<SceneData>>::failure("unknown run " + runId));
    const QString dir = r->dir, title = r->title();
    QJsonObject batch = postBatch(runId);
    writePostBatch(runId, batch);
    QJsonArray arr;
    for (const Plane& p : planes) arr.append(QJsonArray{p.origin.x(), p.origin.y(), p.origin.z(), p.normal.x(), p.normal.y(), p.normal.z()});
    m_cli->exportView(dir, planes, QDir(dir).filePath("app-post.json"),
                      [this, runId, planes, retried, batch, dir, title, arr, key, done](Result<ViewManifest> m) mutable {
        // A CLI from before fields_policy rejects the property, or one
        // unavailable field: fall back once to its own selection of fields.
        if (!m && !retried && (m.error().contains("fields_policy") || (batch.contains("fields") && m.error().contains("unavailable")))) {
            emit consoleLine(ConsoleLevel::Warning, tr("This kutta CLI cannot skip unavailable fields; exporting the fields it can compute."));
            batch.remove("fields_policy");
            batch.remove("fields");
            writePostBatch(runId, batch);
            return exportNow(runId, planes, true, done);
        }
        if (!m) return done(key, Result<std::shared_ptr<SceneData>>::failure(m.error()));
        if (!m.value().unavailableFields.isEmpty())
            emit consoleLine(ConsoleLevel::Info, tr("Not available for this run: %1.").arg(m.value().unavailableFields.join(", ")));
        if (const QJsonObject pc = m.value().postCache; pc.contains("reused") && pc.value("reused").toInt() > 0)
            emit consoleLine(ConsoleLevel::Info, tr("Post-processing reused %1 earlier result(s) and computed %2 in %3 s.")
                                                     .arg(pc.value("reused").toInt()).arg(pc.value("computed").toInt())
                                                     .arg(pc.value("seconds").toDouble(), 0, 'f', 2));
        QFile f(QDir(dir).filePath("app-view.json"));
        if (f.open(QIODevice::WriteOnly)) f.write(QJsonDocument(QJsonObject{{"slices", arr}}).toJson());
        m_views.insert(runId, m.value());
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
                m_views.insert(id, m.value());
                auto scene = kapp::loadScene(m.value(), r->title());
                if (scene) m_scenes.insert(key, scene.value());
                return done(key, scene);
            }
        }
        return exportAndLoad(id, slices(id), done);
    }
    if (key.startsWith("mesh:")) {
        if (!m_case || key != meshSceneKey()) return done(key, Result<std::shared_ptr<SceneData>>::failure("that case is no longer open"));
        // One export per case content, reused whenever the case is reopened.
        const QString dir = QDir(m_dataRoot).filePath("geometry/" + key.mid(5));
        const QString title = caseName();
        auto finish = [this, key, title, done](Result<ViewManifest> m) {
            if (!m) return done(key, Result<std::shared_ptr<SceneData>>::failure(m.error()));
            auto scene = kapp::loadScene(m.value(), tr("Mesh"));
            if (scene) {
                scene.value()->geometryOnly = true;
                m_scenes.insert(key, scene.value());
            }
            done(key, scene);
        };
        const QString manifest = QDir(dir).filePath("geometry/manifest.json");
        if (QFileInfo::exists(manifest))
            if (auto m = loadViewManifest(manifest); m && m.value().geometryOnly) return finish(m);
        QDir().mkpath(dir);
        emit consoleLine(ConsoleLevel::Info, tr("Reading the mesh of %1 for display…").arg(title));
        m_cli->exportGeometry(m_casePath, dir, [this, finish](Result<ViewManifest> m) {
            if (!m) emit consoleLine(ConsoleLevel::Warning, tr("Mesh display unavailable: %1").arg(m.error()));
            finish(m);
        });
        return;
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

void Session::postProcess(const QString& runId, PostKind kind, const QJsonObject& op, std::function<void(Result<PostOutput>)> done)
{
    using R = Result<PostOutput>;
    const RunRecord* r = m_history.find(runId);
    if (!r) return done(R::failure("unknown run"));
    const QString name = op.value("name").toString();
    const QJsonObject before = postBatch(runId);
    QJsonObject batch = before;
    for (const char* k : {"streamlines", "isosurfaces", "probes", "integrals"}) {
        QJsonArray a = batch.value(k).toArray();
        for (qsizetype i = a.size() - 1; i >= 0; --i)
            if (a[i].toObject().value("name").toString() == name) a.removeAt(i);
        if (a.isEmpty()) batch.remove(k);
        else batch[k] = a;
    }
    QJsonArray a = batch.value(batchKey(kind)).toArray();
    a.append(op);
    batch[batchKey(kind)] = a;
    writePostBatch(runId, batch);
    m_scenes.remove("run:" + runId);
    for (const QString& k : m_scenes.keys())
        if (k.startsWith("diff:") && k.contains(runId)) m_scenes.remove(k);
    const QString dir = r->dir;
    exportAndLoad(runId, slices(runId), [this, runId, kind, name, before, dir, done](const QString&, Result<std::shared_ptr<SceneData>> s) {
        if (!s) {
            // A rejected operation must not poison later exports; the CLI kept the old view.
            writePostBatch(runId, before);
            emit consoleLine(ConsoleLevel::Error, "Post-processing failed: " + s.error());
            return done(R::failure(s.error()));
        }
        const ViewManifest& m = m_views[runId];
        PostOutput out;
        out.name = name;
        out.provenance = m.provenance;
        out.synthetic = m.synthetic;
        out.scene = s.value();
        if (kind == PostKind::Integral) {
            QFile f(QDir(m.dir).filePath(m.reportsFile.isEmpty() ? QStringLiteral("reports.json") : m.reportsFile));
            if (!f.open(QIODevice::ReadOnly)) return done(R::failure("kutta wrote no reports.json"));
            for (const QJsonValue& v : QJsonDocument::fromJson(f.readAll()).object().value("reports").toArray())
                if (v.toObject().value("name").toString() == name) out.report = v.toObject();
            if (out.report.isEmpty()) return done(R::failure("reports.json has no report named " + name));
            return done(out);
        }
        const auto& parts = kind == PostKind::Streamlines ? m.streamlines : kind == PostKind::IsoSurface ? m.isosurfaces : m.probes;
        auto it = std::find_if(parts.begin(), parts.end(), [&](const ViewPart& p) { return p.name == name; });
        if (it == parts.end()) return done(R::failure("the exported view has no result named " + name));
        out.part = *it;
        done(out);
    });
}

} // namespace kapp
