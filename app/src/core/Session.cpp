#include "core/Session.h"

#include "core/ChartScale.h"
#include "core/Geometry.h"

#include <QCryptographicHash>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLocale>
#include <QSaveFile>
#include <QRegularExpression>
#include <QSettings>
#include <QTimer>
#include <QTemporaryDir>

#include <cmath>

namespace kapp {

namespace {
// A case meshed from geometry shows each part in the colour it had in the
// geometry view, so the rear wing stays the rear wing's colour throughout.
void keepGeometryColours(SceneData& scene, const SceneData* geometry)
{
    if (!geometry) return;
    for (auto& it : scene.items)
        if (const SceneItem* g = geometry->item(it.id); g && g->color.isValid() && it.color.isValid()) it.color = g->color;
}
} // namespace

Session::Session(const CliInfo& cli, const QString& dataRoot, QObject* parent)
    : QObject(parent), m_cli(new CliClient(cli, this)), m_mesh(new MeshWorkflow(m_cli, this)), m_dataRoot(dataRoot)
{
    connect(m_mesh, &MeshWorkflow::consoleLine, this, &Session::consoleLine);
    // The geometry is built in the app: replace the cached scene on each change.
    connect(m_mesh, &MeshWorkflow::sceneChanged, this, [this] {
        if (auto s = m_mesh->scene()) {
            m_scenes.insert("geometry", s);
        } else {
            m_scenes.remove("geometry");
        }
        emit sceneChanged("geometry");
    });
    connect(m_mesh, &MeshWorkflow::meshStateChanged, this, [this] {
        // A new mesh replaces the cached display of the old one.
        if (!m_mesh->meshing()) m_scenes.remove(generatedMeshSceneKey());
    });
    m_history.setRoot(QDir(dataRoot).filePath("runs"));
    m_notesTimer.setSingleShot(true);
    m_notesTimer.setInterval(300);
    connect(&m_notesTimer, &QTimer::timeout, this, &Session::refreshDraftNotes);
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
    for (const char* ext : {".kutta.json", ".json", ".cas.h5", ".cas", ".msh.h5", ".msh"})
        if (n.endsWith(QLatin1String(ext))) { n.chop(int(qstrlen(ext))); break; }
    return n;
}

void Session::failOpen(const QString& casePath, const QString& message)
{
    m_openError = {casePath, message};
    emit caseFailed(message);
}

QStringList Session::recentCases() const { return QSettings().value("recentCases").toStringList(); }

void Session::openCase(const QString& casePath, const QString& dataPath)
{
    if (running()) {
        failOpen(casePath, tr("Stop the running calculation before opening another case."));
        return;
    }
    for (const QString& p : {casePath, dataPath})
        if (!p.isEmpty() && !QFileInfo(p).isFile()) {
            failOpen(casePath, tr("%1 does not exist (moved, renamed or on a drive that is not connected).").arg(QDir::toNativeSeparators(p)));
            return;
        }
    emit caseOpening(casePath);
    emit consoleLine(ConsoleLevel::Info, QStringLiteral("Reading %1%2").arg(QDir::toNativeSeparators(casePath),
                     dataPath.isEmpty() ? QString() : " with data " + QDir::toNativeSeparators(dataPath)));
    m_cli->inspect(casePath, dataPath, [this, casePath, dataPath](Result<CaseSummary> r) {
        if (!r) {
            emit consoleLine(ConsoleLevel::Error, r.error());
            failOpen(casePath, r.error());
            return;
        }
        m_openError = {};
        m_draft.reset();
        m_draftError.clear();
        m_case = std::move(r).value();
        m_casePath = casePath;
        m_dataPath = dataPath;
        m_approved.clear();
        m_scenes.clear();
        m_settings.reset(m_case->settings);
        m_history.setCase(casePath, m_case->sha256);
        m_studyStore.setDir(m_history.caseDir());
        m_studies = m_studyStore.load();
        m_activeStudy.clear();
        m_aero.clear();
        m_aeroSpec = {};
        m_aeroSpecSaved = false;
        m_runDirAero.clear();
        // A design variant's case starts with its base case's aero settings,
        // copied into the project when the variant was made.
        const QString projectSpec = QDir(QFileInfo(casePath).absolutePath()).filePath("aero-spec.json");
        if (!QFileInfo::exists(QDir(m_history.caseDir()).filePath("aero-spec.json")) && QFileInfo::exists(projectSpec))
            QFile::copy(projectSpec, QDir(m_history.caseDir()).filePath("aero-spec.json"));
        if (QFile f(QDir(m_history.caseDir()).filePath("aero-spec.json")); f.open(QIODevice::ReadOnly)) {
            m_aeroSpec = AeroSpec::fromJson(QJsonDocument::fromJson(f.readAll()).object());
            m_aeroSpecSaved = true;
        } else {
            m_aeroSpec.groups = suggestAeroGroups(carZones());
            // Axes from the case's own force reports. Teams often report downforce
            // as a lift with a downward vector, so a report fixes the up axis, not
            // its sign; a 2D case's up is in its plane.
            const QVector3D up = m_case->planar() ? QVector3D(0, 1, 0) : QVector3D(0, 0, 1);
            m_aeroSpec.lift = up;
            for (const auto& r : m_case->reports) {
                if (!r.direction || r.direction->isNull()) continue;
                const QVector3D d = r.direction->normalized();
                if (r.kind == "drag") m_aeroSpec.drag = d;
                else if (r.kind == "lift" || r.kind == "downforce") m_aeroSpec.lift = QVector3D::dotProduct(d, up) < 0 ? -d : d;
            }
            if (std::abs(QVector3D::dotProduct(m_aeroSpec.drag, m_aeroSpec.lift)) > 1e-6) {
                m_aeroSpec.drag = AeroSpec{}.drag;
                m_aeroSpec.lift = up;
            }
        }
        QSettings s;
        QStringList recent = s.value("recentCases").toStringList();
        recent.removeAll(casePath);
        recent.prepend(casePath);
        s.setValue("recentCases", recent.mid(0, 10));
        const auto& c = *m_case;
        emit consoleLine(ConsoleLevel::Info,
                         QStringLiteral("Case %1: %2 cells, %3 faces, %4 zones, %5").arg(caseName(),
                             QLocale().toString(c.cells), QLocale().toString(c.faces), QString::number(c.zones.size()),
                             c.kuttaCase() ? "Kutta case on " + QFileInfo(c.meshFile).fileName()
                                           : "Fluent " + (c.fluentVersion.isEmpty() ? QStringLiteral("version unknown") : c.fluentVersion)));
        for (const auto& i : c.issues)
            if (i.cls != IssueClass::Informational)
                emit consoleLine(i.cls == IssueClass::Unsupported ? ConsoleLevel::Error : ConsoleLevel::Warning,
                                 QStringLiteral("%1 (%2): %3%4").arg(i.label, i.setting, i.explanationOnly(),
                                     i.measuredEffect.isEmpty() ? QString() : QStringLiteral(" Measured effect: %1.").arg(i.measuredEffect)));
        // A case set up from geometry belongs to its project: record it, and
        // show that project's geometry on the Meshing pages.
        m_draftTarget.clear();
        if (const QString project = MeshWorkflow::projectForCase(casePath); !project.isEmpty()) {
            if (m_mesh->hasProject() && QFileInfo(m_mesh->project().projectFile()) == QFileInfo(project)) {
                if (m_mesh->project().caseFile.isEmpty() || QFileInfo(m_mesh->project().caseFile) == QFileInfo(casePath)) m_mesh->setCaseFile(casePath);
            } else if (!m_mesh->meshing()) {
                m_mesh->openProject(project);
            }
        }
        emit caseOpened();
        emit newCaseChanged();
        emit historyChanged();
        emit studiesChanged();
        emit aeroSpecChanged();
        emit approvalsChanged();
    });
}

void Session::startNewCase(const QString& meshPath, double speed, const QMap<QString, QString>& roles)
{
    m_draftTarget.clear();
    m_pendingWheels.clear();
    const auto revision = ++m_draftRevision;
    m_notesTimer.stop();
    m_notesCurrent = true;
    m_notesError.clear();
    m_draftBusy = false;
    m_draft.reset();
    m_draftError.clear();
    if (!QFileInfo(meshPath).isFile()) {
        m_draftError = tr("%1 does not exist (moved, renamed or on a drive that is not connected).").arg(QDir::toNativeSeparators(meshPath));
        emit newCaseChanged();
        return;
    }
    m_draftBusy = true;
    emit newCaseChanged();
    emit consoleLine(ConsoleLevel::Info, tr("Reading mesh %1").arg(QDir::toNativeSeparators(meshPath)));
    m_cli->newCase({meshPath, speed, roles, {}, {}}, [this, revision, meshPath](Result<NewCaseDraft> r) {
        if (revision != m_draftRevision) return;
        if (!r) {
            m_draftBusy = false;
            m_draftError = r.error();
            emit consoleLine(ConsoleLevel::Error, r.error());
            emit newCaseChanged();
        } else {
            // export-geometry also reads bare meshes. Keep this export local
            // to this request; discarded/replaced drafts cannot receive it.
            auto draft = std::move(r).value();
            auto dir = std::make_shared<QTemporaryDir>();
            m_cli->exportGeometry(meshPath, dir->path(), [this, revision, draft = std::move(draft), dir](Result<ViewManifest> m) mutable {
                if (revision != m_draftRevision) return;
                if (m) {
                    for (const auto& part : m.value().zones) {
                        auto* zone = draft.zone(part.name);
                        if (!zone || zone->kind != "boundary") continue;
                        auto mesh = readMeshFile(QDir(m.value().dir).filePath(part.file), part);
                        if (mesh) zone->wheelEstimate = estimateWheel(mesh.value(), draft.dimension);
                        else draft.geometryError = mesh.error();
                    }
                } else draft.geometryError = m.error();
                // Solids flagged as wheels in the geometry rotate, from the
                // surface estimate (which still needs the user's check).
                bool wheelsSet = false;
                for (const QString& name : std::exchange(m_pendingWheels, {})) {
                    auto* z = draft.zone(name);
                    if (!z || z->kind != "boundary" || z->role != "wall") continue;
                    wheelsSet = true;
                    WheelGeometry manual;
                    manual.axis = draft.dimension == 2 ? std::array<double, 3>{0, 0, 1} : std::array<double, 3>{0, -1, 0};
                    z->wheel = z->wheelEstimate.value_or(manual);
                    z->wheelConfirmed = !z->wheelEstimate;
                }
                m_draft = std::move(draft);
                if (wheelsSet) {
                    // The CLI's notes were for stationary wheels: check again.
                    m_notesCurrent = false;
                    m_draft->notes.clear();
                    m_notesTimer.start();
                }
                m_draftBusy = false;
                emit consoleLine(ConsoleLevel::Info, tr("Mesh %1: %2 cells, %3 boundary zones; roles suggested from the zone names")
                                                     .arg(QFileInfo(m_draft->meshPath).fileName(), QLocale().toString(m_draft->cells))
                                                     .arg(m_draft->boundaries().size()));
                emit newCaseChanged();
            });
        }
    });
}

void Session::setNewCaseRole(const QString& zone, const QString& role)
{
    if (!m_draft || m_draftBusy || role == "interior" || !m_draft->roles.contains(role)) return;
    if (NewCaseZone* z = m_draft->zone(zone); z && z->kind == "boundary" && z->role != role) {
        z->role = role;
        if (role != "wall") z->wheel.reset();
        draftEdited();
    }
}

void Session::setNewCaseSpeed(double speed)
{
    if (!m_draft || m_draftBusy || !(speed > 0) || !std::isfinite(speed) || speed == m_draft->speed) return;
    m_draft->speed = speed;
    draftEdited(false);
}

void Session::setNewCaseWheel(const QString& zone, std::optional<WheelGeometry> wheel)
{
    if (!m_draft || m_draftBusy) return;
    if (auto* z = m_draft->zone(zone); z && z->kind == "boundary" && z->role == "wall") {
        const bool rebuild = z->wheel.has_value() != wheel.has_value();
        z->wheel = std::move(wheel);
        draftEdited(rebuild);
    }
}

void Session::enableNewCaseWheel(const QString& zone, bool enabled)
{
    if (!m_draft || m_draftBusy) return;
    auto* z = m_draft->zone(zone);
    if (!z || z->kind != "boundary" || z->role != "wall") return;
    WheelGeometry manual;
    manual.axis = m_draft->dimension == 2 ? std::array<double, 3>{0, 0, 1} : std::array<double, 3>{0, -1, 0};
    z->wheelConfirmed = !enabled || !z->wheelEstimate;
    setNewCaseWheel(zone, enabled ? std::optional(z->wheelEstimate.value_or(manual)) : std::nullopt);
}

void Session::confirmNewCaseWheel(const QString& zone, bool confirmed)
{
    if (!m_draft || m_draftBusy) return;
    if (auto* z = m_draft->zone(zone); z && z->wheel) {
        z->wheelConfirmed = confirmed;
        draftEdited(false);
    }
}

CliClient::NewCaseRequest Session::draftRequest(const QString& out) const
{
    CliClient::NewCaseRequest req{m_draft->meshPath, m_draft->speed, {}, out, {}};
    for (const auto& z : m_draft->zones) {
        if (z.role != z.suggestedRole) req.roles[z.name] = z.role;
        if (z.role == "wall" && z.wheel) {
            auto wheel = *z.wheel;
            wheel.axis = wheel.rollingAxis(m_draft->dimension);
            req.wheels[z.name] = wheel;
        }
    }
    return req;
}

void Session::draftEdited(bool rebuild)
{
    ++m_draftRevision;
    m_draftError.clear();
    m_notesCurrent = false;
    m_notesError.clear();
    m_draft->notes.clear();  // never present old advice as current
    m_notesTimer.start();
    emit newCaseNotesChanged();
    if (rebuild) emit newCaseChanged();
    else emit newCaseValuesChanged();
}

void Session::refreshDraftNotes()
{
    if (!m_draft || m_draftBusy || m_notesInFlight || !m_draft->wheelError().isEmpty()) return;
    const auto revision = m_draftRevision;
    m_notesInFlight = true;
    m_cli->newCase(draftRequest(), [this, revision](Result<NewCaseDraft> r) {
        m_notesInFlight = false;
        if (!m_draft) return;
        if (revision != m_draftRevision) {
            if (!m_draftBusy && !m_notesCurrent) m_notesTimer.start();
            return;
        }
        if (r) {
            m_draft->notes = r.value().notes;
            m_notesCurrent = true;
        } else {
            m_notesError = tr("Could not refresh the checks: %1").arg(r.error());
        }
        emit newCaseNotesChanged();
    });
}

void Session::discardNewCase()
{
    m_draftTarget.clear();
    ++m_draftRevision;
    m_notesTimer.stop();
    m_draft.reset();
    m_draftError.clear();
    m_draftBusy = false;
    emit newCaseChanged();
}

void Session::createNewCase(const QString& outPath)
{
    if (!m_draft || m_draftBusy) return;
    if (const QString error = m_draft->wheelError(); !error.isEmpty()) {
        m_draftError = error;
        emit newCaseChanged();
        return;
    }
    if (running()) {
        m_draftError = tr("Stop the running calculation before creating a case.");
        emit newCaseChanged();
        return;
    }
    const auto req = draftRequest(outPath);
    const auto revision = ++m_draftRevision;
    m_notesTimer.stop();
    m_draftBusy = true;
    m_draftError.clear();
    emit newCaseChanged();
    m_cli->newCase(req, [this, outPath, revision](Result<NewCaseDraft> r) {
        if (revision != m_draftRevision) return;
        m_draftBusy = false;
        if (!r) {
            m_draftError = r.error();
            emit consoleLine(ConsoleLevel::Error, r.error());
            emit newCaseChanged();
            return;
        }
        emit consoleLine(ConsoleLevel::Success, tr("Wrote %1").arg(QDir::toNativeSeparators(outPath)));
        openCase(outPath);  // clears the draft once the case is open
    });
}

void Session::saveCase(const QString& outPath)
{
    if (!m_case || m_saving) return;
    if (running()) {
        emit saveFailed(tr("Stop the running calculation before saving the case."));
        return;
    }
    // The overrides go to convert as a setup file, as they go to solve.
    QString setupFile;
    const QJsonObject setup = SettingsModel::setupJsonFor(m_settings.caseOverrides(), {});
    if (!setup.value("overrides").toObject().isEmpty()) {
        setupFile = QDir(m_history.caseDir()).filePath("save-setup.json");
        QSaveFile f(setupFile);
        if (!f.open(QIODevice::WriteOnly) || f.write(QJsonDocument(setup).toJson(QJsonDocument::Indented)) < 0 || !f.commit()) {
            emit saveFailed(tr("Kutta cannot write the changed settings in %1.").arg(QDir::toNativeSeparators(m_history.caseDir())));
            return;
        }
    }
    m_saving = true;
    m_cli->convert(m_casePath, setupFile, outPath, [this, outPath, setupFile](Result<ConvertResult> r) {
        m_saving = false;
        if (!setupFile.isEmpty()) QFile::remove(setupFile);
        if (!r) {
            emit consoleLine(ConsoleLevel::Error, tr("Could not save %1: %2").arg(QFileInfo(outPath).fileName(), r.error()));
            emit saveFailed(r.error());
            return;
        }
        for (const auto& o : r->omitted)
            emit consoleLine(ConsoleLevel::Info, tr("Not saved: %1 (%2)").arg(o.path, o.reason));
        emit consoleLine(ConsoleLevel::Success, tr("Saved %1").arg(QDir::toNativeSeparators(outPath)));
        emit caseSaved(outPath);
        // Reopen from the saved file: its settings now include the changes.
        openCase(outPath);
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
    if (!m_case) return no("Open a case first.");
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
    return start(m_settings.caseOverrides(), m_settings.runControls(), m_approved, {});
}

QString Session::start(const QJsonObject& caseOverrides, const QJsonObject& runControls, const QStringList& approximations,
                       const QString& studyId)
{
    RunRecord rec = m_history.create();
    rec.casePath = m_casePath;
    rec.caseSha = m_case->sha256;
    rec.dataPath = m_dataPath;
    rec.overrides = caseOverrides;
    rec.runControls = runControls;
    rec.approximations = approximations;
    for (const QString& a : approximations) rec.approximationLabels << m_case->approximationShort(a);
    rec.studyId = studyId;
    const Setting* iterations = m_settings.find("run/iterations");
    rec.iterationsRequested = runControls.contains("run/iterations") ? runControls.value("run/iterations").toInt()
                              : iterations && iterations->value.isDouble() ? iterations->value.toInt()
                                                                           : 1000;
    rec.liveViewInterval = m_liveView;
    rec.fixture = cliInfo().fixture;
    rec.synthetic = m_case->synthetic;
    rec.cli = cliInfo().program;
    rec.backend = m_backend;
    // A run whose record or setup cannot be written would solve the case
    // without the user's changes, or leave no trace: refuse it instead.
    auto refuse = [this, &rec](const QString& what) {
        m_history.remove(rec.id);
        const QString why = tr("Could not start the run: Kutta cannot write %1 in %2. Check free disk space and that the folder is writable.")
                                .arg(what, QDir::toNativeSeparators(QFileInfo(rec.dir).absolutePath()));
        emit consoleLine(ConsoleLevel::Error, why);
        emit runRefused(why);
        return QString();
    };
    rec.status = RunStatus::Interrupted;  // until the process starts; never left "Running" on disk
    if (!QFileInfo(rec.dir).isDir() || !m_history.save(rec)) return refuse(tr("the run record"));
    rec.status = RunStatus::Running;

    SolveRequest req;
    req.casePath = m_casePath;
    req.dataPath = m_dataPath;
    req.approximations = approximations;
    req.backend = m_backend;
    req.outDir = rec.dir;
    req.liveViewInterval = m_liveView;
    req.iterations = rec.iterationsRequested;
    const QJsonObject setup = SettingsModel::setupJsonFor(caseOverrides, runControls);
    if (!setup.value("overrides").toObject().isEmpty()) {
        req.setupFile = QDir(rec.dir).filePath("setup.json");
        QSaveFile f(req.setupFile);
        if (!f.open(QIODevice::WriteOnly) || f.write(QJsonDocument(setup).toJson(QJsonDocument::Indented)) < 0 || !f.commit())
            return refuse(tr("the changed settings (setup.json)"));
    }
    m_history.save(rec);
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
    // The CLI's own lines, verbatim: history plots and the study table read them back.
    connect(proc, &SolveProcess::progressLine, this, [this](const QByteArray& line) {
        if (!m_current) return;
        QFile log(QDir(m_current->dir).filePath("progress.jsonl"));
        if (log.open(QIODevice::Append)) log.write(line + "\n");
    });
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
    if (m_tableNames.isEmpty()) {
        const bool planar = m_case && m_case->planar();
        for (const QString& n : residualOrder(e.residuals.keys()))
            if (residualShown(n, planar)) m_tableNames << n;
    }
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
    std::visit([&](const auto& e) {
        using T = std::decay_t<decltype(e)>;
        if constexpr (std::is_same_v<T, SetupEvent>) {
            m_live.setup = e;
            m_live.progress.setup = e;
            m_current->backend = e.backend;
            m_current->device = e.device;
            m_current->coupling = e.coupling;
            m_history.save(*m_current);
            m_state = RunState::Running;
            emit runStateChanged(m_state);
            emit consoleLine(ConsoleLevel::Info, QStringLiteral("Solving on %1 (%2), %3 cells, initialization: %4%5")
                                 .arg(e.device, e.backend, QLocale().toString(e.cells), e.init,
                                      e.coupling.isEmpty() ? QString() : QStringLiteral(", %1 pressure-velocity coupling").arg(couplingLabel(e.coupling))));
        } else if constexpr (std::is_same_v<T, IterationEvent>) {
            appendIteration(m_live.progress, e);
            m_live.last = e;
            m_current->iterationsDone = e.iter;
            printResidualRow(e);
            emit iteration(e);
        } else if constexpr (std::is_same_v<T, WarningEvent>) {
            m_live.progress.warnings << e.message;
            emit consoleLine(ConsoleLevel::Warning, e.message);
        } else if constexpr (std::is_same_v<T, DoneEvent>) {
            m_live.progress.done = e;
            m_current->stopReason = e.stopReason;
            m_current->elapsed = e.elapsed;
            m_current->iterationsDone = e.iterations;
            m_current->convergedBy = e.metBy;
            m_current->verdict = e.verdict;
            const ConvergenceSetup conv = m_live.setup ? m_live.setup->convergence : ConvergenceSetup{};
            QString why = stopReasonText(e.verdict, e.metBy, conv);
            // CLIs before stop_reason ran to the iteration limit or nothing.
            if (why.isEmpty()) why = conv.enabled ? QStringLiteral("Iteration limit; convergence criteria not met")
                                                  : QStringLiteral("Iteration limit; convergence checking was off");
            emit consoleLine(e.verdict == Verdict::Converged ? ConsoleLevel::Success
                             : e.verdict == Verdict::ForcesSettledResidualsWaived ? ConsoleLevel::Warning : ConsoleLevel::Info,
                             QStringLiteral("Done: %1 iterations in %2 s. %3.").arg(e.iterations).arg(e.elapsed, 0, 'f', 1).arg(why));
        } else if constexpr (std::is_same_v<T, ErrorEvent>) {
            m_live.progress.error = e;
            m_current->error = e.message;
            if (e.verdict == Verdict::Diverged) m_current->verdict = Verdict::Diverged;
            emit consoleLine(ConsoleLevel::Error, e.message.trimmed());
            for (const auto& i : e.issues)
                if (i.cls != IssueClass::Informational)
                    emit consoleLine(ConsoleLevel::Error, QStringLiteral("  %1: %2").arg(i.setting, i.explanation));
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
            if (m_case && m_case->planar()) makePlanar(*scene.value());
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
        r.status = r.verdict == Verdict::Stopped ? RunStatus::Stopped : RunStatus::Completed;
    }
    m_history.save(r);
    if (r.status != RunStatus::Failed)
        if (const auto summary = m_history.summary(r.id))
            if (const QString poor = poorMeshNumericsText(*summary); !poor.isEmpty()) emit consoleLine(ConsoleLevel::Info, poor + '.');
    const QString id = r.id;
    const RunStatus status = r.status;
    m_current.reset();
    m_state = RunState::Idle;
    if (m_proc) m_proc->deleteLater();
    emit runStateChanged(m_state);
    emit historyChanged();
    emit runFinished(id);
    if (m_activeStudy.isEmpty()) return;
    auto it = std::find_if(m_studies.begin(), m_studies.end(), [this](const Study& x) { return x.id == m_activeStudy; });
    if (it == m_studies.end()) return;
    markFinished(*it, id, status == RunStatus::Completed ? VariantStatus::Done
                          : status == RunStatus::Stopped ? VariantStatus::Stopped
                                                         : VariantStatus::Failed);
    saveStudy(*it);
    // The next variant starts once this run's process has gone.
    QTimer::singleShot(0, this, [this] { startNextVariant(); });
}

// ---- studies ----------------------------------------------------------------------------

const Study* Session::study(const QString& id) const
{
    auto it = std::find_if(m_studies.begin(), m_studies.end(), [&](const Study& s) { return s.id == id; });
    return it == m_studies.end() ? nullptr : &*it;
}

void Session::saveStudy(const Study& s)
{
    if (!m_studyStore.save(s))
        emit consoleLine(ConsoleLevel::Warning, tr("Could not save %1 in %2.").arg(s.title(), QDir::toNativeSeparators(m_history.caseDir())));
    emit studiesChanged();
}

QString Session::addStudy(Study s, QString* error)
{
    if (!m_case) {
        if (error) *error = tr("Open a case first.");
        return {};
    }
    s.number = std::max(m_studyStore.nextNumber(), m_studies.empty() ? 1 : m_studies.back().number + 1);
    s.id = QStringLiteral("study-%1").arg(s.number, 4, 10, QChar('0'));
    s.approximations = m_approved;
    if (!m_studyStore.save(s)) {
        if (error) *error = tr("Kutta cannot write the study in %1.").arg(QDir::toNativeSeparators(m_history.caseDir()));
        return {};
    }
    m_studies.push_back(s);
    emit studiesChanged();
    return s.id;
}

void Session::deleteStudy(const QString& id)
{
    if (id == m_activeStudy) return;
    m_studyStore.remove(id);
    std::erase_if(m_studies, [&](const Study& s) { return s.id == id; });
    emit studiesChanged();
}

bool Session::startStudy(const QString& id, QString* why)
{
    auto no = [&](const QString& r) { if (why) *why = r; return false; };
    auto it = std::find_if(m_studies.begin(), m_studies.end(), [&](const Study& s) { return s.id == id; });
    if (it == m_studies.end()) return no(tr("No such study."));
    if (running()) return no(tr("A calculation is already running."));
    if (!m_case) return no(tr("Open a case first."));
    const auto blocking = m_case->blocking(it->approximations);
    if (!blocking.empty()) return no(tr("This study's runs need approvals the case requires: %1.").arg(blocking.front()->label));
    requeue(*it);
    if (it->nextQueued() < 0) return no(tr("Every variant already ran."));
    m_activeStudy = id;
    saveStudy(*it);
    emit consoleLine(ConsoleLevel::Info, tr("%1: %n variant(s) of %2, one after another", nullptr, int(it->variants.size()))
                                             .arg(it->title(), it->parameterLabel));
    startNextVariant();
    return true;
}

void Session::startNextVariant()
{
    auto it = std::find_if(m_studies.begin(), m_studies.end(), [this](const Study& x) { return x.id == m_activeStudy; });
    if (it == m_studies.end() || running()) return;
    const int next = it->nextQueued();
    if (next < 0) {
        emit consoleLine(ConsoleLevel::Success, tr("%1 finished.").arg(it->title()));
        m_activeStudy.clear();
        emit studiesChanged();
        return;
    }
    const StudyVariant& v = it->variants[std::size_t(next)];
    emit consoleLine(ConsoleLevel::Info, tr("%1: %2 = %3 (%4 of %5)").arg(it->title(), it->parameterLabel, QLocale().toString(v.value, 'g', 8))
                                             .arg(next + 1).arg(it->variants.size()));
    const QString runId = start(it->overridesFor(v), it->runControls, it->approximations, it->id);
    if (runId.isEmpty()) {
        cancelQueued(*it);
        m_activeStudy.clear();
        saveStudy(*it);
        return;
    }
    markStarted(*it, next, runId);
    // The run's label says which variant it is, in Run History and every viewport.
    const QString note = QStringLiteral("%1 · %2").arg(it->title(), QLocale().toString(v.value, 'g', 8));
    m_history.setNote(runId, note);
    if (m_current) m_current->note = note;
    saveStudy(*it);
}

void Session::stopStudy()
{
    auto it = std::find_if(m_studies.begin(), m_studies.end(), [this](const Study& x) { return x.id == m_activeStudy; });
    if (it == m_studies.end()) return;
    cancelQueued(*it);
    saveStudy(*it);
    stop();
    if (!running()) m_activeStudy.clear();
}

// ---- aero -------------------------------------------------------------------------------

QStringList Session::carZones() const
{
    QStringList out;
    if (!m_case) return out;
    for (const auto& r : m_case->reports)
        if (r.kind == "drag" || r.kind == "lift" || r.kind == "downforce" || r.kind == "force")
            for (const QString& z : r.zones)
                if (!out.contains(z)) out << z;
    if (!out.isEmpty()) return out;
    static const QRegularExpression ground(QStringLiteral("ground|road|belt"), QRegularExpression::CaseInsensitiveOption);
    for (const auto& z : m_case->zones)
        if (z.type == "wall" && !ground.match(z.name).hasMatch()) out << z.name;
    return out;
}

void Session::setAeroSpec(const AeroSpec& spec)
{
    m_aeroSpec = spec;
    m_runDirAero.clear();
    QSaveFile f(QDir(m_history.caseDir()).filePath("aero-spec.json"));
    QDir().mkpath(m_history.caseDir());
    m_aeroSpecSaved = f.open(QIODevice::WriteOnly) && f.write(QJsonDocument(spec.toJson()).toJson(QJsonDocument::Indented)) >= 0 && f.commit();
    if (!m_aeroSpecSaved) emit consoleLine(ConsoleLevel::Warning, tr("Could not save the aero settings in %1.").arg(QDir::toNativeSeparators(m_history.caseDir())));
    m_aero.clear();
    emit aeroSpecChanged();
}

std::optional<AeroResult> Session::aeroIfComputed(const QString& runId) const
{
    auto it = m_aero.find(runId);
    return it == m_aero.end() ? std::nullopt : std::optional<AeroResult>(*it);
}

void Session::aero(const QString& runId, std::function<void(Result<AeroResult>)> done)
{
    if (auto it = m_aero.find(runId); it != m_aero.end()) return done(*it);
    const RunRecord* r = m_history.find(runId);
    if (!r) return done(Result<AeroResult>::failure(tr("No such run.")));
    auto& waiting = m_aeroPending[runId];
    waiting.push_back(std::move(done));
    if (waiting.size() > 1) return;
    // The spec the CLI reads is the case's saved one; write it if it never was.
    const QString spec = QDir(m_history.caseDir()).filePath("aero-spec.json");
    if (!QFileInfo::exists(spec)) setAeroSpec(m_aeroSpec);
    m_cli->aero(r->dir, spec, [this, runId](Result<AeroResult> res) {
        if (res) m_aero.insert(runId, res.value());
        const auto callbacks = m_aeroPending.take(runId);
        for (const auto& cb : callbacks) cb(res);
        if (res) emit aeroComputed(runId);
    });
}

void Session::deleteRun(const QString& id)
{
    m_aero.remove(id);
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
    // A 2D case is its own plane: the extrusion's face shows the solution.
    if (!m_case || !m_case->bounds.valid || m_case->planar()) return {};
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
        {
            // Tangled cells the CLI left out of slices and iso-surfaces: drawn without them.
            QStringList skipped;
            for (const auto* parts : {&m.value().slices, &m.value().isosurfaces})
                for (const ViewPart& p : *parts)
                    if (p.skippedCells > 0)
                        skipped << QStringLiteral("%1 (%2)").arg(p.kind == PartKind::Slice ? itemLabelForSlice(p) : p.name, QLocale().toString(p.skippedCells));
            if (!skipped.isEmpty())
                emit consoleLine(ConsoleLevel::Info, tr("Cells kutta could not cut (tangled mesh cells) are left out of %1.").arg(skipped.join(", ")));
        }
        if (const QJsonObject pc = m.value().postCache; pc.contains("reused") && pc.value("reused").toInt() > 0)
            emit consoleLine(ConsoleLevel::Info, tr("Post-processing reused %1 earlier result(s) and computed %2 in %3 s.")
                                                     .arg(pc.value("reused").toInt()).arg(pc.value("computed").toInt())
                                                     .arg(pc.value("seconds").toDouble(), 0, 'f', 2));
        QFile f(QDir(dir).filePath("app-view.json"));
        if (f.open(QIODevice::WriteOnly)) f.write(QJsonDocument(QJsonObject{{"slices", arr}}).toJson());
        m_views.insert(runId, m.value());
        auto scene = kapp::loadScene(m.value(), title);
        if (scene) {
            if (m_case && m_case->planar()) makePlanar(*scene.value());
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
                if (scene) {
                    if (m_case && m_case->planar()) makePlanar(*scene.value());
                    m_scenes.insert(key, scene.value());
                }
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
                if (m_case && m_case->planar()) makePlanar(*scene.value());
                if (m_mesh->hasProject() && m_mesh->project().caseFile == m_casePath) keepGeometryColours(*scene.value(), m_mesh->scene().get());
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
    if (key == "geometry") {
        if (auto g = m_mesh->scene()) {
            m_scenes.insert(key, g);
            return done(key, g);
        }
        return done(key, Result<std::shared_ptr<SceneData>>::failure(tr("no geometry imported")));
    }
    if (key.startsWith("meshfile:")) {
        if (key != generatedMeshSceneKey()) return done(key, Result<std::shared_ptr<SceneData>>::failure(tr("that mesh is no longer the project's")));
        const QString meshFile = m_mesh->meshFile();
        const QString dir = QDir(m_dataRoot).filePath("geometry/" + key.mid(9));
        auto finish = [this, key, done](Result<ViewManifest> m) {
            if (!m) return done(key, Result<std::shared_ptr<SceneData>>::failure(m.error()));
            auto scene = kapp::loadScene(m.value(), tr("Generated mesh"));
            if (scene) {
                scene.value()->geometryOnly = true;
                keepGeometryColours(*scene.value(), m_mesh->scene().get());
                m_scenes.insert(key, scene.value());
            }
            done(key, scene);
        };
        const QString manifest = QDir(dir).filePath("geometry/manifest.json");
        if (QFileInfo::exists(manifest))
            if (auto m = loadViewManifest(manifest); m && m.value().geometryOnly) return finish(m);
        QDir().mkpath(dir);
        emit consoleLine(ConsoleLevel::Info, tr("Reading the generated mesh for display…"));
        m_cli->exportGeometry(meshFile, dir, [this, finish](Result<ViewManifest> m) {
            if (!m) emit consoleLine(ConsoleLevel::Warning, tr("Mesh display unavailable: %1").arg(m.error()));
            finish(m);
        });
        return;
    }
    if (key.startsWith("ext:")) {
        // A run of another case (a design variant), shown beside this case's runs.
        const QString runDir = key.mid(4);
        QFile rec(QDir(runDir).filePath("app-run.json"));
        auto record = rec.open(QIODevice::ReadOnly) ? RunRecord::fromJson(QJsonDocument::fromJson(rec.readAll()).object(), runDir) : std::nullopt;
        if (!record) return done(key, Result<std::shared_ptr<SceneData>>::failure(tr("no run in %1").arg(QDir::toNativeSeparators(runDir))));
        const QString title = QStringLiteral("%1 · %2").arg(QFileInfo(QFileInfo(record->casePath).absolutePath()).fileName(), record->title());
        auto finish = [this, key, title, done](Result<ViewManifest> m) {
            if (!m) return done(key, Result<std::shared_ptr<SceneData>>::failure(m.error()));
            auto scene = kapp::loadScene(m.value(), title);
            if (scene) m_scenes.insert(key, scene.value());
            done(key, scene);
        };
        const QString manifest = QDir(runDir).filePath("view/manifest.json");
        if (QFileInfo::exists(manifest)) return finish(loadViewManifest(manifest));
        m_cli->exportView(runDir, defaultSlices(), {}, finish);
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

// ---- New Case from Geometry --------------------------------------------------------

QString Session::generatedMeshSceneKey() const
{
    if (!m_mesh->hasMesh()) return {};
    const QFileInfo fi(m_mesh->meshFile());
    const QByteArray id = QStringLiteral("%1|%2|%3").arg(fi.absoluteFilePath()).arg(fi.size()).arg(fi.lastModified().toMSecsSinceEpoch()).toUtf8();
    return "meshfile:" + QString::fromLatin1(QCryptographicHash::hash(id, QCryptographicHash::Sha256).toHex().left(16));
}

void Session::setUpCaseFromMesh()
{
    if (!m_mesh->hasMesh()) {
        emit consoleLine(ConsoleLevel::Error, tr("Generate the mesh first."));
        return;
    }
    if (running()) {
        emit consoleLine(ConsoleLevel::Error, tr("Stop the running calculation before setting up another case."));
        return;
    }
    const GeometryProject& p = m_mesh->project();
    const QString target = p.defaultCaseFile();
    if (!p.baseCase.isEmpty()) {
        // A design variant: the base case's settings on this mesh. Zones the
        // new mesh lacks are dropped; new solids become walls in the car's reports.
        QFile in(p.baseCase);
        if (!in.open(QIODevice::ReadOnly)) {
            emit consoleLine(ConsoleLevel::Error, tr("The base case %1 cannot be read; set up the case from the mesh instead.").arg(QDir::toNativeSeparators(p.baseCase)));
            return;
        }
        QJsonObject kase = QJsonDocument::fromJson(in.readAll()).object();
        QStringList zones;
        if (m_mesh->report())
            for (const auto& z : m_mesh->report()->zones) zones << z.name;
        QStringList walls;
        for (const auto& s : p.solids) walls << s.name;
        QJsonObject boundaries = kase.value("boundaries").toObject();
        QStringList baseWalls, dropped, added;
        for (auto it = boundaries.begin(); it != boundaries.end(); ++it)
            if (it.value().toObject().value("type") == "wall" && it.value().toObject().value("motion").toString() != "moving") baseWalls << it.key();
        if (!zones.isEmpty())
            for (const QString& name : boundaries.keys())
                if (!zones.contains(name)) {
                    boundaries.remove(name);
                    dropped << name;
                }
        for (const QString& w : walls)
            if (!boundaries.contains(w)) {
                boundaries[w] = QJsonObject{{"type", "wall"}};
                added << w;
            }
        kase["boundaries"] = boundaries;
        QJsonObject reports = kase.value("reports").toObject();
        for (const QString& name : reports.keys()) {
            QJsonObject r = reports.value(name).toObject();
            QStringList rz;
            for (const QJsonValue& v : r.value("zones").toArray()) rz << v.toString();
            const bool wholeCar = !baseWalls.isEmpty() && std::all_of(baseWalls.begin(), baseWalls.end(), [&](const QString& w) { return rz.contains(w); });
            rz.erase(std::remove_if(rz.begin(), rz.end(), [&](const QString& z) { return dropped.contains(z); }), rz.end());
            if (wholeCar) for (const QString& a : added) if (!rz.contains(a)) rz << a;
            r["zones"] = QJsonArray::fromStringList(rz);
            reports[name] = r;
        }
        kase["reports"] = reports;
        kase["mesh"] = QDir(QFileInfo(target).absolutePath()).relativeFilePath(m_mesh->meshFile());
        QSaveFile out(target);
        if (!out.open(QIODevice::WriteOnly) || out.write(QJsonDocument(kase).toJson(QJsonDocument::Indented)) < 0 || !out.commit()) {
            emit consoleLine(ConsoleLevel::Error, tr("Kutta cannot write %1: %2").arg(QDir::toNativeSeparators(target), out.errorString()));
            return;
        }
        emit consoleLine(ConsoleLevel::Success, tr("Wrote %1 with the settings of %2").arg(QDir::toNativeSeparators(target), QFileInfo(p.baseCase).fileName()));
        if (!added.isEmpty()) emit consoleLine(ConsoleLevel::Info, tr("New walls (stationary, in the whole-car reports): %1").arg(added.join(", ")));
        if (!dropped.isEmpty()) emit consoleLine(ConsoleLevel::Info, tr("Zones not in the new mesh, left out: %1").arg(dropped.join(", ")));
        m_mesh->setCaseFile(target);
        openCase(target);
        return;
    }
    QStringList wheels;
    for (const auto& s : p.solids)
        if (s.wheel) wheels << s.name;
    startNewCase(m_mesh->meshFile(), p.speed, m_mesh->tunnelRoles());
    m_draftTarget = target;
    m_pendingWheels = wheels;
}

std::vector<Session::VariantRun> Session::variantRuns() const
{
    std::vector<VariantRun> out;
    const QString root = QDir(m_dataRoot).filePath("runs");
    for (const auto& v : m_mesh->family()) {
        VariantRun vr{v.name, v.dir, v.caseFile, v.current, v.base, std::nullopt};
        if (!v.caseFile.isEmpty()) {
            std::vector<RunRecord> runs;
            if (m_case && QFileInfo(m_casePath) == QFileInfo(v.caseFile)) {
                runs = m_history.runs();
            } else {
                const QString key = RunHistory::caseKey(v.caseFile);
                for (const QFileInfo& d : QDir(root).entryInfoList({key + "-*", key}, QDir::Dirs | QDir::NoDotAndDotDot)) {
                    RunHistory h;
                    h.setCaseDir(d.absoluteFilePath());
                    for (const auto& r : h.runs())
                        if (QFileInfo(r.casePath) == QFileInfo(v.caseFile)) runs.push_back(r);
                }
            }
            for (const auto& r : runs)
                if ((r.status == RunStatus::Completed || r.status == RunStatus::Stopped) && (!vr.run || r.finished > vr.run->finished)) vr.run = r;
        }
        out.push_back(vr);
    }
    return out;
}

void Session::aeroForRunDir(const QString& runDir, std::function<void(Result<AeroResult>)> done)
{
    if (auto it = m_runDirAero.find(runDir); it != m_runDirAero.end()) return done(*it);
    if (!m_case) return done(Result<AeroResult>::failure(tr("Open a case first.")));
    const QString spec = QDir(m_history.caseDir()).filePath("aero-spec.json");
    if (!QFileInfo::exists(spec)) setAeroSpec(m_aeroSpec);
    m_cli->aero(runDir, spec, [this, runDir, done](Result<AeroResult> r) {
        if (r) m_runDirAero.insert(runDir, r.value());
        done(r);
    });
}

} // namespace kapp
