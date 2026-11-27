#include "core/MeshSensitivity.h"

#include "core/CliClient.h"
#include "core/MeshSetup.h"

#include <QDir>
#include <QFile>
#include <QJsonDocument>
#include <QSaveFile>
#include <QTemporaryDir>

#include <cmath>

namespace kapp {
namespace {
QJsonObject readObject(const QString& path)
{
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) return {};
    return QJsonDocument::fromJson(file.readAll()).object();
}

bool writeObject(const QString& path, const QJsonObject& object)
{
    QSaveFile file(path);
    const auto bytes = QJsonDocument(object).toJson(QJsonDocument::Indented);
    return file.open(QIODevice::WriteOnly) && file.write(bytes) == bytes.size() && file.commit();
}

QString root(const QString& project)
{
    return QDir(project).filePath("mesh-sensitivity");
}

QJsonObject gridSetup(MeshSetup setup, int grid)
{
    const double scale = std::ldexp(1., -grid);
    setup.baseSize *= scale;
    setup.minSize *= scale;
    setup.maxSize *= scale;
    setup.proximityMinSize *= scale;
    setup.firstHeight *= scale;
    for (auto& box : setup.boxes) box.size *= scale;
    return setup.toJson();
}

QStringList zoneSignature(const MeshReport& report)
{
    QStringList names;
    for (const auto& zone : report.zones)
        if (zone.faces > 0) names << zone.type + "/" + zone.name;
    names.sort();
    return names;
}

QJsonValue numberOrNull(const std::optional<double>& value)
{
    return value ? QJsonValue(*value) : QJsonValue(QJsonValue::Null);
}
} // namespace

MeshSensitivity::MeshSensitivity(CliClient* cli, QObject* parent) : QObject(parent), m_cli(cli) {}

QJsonObject MeshSensitivity::latest(const QString& projectDir) const
{
    if (projectDir == m_project && !m_record.isEmpty()) return m_record;
    for (const auto& dir : QDir(root(projectDir)).entryInfoList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::Time)) {
        auto record = readObject(QDir(dir.absoluteFilePath()).filePath("study.json"));
        if (record.value("v").toInt() != 1) continue;
        if (record.value("state") == QStringLiteral("running")) {
            record.insert("state", "interrupted");
            record.insert("message", tr("The study ended before all grids finished."));
        }
        return record;
    }
    return {};
}

bool MeshSensitivity::start(const QString& projectDir, const QStringList& inputs, const MeshSetup& mesh,
                            double speed, const QString& casePath, const QJsonObject& overrides,
                            const QJsonObject& aero, const QString& backend, int iterations, const QStringList& approximations,
                            QString* error)
{
    const auto refuse = [error](const QString& message) {
        if (error) *error = message;
        return false;
    };
    if (m_running) return refuse(tr("A mesh sensitivity study is already running."));
    if (inputs.isEmpty()) return refuse(tr("No geometry is available."));
    auto caseJson = readObject(casePath);
    if (caseJson.value("kutta_case").toInt() != 1) return refuse(tr("Could not read the project's Kutta case."));
    if (!QDir().mkpath(root(projectDir))) return refuse(tr("Could not create the mesh sensitivity folder."));
    QTemporaryDir stage(QDir(root(projectDir)).filePath("study-XXXXXX"));
    if (!stage.isValid()) return refuse(tr("Could not create a study directory."));
    const QDir dir(stage.path());
    QStringList copies;
    for (int i = 0; i < inputs.size(); ++i) {
        const QString folder = QStringLiteral("geometry-%1").arg(i + 1);
        const QString copy = dir.filePath(folder + "/" + QFileInfo(inputs[i]).fileName());
        if (!dir.mkdir(folder) || !QFile::copy(inputs[i], copy)) return refuse(tr("Could not copy %1 into the study.").arg(inputs[i]));
        copies << copy;
    }
    if (!writeObject(dir.filePath("setup.json"), overrides) || !writeObject(dir.filePath("aero-spec.json"), aero)
        || !writeObject(dir.filePath("case-template.json"), caseJson))
        return refuse(tr("Could not save the study inputs."));
    for (int i = 0; i < 3; ++i) {
        const QString grid = QStringLiteral("grid-%1").arg(i + 1);
        if (!dir.mkdir(grid) || !writeObject(dir.filePath(grid + "/mesh.json"), gridSetup(mesh, i)))
            return refuse(tr("Could not save the grid settings."));
    }

    stage.setAutoRemove(false);
    m_dir = dir.path();
    m_project = projectDir;
    m_case = std::move(caseJson);
    m_inputs = std::move(copies);
    m_backend = backend;
    m_iterations = iterations;
    m_approximations = approximations;
    m_speed = speed;
    m_grid = 0;
    m_bounds = {};
    m_rows = {};
    m_error.clear();
    m_zones.clear();
    m_running = true;
    m_stopping = false;
    m_record = {{"v", 1}, {"directory", m_dir}, {"case_file", casePath}};
    nextGrid();
    return true;
}

void MeshSensitivity::nextGrid()
{
    m_gridDir = QDir(m_dir).filePath(QStringLiteral("grid-%1").arg(m_grid + 1));
    const QString setupPath = QDir(m_gridDir).filePath("mesh.json");
    if (m_grid > 0) {
        auto setup = readObject(setupPath);
        if (setup.isEmpty()) return finish("failed", tr("Could not read the grid settings."));
        auto tunnel = setup.value("tunnel").toObject();
        tunnel.insert("fixed_bounds", m_bounds);
        setup.insert("tunnel", tunnel);
        if (!writeObject(setupPath, setup)) return finish("failed", tr("Could not fix the tunnel bounds."));
    }
    m_error.clear();
    m_meshed = false;
    m_mesh = new MeshProcess(m_cli->info(), {m_inputs, setupPath, m_gridDir, m_speed}, this);
    connect(m_mesh, &MeshProcess::finished, this, &MeshSensitivity::meshFinished);
    connect(m_mesh, &MeshProcess::meshEvent, this, [this](const MeshEvent& e) {
        if (std::holds_alternative<MeshDoneEvent>(e)) m_meshed = true;
        if (const auto* error = std::get_if<MeshErrorEvent>(&e)) m_error = error->message;
    });
    connect(m_mesh, &MeshProcess::stderrText, this, [this](const QString& text) { emit m_cli->commandOutput(text); });
    emit m_cli->commandStarted(shellQuote(QStringList{m_cli->info().program} + m_mesh->arguments()));
    update("running", tr("Generating grid %1 of 3…").arg(m_grid + 1));
    m_mesh->start();
}

void MeshSensitivity::meshFinished(int exitCode, bool cancelled)
{
    m_mesh->deleteLater();
    m_mesh = nullptr;
    if (cancelled || m_stopping) return finish("stopped", tr("Study stopped."));
    if (exitCode != 0 || !m_meshed)
        return finish("failed", tr("Grid %1 generation failed. %2").arg(m_grid + 1).arg(m_error));
    const auto report = loadMeshReport(meshReportFileIn(m_gridDir));
    if (!report) return finish("failed", report.error());
    const auto bounds = report->raw.value("bounds").toArray();
    const auto zones = zoneSignature(report.value());
    if (bounds.size() != 2 || zones.isEmpty())
        return finish("failed", tr("Grid %1 has incomplete tunnel or zone metadata.").arg(m_grid + 1));
    if (report->missingSurfaceZones && !report->missingSurfaceZones->isEmpty())
        return finish("failed", tr("Grid %1 is missing wall zones: %2.").arg(m_grid + 1).arg(report->missingSurfaceZones->join(", ")));
    if (m_grid == 0) {
        m_bounds = bounds;
        m_zones = zones;
    } else if (bounds != m_bounds || zones != m_zones) {
        return finish("failed", tr("Grid %1 changed the tunnel or its boundary zones.").arg(m_grid + 1));
    }
    if (m_grid > 0 && report->cells <= m_rows.last().toObject().value("cells").toInteger())
        return finish("failed", tr("Grid %1 did not increase the cell count.").arg(m_grid + 1));

    QJsonObject gridCase = m_case;
    gridCase.insert("mesh", QFileInfo(meshFileIn(m_gridDir)).fileName());
    const auto casePath = QDir(m_gridDir).filePath("grid.kutta.json");
    if (!writeObject(casePath, gridCase)) return finish("failed", tr("Could not write the grid's case."));
    SolveRequest req;
    req.casePath = casePath;
    req.backend = m_backend;
    req.iterations = m_iterations;
    req.approximations = m_approximations;
    req.outDir = m_gridDir;
    req.setupFile = QDir(m_dir).filePath("setup.json");
    m_solve = m_cli->solve(req, this);
    connect(m_solve, &SolveProcess::finished, this, [this](int code, bool crashed) { solveFinished(crashed ? -1 : code); });
    connect(m_solve, &SolveProcess::stderrText, this, [this](const QString& text) { emit m_cli->commandOutput(text); });
    connect(m_solve, &SolveProcess::progress, this, [this](const ProgressEvent& event) {
        if (const auto* error = std::get_if<ErrorEvent>(&event)) m_error = error->message;
    });
    QJsonObject row{{"cells", report->cells}, {"state", "solving"}};
    if (report->gaps) {
        const auto unresolved = std::count_if(report->gaps->begin(), report->gaps->end(),
                                            [](const MeshGap& gap) { return !gap.contact && gap.cellsAcross < 3; });
        row.insert("under_resolved_gaps", qint64(unresolved));
    }
    m_rows.append(row);
    update("running", tr("Solving grid %1 of 3…").arg(m_grid + 1));
    m_solve->start();
}

void MeshSensitivity::solveFinished(int exitCode)
{
    m_solve->deleteLater();
    m_solve = nullptr;
    if (m_stopping) return finish("stopped", tr("Study stopped."));
    if (exitCode != 0)
        return finish("failed", tr("Grid %1 solve failed. %2").arg(m_grid + 1).arg(m_error));
    const auto summary = parseRunSummary(readObject(QDir(m_gridDir).filePath("summary.json")));
    if (!summary || summary->verdict == Verdict::Unknown)
        return finish("failed", tr("Grid %1 did not produce a run summary.").arg(m_grid + 1));
    if (summary->verdict == Verdict::Stopped) return finish("stopped", tr("Study stopped."));
    const auto row = m_rows.last().toObject();
    if (summary->cells != row.value("cells").toInteger())
        return finish("failed", tr("Grid %1 solve has the wrong cell count.").arg(m_grid + 1));
    auto measured = row;
    measured.insert("state", summary->stopReason);
    measured.insert("verdict", verdictKey(summary->verdict));
    measured.insert("iterations", summary->iterations);
    measured.insert("approximate", summary->approximate);
    measured.insert("device", summary->backend + "/" + summary->device);
    if (m_grid > 0 && m_rows.first().toObject().value("device") != measured.value("device"))
        return finish("failed", tr("The solver selected different devices across grids."));
    update("running", tr("Reading loads for grid %1 of 3…").arg(m_grid + 1));

    QPointer<MeshSensitivity> guard(this);
    m_cli->aero(m_gridDir, QDir(m_dir).filePath("aero-spec.json"), [this, guard, measured](Result<AeroResult> result) mutable {
        if (!guard) return;
        if (m_stopping) return finish("stopped", tr("Study stopped."));
        if (!result) return finish("failed", result.error());
        const auto* total = result->total();
        if (!total) return finish("failed", tr("No aerodynamic loads were available for grid %1.").arg(m_grid + 1));
        const QJsonArray zones = QJsonArray::fromStringList(total->zones);
        if (m_grid > 0 && m_rows.first().toObject().value("load_zones") != zones)
            return finish("failed", tr("The load surfaces changed between grids."));
        measured.insert("load_zones", zones);
        measured.insert("drag_n", total->total.drag);
        measured.insert("downforce_n", total->downforce);
        measured.insert("front_percent", numberOrNull(result->balance.frontPercent));
        measured.insert("cd", numberOrNull(total->coefficients ? std::optional(total->coefficients->cd) : std::nullopt));
        measured.insert("cl", numberOrNull(total->coefficients ? std::optional(total->coefficients->cl) : std::nullopt));
        m_rows.replace(m_grid, measured);
        if (++m_grid == 3) return finish("complete", tr("All three grids finished."));
        nextGrid();
    });
}

void MeshSensitivity::stop()
{
    if (!m_running || m_stopping) return;
    m_stopping = true;
    if (m_mesh) m_mesh->cancel();
    if (m_solve) m_solve->stop();
    update("running", tr("Stopping…"));
}

void MeshSensitivity::update(const QString& state, const QString& message)
{
    m_record.insert("state", state);
    m_record.insert("message", message);
    m_record.insert("grids", m_rows);
    if (!writeObject(QDir(m_dir).filePath("study.json"), m_record))
        m_record.insert("message", message + tr(" Could not save the study record."));
    emit changed();
}

void MeshSensitivity::finish(const QString& state, const QString& message)
{
    m_running = false;
    update(state, message);
}

} // namespace kapp
