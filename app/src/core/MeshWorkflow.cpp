#include "core/MeshWorkflow.h"
#include "core/PathDisplay.h"

#include <QCryptographicHash>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLocale>
#include <QRegularExpression>
#include <QSaveFile>
#include <QSet>
#include <QThread>
#include <QUuid>

#include <cmath>

namespace kapp {

// ---- Project file -------------------------------------------------------------

QString GeometryProject::name() const { return QFileInfo(dir).fileName(); }
QString GeometryProject::projectFile() const { return QDir(dir).filePath("geometry.json"); }
QString GeometryProject::meshDir() const { return QDir(dir).filePath("mesh"); }
QString GeometryProject::defaultCaseFile() const { return QDir(dir).filePath(name() + ".kutta.json"); }

QJsonObject GeometryProject::toJson() const
{
    const QDir d(dir);
    auto rel = [&](const QString& p) { return p.isEmpty() ? QString() : d.relativeFilePath(p); };
    QJsonArray f, s;
    for (const QString& p : files) f.append(QJsonObject{{"path", rel(p)}});
    for (const auto& x : solids)
        s.append(QJsonObject{{"file", rel(x.file)}, {"source", x.source}, {"name", x.name}, {"wheel", x.wheel}});
    // v2: "mesh" lengths in metres, as mesh.json (v1 wrote them in "units").
    return {{"v", 2}, {"kind", "kutta-geometry"}, {"family", family}, {"variant_of", rel(variantOf)}, {"base_case", rel(baseCase)},
            {"files", f}, {"solids", s}, {"mesh", setup.toJson()}, {"speed", speed}, {"case", rel(caseFile)}};
}

GeometryProject GeometryProject::fromJson(const QJsonObject& o, const QString& dir)
{
    GeometryProject p;
    p.dir = QDir(dir).absolutePath();
    const QDir d(p.dir);
    auto abs = [&](const QString& r) { return r.isEmpty() ? QString() : QDir::cleanPath(d.absoluteFilePath(r)); };
    p.family = o.value("family").toString();
    p.variantOf = abs(o.value("variant_of").toString());
    p.baseCase = abs(o.value("base_case").toString());
    for (const QJsonValue& v : o.value("files").toArray()) p.files << abs(v.toObject().value("path").toString());
    for (const QJsonValue& v : o.value("solids").toArray()) {
        const QJsonObject x = v.toObject();
        p.solids.push_back({abs(x.value("file").toString()), x.value("source").toString(), x.value("name").toString(), x.value("wheel").toBool()});
    }
    p.setup = MeshSetup::fromJson(o.value("mesh").toObject(), o.value("v").toInt(1) < 2);
    p.speed = o.value("speed").toDouble(p.speed);
    p.caseFile = abs(o.value("case").toString());
    return p;
}

// ---- Workflow ---------------------------------------------------------------------

MeshWorkflow::MeshWorkflow(CliClient* cli, QObject* parent) : QObject(parent), m_cli(cli) {}

MeshWorkflow::~MeshWorkflow()
{
    if (m_proc) {
        m_proc->disconnect(this);
        delete m_proc.data();
    }
}

void MeshWorkflow::fail(const QString& message)
{
    m_error = message;
    emit consoleLine(ConsoleLevel::Error, message);
    emit projectChanged();
}

void MeshWorkflow::readAsync(const QStringList& paths, LengthUnit unit, std::function<void(std::vector<Loaded>)> done)
{
    ++m_loading;
    const quint64 generation = m_generation;
    // Large CAD exports take seconds to parse: off the UI thread.
    QThread* t = QThread::create([this, paths, unit, generation, done] {
        auto out = std::make_shared<std::vector<Loaded>>();
        for (const QString& p : paths) {
            Loaded l;
            l.path = QFileInfo(p).absoluteFilePath();
            auto r = readStl(p, unit);
            if (!r) l.error = r.error();
            else {
                l.binary = r.value().binary;
                for (auto& s : r.value().solids)
                    l.solids.push_back({s.name, s.triangles, s.area, std::make_shared<const PreparedMesh>(std::move(s.mesh))});
            }
            out->push_back(std::move(l));
        }
        QMetaObject::invokeMethod(this, [this, out, generation, done] {
            --m_loading;
            if (generation != m_generation) return;  // the project closed or changed meanwhile
            done(std::move(*out));
        }, Qt::QueuedConnection);
    });
    connect(t, &QThread::finished, t, &QObject::deleteLater);
    t->start();
    emit projectChanged();
}

void MeshWorkflow::importGeometry(const QStringList& stls, const QString& projectDir)
{
    if (meshing()) return fail(tr("Cancel the mesh being generated before importing new geometry."));
    if (stls.isEmpty()) return;
    ++m_generation;
    m_error.clear();
    // The unit is guessed from the first file's size: CAD exports for FSAE
    // are usually millimetres; the project can switch it afterwards.
    QString dir = projectDir;
    if (dir.isEmpty()) {
        const QFileInfo first(stls.first());
        QString base = QDir(first.absolutePath()).filePath(first.completeBaseName() + "-kutta");
        dir = base;
        for (int n = 2; QFileInfo::exists(QDir(dir).filePath("geometry.json")); ++n) dir = base + "-" + QString::number(n);
    }
    if (!QDir().mkpath(dir)) return fail(tr("Kutta cannot create the project folder %1.").arg(QDir::toNativeSeparators(dir)));
    GeometryProject p;
    p.dir = QDir(dir).absolutePath();
    p.family = QUuid::createUuid().toString(QUuid::WithoutBraces);
    m_project = p;
    m_loaded.clear();
    m_files.clear();
    m_report.reset();
    m_reportError.clear();
    m_outcome = Outcome::None;
    m_log.clear();
    emit consoleLine(ConsoleLevel::Info, tr("Importing %1 into %2").arg(stls.join(", "), QDir::toNativeSeparators(p.dir)));
    readAsync(stls, LengthUnit::Metre, [this](std::vector<Loaded> loaded) {
        Bounds raw;
        for (const auto& l : loaded)
            for (const auto& s : l.solids) raw.expand(s.mesh->bounds());
        const LengthUnit unit = guessUnit(raw);
        if (unit != LengthUnit::Metre) {
            // Re-scale rather than re-read.
            for (auto& l : loaded)
                for (auto& s : l.solids) {
                    MeshData m = s.mesh->data();
                    for (float& x : m.positions) x = float(x * unitScale(unit));
                    s.area *= unitScale(unit) * unitScale(unit);
                    s.mesh = std::make_shared<const PreparedMesh>(std::move(m));
                }
            emit consoleLine(ConsoleLevel::Info, tr("Geometry is larger than 20 m in metres: read as millimetres. Change Units on the Geometry page if that is wrong."));
        }
        m_project->setup = defaultMeshSetup(Bounds{}, unit);
        adopt(std::move(loaded), true);
        m_project->setup = defaultMeshSetup(geometryBounds(), unit);
        save();
        rebuildScene();  // the sized setup brings the tunnel and refinement regions
        emit setupChanged();
    });
}

void MeshWorkflow::adopt(std::vector<Loaded> loaded, bool replaceAll)
{
    if (replaceAll) {
        m_loaded.clear();
        m_project->files.clear();
        m_project->solids.clear();
    }
    QStringList errors;
    for (auto& l : loaded) {
        if (!l.error.isEmpty()) {
            errors << l.error;
            continue;
        }
        const qsizetype existing = m_project->files.indexOf(l.path);
        if (existing >= 0) {
            m_loaded[std::size_t(existing)] = std::move(l);
            continue;
        }
        m_project->files << l.path;
        for (const auto& s : l.solids) {
            // Zone names must be unique across files: suffix repeats.
            QString name = s.name;
            QSet<QString> taken;
            for (const auto& x : m_project->solids) taken.insert(x.name);
            for (int n = 2; taken.contains(name); ++n) name = s.name + "-" + QString::number(n);
            const QString lower = name.toLower();
            const bool wheel = lower.contains("wheel") || lower.contains("tyre") || lower.contains("tire");
            m_project->solids.push_back({l.path, s.name, name, wheel});
        }
        m_loaded.push_back(std::move(l));
    }
    m_files.clear();
    for (const auto& l : m_loaded) {
        FileInfo f{l.path, l.binary, l.error, 0, 0, {}};
        for (const auto& s : l.solids) {
            f.triangles += s.triangles;
            f.area += s.area;
            f.bounds.expand(s.mesh->bounds());
        }
        m_files.push_back(f);
    }
    m_error = errors.join("\n");
    for (const QString& e : errors) emit consoleLine(ConsoleLevel::Error, e);
    if (!m_loaded.empty())
        emit consoleLine(ConsoleLevel::Info, tr("Geometry: %1 solid(s), %2 triangles, %3 m long")
                                                 .arg(m_project->solids.size()).arg(QLocale().toString(qulonglong(triangles())))
                                                 .arg(QLocale().toString(geometryBounds().valid ? double(geometryBounds().max.x() - geometryBounds().min.x()) : 0.0, 'f', 3)));
    save();
    rebuildScene();
    emit projectChanged();
}

void MeshWorkflow::openProject(const QString& pathOrDir)
{
    if (meshing()) return fail(tr("Cancel the mesh being generated before opening another project."));
    const QFileInfo fi(pathOrDir);
    const QString file = fi.isDir() ? QDir(pathOrDir).filePath("geometry.json") : pathOrDir;
    QFile f(file);
    if (!f.open(QIODevice::ReadOnly)) return fail(tr("No geometry project at %1.").arg(QDir::toNativeSeparators(pathOrDir)));
    const QJsonObject o = QJsonDocument::fromJson(f.readAll()).object();
    if (o.value("kind").toString() != "kutta-geometry") return fail(tr("%1 is not a Kutta geometry project.").arg(QDir::toNativeSeparators(file)));
    ++m_generation;
    m_error.clear();
    m_project = GeometryProject::fromJson(o, QFileInfo(file).absolutePath());
    m_loaded.clear();
    m_files.clear();
    m_outcome = Outcome::None;
    m_log.clear();
    loadReport();
    const GeometryProject saved = *m_project;
    readAsync(saved.files, saved.setup.units, [this, saved](std::vector<Loaded> loaded) {
        adopt(std::move(loaded), true);
        // Keep the saved names and roles over the files' own.
        for (auto& s : m_project->solids)
            for (const auto& x : saved.solids)
                if (x.file == s.file && x.source == s.source) {
                    s.name = x.name;
                    s.wheel = x.wheel;
                }
        m_project->setup = saved.setup;
        m_project->caseFile = saved.caseFile;
        save();
        rebuildScene();
        emit projectChanged();
        emit setupChanged();
    });
}

void MeshWorkflow::closeProject()
{
    if (meshing()) cancel();
    ++m_generation;
    m_project.reset();
    m_loaded.clear();
    m_files.clear();
    m_scene.reset();
    m_report.reset();
    m_error.clear();
    emit projectChanged();
    emit sceneChanged();
}

void MeshWorkflow::addGeometry(const QStringList& stls)
{
    if (!m_project) return importGeometry(stls);
    if (meshing()) return fail(tr("Cancel the mesh being generated before changing the geometry."));
    m_error.clear();
    readAsync(stls, m_project->setup.units, [this](std::vector<Loaded> loaded) { adopt(std::move(loaded), false); });
}

void MeshWorkflow::replaceGeometry(int index, const QString& stl)
{
    if (!m_project || index < 0 || index >= m_project->files.size()) return;
    if (meshing()) return fail(tr("Cancel the mesh being generated before changing the geometry."));
    m_error.clear();
    const QString old = m_project->files[index];
    // Names and roles carry over to solids of the same name in the new file:
    // swapping rear-wing-v1.stl for rear-wing-v2.stl keeps the zone "rear-wing".
    std::vector<GeometrySolid> previous;
    for (const auto& s : m_project->solids)
        if (s.file == old) previous.push_back(s);
    readAsync({stl}, m_project->setup.units, [this, old, previous](std::vector<Loaded> loaded) {
        if (loaded.empty() || !loaded.front().error.isEmpty()) {
            m_error = loaded.empty() ? tr("Nothing was read.") : loaded.front().error;
            emit consoleLine(ConsoleLevel::Error, m_error);
            emit projectChanged();
            return;
        }
        Loaded l = std::move(loaded.front());
        const int at = int(m_project->files.indexOf(old));
        if (at < 0) return;
        std::vector<GeometrySolid> kept;
        for (const auto& s : m_project->solids)
            if (s.file != old) kept.push_back(s);
        std::vector<GeometrySolid> added;
        for (std::size_t i = 0; i < l.solids.size(); ++i) {
            const auto& s = l.solids[i];
            GeometrySolid g{l.path, s.name, s.name, false};
            // Same source name, else the same position in the file (a single-solid swap).
            const GeometrySolid* match = nullptr;
            for (const auto& p : previous)
                if (p.source == s.name) match = &p;
            if (!match && previous.size() == l.solids.size()) match = &previous[i];
            if (match) {
                g.name = match->name;
                g.wheel = match->wheel;
            }
            added.push_back(g);
        }
        m_project->files[at] = l.path;
        m_loaded[std::size_t(at)] = std::move(l);
        // Keep the solids in file order.
        std::vector<GeometrySolid> solids;
        for (const QString& f : m_project->files) {
            const auto& from = f == m_project->files[at] ? added : kept;
            for (const auto& s : from)
                if (s.file == f) solids.push_back(s);
        }
        m_project->solids = solids;
        emit consoleLine(ConsoleLevel::Info, tr("Replaced %1 with %2").arg(QFileInfo(old).fileName(), QFileInfo(m_project->files[at]).fileName()));
        adopt({}, false);
    });
}

void MeshWorkflow::removeGeometry(int index)
{
    if (!m_project || index < 0 || index >= m_project->files.size() || meshing()) return;
    const QString file = m_project->files.takeAt(index);
    m_loaded.erase(m_loaded.begin() + index);
    std::erase_if(m_project->solids, [&](const GeometrySolid& s) { return s.file == file; });
    emit consoleLine(ConsoleLevel::Info, tr("Removed %1 from the geometry").arg(QFileInfo(file).fileName()));
    adopt({}, false);
}

void MeshWorkflow::setUnits(LengthUnit u)
{
    if (!m_project || u == m_project->setup.units || meshing()) return;
    const double k = unitScale(u) / unitScale(m_project->setup.units);
    for (auto& l : m_loaded)
        for (auto& s : l.solids) {
            MeshData m = s.mesh->data();
            for (float& x : m.positions) x = float(x * k);
            s.area *= k * k;
            s.mesh = std::make_shared<const PreparedMesh>(std::move(m));
        }
    // Settings follow the geometry's new size.
    m_project->setup = defaultMeshSetup({}, u);
    adopt({}, false);
    m_project->setup = defaultMeshSetup(geometryBounds(), u);
    save();
    rebuildScene();
    emit setupChanged();
}

void MeshWorkflow::renameSolid(int index, const QString& name)
{
    if (!m_project || index < 0 || index >= int(m_project->solids.size())) return;
    const QString n = name.trimmed();
    if (n.isEmpty() || n == m_project->solids[std::size_t(index)].name) return;
    for (const auto& s : m_project->solids)
        if (s.name == n) {
            m_error = tr("Another solid is already called %1; zone names must be unique.").arg(n);
            emit projectChanged();
            return;
        }
    if (n.contains(QRegularExpression("\\s"))) {
        m_error = tr("Zone names cannot contain spaces (Fluent mesh zone names): %1").arg(n);
        emit projectChanged();
        return;
    }
    m_error.clear();
    m_project->solids[std::size_t(index)].name = n;
    save();
    rebuildScene();
    emit projectChanged();
}

void MeshWorkflow::setSolidWheel(int index, bool wheel)
{
    if (!m_project || index < 0 || index >= int(m_project->solids.size())) return;
    m_project->solids[std::size_t(index)].wheel = wheel;
    save();
    emit projectChanged();
}

void MeshWorkflow::setSetup(const MeshSetup& s)
{
    if (!m_project || meshing()) return;
    if (s == m_project->setup) return;
    m_project->setup = s;
    save();
    rebuildScene();
    emit setupChanged();
}

void MeshWorkflow::resetSetup()
{
    if (!m_project) return;
    setSetup(defaultMeshSetup(geometryBounds(), m_project->setup.units));
}

void MeshWorkflow::setSpeed(double v)
{
    if (!m_project || meshing() || !(v > 0) || !std::isfinite(v) || v == m_project->speed) return;
    m_project->speed = v;
    save();
    emit setupChanged();
}

bool MeshWorkflow::loadSetupFile(const QString& path)
{
    if (!m_project || meshing()) return false;
    QFile f(path);
    QJsonParseError pe;
    const QJsonDocument doc = f.open(QIODevice::ReadOnly) ? QJsonDocument::fromJson(f.readAll(), &pe) : QJsonDocument();
    if (!doc.isObject() || doc.object().value("v").toInt() != 1) {
        fail(tr("%1 is not a kutta mesh setup (mesh.json, \"v\": 1).").arg(QDir::toNativeSeparators(path)));
        return false;
    }
    MeshSetup s = MeshSetup::fromJson(doc.object());
    s.units = m_project->setup.units;  // the file's lengths are metres; the unit is the STL's
    m_error.clear();
    emit consoleLine(ConsoleLevel::Info, tr("Mesh settings from %1").arg(QDir::toNativeSeparators(path)));
    setSetup(s);
    return true;
}

Bounds MeshWorkflow::geometryBounds() const
{
    Bounds b;
    for (const auto& f : m_files) b.expand(f.bounds);
    return b;
}

double MeshWorkflow::surfaceArea() const
{
    double a = 0;
    for (const auto& f : m_files) a += f.area;
    return a;
}

std::size_t MeshWorkflow::triangles() const
{
    std::size_t n = 0;
    for (const auto& f : m_files) n += f.triangles;
    return n;
}

CellEstimate MeshWorkflow::estimate() const
{
    if (!m_project) return {};
    return estimateCells(m_project->setup, geometryBounds(), surfaceArea());
}

QStringList MeshWorkflow::problems() const
{
    if (!m_project) return {tr("Import geometry first.")};
    QStringList p;
    if (loading()) p << tr("The geometry is still being read.");
    if (m_project->solids.empty()) p << tr("There is no geometry to mesh.");
    p += m_project->setup.problems();
    for (const auto& s : m_project->solids)
        if (reservedZoneNames().contains(s.name.toLower())) p << tr("The solid %1 has the name of a tunnel boundary; rename it.").arg(s.name);
    const Bounds g = geometryBounds();
    // kutta mesh puts the ground at z = 0 and rejects geometry below it.
    if (g.valid && m_project->setup.tunnel.ground && g.min.z() < -1e-7f * std::max(1e-3f, g.max.x() - g.min.x()))
        p << tr("The ground is the plane z = 0 and the geometry goes %1 m below it: move the car up in CAD so the tyres stand on z = 0, or turn off the ground.")
                 .arg(QLocale().toString(-double(g.min.z()), 'g', 3));
    if (g.valid && m_project->setup.tunnel.symmetry == "y" && !(g.max.y() > 0))
        p << tr("A half car keeps y ≥ 0, but the geometry is entirely below y = 0.");
    return p;
}

void MeshWorkflow::save()
{
    if (!m_project) return;
    QSaveFile f(m_project->projectFile());
    if (!f.open(QIODevice::WriteOnly) || f.write(QJsonDocument(m_project->toJson()).toJson(QJsonDocument::Indented)) < 0 || !f.commit())
        m_error = tr("Kutta cannot save the project in %1: %2").arg(QDir::toNativeSeparators(m_project->dir), f.errorString());
}

void MeshWorkflow::setCaseFile(const QString& path)
{
    if (!m_project) return;
    m_project->caseFile = path;
    save();
    emit projectChanged();
}

// ---- Scene ----------------------------------------------------------------------------

namespace {

MeshData boxOutline(const Bounds& b)
{
    MeshData m;
    m.polylines = true;
    const QVector3D lo = b.min, hi = b.max;
    auto corner = [&](int i) { return QVector3D(i & 1 ? hi.x() : lo.x(), i & 2 ? hi.y() : lo.y(), i & 4 ? hi.z() : lo.z()); };
    // Two closed rectangles and four uprights.
    const std::vector<std::vector<int>> lines{{0, 1, 3, 2, 0}, {4, 5, 7, 6, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    m.indices.push_back(0);
    for (const auto& line : lines) {
        for (int c : line) {
            const QVector3D p = corner(c);
            m.positions.insert(m.positions.end(), {p.x(), p.y(), p.z()});
        }
        m.indices.push_back(std::uint32_t(m.positions.size() / 3));
    }
    return m;
}

MeshData quad(const Bounds& b, int axis, bool atMax)
{
    MeshData m;
    const int u = (axis + 1) % 3, v = (axis + 2) % 3;
    const float w = atMax ? b.max[axis] : b.min[axis];
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i) {
            QVector3D p;
            p[axis] = w;
            p[u] = i ? b.max[u] : b.min[u];
            p[v] = j ? b.max[v] : b.min[v];
            m.positions.insert(m.positions.end(), {p.x(), p.y(), p.z()});
        }
    m.indices = atMax ? std::vector<std::uint32_t>{0, 1, 3, 0, 3, 2} : std::vector<std::uint32_t>{0, 3, 1, 0, 2, 3};
    return m;
}

} // namespace

void MeshWorkflow::rebuildScene()
{
    if (!m_project) {
        m_scene.reset();
        emit sceneChanged();
        return;
    }
    auto s = std::make_shared<SceneData>();
    s->label = tr("Geometry");
    s->geometryOnly = true;
    s->geometrySetup = true;
    int colour = 0;
    for (const auto& l : m_loaded)
        for (const auto& solid : l.solids) {
            QString name = solid.name;
            for (const auto& g : m_project->solids)
                if (g.file == l.path && g.source == solid.name) name = g.name;
            SceneItem it;
            it.id = "zone:" + name;
            it.label = name;
            it.kind = ItemKind::Zone;
            it.zoneType = "wall";
            it.mesh = solid.mesh;
            it.color = zoneColor(colour++);
            s->items.push_back(it);
            s->bounds.expand(solid.mesh->bounds());
        }
    const Bounds g = geometryBounds();
    const MeshSetup& setup = m_project->setup;
    if (g.valid && setup.problems().filter(tr("Tunnel")).isEmpty() && parseAxis(setup.flow).first != parseAxis(setup.up).first) {
        const Bounds t = tunnelBounds(setup, g);
        s->bounds.expand(t);
        QSet<QString> seen;
        for (const auto& f : tunnelFaces(setup)) {
            SceneItem it;
            it.id = QStringLiteral("tunnel:%1:%2:%3").arg(f.name).arg(f.axis).arg(f.atMax);
            it.label = f.name;
            it.kind = ItemKind::Zone;
            it.zoneType = f.role == "moving-ground" ? "wall" : f.role;
            if (f.role == "moving-ground") it.id = "zone:ground";
            it.mesh = std::make_shared<const PreparedMesh>(quad(t, f.axis, f.atMax));
            if (seen.contains(it.id)) continue;
            seen.insert(it.id);
            s->items.push_back(it);
        }
        SceneItem outline;
        outline.id = "outline:tunnel";
        outline.label = tr("Wind tunnel");
        outline.kind = ItemKind::Probe;
        outline.mesh = std::make_shared<const PreparedMesh>(boxOutline(t));
        s->items.push_back(outline);
    }
    for (const auto& b : setup.boxes) {
        if (!(b.min.x() < b.max.x() && b.min.y() < b.max.y() && b.min.z() < b.max.z())) continue;
        SceneItem it;
        it.id = "box:" + b.name;
        it.label = tr("%1 (refinement, %2 m)").arg(b.name, QLocale().toString(b.size, 'g', 3));
        it.kind = ItemKind::Probe;
        Bounds bb;
        bb.min = b.min;
        bb.max = b.max;
        bb.valid = true;
        it.mesh = std::make_shared<const PreparedMesh>(boxOutline(bb));
        s->items.push_back(it);
    }
    m_scene = s;
    emit sceneChanged();
}

// ---- Meshing ------------------------------------------------------------------------------

QString MeshWorkflow::meshFile() const { return m_project ? meshFileIn(m_project->meshDir()) : QString(); }
bool MeshWorkflow::hasMesh() const { return m_project && QFileInfo::exists(meshFile()) && m_report.has_value(); }

QString MeshWorkflow::inputsHash() const
{
    if (!m_project) return {};
    QCryptographicHash h(QCryptographicHash::Sha256);
    h.addData(QJsonDocument(m_project->setup.toJson()).toJson(QJsonDocument::Compact));
    h.addData(QByteArray::number(m_project->speed, 'g', 10));
    for (const QString& f : m_project->files) {
        const QFileInfo fi(f);
        h.addData(QStringLiteral("%1|%2|%3").arg(fi.absoluteFilePath()).arg(fi.size()).arg(fi.lastModified().toMSecsSinceEpoch()).toUtf8());
    }
    for (const auto& s : m_project->solids) h.addData(QStringLiteral("%1|%2|%3").arg(s.file, s.source, s.name).toUtf8());
    return QString::fromLatin1(h.result().toHex());
}

bool MeshWorkflow::meshSynthetic() const
{
    if (!hasMesh()) return false;
    QFile f(QDir(m_project->meshDir()).filePath("app-mesh.json"));
    return f.open(QIODevice::ReadOnly) && QJsonDocument::fromJson(f.readAll()).object().value("fixture").toBool();
}

bool MeshWorkflow::meshStale() const
{
    if (!hasMesh()) return false;
    QFile f(QDir(m_project->meshDir()).filePath("app-mesh.json"));
    if (!f.open(QIODevice::ReadOnly)) return true;
    return QJsonDocument::fromJson(f.readAll()).object().value("inputs").toString() != inputsHash();
}

void MeshWorkflow::loadReport()
{
    m_report.reset();
    m_reportError.clear();
    if (!m_project) return;
    const QString path = meshReportFileIn(m_project->meshDir());
    if (!QFileInfo::exists(path)) return;
    auto r = loadMeshReport(path);
    if (r) m_report = r.value();
    else m_reportError = r.error();
}

QStringList MeshWorkflow::meshInputs(QString* error)
{
    QStringList out;
    for (std::size_t i = 0; i < m_loaded.size(); ++i) {
        const auto& l = m_loaded[i];
        bool renamed = false;
        for (const auto& s : m_project->solids)
            if (s.file == l.path && s.name != s.source) renamed = true;
        if (!renamed) {
            out << l.path;
            continue;
        }
        // The mesher names wall zones after the solids: hand it a copy with the new names.
        std::vector<std::pair<QString, const MeshData*>> solids;
        for (const auto& solid : l.solids) {
            QString name = solid.name;
            for (const auto& s : m_project->solids)
                if (s.file == l.path && s.source == solid.name) name = s.name;
            solids.push_back({name, &solid.mesh->data()});
        }
        const QString dir = QDir(m_project->dir).filePath("geometry");
        QDir().mkpath(dir);
        const QString copy = QDir(dir).filePath(QFileInfo(l.path).completeBaseName() + ".stl");
        QString err;
        if (!writeAsciiStl(copy, solids, m_project->setup.units, &err)) {
            if (error) *error = tr("Kutta cannot write the renamed geometry %1: %2").arg(QDir::toNativeSeparators(copy), err);
            return {};
        }
        out << copy;
    }
    return out;
}

bool MeshWorkflow::generate()
{
    if (!m_project || meshing()) return false;
    const QStringList problems = this->problems();
    if (!problems.isEmpty()) {
        m_error = problems.join("\n");
        emit projectChanged();
        return false;
    }
    m_error.clear();
    QSaveFile setupFile(QDir(m_project->dir).filePath("mesh.json"));
    if (!setupFile.open(QIODevice::WriteOnly) || setupFile.write(QJsonDocument(m_project->setup.toJson()).toJson(QJsonDocument::Indented)) < 0 ||
        !setupFile.commit()) {
        fail(tr("Kutta cannot write %1: %2").arg(QDir::toNativeSeparators(setupFile.fileName()), setupFile.errorString()));
        return false;
    }
    QString err;
    const QStringList inputs = meshInputs(&err);
    if (inputs.isEmpty()) {
        fail(err.isEmpty() ? tr("There is no geometry to mesh.") : err);
        return false;
    }
    QDir().mkpath(m_project->meshDir());
    MeshRequest req{inputs, QDir(m_project->dir).filePath("mesh.json"), m_project->meshDir(), m_project->speed};
    m_meshFromFixture = m_cli->info().fixture;
    m_proc = new MeshProcess(m_cli->info(), req, this);
    m_progress = {};
    m_outcome = Outcome::None;
    m_failure.clear();
    m_advice.clear();
    m_log.clear();
    m_lastStage.clear();
    m_sawDone = false;
    emit m_cli->commandStarted(shellQuote(QStringList{QFileInfo(m_cli->info().program).fileName()} + m_proc->arguments()));
    connect(m_proc, &MeshProcess::meshEvent, this, [this](const MeshEvent& e) {
        if (const auto* p = std::get_if<MeshProgressEvent>(&e)) {
            // The mesher's fraction can step back within a stage (layer
            // retreat rounds); the bar only moves forward.
            m_progress = {p->stage, std::max(m_progress.fraction, p->fraction), p->message, m_proc ? m_proc->elapsedSeconds() : 0};
            if (p->stage != m_lastStage) {
                m_lastStage = p->stage;
                const QString line = tr("%1%2").arg(meshStageLabel(p->stage), p->message.isEmpty() ? QString() : ": " + p->message);
                m_log << line;
                emit consoleLine(ConsoleLevel::Info, line);
            }
        } else if (const auto* d = std::get_if<MeshDoneEvent>(&e)) {
            m_sawDone = true;
            m_progress.fraction = 1;
            Q_UNUSED(d);
        } else if (const auto* x = std::get_if<MeshErrorEvent>(&e)) {
            m_failure = x->message;
            m_advice = x->advice;
            emit consoleLine(ConsoleLevel::Error, x->message + (x->advice.isEmpty() ? QString() : "  " + x->advice));
        }
        emit meshStateChanged();
    });
    connect(m_proc, &MeshProcess::stderrText, this, [this](const QString& t) {
        for (const QString& l : t.split('\n', Qt::SkipEmptyParts)) emit consoleLine(ConsoleLevel::Info, l);
    });
    connect(m_proc, &MeshProcess::protocolError, this, [this](const QString& m) {
        if (m_failure.isEmpty() && m.startsWith("Could not start")) m_failure = m;
        emit consoleLine(ConsoleLevel::Warning, m);
    });
    connect(m_proc, &MeshProcess::finished, this, &MeshWorkflow::finishMesh);
    m_proc->start();
    emit meshStateChanged();
    return true;
}

void MeshWorkflow::cancel()
{
    if (!m_proc) return;
    emit consoleLine(ConsoleLevel::Warning, tr("Cancelling mesh generation…"));
    m_proc->cancel();
    emit meshStateChanged();
}

void MeshWorkflow::finishMesh(int code, bool cancelled)
{
    const double seconds = m_proc ? m_proc->elapsedSeconds() : 0;
    if (m_proc) m_proc->deleteLater();
    m_proc = nullptr;
    if (cancelled) {
        m_outcome = Outcome::Cancelled;
        emit consoleLine(ConsoleLevel::Warning, tr("Mesh generation cancelled after %1 s; the previous mesh, if any, is unchanged.").arg(seconds, 0, 'f', 1));
    } else if (code != 0) {
        m_outcome = Outcome::Failed;
        if (m_failure.isEmpty()) {
            // Exit 2 without an event is a usage error: a CLI from before kutta mesh.
            m_failure = code == 2 ? tr("This kutta CLI cannot generate meshes (no kutta mesh). Update Kutta, or mesh in your own mesher and use File ▸ New Case from Mesh.")
                        : code < 0 ? tr("kutta mesh crashed.") : tr("kutta mesh exited with code %1.").arg(code);
        }
        emit consoleLine(ConsoleLevel::Error, tr("Mesh generation failed: %1").arg(m_failure));
    } else {
        loadReport();
        if (!m_report) {
            m_outcome = Outcome::Failed;
            m_failure = m_reportError.isEmpty() ? tr("kutta mesh finished without writing its report.") : m_reportError;
            emit consoleLine(ConsoleLevel::Error, m_failure);
        } else {
            m_outcome = Outcome::Succeeded;
            QSaveFile f(QDir(m_project->meshDir()).filePath("app-mesh.json"));
            if (f.open(QIODevice::WriteOnly)) {
                f.write(QJsonDocument(QJsonObject{{"inputs", inputsHash()}, {"seconds", seconds}, {"fixture", m_meshFromFixture}}).toJson());
                f.commit();
            }
            emit consoleLine(ConsoleLevel::Success, tr("Mesh: %1 cells, %2 faces in %3 s%4")
                                                        .arg(QLocale().toString(m_report->cells), QLocale().toString(m_report->faces))
                                                        .arg(m_report->wallSeconds.value_or(seconds), 0, 'f', 1)
                                                        .arg(meshSynthetic() ? tr(" (SYNTHETIC fixture mesh)") : QString()));
        }
    }
    emit meshStateChanged();
}

// ---- Variants ---------------------------------------------------------------------

QString MeshWorkflow::projectForCase(const QString& casePath)
{
    if (casePath.isEmpty()) return {};
    const QString file = QDir(QFileInfo(casePath).absolutePath()).filePath("geometry.json");
    return QFileInfo::exists(file) ? file : QString();
}

QString MeshWorkflow::duplicateAsVariant(const QJsonObject& aeroSpec)
{
    if (!m_project) return {};
    if (meshing()) {
        fail(tr("Cancel the mesh being generated before duplicating the project."));
        return {};
    }
    const QDir parent = QFileInfo(m_project->dir).dir();
    QString stem = m_project->name();
    stem.remove(QRegularExpression("-v\\d+$"));
    QString dir;
    for (int n = 2;; ++n) {
        dir = parent.filePath(QStringLiteral("%1-v%2").arg(stem).arg(n));
        if (!QFileInfo::exists(dir)) break;
    }
    if (!QDir().mkpath(dir)) {
        fail(tr("Kutta cannot create %1.").arg(QDir::toNativeSeparators(dir)));
        return {};
    }
    GeometryProject v = *m_project;
    v.dir = QDir(dir).absolutePath();
    v.variantOf = m_project->dir;
    v.baseCase = m_project->caseFile.isEmpty() ? m_project->baseCase : m_project->caseFile;
    v.caseFile.clear();
    QSaveFile f(v.projectFile());
    if (!f.open(QIODevice::WriteOnly) || f.write(QJsonDocument(v.toJson()).toJson(QJsonDocument::Indented)) < 0 || !f.commit()) {
        fail(tr("Kutta cannot write %1.").arg(QDir::toNativeSeparators(v.projectFile())));
        return {};
    }
    if (!aeroSpec.isEmpty()) {
        QSaveFile a(QDir(dir).filePath("aero-spec.json"));
        if (a.open(QIODevice::WriteOnly)) {
            a.write(QJsonDocument(aeroSpec).toJson(QJsonDocument::Indented));
            a.commit();
        }
    }
    emit consoleLine(ConsoleLevel::Success, tr("Design variant %1 created from %2: same geometry, mesh settings and case settings. Replace or add parts, then generate its mesh.")
                                                .arg(v.name(), m_project->name()));
    // The STLs are already in memory: the variant starts with them.
    ++m_generation;
    m_project = v;
    m_report.reset();
    m_reportError.clear();
    m_outcome = Outcome::None;
    m_log.clear();
    save();
    rebuildScene();
    emit projectChanged();
    emit setupChanged();
    return dir;
}

std::vector<MeshWorkflow::Variant> MeshWorkflow::family() const
{
    std::vector<Variant> out;
    if (!m_project || m_project->family.isEmpty()) return out;
    const QDir parent = QFileInfo(m_project->dir).dir();
    for (const QFileInfo& d : parent.entryInfoList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name)) {
        QFile f(QDir(d.absoluteFilePath()).filePath("geometry.json"));
        if (!f.open(QIODevice::ReadOnly)) continue;
        const QJsonObject o = QJsonDocument::fromJson(f.readAll()).object();
        if (o.value("family").toString() != m_project->family) continue;
        const GeometryProject p = GeometryProject::fromJson(o, d.absoluteFilePath());
        out.push_back({p.dir, p.name(), p.caseFile, p.dir == m_project->dir, p.variantOf.isEmpty()});
    }
    std::stable_sort(out.begin(), out.end(), [](const Variant& a, const Variant& b) { return a.base > b.base; });
    return out;
}

} // namespace kapp
