#pragma once

// New Case from Geometry: a project folder holding STL references, the
// meshing settings and, once meshed, the volume mesh and its case.
//
//   <project>/geometry.json        this app's project file (STLs, names, roles, mesh.json)
//   <project>/mesh.json            what kutta mesh is asked for (written on Generate)
//   <project>/geometry/*.stl       copies with renamed solids, only when names changed
//   <project>/mesh/                kutta mesh --out: mesh.msh, mesh-report.json
//   <project>/<project>.kutta.json the case, from kutta new-case or a variant's base case
//
// STL surfaces are read in the app and shown at once; the CLI runs only to mesh.

#include "core/CliClient.h"
#include "core/Console.h"
#include "core/MeshCli.h"
#include "core/MeshSetup.h"
#include "core/SceneData.h"

#include <QObject>
#include <QPointer>

#include <memory>
#include <optional>

namespace kapp {

struct GeometrySolid {
    QString file;    // absolute STL path
    QString source;  // the solid's name in the file
    QString name;    // the wall zone it becomes (editable)
    bool wheel = false;  // rotates as a wheel in the case
};

struct GeometryProject {
    QString dir;
    QString family;     // shared by a design and its variants
    QString variantOf;  // absolute dir of the project this one was duplicated from
    QString baseCase;   // variants: the case whose settings the new case reuses
    QStringList files;  // absolute STL paths
    std::vector<GeometrySolid> solids;
    MeshSetup setup;
    double speed = 20;  // m/s: the y+ estimate's speed, and the case's freestream
    GroundMotion ground = GroundMotion::Moving;  // the case's ground: moving (cars) or stationary (wind-tunnel models)
    QString caseFile;   // absolute; empty until a case is set up
    QString name() const;
    QString projectFile() const;
    QString meshDir() const;
    QString defaultCaseFile() const;
    QJsonObject toJson() const;
    static GeometryProject fromJson(const QJsonObject& o, const QString& dir);
};

class MeshWorkflow : public QObject {
    Q_OBJECT
public:
    MeshWorkflow(CliClient* cli, QObject* parent = nullptr);
    ~MeshWorkflow() override;

    bool hasProject() const { return m_project.has_value(); }
    const GeometryProject& project() const { return *m_project; }
    bool loading() const { return m_loading > 0; }
    // The last problem, shown on the geometry pages; cleared by the next action.
    QString error() const { return m_error; }

    // A new project from STL files, in `projectDir` (default: "<first STL>-kutta" beside it).
    void importGeometry(const QStringList& stls, const QString& projectDir = {});
    void openProject(const QString& pathOrDir);
    void closeProject();
    void addGeometry(const QStringList& stls);
    void replaceGeometry(int fileIndex, const QString& stl);
    void removeGeometry(int fileIndex);
    void setUnits(LengthUnit u);
    void renameSolid(int index, const QString& name);
    void setSolidWheel(int index, bool wheel);
    void setSetup(const MeshSetup& s);
    void resetSetup();
    void setSpeed(double metresPerSecond);
    void setGroundMotion(GroundMotion g);
    // Roles the case's new-case call gives the tunnel's faces where they
    // differ from new-case's suggestion (a stationary ground).
    QMap<QString, QString> tunnelRoles() const;
    // Takes a mesh.json's settings (lengths in metres; keys it omits get
    // kutta mesh's defaults); the project keeps its geometry units.
    bool loadSetupFile(const QString& path);

    // The geometry as read (metres).
    struct FileInfo { QString path; bool binary = false; QString error; std::size_t triangles = 0; double area = 0; Bounds bounds; };
    const std::vector<FileInfo>& files() const { return m_files; }
    Bounds geometryBounds() const;
    double surfaceArea() const;
    std::size_t triangles() const;
    CellEstimate estimate() const;
    // Why Generate is refused; empty when it can run.
    QStringList problems() const;

    // Meshing
    enum class Outcome { None, Succeeded, Failed, Cancelled };
    bool meshing() const { return m_proc != nullptr; }
    bool generate();
    void cancel();
    struct Progress { QString stage; double fraction = 0; QString message; double seconds = 0; };
    const Progress& progress() const { return m_progress; }
    Outcome outcome() const { return m_outcome; }
    QString failure() const { return m_failure; }
    QString advice() const { return m_advice; }
    QStringList log() const { return m_log; }
    bool hasMesh() const;
    QString meshFile() const;
    // The geometry or settings changed since the mesh was generated.
    bool meshStale() const;
    const std::optional<MeshReport>& report() const { return m_report; }
    QString reportError() const { return m_reportError; }
    // The mesh came from the fixture CLI (app development), not kutta mesh.
    bool meshSynthetic() const;

    // Records the case written for this project's mesh.
    void setCaseFile(const QString& path);

    // Design variants: duplicates this project (same STLs, settings and case
    // settings) beside it and opens the copy, ready for a swapped part.
    // `aeroSpec` (optional) is copied for the variant's case. Returns the new
    // project folder, or empty with error().
    QString duplicateAsVariant(const QJsonObject& aeroSpec = {});
    struct Variant { QString dir, name, caseFile; bool current = false; bool base = false; };
    // The family's projects in the parent folder, base first.
    std::vector<Variant> family() const;
    // The project that a case file belongs to, if any (its folder's geometry.json).
    static QString projectForCase(const QString& casePath);

    // The geometry for the viewports: solids, tunnel faces and outlines, and
    // refinement regions. Rebuilt after each change.
    std::shared_ptr<SceneData> scene() const { return m_scene; }

signals:
    void projectChanged();   // files, solids or the project itself
    void setupChanged();     // mesh settings
    void sceneChanged();
    void meshStateChanged(); // started, progress, finished
    void consoleLine(kapp::ConsoleLevel level, const QString& text);

private:
    struct LoadedSolid { QString name; std::size_t triangles = 0; double area = 0; std::shared_ptr<const PreparedMesh> mesh; };
    struct Loaded { QString path; bool binary = false; QString error; std::vector<LoadedSolid> solids; };
    void readAsync(const QStringList& paths, LengthUnit unit, std::function<void(std::vector<Loaded>)> done);
    void adopt(std::vector<Loaded> loaded, bool replaceAll);
    void rebuildScene();
    void save();
    void fail(const QString& message);
    QString inputsHash() const;
    void loadReport();
    QStringList meshInputs(QString* error);
    void finishMesh(int code, bool cancelled);

    CliClient* m_cli;
    std::optional<GeometryProject> m_project;
    std::vector<Loaded> m_loaded;  // parallel to project files
    std::vector<FileInfo> m_files;
    std::shared_ptr<SceneData> m_scene;
    QString m_error;
    int m_loading = 0;
    quint64 m_generation = 0;
    QPointer<MeshProcess> m_proc;
    Progress m_progress;
    Outcome m_outcome = Outcome::None;
    QString m_failure, m_advice;
    QStringList m_log;
    std::optional<MeshReport> m_report;
    QString m_reportError;
    QString m_lastStage;
    bool m_sawDone = false;
    bool m_meshFromFixture = false;
};

} // namespace kapp
