#pragma once

// CPU-side scene content: loaded view parts plus the derived arrays the GPU
// needs (normals, scalar channels, edges). Shared read-only between viewports.

#include "core/Colormap.h"
#include "core/Fields.h"
#include "core/ViewFormat.h"

#include <QColor>
#include <QJsonObject>
#include <QString>

#include <memory>
#include <mutex>
#include <unordered_map>

namespace kapp {

enum class ItemKind { Zone, Slice, IsoSurface, Streamlines, Probe };

// Value written for missing samples; shaders treat anything above 1e37 as missing.
constexpr float kMissingSentinel = 3.0e38f;

class PreparedMesh {
public:
    explicit PreparedMesh(MeshData data, bool flatNormal = false, QVector3D planeNormal = {});
    const MeshData& data() const { return m_data; }
    const std::vector<float>& positionsNormals() const { return m_pn; }  // 6 floats per vertex
    const Bounds& bounds() const { return m_bounds; }
    std::size_t triangleCount() const { return m_data.polylines ? 0 : m_data.indices.size() / 3; }

    // Scalar channel for a field (component, magnitude or plain scalar) with
    // missing samples replaced by the sentinel. Cached; thread-safe.
    std::shared_ptr<const std::vector<float>> scalar(const FieldRef& f) const;
    std::optional<std::pair<double, double>> scalarRange(const FieldRef& f) const;
    bool hasField(const FieldRef& f) const;
    QStringList fieldNames() const;
    // Mesh edges as a line list: the zone's polygon edges when the view
    // carries them, else unique triangle edges (lazily built).
    const std::vector<std::uint32_t>& edges() const;
    void setFaceEdges(std::vector<std::uint32_t> lines);
    bool hasFaceEdges() const { return m_faceEdges; }
    // Interpolated value at a position on a triangle near `p` (CPU fallback/readout).
    std::optional<double> sampleNear(const QVector3D& p, const FieldRef& f) const;

private:
    MeshData m_data;
    std::vector<float> m_pn;
    Bounds m_bounds;
    mutable std::mutex m_mutex;
    mutable std::unordered_map<std::string, std::shared_ptr<const std::vector<float>>> m_scalars;
    mutable std::vector<std::uint32_t> m_edges;
    mutable bool m_edgesBuilt = false;
    bool m_faceEdges = false;
};

struct SceneItem {
    QString id;         // stable: "zone:<name>", "slice:<n>", "post:<operation name>"
    QString label;      // shown in lists: Fluent zone name, "z = 0.025 m", ...
    ItemKind kind = ItemKind::Zone;
    QString zoneType;   // wall, velocity-inlet, ... for zones
    std::shared_ptr<const PreparedMesh> mesh;
    bool defaultVisible = true;
    QColor color;       // zones: a distinct colour for "colour by zone" displays
};

// One dataset: a run's exported view (+ post-processing), or a derived comparison.
struct SceneData {
    QString label;          // "Run 3 · 16 m/s", "Run 3 − Run 2"
    QString provenance;     // approximation/modified notice from the CLI; empty when faithful
    bool synthetic = false; // fixture data invented for demonstration
    bool reference = false; // Fluent fields present
    bool runDelta = false;  // fields are Run A − Run B
    bool geometryOnly = false;  // a case's boundary mesh, before any run
    // A 2D case, solved as one cell between two symmetry caps: shown as the
    // plane (one cap), seen face on, panned and zoomed but not orbited.
    bool planar = false;
    Bounds bounds;
    std::vector<SceneItem> items;

    const SceneItem* item(const QString& id) const;
    // Fields offered for colouring, Fluent order, per role available.
    std::vector<FieldRef> availableFields() const;
    std::optional<std::pair<double, double>> range(const FieldRef& f, const QStringList& visibleIds) const;
    QString unitOf(const QString& quantity) const;
};

// Load a view manifest's zones, slices, iso-surfaces, traces and probes into a scene.
Result<std::shared_ptr<SceneData>> loadScene(const ViewManifest& m, const QString& label);
// Run A − Run B on parts with identical topology; parts that differ are reported, not guessed.
std::shared_ptr<SceneData> runDifference(const SceneData& a, const SceneData& b, QStringList* skipped);

QString itemLabelForSlice(const ViewPart& s);
// Distinct, colour-vision-friendly zone colours, stable by zone order.
QColor zoneColor(int index);
// Zones a mesh view shows first: bodies and walls, not the planar faces of an
// enclosing box (ground, tunnel sides, far-field inlet/outlet) nor symmetry
// planes that would hide the inside of the domain.
QStringList defaultMeshItems(const SceneData& s);
// A face of the enclosing box (a ground plane, tunnel sides, far-field
// planes): flat in one axis and spanning the domain in the other two.
bool enclosesDomain(const SceneItem& it, const Bounds& domain);
// External aero: walls that enclose the domain (or are named ground) beside
// walls that do not, the body. Pipes and ducts are internal.
bool isExternalFlow(const SceneData& s);
bool isEnclosingWall(const SceneItem& it, const SceneData& s);
// Marks a 2D case's scene: its plane is the +z extrusion cap, labelled as the
// domain; the opposite cap (the same cells, facing away) is dropped.
void makePlanar(SceneData& s);
// A 2D case's plane: the extrusion cap facing +z (the camera's side in the
// Top view); empty when the scene is not planar or has no such cap.
QString planarFace(const SceneData& s);
// What a results view shows first, and the field it is coloured by.
QStringList defaultVisibleItems(const SceneData& s);
FieldRef defaultField(const SceneData& s);

} // namespace kapp
