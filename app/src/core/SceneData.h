#pragma once

// CPU-side scene content: loaded view parts plus the derived arrays the GPU
// needs (normals, scalar channels, edges). Shared read-only between viewports.

#include "core/Colormap.h"
#include "core/Fields.h"
#include "core/ViewFormat.h"

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
    // Unique triangle edges as a line list (lazily built).
    const std::vector<std::uint32_t>& edges() const;
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
};

struct SceneItem {
    QString id;         // stable: "zone:<name>", "slice:<n>", "post:<operation name>"
    QString label;      // shown in lists: Fluent zone name, "z = 0.025 m", ...
    ItemKind kind = ItemKind::Zone;
    QString zoneType;   // wall, velocity-inlet, ... for zones
    std::shared_ptr<const PreparedMesh> mesh;
    bool defaultVisible = true;
};

// One dataset: a run's exported view (+ post-processing), or a derived comparison.
struct SceneData {
    QString label;          // "Run 3 · 16 m/s", "Run 3 − Run 2"
    QString provenance;     // approximation/modified notice from the CLI; empty when faithful
    bool synthetic = false; // fixture data invented for demonstration
    bool reference = false; // Fluent fields present
    bool runDelta = false;  // fields are Run A − Run B
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

} // namespace kapp
