#pragma once

// Kutta view data (src/kutta/io/VIEW_FORMAT.md, version 1): KUTV triangles,
// KUTL polylines and the manifest, including export-view --post artifacts.
// Reading is the trust boundary: every count, offset and size is checked
// before any array is touched.

#include "core/Result.h"

#include <QJsonObject>
#include <QString>
#include <QVector3D>

#include <array>
#include <cstdint>
#include <optional>
#include <vector>

namespace kapp {

struct Bounds {
    QVector3D min{0, 0, 0}, max{0, 0, 0};
    bool valid = false;
    QVector3D center() const { return (min + max) * 0.5f; }
    float radius() const { return (max - min).length() * 0.5f; }
    void expand(const QVector3D& p);
    void expand(const Bounds& b);
};

struct FieldDesc {
    QString name;      // e.g. "pressure", "pressure@fluent", "pressure@delta"
    QString unit;      // e.g. "Pa", "1"
    int components = 1;
    qint64 offsetBytes = -1;  // per-file descriptors only
    std::optional<double> min, max;  // null when every sample is missing
};

enum class PartKind { Zone, Slice, IsoSurface, Streamlines, Probe };

inline bool isPolyline(PartKind k) { return k == PartKind::Streamlines || k == PartKind::Probe; }

// One binary file listed in a manifest.
struct ViewPart {
    PartKind kind = PartKind::Zone;
    int id = -1;          // Fluent zone id for zones
    QString name;         // zone name, operation name; slices get a generated label
    QString type;         // Fluent zone type for zones; "line"/"section" for probes
    QString file;
    quint32 vertices = 0, indices = 0, lineCount = 0;  // lineCount: KUTL only
    QVector3D origin, normal;           // slices
    std::vector<FieldDesc> fields;
    QJsonObject meta;     // the whole descriptor (iso value, trace lines, probe csv, ...)
    const FieldDesc* field(const QString& name) const;
};

// A field export-view can sample (manifest available_fields).
struct AvailableField {
    QString name, unit;
    int components = 1;
    bool volume = false;  // false: boundary only (no iso-surfaces or volume probes)
};

struct ViewManifest {
    QString dir;           // directory holding the manifest
    bool reference = false;
    QString provenance;    // approximation/modified notice; empty when faithful
    bool synthetic = false;
    Bounds bounds;
    std::vector<ViewPart> zones, slices;
    std::vector<ViewPart> streamlines, isosurfaces, probes;  // export-view --post
    std::vector<AvailableField> availableFields;
    QString reportsFile;            // reports.json beside the manifest, when integrals ran
    QJsonObject referenceValues;
    std::vector<FieldDesc> fields;  // combined ranges across files
    const FieldDesc* field(const QString& name) const;
};

Result<ViewManifest> parseViewManifest(const QJsonObject& json, const QString& dir);
Result<ViewManifest> loadViewManifest(const QString& manifestPath);

// Triangles (zones, slices, iso-surfaces) or polylines (traces, probes).
struct MeshData {
    std::vector<float> positions;       // xyz interleaved
    // Triangle list, or for polylines L+1 offsets into the (sequential) vertices.
    // KUTL index indirection is resolved on read so renderers see plain strips.
    std::vector<std::uint32_t> indices;
    struct Field { FieldDesc desc; std::vector<float> values; };  // V * components
    std::vector<Field> fields;
    bool polylines = false;

    std::size_t vertexCount() const { return positions.size() / 3; }
    const Field* field(const QString& name) const;
    Bounds bounds() const;
};

// The magic decides triangles (KUTV) versus polylines (KUTL).
Result<MeshData> readMeshFile(const QString& path, const ViewPart& part);
Result<MeshData> readMeshBytes(const QByteArray& bytes, const ViewPart& part);

// Writers for the fixture CLI and tests. Field offsets in `part` are filled in.
QByteArray writeMeshBytes(const MeshData& mesh, ViewPart* part);

QJsonObject toJson(const FieldDesc& f, bool withOffset);
QJsonObject toJson(const ViewPart& p);

// Range of a field's samples (magnitude for vectors), ignoring non-finite.
std::optional<std::pair<double, double>> finiteRange(const MeshData::Field& f);

} // namespace kapp
