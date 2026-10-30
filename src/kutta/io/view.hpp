#pragma once

#include "kutta/core/json_value.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/topology.hpp"
#include <optional>

namespace kutta::view {
// Snapshot fields use solver mesh order. Boundary fields index f-n_internal.
// NaN means unavailable reference coverage; it is never turned into zero.
struct Field {
  std::string name, unit;
  std::uint32_t components=1;
  std::vector<float> values;
};
struct Reference {
  std::optional<double> density,viscosity,reference_density,velocity,pressure,area,length;
};
struct Fields {
  std::vector<Field> cells, boundary;
  bool reference=false;
  std::string provenance;
  Reference properties;
};
struct Snapshot : Fields { PolyMesh mesh; };
struct Plane { Vec3d origin, normal; };
Plane parse_plane(std::string_view);
// Private, versioned run snapshot. View files are documented in VIEW_FORMAT.md.
void write_snapshot(const std::filesystem::path&,const PolyMesh&,const Fields&);
inline void write_snapshot(const std::filesystem::path& path,const Snapshot& s) {write_snapshot(path,s.mesh,s);}
Snapshot read_snapshot(const std::filesystem::path&);
// Write a surface-only view into a new, caller-owned staging directory.
// Borrows the mesh and never constructs a volume snapshot or copies topology.
// Live views omit face edges; exports and case geometry include them.
// Returns the written manifest.
Json::Object write_surface_view(const std::filesystem::path&,const PolyMesh&,const Fields&,bool edges=false,Json::Object extra={});
// Ordinary surfaces/slices keep storage proportional to visible geometry.
// Post-processing publishes through post::export_view instead.
void export_view(const std::filesystem::path& directory,const Snapshot&,std::span<const Plane> planes={});

// Building blocks shared by export_view, live view, case geometry and the
// incremental post exporter. `local` is node-sized scratch, all -1 on entry
// and on return.
Json::Object manifest_header(bool reference,const std::string& provenance,std::span<const Vec3d> points);
// Writes zone-N.surface.bin (and zone-N.edges.bin) into dir; returns the
// manifest zone entry.
Json::Object write_zone(const std::filesystem::path& dir,std::size_t ordinal,const PolyMesh&,const Fields&,std::vector<Index>& local,bool edges);
double slice_tolerance(const PolyMesh&);
// Writes one slice file; returns its descriptor including origin/normal.
Json::Object write_slice(const std::filesystem::path& file,const Snapshot&,const Plane&,const kutta::CellFaces&,double eps);
// Adds the top-level "fields" catalog, combined from every geometry
// descriptor, and writes the closed manifest.
void write_manifest(const std::filesystem::path&,Json::Object manifest);
// A fresh staging directory next to the published DIR/view.
std::filesystem::path staging(const std::filesystem::path& directory);
// Replace DIR/view by stage. Readers see the old or new complete view.
void publish(const std::filesystem::path& stage,const std::filesystem::path& directory);
} // namespace kutta::view
