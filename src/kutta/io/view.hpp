#pragma once

#include "kutta/io/mesh_io.hpp"
#include <optional>

namespace kutta::post { struct Request; }

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
// Ordinary surfaces/slices keep storage proportional to visible geometry.
// Optional post requests may construct volume gradients and nodal caches.
void export_view(const std::filesystem::path& directory,const Snapshot&,std::span<const Plane> planes={},const post::Request* request=nullptr);
} // namespace kutta::view
