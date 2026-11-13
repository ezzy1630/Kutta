#pragma once

#include "kutta/core/json_value.hpp"
#include "kutta/mesh/mesh.hpp"
#include <array>
#include <cmath>
#include <filesystem>
#include <functional>
#include <memory>
#include <limits>
#include <span>

namespace kutta::mesher {
struct Box { Vec3d min, max; };
struct RefinementBox { std::string name; Box bounds; double size; };
struct Settings {
  std::string units = "m";
  double upstream = 3, downstream = 6, side = 2, top = 2;
  bool ground = true, symmetry_y = false;
  double base_size = 0.25, min_size = 0.006, max_size = 0.03;
  double curvature_degrees = 20;
  int proximity_cells = 3;
  // Gap (slot, ground clearance) sizing may refine below min_size down to
  // this floor; 0 means min_size/4.
  double proximity_min_size = 0;
  double proximity_floor() const { return proximity_min_size > 0 ? proximity_min_size : min_size / 4; }
  std::vector<RefinementBox> boxes;
  // Default stack: first-cell centroid near y+ 50 at 20 m/s (wall-function
  // range of Kutta's y+-insensitive SST wall treatment).
  int layer_count = 4;
  double first_height = 0.002, growth = 1.3;
  double last_layer() const { return first_height * std::pow(growth, layer_count - 1); }
  std::vector<std::string> layer_zones{"*"};
  Index max_cells = 10000000;
};
Settings parse_settings(std::string_view json);

// radius: smallest curvature radius estimated across smooth (non-feature)
// edges, centroid distance over normal angle; infinite on flat regions.
struct Triangle { std::array<Index, 3> nodes; Index zone = 0, component = 0; std::array<bool,3> feature_edges{}; double radius = std::numeric_limits<double>::infinity(); };
struct Surface {
  std::vector<Vec3d> points;
  std::vector<Triangle> triangles;
  std::vector<std::string> zones;
  Box bounds;
  Json::Object repairs;
};
// Binary detection uses the exact byte count, including headers starting "solid".
// Weld tolerance is 1e-7 of the body diagonal; no unbounded hole filling.
Surface read_stl(std::span<const std::filesystem::path>, const Settings&);
struct Sample { double distance; Index triangle; bool inside; };
// The Surface must outlive its query; queries never mutate the repaired input.
class SurfaceQuery {
 public:
  explicit SurfaceQuery(const Surface&);
  ~SurfaceQuery();
  SurfaceQuery(SurfaceQuery&&) noexcept;
  SurfaceQuery& operator=(SurfaceQuery&&) noexcept;
  Sample sample(Vec3d) const;
  // Distance and nearest facet only; inside is false, no ray classification.
  Sample nearest(Vec3d) const;
  // Closest point on a different connected component, for slot refinement.
  double other_distance(Vec3d, Index component) const;
  // Distance to the first facet hit along unit direction d, excluding facet
  // `skip`; infinity when nothing is hit. `hit` receives the facet.
  double ray(Vec3d p, Vec3d d, Index skip = -1, Index* hit = nullptr) const;
  // True when some facet crosses the box interior (touching does not count).
  bool crosses(const Box&) const;
 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
using Progress = std::function<void(std::string_view stage, double fraction, std::string_view message)>;
struct Result { PolyMesh mesh; Json::Object report; };
Result generate(const Surface&, const Settings&, double speed = 20, const Progress& = {});
// Fluent ASCII, mixed polyhedral cell zone and outward owner winding.
void write_msh(const std::filesystem::path&, const PolyMesh&);
}  // namespace kutta::mesher
