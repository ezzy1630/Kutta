#pragma once

#include "kutta/core/json_value.hpp"
#include "kutta/mesh/mesh.hpp"
#include <array>
#include <filesystem>
#include <functional>
#include <memory>
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
  std::vector<RefinementBox> boxes;
  int layer_count = 0;
  double first_height = 0.00005, growth = 1.3;
  std::vector<std::string> layer_zones{"*"};
  Index max_cells = 10000000;
};
Settings parse_settings(std::string_view json);

struct Triangle { std::array<Index, 3> nodes; Index zone = 0, component = 0; std::array<bool,3> feature_edges{}; };
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
  // Closest point on a different connected component, for slot refinement.
  double other_distance(Vec3d, Index component) const;
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
