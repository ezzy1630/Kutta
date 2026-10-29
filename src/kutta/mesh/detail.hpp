#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <span>

#include "kutta/mesh/mesh.hpp"

namespace kutta::detail {
inline Vec3d operator+(Vec3d a, Vec3d b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline Vec3d operator-(Vec3d a, Vec3d b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline Vec3d operator*(Vec3d a, double s) {
  return {a.x * s, a.y * s, a.z * s};
}
inline Vec3d operator/(Vec3d a, double s) { return a * (1 / s); }
inline double dot(Vec3d a, Vec3d b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline Vec3d cross(Vec3d a, Vec3d b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(Vec3d a) { return std::sqrt(dot(a, a)); }
inline bool finite(Vec3d a) {
  return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}
inline Index index(std::size_t n) {
  KUTTA_REQUIRE(
      n <= static_cast<std::size_t>(std::numeric_limits<Index>::max()),
      "mesh exceeds 32-bit index capacity");
  return static_cast<Index>(n);
}
inline std::span<const Index> nodes(const PolyMesh& m, Index f) {
  return std::span<const Index>(m.face_nodes)
      .subspan(m.face_offsets[f], m.face_offsets[f + 1] - m.face_offsets[f]);
}
inline Vec3d face_center(const PolyMesh& m, Index f) {
  const auto ns = nodes(m, f);
  const auto origin = m.points[ns[0]];
  Vec3d d{};
  for (auto n : ns) d = d + (m.points[n] - origin);
  return origin + d / static_cast<double>(ns.size());
}
// Records are an internal construction format; pack restores public ordering.
struct Face {
  std::vector<Index> nodes;
  Index owner;
  Index neighbour = -1;
  Index patch = -1;
};
PolyMesh pack(std::vector<Vec3d> points, Index cells, std::vector<Face> faces,
              std::vector<Patch> patches, std::vector<CellZone> zones = {});
using CellFaces = std::vector<std::vector<Index>>;
PolyMesh from_cells(
    std::vector<Vec3d> points, const std::vector<CellFaces>& cells,
    std::vector<Patch> patches,
    const std::function<Index(const std::vector<Index>&)>& boundary_patch);
}  // namespace kutta::detail
