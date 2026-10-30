#include "kutta/mesh/detail.hpp"
#include "view_geometry.hpp"
#include <algorithm>
#include <numeric>

namespace kutta::view {
using namespace detail;
namespace {
Vec3d polygon_area(std::span<const Vec3d> p) {
  Vec3d area;
  for (std::size_t i = 1; i + 1 < p.size(); ++i)
    area = area + cross(p[i] - p[0], p[i + 1] - p[0]) * .5;
  return area;
}
} // namespace
// Ear clipping in the dominant projection preserves concave Fluent polygons.
// Convex faces take a linear fan path. Warped faces use the same projection;
// triangles retain original 3D vertices and outward winding.
std::vector<std::uint32_t> triangulate(std::span<const Vec3d> p) {
  KUTTA_REQUIRE(p.size() >= 3, "view polygon has fewer than three vertices");
  const auto area = polygon_area(p);
  KUTTA_REQUIRE(norm(area) > 0, "view polygon has zero area");
  int axis = 0;
  if (std::abs(area.y) > std::abs(area.x))
    axis = 1;
  if (std::abs(area.z) > (axis == 0 ? std::abs(area.x) : std::abs(area.y)))
    axis = 2;
  const double direction = (axis == 0   ? area.x
                            : axis == 1 ? area.y
                                        : area.z) > 0
                               ? 1
                               : -1;
  const auto turn = [&](Vec3d a, Vec3d b, Vec3d c) {
    const auto v = cross(b - a, c - a);
    return direction * (axis == 0 ? v.x : axis == 1 ? v.y : v.z);
  };
  double scale = 0;
  for (auto v : p)
    scale = std::max(scale, norm(v - p[0]));
  const double eps = scale * scale * 1e-13;
  bool convex = true;
  for (std::size_t i = 0; i < p.size(); ++i)
    convex &= turn(p[i], p[(i + 1) % p.size()], p[(i + 2) % p.size()]) > eps;
  std::vector<std::uint32_t> out;
  out.reserve((p.size() - 2) * 3);
  if (convex) {
    for (std::uint32_t i = 1; i + 1 < p.size(); ++i)
      out.insert(out.end(), {0, i, i + 1});
    return out;
  }
  std::vector<std::uint32_t> ring(p.size());
  std::iota(ring.begin(), ring.end(), 0);
  // Slice triangles can introduce collinear points along face diagonals.
  for (std::size_t i = 0; ring.size() > 3 && i < ring.size();) {
    const auto a = ring[(i + ring.size() - 1) % ring.size()], b = ring[i],
               c = ring[(i + 1) % ring.size()];
    if (std::abs(turn(p[a], p[b], p[c])) <= eps) {
      ring.erase(ring.begin() + std::ptrdiff_t(i));
      i = 0;
    } else
      ++i;
  }
  while (ring.size() > 3) {
    bool found = false;
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const auto a = ring[(i + ring.size() - 1) % ring.size()], b = ring[i],
                 c = ring[(i + 1) % ring.size()];
      if (turn(p[a], p[b], p[c]) <= eps)
        continue;
      bool inside = false;
      for (auto j : ring)
        if (j != a && j != b && j != c && turn(p[a], p[b], p[j]) >= -eps &&
            turn(p[b], p[c], p[j]) >= -eps && turn(p[c], p[a], p[j]) >= -eps) {
          inside = true;
          break;
        }
      if (inside)
        continue;
      out.insert(out.end(), {a, b, c});
      ring.erase(ring.begin() + std::ptrdiff_t(i));
      found = true;
      break;
    }
    KUTTA_REQUIRE(
        found,
        "cannot triangulate self-intersecting or degenerate view polygon");
  }
  out.insert(out.end(), ring.begin(), ring.end());
  return out;
}
void append_deltas(std::vector<Field> &fields) {
  const auto count = fields.size();
  for (std::size_t i = 0; i < count; ++i) {
    const auto &f = fields[i];
    if (!f.name.ends_with("@fluent"))
      continue;
    const auto name = f.name.substr(0, f.name.size() - 7);
    auto base =
        std::find_if(fields.begin(), fields.begin() + std::ptrdiff_t(count),
                     [&](const auto &b) { return b.name == name; });
    if (base == fields.begin() + std::ptrdiff_t(count) ||
        std::any_of(fields.begin(), fields.end(),
                    [&](const auto &v) { return v.name == name + "@delta"; }))
      continue;
    Field d{name + "@delta", f.unit, f.components, {}};
    d.values.reserve(f.values.size());
    for (std::size_t j = 0; j < f.values.size(); ++j)
      d.values.push_back(base->values[j] - f.values[j]);
    fields.push_back(std::move(d));
  }
}

} // namespace kutta::view
