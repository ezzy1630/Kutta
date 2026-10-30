#pragma once
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/post/post.hpp"
#include <functional>

namespace post_fixture {
using namespace kutta;
using namespace kutta::detail;
inline view::Snapshot snapshot(int nx, int ny, int nz, Vec3d lo, Vec3d hi,
                               const std::function<Vec3d(Vec3d)> &velocity,
                               const std::function<double(Vec3d)> &pressure) {
  view::Snapshot s;
  s.mesh = make_box_hex(nx, ny, nz, lo, hi);
  s.properties = {2, .01, 2, 2, 1, 3, 2};
  const auto geometry = compute_geometry(s.mesh);
  view::Field u{"velocity", "m/s", 3, {}}, p{"pressure", "Pa", 1, {}},
      phi{"phi", "m", 1, {}};
  for (auto x : geometry.cell_centroid) {
    const auto v = velocity(x);
    u.values.insert(u.values.end(), {float(v.x), float(v.y), float(v.z)});
    p.values.push_back(float(pressure(x)));
    phi.values.push_back(float(x.x));
  }
  s.cells = {
      u,
      p,
      phi,
      {"k", "m^2/s^2", 1, std::vector<float>(s.mesh.n_cells, 1)},
      {"omega", "1/s", 1, std::vector<float>(s.mesh.n_cells, 10)},
      {"viscosity_ratio", "1", 1, std::vector<float>(s.mesh.n_cells, 0)}};
  view::Field ub{"velocity", "m/s", 3, {}}, pb{"pressure", "Pa", 1, {}},
      shear{"wall_shear", "Pa", 3, {}}, flux{"mass_flux", "kg/s", 1, {}};
  for (Index f = s.mesh.n_internal_faces(); f < s.mesh.n_faces(); ++f) {
    const auto x = geometry.face_centroid[f], v = velocity(x);
    ub.values.insert(ub.values.end(), {float(v.x), float(v.y), float(v.z)});
    pb.values.push_back(float(pressure(x)));
    shear.values.insert(shear.values.end(), {0, -3, 0});
    flux.values.push_back(float(2 * dot(v, geometry.face_area[f])));
  }
  s.boundary = {
      ub,
      pb,
      shear,
      flux,
      {"y_plus", "1", 1,
       std::vector<float>(s.mesh.n_faces() - s.mesh.n_internal_faces(), 7)}};
  return s;
}
inline view::Snapshot uniform(int n = 3) {
  return snapshot(
      n, n, n, {0, 0, 0}, {1, 1, 1}, [](Vec3d) { return Vec3d{1, 0, 0}; },
      [](Vec3d p) { return 2 + 3 * p.x - 2 * p.y + p.z; });
}
inline const Json::Object &object(const Json &j) {
  return std::get<Json::Object>(j.value);
}
inline const Json::Array &array(const Json &j) {
  return std::get<Json::Array>(j.value);
}
inline double number(const Json &j) { return std::get<double>(j.value); }
inline Vec3d point(const view::Geometry &g, std::size_t n) {
  return {g.positions[n * 3], g.positions[n * 3 + 1], g.positions[n * 3 + 2]};
}
} // namespace post_fixture
