#pragma once

#include <array>
#include <limits>

#include "kutta/fv/turbulence/sst.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace sst_test {
using namespace kutta;
using namespace kutta::detail;

inline Buffer<float> filled(Device& dev, std::size_t n, float value) {
  auto b = dev.alloc<float>(n); dev.fill(b, value); return b;
}
inline DeviceVec3 vector(Device& dev, const std::vector<Vec3d>& values) {
  std::vector<float> x, y, z;
  for (auto p : values) { x.push_back(float(p.x)); y.push_back(float(p.y)); z.push_back(float(p.z)); }
  return {dev.upload<float>(x), dev.upload<float>(y), dev.upload<float>(z)};
}
inline DeviceVec3 constant_vector(Device& dev, std::size_t n, Vec3d value = {}) {
  return vector(dev, std::vector<Vec3d>(n, value));
}
inline Buffer<float> scalar(Device& dev, const std::vector<double>& values) {
  return dev.upload<float>(std::vector<float>(values.begin(), values.end()));
}

// TEST ONLY: populate DeviceMesh directly; deliberately independent of lane D's
// upload_mesh(). wall_distances is supplied by the test (saved SV_WALL_DIST for
// Fluent; exact box distances for unit tests), not a production distance search.
inline DeviceMesh mesh(Device& dev, const PolyMesh& m, const MeshGeometry& g,
                       const std::vector<double>& wall_distances) {
  DeviceMesh d;
  d.n_cells = m.n_cells; d.n_faces = m.n_faces(); d.n_internal_faces = m.n_internal_faces();
  d.owner = dev.upload<Index>(m.owner); d.neighbour = dev.upload<Index>(m.neighbour);
  d.area = vector(dev, g.face_area);
  std::vector<Vec3d> rp(m.n_faces()), rn(m.n_internal_faces()), pn(m.n_internal_faces());
  std::vector<float> weights(m.n_internal_faces());
  std::vector<double> normal_distances(m.n_faces()-m.n_internal_faces());
  std::vector<bool> wall(normal_distances.size());
  for (const auto& p : m.patches) {
    d.patches.push_back({p.name, p.kind, p.start, p.size});
    if (p.kind == PatchKind::Wall)
      for (Index f = p.start; f < p.start+p.size; ++f) wall[f-m.n_internal_faces()] = true;
  }
  for (Index f = 0; f < m.n_faces(); ++f) {
    rp[f] = g.face_centroid[f] - g.cell_centroid[m.owner[f]];
    if (f < m.n_internal_faces()) {
      rn[f] = g.face_centroid[f] - g.cell_centroid[m.neighbour[f]];
      pn[f] = g.cell_centroid[m.neighbour[f]] - g.cell_centroid[m.owner[f]];
      weights[f] = float(dot(g.face_area[f], rn[f] * -1.0) / dot(g.face_area[f], pn[f]));
    } else if (wall[f-m.n_internal_faces()]) {
      normal_distances[f-m.n_internal_faces()] = dot(rp[f], g.face_area[f]) / norm(g.face_area[f]);
    }
  }
  d.owner_to_face = vector(dev, rp); d.neighbour_to_face = vector(dev, rn);
  d.owner_to_neighbour = vector(dev, pn); d.weight = dev.upload<float>(weights);
  d.volume = scalar(dev, g.cell_volume);
  auto cf = build_cell_faces(m);
  d.cell_face_offsets = dev.upload<Index>(cf.offsets); d.cell_faces = dev.upload<Index>(cf.faces);
  std::vector<float> signs(cf.signs.begin(), cf.signs.end());
  d.cell_face_sign = dev.upload<float>(signs);
  d.wall_distance = scalar(dev, wall_distances); d.wall_normal_distance = scalar(dev, normal_distances);
  return d;
}
inline DeviceMesh box(Device& dev, int nx = 1) {
  auto m = make_box_hex(nx, 1, 1, {0,0,0}, {1,1,1});
  auto g = compute_geometry(m);
  std::vector<double> distance;
  for (auto c : g.cell_centroid) distance.push_back(std::min({c.x,1-c.x,c.y,1-c.y,c.z,1-c.z}));
  return mesh(dev, m, g, distance);
}

struct State {
  DeviceVec3 velocity, gu, gv, gw, gk, go, wall_velocity;
  Buffer<float> rho, mu, k, omega;
  State(Device& dev, const DeviceMesh& m)
      : velocity(constant_vector(dev,m.n_cells)), gu(constant_vector(dev,m.n_cells)),
        gv(constant_vector(dev,m.n_cells)), gw(constant_vector(dev,m.n_cells)),
        gk(constant_vector(dev,m.n_cells)), go(constant_vector(dev,m.n_cells)),
        wall_velocity(constant_vector(dev,m.n_boundary_faces())),
        rho(filled(dev,m.n_cells,1)), mu(filled(dev,m.n_cells,1e-3f)),
        k(filled(dev,m.n_cells,1)), omega(filled(dev,m.n_cells,2)) {}
  SstInputs inputs() const { return {velocity,{gu,gv,gw},rho,mu,k,omega,gk,go,wall_velocity}; }
};

inline std::vector<std::vector<float>> snapshot(Device& d, const SstFields& f) {
  return {d.download(f.mu_t), d.download(f.k_source), d.download(f.k_sink),
    d.download(f.omega_source), d.download(f.omega_sink),
    d.download(f.face_diffusivity_k), d.download(f.face_diffusivity_omega),
    d.download(f.fixed_omega), d.download(f.wall_production),
    d.download(f.wall_tau), d.download(f.wall_u_star), d.download(f.wall_u_tau),
    d.download(f.wall_y_plus), d.download(f.wall_y_plus_utau),
    d.download(f.wall_shear.x), d.download(f.wall_shear.y), d.download(f.wall_shear.z),
    d.download(f.strain), d.download(f.vorticity), d.download(f.f1), d.download(f.f2),
    d.download(f.cdkw), d.download(f.diffusivity_k), d.download(f.diffusivity_omega),
    d.download(f.production_k)};
}
}  // namespace sst_test
