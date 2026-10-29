#include "kutta/fv/turbulence/sst.hpp"

#include <cmath>

#include "kernels.hpp"

namespace kutta {
namespace {
void positive(float value, const char* name) {
  KUTTA_REQUIRE(std::isfinite(value) && value > 0, std::string("SST: invalid ") + name);
}
void validate(const SstOptions& o) {
  positive(o.production_clip, "production clip factor");
  positive(o.cdkw_floor, "CDkw floor");
  positive(o.omega_floor, "omega floor");
  positive(o.k_floor, "k floor");
  positive(o.max_viscosity_ratio, "maximum viscosity ratio");
  KUTTA_REQUIRE(!o.curvature_correction, "SST: curvature correction is unsupported (requires strain material derivatives)");
  KUTTA_REQUIRE(!o.compressibility, "SST: compressibility is unsupported by the incompressible model");
  KUTTA_REQUIRE(!(o.kato_launder && o.vorticity_production), "SST: Kato-Launder and vorticity production are mutually exclusive");
  switch (o.gamma_set) {
    case SstGammaSet::Computed1994: case SstGammaSet::Menter2003: break;
    default: throw Error("SST: invalid gamma set");
  }
  switch (o.viscosity_limiter) {
    case SstViscosityLimiter::Strain: case SstViscosityLimiter::Vorticity: break;
    default: throw Error("SST: invalid viscosity limiter");
  }
  switch (o.wall_omega) {
    case SstWallOmega::Correlation: break;
    case SstWallOmega::Tabulated:
      throw Error("SST: tabulated wall omega is unsupported; Fluent's table is unpublished");
    default: throw Error("SST: invalid wall omega treatment");
  }
}
template<class T>
void size(const Buffer<T>& b, std::size_t n, const char* name) {
  KUTTA_REQUIRE(b.size() == n && b.handle(), std::string("SST: wrong size/unallocated ") + name);
}
void size(const DeviceVec3& v, std::size_t n, const char* name) {
  size(v.x, n, name); size(v.y, n, name); size(v.z, n, name);
}
}  // namespace

SstModel::SstModel(Device& dev, const DeviceMesh& mesh, SstOptions options)
    : device_(dev), mesh_(mesh), options_(options) {
  validate(options_);
  KUTTA_REQUIRE(mesh.n_cells > 0 && mesh.n_internal_faces >= 0 &&
                    mesh.n_faces >= mesh.n_internal_faces, "SST: invalid mesh counts");
  const auto nc = mesh.n_cells, nf = mesh.n_faces, nb = mesh.n_boundary_faces();
  size(mesh.owner, nf, "mesh owner"); size(mesh.neighbour, mesh.n_internal_faces, "mesh neighbour");
  size(mesh.weight, mesh.n_internal_faces, "mesh weights");
  size(mesh.area, nf, "mesh area");
  size(mesh.wall_distance, nc, "wall distance");
  size(mesh.wall_normal_distance, nb, "wall normal distance");
  size(mesh.cell_face_offsets, std::size_t(nc)+1, "cell face offsets");
  KUTTA_REQUIRE(mesh.cell_faces.handle(), "SST: missing cell face CSR");
  std::vector<std::int32_t> wall(nb, 0);
  Index next = mesh.n_internal_faces;
  for (const auto& patch : mesh.patches) {
    KUTTA_REQUIRE(patch.start == next && patch.size >= 0 && patch.size <= nf-next,
                  "SST: boundary patches must cover faces contiguously");
    if (patch.kind == PatchKind::Wall)
      for (Index f = patch.start; f < patch.start+patch.size; ++f)
        wall[f-mesh.n_internal_faces] = 1;
    next += patch.size;
  }
  KUTTA_REQUIRE(next == nf, "SST: boundary patches do not cover faces");
  wall_mask_ = dev.upload<std::int32_t>(wall);
  auto cell = [&] { return dev.alloc<float>(nc); };
  auto boundary = [&] { return dev.alloc<float>(nb); };
  auto& f = fields_;
  f.mu_t = cell(); f.k_source = cell(); f.k_sink = cell();
  f.omega_source = cell(); f.omega_sink = cell();
  f.strain = cell(); f.vorticity = cell(); f.f1 = cell(); f.f2 = cell(); f.cdkw = cell();
  f.diffusivity_k = cell(); f.diffusivity_omega = cell(); f.production_k = cell();
  f.fixed_omega_mask = dev.alloc<std::int32_t>(nc);
  f.selected_wall_face = dev.alloc<std::int32_t>(nc);
  f.fixed_omega = cell(); f.wall_production = cell();
  f.face_diffusivity_k = dev.alloc<float>(nf); f.face_diffusivity_omega = dev.alloc<float>(nf);
  f.wall_tau = boundary(); f.wall_u_star = boundary(); f.wall_u_tau = boundary();
  f.wall_y_plus = boundary(); f.wall_y_plus_utau = boundary();
  f.wall_shear = {boundary(), boundary(), boundary()};
  face_omega_ = boundary(); face_production_ = boundary();
}

void SstModel::update(const SstInputs& in) {
  const auto nc = mesh_.n_cells, nf = mesh_.n_faces, nb = mesh_.n_boundary_faces();
  size(in.velocity, nc, "velocity");
  const auto& g = in.velocity_gradient;
  size(g.u, nc, "grad u"); size(g.v, nc, "grad v"); size(g.w, nc, "grad w");
  size(in.density, nc, "density"); size(in.viscosity, nc, "viscosity");
  size(in.k, nc, "k"); size(in.omega, nc, "omega");
  size(in.k_gradient, nc, "grad k"); size(in.omega_gradient, nc, "grad omega");
  size(in.wall_velocity, nb, "wall velocity");
  auto& d = device_; auto& f = fields_; const auto& o = options_;
  kern::sst_rates(d, nc, g.u.x, g.u.y, g.u.z, g.v.x, g.v.y, g.v.z, g.w.x, g.w.y, g.w.z,
                  f.strain, f.vorticity);
  kern::sst_blending(d, nc, in.density, in.viscosity, in.k, in.omega, mesh_.wall_distance,
      in.k_gradient.x, in.k_gradient.y, in.k_gradient.z,
      in.omega_gradient.x, in.omega_gradient.y, in.omega_gradient.z,
      f.f1, f.f2, f.cdkw, o.cdkw_floor, o.k_floor, o.omega_floor);
  kern::sst_viscosity(d, nc, in.density, in.viscosity, in.k, in.omega,
      f.strain, f.vorticity, f.f1, f.f2, f.mu_t, f.diffusivity_k, f.diffusivity_omega,
      o.low_re, o.viscosity_limiter == SstViscosityLimiter::Vorticity,
      o.max_viscosity_ratio, o.k_floor, o.omega_floor);
  kern::sst_wall(d, nb, mesh_.owner, wall_mask_, mesh_.n_internal_faces,
      mesh_.area.x, mesh_.area.y, mesh_.area.z,
      in.velocity.x, in.velocity.y, in.velocity.z,
      in.wall_velocity.x, in.wall_velocity.y, in.wall_velocity.z,
      mesh_.wall_normal_distance, in.density, in.viscosity, in.k, in.omega,
      f.wall_tau, f.wall_shear.x, f.wall_shear.y, f.wall_shear.z,
      f.wall_u_star, f.wall_u_tau, f.wall_y_plus, f.wall_y_plus_utau,
      face_omega_, face_production_, o.k_floor, o.omega_floor);
  kern::sst_wall_cells(d, nc, mesh_.cell_face_offsets, mesh_.cell_faces, wall_mask_,
      mesh_.wall_normal_distance, face_omega_, face_production_, mesh_.n_internal_faces,
      f.fixed_omega_mask, f.selected_wall_face, f.fixed_omega, f.wall_production);
  const int production_mode = o.kato_launder ? 1 : (o.vorticity_production ? 2 : 0);
  kern::sst_sources(d, nc, in.density, in.viscosity, in.k, in.omega,
      f.strain, f.vorticity, f.f1, f.cdkw, f.mu_t, f.fixed_omega_mask, f.wall_production,
      f.k_source, f.k_sink, f.omega_source, f.omega_sink, f.production_k,
      o.low_re, o.gamma_set == SstGammaSet::Menter2003, production_mode,
      o.production_limiter, o.limit_omega_production, o.wall_production,
      o.production_clip, o.k_floor, o.omega_floor);
  kern::sst_face_diffusion(d, nf, mesh_.owner, mesh_.neighbour, mesh_.weight, wall_mask_,
      in.viscosity, f.diffusivity_k, f.diffusivity_omega, mesh_.n_internal_faces,
      f.face_diffusivity_k, f.face_diffusivity_omega);
}
}  // namespace kutta
