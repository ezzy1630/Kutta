#pragma once

// One model per device/mesh; neither may move or die before the model.
//
//   SstModel sst(dev, mesh, SstOptions{});
//   // Each outer iteration, after computing UNLIMITED cell gradients:
//   sst.update({velocity, {grad_u, grad_v, grad_w}, rho, mu, k, omega,
//               grad_k, grad_omega, wall_velocity});
//   const auto& t = sst.fields();
//   // Momentum: t.mu_t; wall-function stress REPLACES wall-face viscous flux.
//   // Freeze C=t.wall_shear_coefficient; implicitly apply area*C*(I-n*n^T)
//   // to cell velocity, with area*C*(I-n*n^T)*U_wall on the RHS.
//   // t.wall_shear is the equivalent explicit stress on the fluid.
//   // Use t.wall_face_omega for wall gradient/reconstruction values.
//   // Transport: add V*source to RHS and V*sink to diagonal; source - sink*phi.
//   // Use t.face_diffusivity_{k,omega}; wall k diffusivity is zero.
//   // Replace omega rows where t.fixed_omega_mask != 0 by identity with
//   // RHS t.fixed_omega. wall_face_omega is NOT a transport Dirichlet BC.
//
// Fields remain on device; update allocates nothing and never downloads data.

#include <optional>

#include "kutta/fv/device_mesh.hpp"

namespace kutta {

enum class SstGammaSet { Computed1994, Menter2003 };
enum class SstViscosityLimiter { Strain, Vorticity };
enum class SstWallOmega { Correlation, Tabulated };

// k-equation side of Fluent's Intermittency Transition Model (TG 4-205..207
// and Menter et al. 2015, see transition.md): k production gamma*G~_k + P_k^lim,
// destruction min(max(gamma,.1),1)*Y_k and F1 = max(F1, F3). The gamma
// equation itself lives in IntermittencyTransport.
struct SstIntermittency {
  float ck = 1.0f;    // C_k in P_k^lim; Fluent rpvar trans-one-eqn-tke-term-coef
  float csep = 1.0f;  // C_SEP in P_k^lim; Fluent rpvar trans-one-eqn-csep
};

struct SstOptions {
  bool production_limiter = true;
  float production_clip = 10.0f;  // C_lim
  bool limit_omega_production = true;  // documented inference; data does not settle this
  bool low_re = false;
  SstGammaSet gamma_set = SstGammaSet::Computed1994;
  SstViscosityLimiter viscosity_limiter = SstViscosityLimiter::Strain;
  bool kato_launder = false;
  bool vorticity_production = false;
  bool curvature_correction = false;  // unsupported: needs strain material derivative
  bool compressibility = false;      // unsupported: incompressible interface
  SstWallOmega wall_omega = SstWallOmega::Correlation;  // Tabulated throws
  bool wall_production = true;  // TG 4-419..422, not verified from saved fields
  float cdkw_floor = 1e-10f;   // floor for F1 ONLY; transport CDkw keeps its sign
  // Inferred from Fluent's sst-limit-cross-diff-to-n-times-w-destruct rpvar;
  // clip the signed transport D_omega to +/- n*rho*beta*omega^2, not F1's CDkw.
  bool limit_cross_diffusion = true;
  float cross_diffusion_clip = 1e6f;
  float omega_floor = 1e-20f;
  float k_floor = 1e-14f;
  float max_viscosity_ratio = 1e5f;
  // Absent: plain SST, with no extra kernel or buffer.
  std::optional<SstIntermittency> intermittency;
};

// Tensor rows are velocity components: u.x = du/dx, u.y = du/dy, etc.
struct SstVelocityGradients {
  const DeviceVec3& u;
  const DeviceVec3& v;
  const DeviceVec3& w;
};

struct SstInputs {
  const DeviceVec3& velocity;
  // Derivative for rates, production and mu_t, independently of SOU.
  SstVelocityGradients velocity_gradient;
  const Buffer<float>& density;    // per cell, finite and > 0
  const Buffer<float>& viscosity;  // dynamic molecular viscosity, finite and > 0
  const Buffer<float>& k;          // finite; clipped to the explicit k_floor
  const Buffer<float>& omega;      // finite; clipped to the explicit omega_floor
  const DeviceVec3& k_gradient;
  const DeviceVec3& omega_gradient;
  const DeviceVec3& wall_velocity;  // per BOUNDARY face; stationary walls require zeros
  // Required exactly when SstOptions::intermittency is set; per cell, in [0, 1].
  const Buffer<float>* intermittency = nullptr;
};

struct SstFields {
  Buffer<float> mu_t;
  Buffer<float> k_source, k_sink, omega_source, omega_sink;
  Buffer<float> face_diffusivity_k, face_diffusivity_omega;

  // Constraints/overrides are per cell. Mask == 0 means no wall constraint.
  Buffer<std::int32_t> fixed_omega_mask;
  Buffer<float> fixed_omega, wall_production;
  // Nearest wall face's BOUNDARY index, or -1; tie -> smallest global face id.
  Buffer<std::int32_t> selected_wall_face;

  // All wall diagnostics use boundary-face indexing; non-wall entries are zero.
  Buffer<float> wall_tau, wall_u_star, wall_u_tau, wall_y_plus, wall_y_plus_utau;
  DeviceVec3 wall_shear;  // Pa, antiparallel to relative tangential velocity
  // C = tau_w/|U_t|, with its finite analytic limit at |U_t|=0. See the
  // implicit momentum contract above; lag C, not the tangential cell velocity.
  Buffer<float> wall_shear_coefficient;
  // Gradient/reconstruction ONLY: blend(6*nu/(.072*dy^2),
  // 2*u_star/(.3*.4187*dy), 1.3). Forensics 3.5 suggests this form but does
  // not establish its role; fixed constants miss saved face SV_O by ~0.42%.
  // Independent of input gradients; a bootstrap update can supply the first
  // reconstruction, then refresh from the current nonlinear iterate.
  Buffer<float> wall_face_omega;

  // Diagnostics and cell diffusivities, also useful to the transport assembler.
  Buffer<float> strain, vorticity, f1, f2, cdkw;
  Buffer<float> diffusivity_k, diffusivity_omega, production_k;
  // P_k^lim, allocated only with SstOptions::intermittency. production_k stays
  // SST's G~_k; k_source/k_sink carry the intermittency-coupled terms.
  Buffer<float> production_limit;
};

class SstModel {
 public:
  // Mesh connectivity/geometry obey DeviceMesh's contract. Wall distances are
  // positive; non-wall wall_normal_distance entries are zero. PatchKind::Wall
  // means a smooth no-slip wall. Rough/zero-shear walls need a different BC.
  SstModel(Device& device, const DeviceMesh& mesh, SstOptions options = {});
  void update(const SstInputs& inputs);
  const SstFields& fields() const { return fields_; }
  const SstOptions& options() const { return options_; }

 private:
  Device& device_;
  const DeviceMesh& mesh_;
  SstOptions options_;
  SstFields fields_;
  Buffer<std::int32_t> wall_mask_;
  Buffer<float> wall_cell_omega_, face_production_;
};

}  // namespace kutta
