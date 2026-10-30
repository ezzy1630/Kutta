#pragma once
#include "kutta/fv/turbulence/sst.hpp"

namespace kutta {
enum class KeWallTreatment { Scalable, Enhanced, Standard };
// The k-epsilon family shares transport, boundaries and wall-cell coupling.
// Standard (TG 4.3.1) uses constant Cmu, c1/c2_standard and sigma_epsilon=1.3;
// select it with standard_ke_options().
enum class KeVariant { Realizable, Standard };
struct RealizableKeOptions {
  KeVariant variant=KeVariant::Realizable;
  float a0=4.04f,c2=1.9f,sigma_k=1,sigma_epsilon=1.2f;
  float cmu_standard=.09f,c1_standard=1.44f,c2_standard=1.92f;
  float k_floor=1e-14f,epsilon_floor=1e-20f,max_viscosity_ratio=1e5f;
  // TG 4-61: high-Re k-epsilon production uses mu_t+C*mu, default C=1.
  float molecular_production=1;
  KeWallTreatment wall=KeWallTreatment::Scalable;
};
// Standard k-epsilon with standard wall functions, Fluent's defaults.
inline RealizableKeOptions standard_ke_options() {
  RealizableKeOptions o;o.variant=KeVariant::Standard;o.sigma_epsilon=1.3f;o.wall=KeWallTreatment::Standard;return o;
}
struct RealizableKeInputs {
  const DeviceVec3& velocity;SstVelocityGradients velocity_gradient;
  const Buffer<float>& density;const Buffer<float>& viscosity;
  const Buffer<float>& k;const Buffer<float>& epsilon;
  const DeviceVec3& wall_velocity;
};
struct RealizableKeFields {
  Buffer<float> mu_t,cmu,strain,k_source,k_sink,epsilon_source,epsilon_sink;
  Buffer<float> diffusivity_k,diffusivity_epsilon,face_diffusivity_k,face_diffusivity_epsilon;
  Buffer<float> wall_coefficient,wall_y_star,wall_tau,wall_epsilon,wall_production;
  DeviceVec3 wall_shear;
  Buffer<Index> fixed_epsilon_mask;Buffer<float> fixed_epsilon;
};
// Inertial-frame, incompressible, smooth-wall realizable k-epsilon. Device-only
// updates; source-sink*phi per volume. Wall shear replaces wall diffusion, with
// the same implicit tangential momentum contract as SST. k has zero wall flux;
// replace wall-cell epsilon rows with fixed_epsilon. No rotating-frame term.
class RealizableKeModel {
 public:
  RealizableKeModel(Device&,const DeviceMesh&,RealizableKeOptions={});
  void update(const RealizableKeInputs&);
  const RealizableKeFields& fields() const {return fields_;}
 private:
  Device& dev_;const DeviceMesh& mesh_;RealizableKeOptions options_;RealizableKeFields fields_;
  Buffer<Index> wall_,selected_;
  Buffer<float> wall_production_cell_;
};
// Documented EWT algebra only, NOT a complete wall treatment. Fluent does not
// specify the blend width or the discrete epsilon coupling. Requiring them
// here prevents a guessed choice from being represented as a Fluent default.
struct EwtBlend {float delta_re_y;float cmu;};
struct EwtFields {Buffer<float> viscosity,epsilon,outer_weight;};
EwtFields allocate_ewt_fields(Device&,std::size_t);
void evaluate_ewt(Device&,const DeviceMesh&,const Buffer<float>& rho,const Buffer<float>& mu,
  const Buffer<float>& k,const Buffer<float>& outer_viscosity,EwtBlend,EwtFields&);
} // namespace kutta
