#pragma once

#include "kutta/fv/device_mesh.hpp"

namespace kutta {
enum class SaProduction { Vorticity, StrainVorticity };
enum class SaFormulation { Fluent, TmrStandard };
enum class SaWallTreatment { Resolved, FluentAutomatic };
struct SaOptions {
  SaProduction production=SaProduction::Vorticity;
  // Fluent TG 4-15..29 omits ft2; TMR standard SA includes it.
  SaFormulation formulation=SaFormulation::Fluent;
  SaWallTreatment wall=SaWallTreatment::Resolved;
  float cb1=.1355f,cb2=.622f,sigma=2.f/3,cv1=7.1f,cw2=.3f,cw3=2,kappa=.41f;
  float max_viscosity_ratio=1e5f;
};
struct SaInputs {
  const Buffer<float>& density;
  const Buffer<float>& viscosity;
  const Buffer<float>& modified_viscosity; // nu-tilde, NOT mu_t/rho
  const DeviceVec3& gradient;
  const Buffer<float>& strain;
  const Buffer<float>& vorticity;
};
struct SaFields {
  Buffer<float> mu_t,source,sink,diffusivity,face_diffusivity;
  Buffer<float> fv1,fv2,modified_rate,fw;
};
struct SaWallFields {DeviceVec3 shear;Buffer<float> y_plus;};
// All fields remain on Device. Sources are per volume: source-sink*nu_tilde.
// Resolved smooth wall: nu_tilde=0; molecular wall stress. No automatic wall
// blend is claimed: Fluent publishes its asymptotes but not the calibrated blend.
class SaModel {
 public:
  SaModel(Device&,const DeviceMesh&,SaOptions={});
  void update(const SaInputs&);
  const SaFields& fields() const {return fields_;}
  const SaOptions& options() const {return options_;}
 private:
  Device& device_;const DeviceMesh& mesh_;SaOptions options_;SaFields fields_;
  Buffer<Index> wall_;
};
} // namespace kutta
