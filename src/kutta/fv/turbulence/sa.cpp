#include "kutta/fv/turbulence/sa.hpp"
#include <cmath>
#include "kernels.hpp"

namespace kutta {
namespace {
template<class T> void size(const Buffer<T>& b,std::size_t n) {
  KUTTA_REQUIRE(b.size()==n && b.handle(),"SA: wrong size/unallocated input");
}
}
SaModel::SaModel(Device& d,const DeviceMesh& m,SaOptions o):device_(d),mesh_(m),options_(o) {
  for(float x:{o.cb1,o.cb2,o.sigma,o.cv1,o.cw2,o.cw3,o.kappa,o.max_viscosity_ratio})
    KUTTA_REQUIRE(std::isfinite(x) && x>0,"SA: invalid constant");
  KUTTA_REQUIRE(o.production==SaProduction::Vorticity || o.production==SaProduction::StrainVorticity,"SA: unknown production");
  KUTTA_REQUIRE(o.formulation==SaFormulation::Fluent || o.formulation==SaFormulation::TmrStandard,"SA: unknown formulation");
  KUTTA_REQUIRE(o.wall==SaWallTreatment::Resolved,"SA: Fluent automatic wall blending is unverified and unavailable");
  KUTTA_REQUIRE(m.n_cells>0,"SA: empty mesh");
  size(m.wall_distance,m.n_cells);size(m.owner,m.n_faces);
  size(m.neighbour,m.n_internal_faces);size(m.weight,m.n_internal_faces);
  std::vector<Index> wall(m.n_boundary_faces(),0);Index next=m.n_internal_faces;
  for(const auto& p:m.patches) {
    KUTTA_REQUIRE(p.start==next && p.size>=0 && p.size<=m.n_faces-next,"SA: invalid patch coverage");
    if(p.kind==PatchKind::Wall) for(Index f=p.start;f<p.start+p.size;++f) wall[f-m.n_internal_faces]=1;
    next+=p.size;
  }
  KUTTA_REQUIRE(next==m.n_faces,"SA: incomplete patch coverage");wall_=d.upload<Index>(wall);
  auto& f=fields_;
  for(auto* b:{&f.mu_t,&f.source,&f.sink,&f.diffusivity,&f.fv1,&f.fv2,&f.modified_rate,&f.fw}) *b=d.alloc<float>(m.n_cells);
  f.face_diffusivity=d.alloc<float>(m.n_faces);
  rate_derivative_=d.alloc<float>(m.n_cells);fw_derivative_=d.alloc<float>(m.n_cells);
}
void SaModel::update(const SaInputs& in) {
  const auto n=mesh_.n_cells;auto& f=fields_;const auto& o=options_;auto& d=device_;
  for(const auto* b:{&in.density,&in.viscosity,&in.modified_viscosity,&in.gradient.x,&in.gradient.y,&in.gradient.z,&in.strain,&in.vorticity}) size(*b,n);
  kern::sa_coefficients(d,n,in.density,in.viscosity,in.modified_viscosity,mesh_.wall_distance,in.strain,in.vorticity,
    f.mu_t,f.diffusivity,f.fv1,f.fv2,f.modified_rate,f.fw,rate_derivative_,fw_derivative_,o.production==SaProduction::StrainVorticity,
    o.formulation==SaFormulation::TmrStandard,o.sigma,o.cv1,o.cw2,o.cw3,o.kappa,o.max_viscosity_ratio);
  kern::sa_sources(d,n,in.density,in.viscosity,in.modified_viscosity,mesh_.wall_distance,
    in.gradient.x,in.gradient.y,in.gradient.z,f.modified_rate,f.fw,rate_derivative_,fw_derivative_,f.source,f.sink,
    o.formulation==SaFormulation::TmrStandard,o.cb1,o.cb2,o.sigma,o.kappa);
  kern::sa_face_diffusion(d,mesh_.n_faces,mesh_.owner,mesh_.neighbour,mesh_.weight,wall_,in.viscosity,
    f.diffusivity,f.face_diffusivity,mesh_.n_internal_faces,o.sigma);
}
} // namespace kutta
