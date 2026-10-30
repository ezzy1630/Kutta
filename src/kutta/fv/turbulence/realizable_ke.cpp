#include "kutta/fv/turbulence/realizable_ke.hpp"
#include <cmath>
#include "kernels.hpp"
namespace kutta {
namespace {
template<class T> void size(const Buffer<T>& b,std::size_t n) {KUTTA_REQUIRE(b.handle() && b.size()==n,"realizable k-epsilon: invalid buffer");}
}
RealizableKeModel::RealizableKeModel(Device& d,const DeviceMesh& m,RealizableKeOptions o):dev_(d),mesh_(m),options_(o) {
  for(float x:{o.a0,o.c2,o.sigma_k,o.sigma_epsilon,o.k_floor,o.epsilon_floor,o.max_viscosity_ratio})
    KUTTA_REQUIRE(std::isfinite(x) && x>0,"realizable k-epsilon: invalid constant");
  KUTTA_REQUIRE(std::isfinite(o.molecular_production) && o.molecular_production>=0,"realizable k-epsilon: invalid molecular production");
  KUTTA_REQUIRE(o.wall==KeWallTreatment::Scalable,"realizable k-epsilon: full EWT transport is unavailable; blend width and discrete epsilon coupling are not specified by Fluent TG");
  KUTTA_REQUIRE(m.n_cells>0,"realizable k-epsilon: empty mesh");
  size(m.owner,m.n_faces);size(m.neighbour,m.n_internal_faces);size(m.weight,m.n_internal_faces);
  size(m.area.x,m.n_faces);size(m.area.y,m.n_faces);size(m.area.z,m.n_faces);
  size(m.wall_normal_distance,m.n_boundary_faces());size(m.cell_face_offsets,std::size_t(m.n_cells)+1);
  KUTTA_REQUIRE(m.cell_faces.handle(),"realizable k-epsilon: missing incidence CSR");
  std::vector<Index> wall(m.n_boundary_faces(),0);Index next=m.n_internal_faces;
  for(const auto& p:m.patches) {
    KUTTA_REQUIRE(p.start==next && p.size>=0 && p.size<=m.n_faces-next,"realizable k-epsilon: invalid patches");
    if(p.kind==PatchKind::Wall) for(Index f=p.start;f<p.start+p.size;++f) wall[f-m.n_internal_faces]=1;
    next+=p.size;
  }
  KUTTA_REQUIRE(next==m.n_faces,"realizable k-epsilon: incomplete patches");wall_=d.upload<Index>(wall);
  auto& f=fields_;const auto nc=m.n_cells,nb=m.n_boundary_faces();
  for(auto* b:{&f.mu_t,&f.cmu,&f.strain,&f.k_source,&f.k_sink,&f.epsilon_source,&f.epsilon_sink,&f.diffusivity_k,&f.diffusivity_epsilon,&f.fixed_epsilon}) *b=d.alloc<float>(nc);
  for(auto* b:{&f.wall_coefficient,&f.wall_y_star,&f.wall_tau,&f.wall_epsilon,&f.wall_production,&f.wall_shear.x,&f.wall_shear.y,&f.wall_shear.z}) *b=d.alloc<float>(nb);
  f.face_diffusivity_k=d.alloc<float>(m.n_faces);f.face_diffusivity_epsilon=d.alloc<float>(m.n_faces);
  f.fixed_epsilon_mask=d.alloc<Index>(nc);selected_=d.alloc<Index>(nc);
  wall_production_cell_=d.alloc<float>(nc);
}
void RealizableKeModel::update(const RealizableKeInputs& in) {
  const auto n=mesh_.n_cells,nb=mesh_.n_boundary_faces();const auto& g=in.velocity_gradient;auto& f=fields_;auto& d=dev_;const auto& o=options_;
  for(const auto* b:{&in.density,&in.viscosity,&in.k,&in.epsilon,&in.velocity.x,&in.velocity.y,&in.velocity.z,
      &g.u.x,&g.u.y,&g.u.z,&g.v.x,&g.v.y,&g.v.z,&g.w.x,&g.w.y,&g.w.z}) size(*b,n);
  for(const auto* b:{&in.wall_velocity.x,&in.wall_velocity.y,&in.wall_velocity.z}) size(*b,nb);
  kern::rke_viscosity(d,n,in.density,in.viscosity,in.k,in.epsilon,g.u.x,g.u.y,g.u.z,g.v.x,g.v.y,g.v.z,g.w.x,g.w.y,g.w.z,
    f.mu_t,f.cmu,f.strain,o.a0,o.k_floor,o.epsilon_floor,o.max_viscosity_ratio);
  kern::rke_scalable_wall(d,nb,mesh_.owner,wall_,mesh_.wall_normal_distance,mesh_.area.x,mesh_.area.y,mesh_.area.z,
    in.velocity.x,in.velocity.y,in.velocity.z,in.wall_velocity.x,in.wall_velocity.y,in.wall_velocity.z,in.density,in.viscosity,in.k,
    f.wall_coefficient,f.wall_y_star,f.wall_tau,f.wall_epsilon,f.wall_production,f.wall_shear.x,f.wall_shear.y,f.wall_shear.z,mesh_.n_internal_faces,o.k_floor);
  // Cell source production is separate from the per-face diagnostics.
  auto& prod=wall_production_cell_;
  kern::rke_wall_cells(d,n,mesh_.cell_face_offsets,mesh_.cell_faces,wall_,mesh_.wall_normal_distance,
    f.wall_epsilon,f.wall_production,f.fixed_epsilon_mask,selected_,f.fixed_epsilon,prod,mesh_.n_internal_faces);
  kern::rke_sources(d,n,in.density,in.viscosity,in.k,in.epsilon,f.mu_t,f.strain,f.fixed_epsilon_mask,prod,f.fixed_epsilon,
    f.k_source,f.k_sink,f.epsilon_source,f.epsilon_sink,f.diffusivity_k,f.diffusivity_epsilon,o.c2,o.sigma_k,o.sigma_epsilon,
    o.molecular_production,o.k_floor,o.epsilon_floor);
  kern::rke_face_diffusion(d,mesh_.n_faces,mesh_.owner,mesh_.neighbour,mesh_.weight,wall_,in.viscosity,
    f.diffusivity_k,f.diffusivity_epsilon,f.face_diffusivity_k,f.face_diffusivity_epsilon,mesh_.n_internal_faces);
}
EwtFields allocate_ewt_fields(Device& d,std::size_t n) {return {d.alloc<float>(n),d.alloc<float>(n),d.alloc<float>(n)};}
void evaluate_ewt(Device& d,const DeviceMesh& m,const Buffer<float>& rho,const Buffer<float>& mu,
    const Buffer<float>& k,const Buffer<float>& outer,EwtBlend b,EwtFields& f) {
  KUTTA_REQUIRE(std::isfinite(b.delta_re_y) && b.delta_re_y>0 && std::isfinite(b.cmu) && b.cmu>0,"EWT: invalid explicit blend parameters");
  for(const auto* p:{&rho,&mu,&k,&outer}) size(*p,m.n_cells);
  for(const auto* p:{&f.viscosity,&f.epsilon,&f.outer_weight}) size(*p,m.n_cells);
  kern::rke_ewt_algebra(d,m.n_cells,rho,mu,k,m.wall_distance,outer,f.viscosity,f.epsilon,f.outer_weight,b.delta_re_y,b.cmu);
}
} // namespace kutta
