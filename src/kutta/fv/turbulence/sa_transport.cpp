#include "kutta/fv/turbulence/sa_transport.hpp"
#include <cmath>
#include <limits>
#include "kernels.hpp"

namespace kutta::fv {
SaTransport::SaTransport(Device& d,const DeviceMesh& m,const Stencil& s,const PolyMesh& h,
    const MeshGeometry& g,const std::map<std::string,BoundaryCondition>& bc,const SolverSettings& o,
    std::shared_ptr<CsrMatrix> matrix)
    :dev_(d),mesh_(m),stencil_(s),settings_(*o.sa),rho_(float(o.density)),mu_(float(o.viscosity)),limit_(o.limit),
    model_(d,m,settings_.model),equation_(d,m,s,h,g,std::move(matrix)) {
  KUTTA_REQUIRE(std::isfinite(settings_.relaxation) && settings_.relaxation>0 && settings_.relaxation<=1,"SA: invalid relaxation");
  const auto n=m.n_cells,nb=m.n_boundary_faces();
  for(auto* b:{&q_,&rho_cell_,&mu_cell_,&strain_,&rotation_,&residual_,&absolute_,&scale_,&fluent_scale_}) *b=d.alloc<float>(n);
  d.fill(q_,0);d.fill(rho_cell_,rho_);d.fill(mu_cell_,mu_);
  boundary_=d.alloc<float>(nb);kind_=d.alloc<Index>(nb);zero_=d.alloc<Index>(n);d.fill(zero_,0);
  std::vector<Index> spec(nb,0);std::vector<float> values(nb,0),lengths(nb,0),wall(nb,0);
  std::vector<BoundaryObservation> observations(nb,BoundaryObservation::Omitted);
  for(const auto& p:m.patches) {
    const auto& condition=bc.at(p.name);const std::optional<TurbulenceCondition>* t=nullptr;
    if(const auto* inlet=std::get_if<VelocityInlet>(&condition)) t=&inlet->turbulence;
    if(const auto* outlet=std::get_if<PressureOutlet>(&condition)) t=&outlet->turbulence;
    float value=0,length=0;Index method=0;
    if(t) {
      KUTTA_REQUIRE(t->has_value(),"SA: missing inlet/backflow data: "+p.name);
      if(const auto* direct=std::get_if<ModifiedViscosity>(&**t)) value=float(direct->value);
      else if(const auto* ratio=std::get_if<SaViscosityRatio>(&**t)) {method=1;value=float(ratio->ratio);}
      else if(const auto* v=std::get_if<IntensityLength>(&**t)) {method=2;value=float(v->intensity);length=float(v->length);}
      else if(const auto* v=std::get_if<IntensityDiameter>(&**t)) {method=2;value=float(v->intensity);length=float(.07*v->diameter);}
      else throw Error("SA: incompatible boundary turbulence data: "+p.name);
      if(method==2) KUTTA_REQUIRE(std::isfinite(length) && length>0,"SA: invalid turbulence length: "+p.name);
      KUTTA_REQUIRE(std::isfinite(value) && value>=0,"SA: invalid boundary data: "+p.name);
    }
    for(Index f=p.start;f<p.start+p.size;++f) {
      auto b=f-m.n_internal_faces;spec[b]=method;values[b]=value;lengths[b]=length;
      wall[b]=p.kind==PatchKind::Wall?1:0;
      if(p.kind==PatchKind::Wall || p.kind==PatchKind::VelocityInlet) observations[b]=BoundaryObservation::FaceValue;
    }
  }
  spec_=d.upload<Index>(spec);values_=d.upload<float>(values);length_=d.upload<float>(lengths);lsq_=boundary_least_squares(d,h,g,observations);
  wall_mask_=d.upload<float>(wall);wall_fields_.shear=vector_field(d,nb);wall_fields_.y_plus=d.alloc<float>(nb);
}
void SaTransport::initialize(std::span<const double> values,const MeshPermutation& permutation) {
  KUTTA_REQUIRE(values.size()==std::size_t(mesh_.n_cells),"SA: initialization size mismatch");
  std::vector<float> q(values.size());
  for(std::size_t c=0;c<q.size();++c) {
    const double x=values[permutation.cells_new_to_old[c]];q[c]=float(x);
    KUTTA_REQUIRE(std::isfinite(x) && x>=0 && std::isfinite(q[c]),"SA: invalid initial modified viscosity");
  }
  dev_.upload<float>(q_,q);
}
void SaTransport::update(const DeviceVec3& velocity,const DeviceVec3& wall_velocity,const DeviceVec3& u,const DeviceVec3& v,const DeviceVec3& w,
    const Buffer<Index>& bc,const Buffer<float>& mass,Buffer<float>& face_viscosity) {
  mass_=&mass;
  kern::sa_boundary(dev_,mesh_.n_boundary_faces(),mesh_.owner,bc,spec_,values_,length_,wall_velocity.x,wall_velocity.y,wall_velocity.z,mass,q_,boundary_,kind_,mesh_.n_internal_faces,mu_/rho_,settings_.model.cv1);
  gradient(dev_,mesh_,stencil_,lsq_,q_,boundary_,equation_.grad);
  kern::sst_rates(dev_,mesh_.n_cells,u.x,u.y,u.z,v.x,v.y,v.z,w.x,w.y,w.z,strain_,rotation_);
  model_.update({rho_cell_,mu_cell_,q_,equation_.grad,strain_,rotation_});
  kern::sa_momentum_viscosity(dev_,mesh_.n_faces,mesh_.owner,mesh_.neighbour,bc,mesh_.weight,
    model_.fields().mu_t,face_viscosity,mesh_.n_internal_faces,mu_);
  auto& t=wall_fields_.shear;
  kern::fv_wall_viscous_force(dev_,mesh_.n_boundary_faces(),mesh_.owner,wall_mask_,mesh_.area.x,mesh_.area.y,mesh_.area.z,
    mesh_.owner_to_face.x,mesh_.owner_to_face.y,mesh_.owner_to_face.z,velocity.x,velocity.y,velocity.z,
    u.x,u.y,u.z,v.x,v.y,v.z,w.x,w.y,w.z,wall_velocity.x,wall_velocity.y,wall_velocity.z,
    t.x,t.y,t.z,mesh_.n_internal_faces,mu_);
  kern::sa_wall_diagnostics(dev_,mesh_.n_boundary_faces(),mesh_.area.x,mesh_.area.y,mesh_.area.z,
    mesh_.wall_normal_distance,t.x,t.y,t.z,wall_fields_.y_plus,mesh_.n_internal_faces,rho_,mu_);
}
void SaTransport::assemble(float time_step) {
  KUTTA_REQUIRE(mass_,"SA: update before assembly");const auto& f=model_.fields();
  TransportSettings t;t.density=rho_;t.pseudo_time_step=time_step;t.relaxation=float(settings_.relaxation);
  t.limit=limit_;t.convection=settings_.convection;
  equation_.assemble_reconstructed(q_,boundary_,kind_,*mass_,f.face_diffusivity,f.source,t,{&f.sink});
}
SolveStats SaTransport::solve(float time_step) {
  assemble(time_step);auto result=equation_.solve(q_,settings_.linear);
  kern::fv_floor(dev_,mesh_.n_cells,q_,0);return result;
}
ResidualValue SaTransport::monitor(float time_step) {
  assemble(time_step);residual(dev_,equation_.matrix,q_,equation_.rhs,residual_);
  kern::fv_scalar_monitor(dev_,mesh_.n_cells,q_,equation_.diagonal,residual_,zero_,mesh_.volume,
    absolute_,scale_,fluent_scale_,rho_,time_step);
  auto raw=dev_.sum(absolute_,mesh_.n_cells),den=dev_.sum(scale_,mesh_.n_cells),fd=dev_.sum(fluent_scale_,mesh_.n_cells);
  auto scaled=[&](double d){return d>0?raw/d:(raw==0?0:std::numeric_limits<double>::infinity());};
  return {raw,scaled(den),scaled(fd),fd};
}
} // namespace kutta::fv
