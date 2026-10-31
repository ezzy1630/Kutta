#include "kutta/fv/turbulence/intermittency.hpp"
#include <cmath>
#include <limits>
#include "kernels.hpp"

namespace kutta::fv {
namespace {
void finite(float value,const char* name) {
  KUTTA_REQUIRE(std::isfinite(value),std::string("transition: non-finite ")+name);
}
}
IntermittencyTransport::IntermittencyTransport(Device& d,const DeviceMesh& m,const Stencil& s,const PolyMesh& h,
    const MeshGeometry& g,const std::map<std::string,BoundaryCondition>& bc,const SolverSettings& o,
    std::shared_ptr<CsrMatrix> matrix,const BoundaryLeastSquares& k_lsq,std::shared_ptr<KrylovWorkspace> workspace)
    :dev_(d),mesh_(m),stencil_(s),lsq_(k_lsq),settings_(*o.sst->transition),rho_(float(o.density)),
    k_floor_(o.sst->model.k_floor),omega_floor_(o.sst->model.omega_floor),limit_(o.limit),
    equation_(d,m,s,h,g,std::move(matrix),nullptr,std::move(workspace)) {
  const auto& t=settings_.model;
  for(auto [v,name]:{std::pair{t.ctu1,"C_TU1"},{t.ctu2,"C_TU2"},{t.ctu3,"C_TU3"},{t.cpg1,"C_PG1"},{t.cpg2,"C_PG2"},
      {t.cpg3,"C_PG3"},{t.cpg1_lim,"C_PG1_lim"},{t.cpg2_lim,"C_PG2_lim"}}) finite(v,name);
  KUTTA_REQUIRE(std::isfinite(t.sigma) && t.sigma>0,"transition: sigma_gamma must be positive");
  KUTTA_REQUIRE(std::isfinite(t.minimum) && std::isfinite(t.maximum) && t.minimum>0 && t.minimum<t.maximum,
                "transition: intermittency limits must satisfy 0 < min < max");
  KUTTA_REQUIRE(std::isfinite(settings_.relaxation) && settings_.relaxation>0 && settings_.relaxation<=1,
                "transition: invalid intermittency relaxation");
  const auto n=m.n_cells,nb=m.n_boundary_faces();
  for(auto* b:{&gamma_,&residual_,&absolute_,&scale_,&fluent_scale_,&fields_.source,&fields_.sink,
               &fields_.normal_derivative,&fields_.re_theta_c,&fields_.onset}) *b=d.alloc<float>(n);
  fields_.face_diffusivity=d.alloc<float>(m.n_faces);
  boundary_=d.alloc<float>(nb);kind_=d.alloc<Index>(nb);zero_=d.alloc<Index>(n);d.fill(zero_,0);
  std::vector<float> inflow(nb,1);std::vector<Index> wall(nb,0);
  for(const auto& p:m.patches) {
    const auto& condition=bc.at(p.name);double value=1;
    if(const auto* v=std::get_if<VelocityInlet>(&condition)) value=v->intermittency;
    if(const auto* v=std::get_if<PressureOutlet>(&condition)) value=v->intermittency;
    if(const auto* v=std::get_if<PressureInlet>(&condition)) value=v->intermittency;
    if(const auto* v=std::get_if<MassFlowInlet>(&condition)) value=v->intermittency;
    KUTTA_REQUIRE(std::isfinite(value) && value>=0 && value<=1,"transition: inflow intermittency must lie in [0, 1]: "+p.name);
    for(Index f=p.start;f<p.start+p.size;++f) {
      const auto b=f-m.n_internal_faces;inflow[b]=float(value);wall[b]=p.kind==PatchKind::Wall?1:0;
    }
  }
  inflow_=d.upload<float>(inflow);wall_=d.upload<Index>(wall);
  fields_.wall_direction=wall_direction(d,h,g);
  d.fill(gamma_,1);
}
void IntermittencyTransport::initialize(std::span<const double> values,const MeshPermutation& permutation) {
  KUTTA_REQUIRE(values.size()==std::size_t(mesh_.n_cells),"transition: initialization size mismatch");
  std::vector<float> q(values.size());
  for(std::size_t c=0;c<q.size();++c) {
    const double x=values[permutation.cells_new_to_old[c]];q[c]=float(x);
    KUTTA_REQUIRE(std::isfinite(x) && x>=0 && x<=1,"transition: initial intermittency must lie in [0, 1]");
  }
  dev_.upload<float>(gamma_,q);
}
void IntermittencyTransport::fill(double gamma) {
  KUTTA_REQUIRE(std::isfinite(gamma) && gamma>=0 && gamma<=1,"transition: initial intermittency must lie in [0, 1]");
  dev_.fill(gamma_,float(gamma));
}
void IntermittencyTransport::update(const SstFields& sst,const Buffer<float>& rho,const Buffer<float>& mu,
    const Buffer<float>& k,const Buffer<float>& omega,const DeviceVec3 (&g)[3],
    const Buffer<Index>& bc,const Buffer<float>& mass) {
  mass_=&mass;const auto n=mesh_.n_cells;const auto& t=settings_.model;auto& f=fields_;
  const auto& w=f.wall_direction;
  kern::transition_normal_derivative(dev_,n,w.x,w.y,w.z,g[0].x,g[0].y,g[0].z,g[1].x,g[1].y,g[1].z,
    g[2].x,g[2].y,g[2].z,f.normal_derivative);
  // SST's floors keep Tu_L and R_T finite where k or omega reach them.
  kern::transition_sources(dev_,n,rho,mu,k,omega,mesh_.wall_distance,f.normal_derivative,sst.strain,sst.vorticity,
    gamma_,f.source,f.sink,f.re_theta_c,f.onset,t.ctu1,t.ctu2,t.ctu3,t.cpg1,t.cpg2,t.cpg3,t.cpg1_lim,t.cpg2_lim,
    k_floor_,omega_floor_);
  kern::transition_face_diffusion(dev_,mesh_.n_faces,mesh_.owner,mesh_.neighbour,mesh_.weight,wall_,mu,sst.mu_t,
    mesh_.n_internal_faces,f.face_diffusivity,1.0f/t.sigma);
  kern::transition_boundary(dev_,mesh_.n_boundary_faces(),mesh_.owner,bc,inflow_,mass,gamma_,boundary_,kind_,mesh_.n_internal_faces);
  gradient(dev_,mesh_,stencil_,lsq_,gamma_,boundary_,equation_.grad);
}
void IntermittencyTransport::assemble(float time_step) {
  KUTTA_REQUIRE(mass_,"transition: update before assembly");
  TransportSettings t;t.density=rho_;t.pseudo_time_step=time_step;t.relaxation=float(settings_.relaxation);
  t.limit=limit_;t.convection=settings_.convection;
  equation_.assemble_reconstructed(gamma_,boundary_,kind_,*mass_,fields_.face_diffusivity,fields_.source,t,{&fields_.sink});
}
SolveStats IntermittencyTransport::solve(float time_step) {
  assemble(time_step);auto result=equation_.solve(gamma_,settings_.linear);
  // Fluent's trans-intermit/min-limit and max-limit.
  kern::transition_clamp(dev_,mesh_.n_cells,gamma_,settings_.model.minimum,settings_.model.maximum);
  return result;
}
ResidualValue IntermittencyTransport::monitor(float time_step) {
  assemble(time_step);residual(dev_,equation_.matrix,gamma_,equation_.rhs,residual_);
  kern::fv_scalar_monitor(dev_,mesh_.n_cells,gamma_,equation_.diagonal,residual_,zero_,mesh_.volume,
    absolute_,scale_,fluent_scale_,rho_,time_step);
  const auto raw=dev_.sum(absolute_,mesh_.n_cells),den=dev_.sum(scale_,mesh_.n_cells),fd=dev_.sum(fluent_scale_,mesh_.n_cells);
  const auto scaled=[&](double d){return d>0?raw/d:(raw==0?0:std::numeric_limits<double>::infinity());};
  return {raw,scaled(den),scaled(fd),fd};
}
} // namespace kutta::fv
