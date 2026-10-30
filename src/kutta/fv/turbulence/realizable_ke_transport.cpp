#include "kutta/fv/turbulence/realizable_ke_transport.hpp"
#include <cmath>
#include <limits>
#include "kernels.hpp"
namespace kutta::fv {
RealizableKeTransport::RealizableKeTransport(Device& d,const DeviceMesh& m,const Stencil& s,const PolyMesh& h,
    const MeshGeometry& g,const std::map<std::string,BoundaryCondition>& bc,const SolverSettings& o,std::shared_ptr<CsrMatrix> matrix,
    std::shared_ptr<const WalkFaceReconstruction> shared_faces)
    :d_(d),mesh_(m),stencil_(s),settings_(*o.realizable_ke),rho_(float(o.density)),mu_(float(o.viscosity)),limit_(o.limit),
    model_(d,m,settings_.model),k_equation_(d,m,s,h,g,matrix),e_equation_(d,m,s,h,g,std::move(matrix)) {
  for(double relaxation:{settings_.k_relaxation,settings_.epsilon_relaxation})
    KUTTA_REQUIRE(std::isfinite(relaxation) && relaxation>0 && relaxation<=1,"realizable k-epsilon: invalid relaxation");
  const auto n=m.n_cells,nb=m.n_boundary_faces();
  for(auto* p:{&k_,&epsilon_,&rho_cell_,&mu_cell_,&residual_,&absolute_,&scale_,&fluent_scale_}) *p=d.alloc<float>(n);
  d.fill(k_,settings_.model.k_floor);d.fill(epsilon_,settings_.model.epsilon_floor);d.fill(rho_cell_,rho_);d.fill(mu_cell_,mu_);
  bk_=d.alloc<float>(nb);be_=d.alloc<float>(nb);kind_=d.alloc<Index>(nb);zero_=d.alloc<Index>(n);d.fill(zero_,0);
  std::vector<Index> methods(nb,0);std::vector<float> first(nb,0),second(nb,0);
  // Wall-function k has dk/dn=0 at walls (TG "Wall Treatment for
  // epsilon-based Models"): a zero-difference normal row, which also closes
  // wall-cell fits whose interior neighbours lie in one plane. Wall epsilon
  // is fixed in the wall cell, so its face data stay omitted.
  std::vector<BoundaryObservation> observations(nb,BoundaryObservation::Omitted),k_observations(nb,BoundaryObservation::Omitted);
  for(const auto& p:m.patches) {
    const auto& condition=bc.at(p.name);const std::optional<TurbulenceCondition>* spec=nullptr;
    if(const auto* inlet=std::get_if<VelocityInlet>(&condition)) spec=&inlet->turbulence;
    if(const auto* outlet=std::get_if<PressureOutlet>(&condition)) spec=&outlet->turbulence;
    Index method=0;float a=0,b=0;
    if(spec) {
      KUTTA_REQUIRE(spec->has_value(),"realizable k-epsilon: missing inlet/backflow data: "+p.name);
      if(const auto* v=std::get_if<KeValues>(&**spec)) {a=float(v->k);b=float(v->epsilon);}
      else if(const auto* v=std::get_if<IntensityLength>(&**spec)) {method=1;a=float(v->intensity);b=float(v->length);}
      else if(const auto* v=std::get_if<IntensityViscosityRatio>(&**spec)) {method=2;a=float(v->intensity);b=float(v->ratio);}
      else if(const auto* v=std::get_if<IntensityDiameter>(&**spec)) {method=3;a=float(v->intensity);b=float(v->diameter);}
      else throw Error("realizable k-epsilon: wrong boundary turbulence model: "+p.name);
      KUTTA_REQUIRE(std::isfinite(a) && a>=0 && std::isfinite(b) && b>0,"realizable k-epsilon: invalid boundary data: "+p.name);
    }
    for(Index f=p.start;f<p.start+p.size;++f) {
      auto j=f-m.n_internal_faces;methods[j]=method;first[j]=a;second[j]=b;
      if(p.kind==PatchKind::VelocityInlet) observations[j]=k_observations[j]=BoundaryObservation::FaceValue;
      if(p.kind==PatchKind::Wall) k_observations[j]=BoundaryObservation::NormalFootValue;
    }
  }
  spec_=d.upload<Index>(methods);a_=d.upload<float>(first);b_=d.upload<float>(second);
  // Green-Gauss observes every boundary face value (operators.hpp).
  if(o.gradient==GradientScheme::GreenGaussCell) {k_lsq_=green_gauss_cell(d,h,g);e_lsq_=green_gauss_cell(d,h,g);}
  else {k_lsq_=boundary_least_squares(d,h,g,k_observations);e_lsq_=boundary_least_squares(d,h,g,observations);}
  // A Green-Gauss case keeps Fluent's documented face split (steady_solver.cpp).
  if(o.gradient!=GradientScheme::GreenGaussCell) {
    const auto faces=walk_reconstruction(std::move(shared_faces),d,h,g,m,s,observations);
    k_equation_.set_reconstruction(faces);e_equation_.set_reconstruction(faces);
  }
}
void RealizableKeTransport::initialize(std::span<const double> kinetic,std::span<const double> dissipation,const MeshPermutation& p) {
  KUTTA_REQUIRE(kinetic.size()==std::size_t(mesh_.n_cells) && dissipation.size()==kinetic.size(),"realizable k-epsilon: initialization size mismatch");
  std::vector<float> k(kinetic.size()),e(k.size());
  for(std::size_t c=0;c<k.size();++c) {
    k[c]=float(kinetic[p.cells_new_to_old[c]]);e[c]=float(dissipation[p.cells_new_to_old[c]]);
    KUTTA_REQUIRE(std::isfinite(k[c]) && k[c]>=0 && std::isfinite(e[c]) && e[c]>0,"realizable k-epsilon: invalid initial fields");
  }
  d_.upload<float>(k_,k);d_.upload<float>(epsilon_,e);
}
void RealizableKeTransport::update(const DeviceVec3& velocity,const DeviceVec3& boundary,SstVelocityGradients g,
    const Buffer<Index>& kind,const Buffer<float>& mass,Buffer<float>& viscosity) {
  mass_=&mass;model_.update({velocity,g,rho_cell_,mu_cell_,k_,epsilon_,boundary});
  kern::rke_boundary(d_,mesh_.n_boundary_faces(),mesh_.owner,kind,spec_,a_,b_,mass,boundary.x,boundary.y,boundary.z,
    k_,epsilon_,bk_,be_,kind_,mesh_.n_internal_faces,rho_,mu_);
  // Outlet observations are extrapolated and omitted, like SST k.
  // model_.update does not depend on these scalar gradients.
  kern::fv_effective_viscosity(d_,mesh_.n_faces,mesh_.owner,mesh_.neighbour,kind,mesh_.weight,
    model_.fields().mu_t,viscosity,mesh_.n_internal_faces,mu_);
}
void RealizableKeTransport::assemble(bool epsilon,float time_step) {
  KUTTA_REQUIRE(mass_,"realizable k-epsilon: update before assembly");
  const auto& f=model_.fields();auto& eq=epsilon?e_equation_:k_equation_;
  TransportSettings t;t.density=rho_;t.pseudo_time_step=time_step;t.limit=limit_;
  t.convection=epsilon?settings_.epsilon_convection:settings_.k_convection;
  t.relaxation=float(epsilon?settings_.epsilon_relaxation:settings_.k_relaxation);
  gradient(d_,mesh_,stencil_,epsilon?e_lsq_:k_lsq_,epsilon?epsilon_:k_,epsilon?be_:bk_,eq.grad);
  eq.assemble_reconstructed(epsilon?epsilon_:k_,epsilon?be_:bk_,kind_,*mass_,epsilon?f.face_diffusivity_epsilon:f.face_diffusivity_k,
    epsilon?f.epsilon_source:f.k_source,t,{epsilon?&f.epsilon_sink:&f.k_sink});
  if(epsilon) kern::fv_fixed_rows(d_,mesh_.n_cells,eq.matrix.row_ptr,eq.matrix.diag,f.fixed_epsilon_mask,f.fixed_epsilon,eq.matrix.val,eq.rhs);
}
SolveStats RealizableKeTransport::solve(bool epsilon,float time_step) {
  assemble(epsilon,time_step);auto& q=epsilon?epsilon_:k_;auto& eq=epsilon?e_equation_:k_equation_;
  auto stats=eq.solve(q,settings_.linear);kern::fv_floor(d_,mesh_.n_cells,q,epsilon?settings_.model.epsilon_floor:settings_.model.k_floor);return stats;
}
std::array<ResidualValue,2> RealizableKeTransport::monitor(float time_step) {
  std::array<ResidualValue,2> result;
  for(int j=0;j<2;++j) {
    const bool epsilon=j==1;assemble(epsilon,time_step);auto& eq=epsilon?e_equation_:k_equation_;const auto& q=epsilon?epsilon_:k_;
    residual(d_,eq.matrix,q,eq.rhs,residual_);
    kern::fv_scalar_monitor(d_,mesh_.n_cells,q,eq.diagonal,residual_,epsilon?model_.fields().fixed_epsilon_mask:zero_,mesh_.volume,
      absolute_,scale_,fluent_scale_,rho_,time_step);
    auto raw=d_.sum(absolute_,mesh_.n_cells),den=d_.sum(scale_,mesh_.n_cells),fd=d_.sum(fluent_scale_,mesh_.n_cells);
    auto scaled=[&](double x){return x>0?raw/x:(raw==0?0:std::numeric_limits<double>::infinity());};result[j]={raw,scaled(den),scaled(fd),fd};
  }
  return result;
}
} // namespace kutta::fv
