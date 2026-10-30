#include "kutta/fv/steady_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>

#include "kernels.hpp"
#include "kutta/kernels/prelude.h"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta::fv {
using namespace detail;
namespace {
template<class Profile>
Vec3d evaluate(const Profile& profile,Vec3d point) {
  if(const auto* uniform=std::get_if<Vec3d>(&profile)) return *uniform;
  return std::get<std::function<Vec3d(Vec3d)>>(profile)(point);
}
void upload_vectors(Device& dev,DeviceVec3& to,const std::vector<Vec3d>& values) {
  std::vector<float> x,y,z;
  for(auto p:values) {
    KUTTA_REQUIRE(finite(p),"FV nonfinite field value");
    x.push_back(float(p.x));y.push_back(float(p.y));z.push_back(float(p.z));
  }
  dev.upload<float>(to.x,x);dev.upload<float>(to.y,y);dev.upload<float>(to.z,z);
}
// The cell whose face planes all enclose the point, else the nearest centroid
// (points outside the domain, or nonconvex cells).
Index locate_cell(const PolyMesh& mesh,const MeshGeometry& geom,Vec3d point) {
  KUTTA_REQUIRE(finite(point),"FV nonfinite pressure reference point");
  const auto cf=build_cell_faces(mesh);
  Index nearest=0;double best=std::numeric_limits<double>::infinity();
  for(Index c=0;c<mesh.n_cells;++c) {
    bool inside=true;
    for(Index j=cf.offsets[c];j<cf.offsets[c+1] && inside;++j) {
      const auto f=cf.faces[j];const auto s=geom.face_area[f]*double(cf.signs[j]);
      inside=dot(s,point-geom.face_centroid[f])<=1e-12*norm(s)*std::cbrt(geom.cell_volume[c]);
    }
    if(inside) return c;
    const auto d=point-geom.cell_centroid[c];
    if(dot(d,d)<best) {best=dot(d,d);nearest=c;}
  }
  return nearest;
}
}

struct SteadySolver::Impl {
  Device& dev;
  PolyMesh host;
  MeshGeometry geom;
  SolverSettings settings;
  DeviceMesh mesh;
  Stencil stencil;
  BoundaryLeastSquares velocity_lsq,pressure_lsq,k_lsq,omega_lsq,strain_lsq,jump_fit;
  DeviceVec3 strain_base[3],strain_gradient[3];
  MeshPermutation permutation;
  std::unique_ptr<SstModel> sst;
  std::unique_ptr<TransportEquation> k_equation,omega_equation;
  Buffer<float> kinetic,omega,rho_cell,mu_cell,boundary_k,boundary_omega,turbulence_a,turbulence_b;
  Buffer<Index> turbulence_spec,turbulence_kind,zero_mask,total_pressure;
  Buffer<float> wall_diagonal,wall_rhs,physical_diagonal,consistent_inverse;
  DeviceVec3 transpose_stress;
  DeviceVec3 velocity,body,boundary_velocity,prescribed_velocity,source,face_velocity;
  Buffer<float> pressure,mass,viscosity,boundary_pressure,prescribed_pressure;
  Buffer<Index> boundary_kind,scalar_kind;
  DeviceVec3 grad_p,limited_gradient[3];
  std::unique_ptr<TransportEquation> momentum[3];
  std::shared_ptr<CsrMatrix> pressure_matrix;
  std::unique_ptr<LinearSolver> pressure_solver;
  Buffer<float> pressure_rhs,correction,correction_boundary,coefficient,face_pressure,pressure_jump;
  DeviceVec3 grad_correction;
  Buffer<float> residual_u,residual_v,residual_w,continuity,provisional_continuity,scale,speed;
  DeviceVec3 absolute_residual,fluent_scale;
  IterationReport report;
  std::map<std::string,BoundaryCondition> boundaries;
  Index reference=-1;
  double length=0,bc_speed=0,pressure_speed=0,continuity_scale=0;
  double first_five_continuity[5]{};
  float time_step=0;

  Impl(Device& device,PolyMesh input,std::vector<PatchCondition> conditions,SolverSettings options)
      :dev(device),host(std::move(input)),settings(options) {
    KUTTA_REQUIRE(settings.density>0 && std::isfinite(settings.density) && settings.viscosity>0 &&
      std::isfinite(settings.viscosity) && settings.momentum_relaxation>0 && settings.momentum_relaxation<=1 &&
      settings.pressure_relaxation>0 && settings.pressure_relaxation<=1 &&
      settings.pseudo_time_scale>0 && settings.pseudo_time_step>=0 && settings.monitor_interval>0,
      "FV invalid steady solver settings");
    for(auto& bc:conditions)
      KUTTA_REQUIRE(boundaries.emplace(bc.patch,std::move(bc.condition)).second,"FV duplicate patch condition");
    KUTTA_REQUIRE(boundaries.size()==host.patches.size(),"FV supply exactly one condition per patch");
    bool outlet=false;
    for(auto& patch:host.patches) {
      const auto it=boundaries.find(patch.name);
      KUTTA_REQUIRE(it!=boundaries.end(),"FV missing boundary: "+patch.name);
      if(std::holds_alternative<Wall>(it->second)) patch.kind=PatchKind::Wall;
      if(std::holds_alternative<VelocityInlet>(it->second)) patch.kind=PatchKind::VelocityInlet;
      if(std::holds_alternative<PressureOutlet>(it->second)) {patch.kind=PatchKind::PressureOutlet;outlet=true;}
      if(std::holds_alternative<Symmetry>(it->second)) patch.kind=PatchKind::Symmetry;
    }
    permutation=renumber_cells_rcm(host);geom=compute_geometry(host);
    KUTTA_REQUIRE(host.n_cells>0,"FV steady solver needs a nonempty mesh");
    if(!outlet) reference=locate_cell(host,geom,settings.reference_point);
    // A single pin cannot constrain disconnected Neumann components.
    const auto neighbours=build_cell_neighbours(host);
    std::vector<bool> seen(host.n_cells);std::vector<Index> stack{0};seen[0]=true;
    for(std::size_t k=0;k<stack.size();++k) for(Index j=neighbours.offsets[stack[k]];j<neighbours.offsets[stack[k]+1];++j)
      if(!seen[neighbours.cells[j]]) {seen[neighbours.cells[j]]=true;stack.push_back(neighbours.cells[j]);}
    KUTTA_REQUIRE(stack.size()==std::size_t(host.n_cells),"FV disconnected domains require separate solvers");
    mesh=upload_mesh(dev,host,geom);stencil=make_stencil(dev,host,geom);
    const auto n=mesh.n_cells,nf=mesh.n_faces,nb=mesh.n_boundary_faces();
    velocity=vector_field(dev,n);face_velocity=vector_field(dev,nf);body=vector_field(dev,n);source=vector_field(dev,n);
    boundary_velocity=vector_field(dev,nb);prescribed_velocity=vector_field(dev,nb);
    for(auto& g:limited_gradient) g=vector_field(dev,n);
    grad_p=vector_field(dev,n);grad_correction=vector_field(dev,n);absolute_residual=vector_field(dev,n);
    fluent_scale=vector_field(dev,n);
    pressure=dev.alloc<float>(n);dev.fill(pressure,0);
    mass=dev.alloc<float>(nf);dev.fill(mass,0);viscosity=dev.alloc<float>(nf);dev.fill(viscosity,float(settings.viscosity));
    boundary_pressure=dev.alloc<float>(nb);dev.fill(boundary_pressure,0);
    prescribed_pressure=dev.alloc<float>(nb);boundary_kind=dev.alloc<Index>(nb);scalar_kind=dev.alloc<Index>(nb);
    pressure_rhs=dev.alloc<float>(n);correction=dev.alloc<float>(n);dev.fill(correction,0);
    correction_boundary=dev.alloc<float>(nb);dev.fill(correction_boundary,0);
    coefficient=dev.alloc<float>(nf);face_pressure=dev.alloc<float>(nf);pressure_jump=dev.alloc<float>(nf);
    residual_u=dev.alloc<float>(n);residual_v=dev.alloc<float>(n);residual_w=dev.alloc<float>(n);
    continuity=dev.alloc<float>(n);provisional_continuity=dev.alloc<float>(n);dev.fill(provisional_continuity,0);scale=dev.alloc<float>(n);speed=dev.alloc<float>(n);
    wall_diagonal=dev.alloc<float>(n);wall_rhs=dev.alloc<float>(n);
    total_pressure=dev.alloc<Index>(nb);
    std::vector<Index> totals(nb,0);
    physical_diagonal=dev.alloc<float>(n);consistent_inverse=dev.alloc<float>(n);
    dev.fill(physical_diagonal,1);dev.fill(consistent_inverse,0);
    std::vector<Index> kinds(nb);std::vector<float> pressures(nb,0);std::vector<Vec3d> velocities(nb);
    double pmin=std::numeric_limits<double>::infinity(),pmax=-pmin;
    for(const auto& patch:host.patches) for(Index f=patch.start;f<patch.start+patch.size;++f) {
      const auto b=f-mesh.n_internal_faces;const auto& bc=boundaries.at(patch.name);
      if(const auto* wall=std::get_if<Wall>(&bc)) {kinds[b]=0;velocities[b]=wall->velocity;}
      if(const auto* inlet=std::get_if<VelocityInlet>(&bc)) {kinds[b]=1;velocities[b]=evaluate(inlet->velocity,geom.face_centroid[f]);}
      if(const auto* out=std::get_if<PressureOutlet>(&bc)) {
        const bool normal=std::holds_alternative<NormalBackflow>(out->backflow_velocity);
        kinds[b]=normal?4:2;
        velocities[b]=normal?Vec3d{}:evaluate(out->backflow_velocity,geom.face_centroid[f]);
        KUTTA_REQUIRE(std::isfinite(out->pressure),"FV nonfinite outlet pressure");
        totals[b]=out->backflow_pressure==BackflowPressure::Total?1:0;
        pressures[b]=float(out->pressure);pmin=std::min(pmin,out->pressure);pmax=std::max(pmax,out->pressure);
      }
      if(std::holds_alternative<Symmetry>(bc)) kinds[b]=3;
      bc_speed=std::max(bc_speed,norm(velocities[b]));
    }
    pressure_speed=outlet?.1*std::sqrt((pmax-pmin)/settings.density):0;
    dev.upload<Index>(total_pressure,totals);
    using enum BoundaryObservation;
    std::vector<BoundaryObservation> independent_velocity(nb),independent_pressure(nb),independent_k(nb),independent_omega(nb);
    for(Index b=0;b<nb;++b) {
      independent_velocity[b]=kinds[b]==3?NormalFootValue:((kinds[b]==0 || kinds[b]==1)?FaceValue:Omitted);
      independent_pressure[b]=(kinds[b]==2 || kinds[b]==4)?FaceValue:Omitted;
      independent_k[b]=kinds[b]==1?FaceValue:Omitted;independent_omega[b]=(kinds[b]==0 || kinds[b]==1)?FaceValue:Omitted;
    }
    dev.upload<Index>(boundary_kind,kinds);dev.upload<float>(prescribed_pressure,pressures);
    upload_vectors(dev,prescribed_velocity,velocities);
    pressure_matrix=std::make_shared<CsrMatrix>(upload(dev,stencil.pattern));
    for(auto& eq:momentum) eq=std::make_unique<TransportEquation>(dev,mesh,stencil,host,geom,pressure_matrix);
    velocity_lsq=boundary_least_squares(dev,host,geom,independent_velocity);
    pressure_lsq=boundary_least_squares(dev,host,geom,independent_pressure);
    if(settings.sst) {
      std::vector<BoundaryObservation> strain_rows(nb);
      for(Index b=0;b<nb;++b) strain_rows[b]=kinds[b]!=0?FaceValue:Omitted;
      strain_lsq=boundary_least_squares(dev,host,geom,strain_rows);
      jump_fit=face_jump_least_squares(dev,host,geom);
      for(int k=0;k<3;++k) {strain_base[k]=vector_field(dev,n);strain_gradient[k]=vector_field(dev,n);}
      const auto& t=*settings.sst;
      KUTTA_REQUIRE(t.k_relaxation>0 && t.k_relaxation<=1 && t.omega_relaxation>0 && t.omega_relaxation<=1,"FV invalid SST relaxation");
      k_lsq=boundary_least_squares(dev,host,geom,independent_k);
      omega_lsq=boundary_least_squares(dev,host,geom,independent_omega);
      kinetic=dev.alloc<float>(n);omega=dev.alloc<float>(n);
      dev.fill(kinetic,t.model.k_floor);dev.fill(omega,t.model.omega_floor);
      rho_cell=dev.alloc<float>(n);mu_cell=dev.alloc<float>(n);
      dev.fill(rho_cell,float(settings.density));dev.fill(mu_cell,float(settings.viscosity));
      boundary_k=dev.alloc<float>(nb);boundary_omega=dev.alloc<float>(nb);
      turbulence_kind=dev.alloc<Index>(nb);zero_mask=dev.alloc<Index>(n);dev.fill(zero_mask,0);
      transpose_stress=vector_field(dev,nf);
      std::vector<Index> specs(nb,0);std::vector<float> a(nb,0),b(nb,0);
      for(const auto& patch:host.patches) {
        const auto& bc=boundaries.at(patch.name);
        const std::optional<TurbulenceCondition>* spec=nullptr;
        if(const auto* in=std::get_if<VelocityInlet>(&bc)) spec=&in->turbulence;
        if(const auto* out=std::get_if<PressureOutlet>(&bc)) spec=&out->turbulence;
        if(!spec) continue;
        KUTTA_REQUIRE(spec->has_value(),"FV missing turbulence specification: "+patch.name);
        for(Index f=patch.start;f<patch.start+patch.size;++f) {
          const auto bf=f-mesh.n_internal_faces;
          if(const auto* v=std::get_if<TurbulenceValues>(&**spec)) {a[bf]=float(v->k);b[bf]=float(v->omega);}
          if(const auto* v=std::get_if<IntensityLength>(&**spec)) {specs[bf]=1;a[bf]=float(v->intensity);b[bf]=float(v->length);}
          if(const auto* v=std::get_if<IntensityViscosityRatio>(&**spec)) {specs[bf]=2;a[bf]=float(v->intensity);b[bf]=float(v->ratio);}
          KUTTA_REQUIRE(a[bf]>=0 && b[bf]>0 && std::isfinite(a[bf]) && std::isfinite(b[bf]),"FV invalid turbulence specification: "+patch.name);
        }
      }
      turbulence_spec=dev.upload<Index>(specs);turbulence_a=dev.upload<float>(a);turbulence_b=dev.upload<float>(b);
      k_equation=std::make_unique<TransportEquation>(dev,mesh,stencil,host,geom,pressure_matrix);
      omega_equation=std::make_unique<TransportEquation>(dev,mesh,stencil,host,geom,pressure_matrix);
      sst=std::make_unique<SstModel>(dev,mesh,t.model);
    }
    CsrMatrixHost initial{stencil.pattern,std::vector<double>(stencil.pattern.nnz(),0)};
    const auto cf=build_cell_faces(host);
    for(Index c=0;c<n;++c) {
      kupper_sum ap{0,0};
      for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
        auto f=cf.faces[j];auto other=f<mesh.n_internal_faces?(host.owner[f]==c?host.neighbour[f]:host.owner[f]):-1;
        if(other<0 && kinds[f-mesh.n_internal_faces]!=2 && kinds[f-mesh.n_internal_faces]!=4) continue;
        auto d=(other>=0?geom.cell_centroid[other]:geom.face_centroid[f])-geom.cell_centroid[c];
        const float a=float(dot(geom.face_area[f],geom.face_area[f])/std::abs(dot(geom.face_area[f],d)));
        if(other<0 || other==reference) ap=kupper_add(ap,a);
        if(other>=0 && other!=reference && c!=reference) {
          auto index=std::lower_bound(initial.pattern.col.begin()+initial.pattern.row_ptr[c],
            initial.pattern.col.begin()+initial.pattern.row_ptr[c+1],other)-initial.pattern.col.begin();
          initial.val[index]=float(float(initial.val[index])-a);
        }
      }
      // The hierarchy is built from the same stored-float M-matrix contract as
      // runtime assembly, rather than independently rounded double diagonals.
      for(Index j=initial.pattern.row_ptr[c];j<initial.pattern.row_ptr[c+1];++j) {
        if(j!=initial.pattern.diag[c]) ap=kupper_add(ap,-float(initial.val[j]));
      }
      initial.val[initial.pattern.diag[c]]=c==reference?1:kupper_value(ap);
    }
    pressure_solver=make_solver(pressure_solver_kind,dev,initial,
      {.structure=OperatorStructure::SymmetricMMatrix});
    const double volume=std::accumulate(geom.cell_volume.begin(),geom.cell_volume.end(),0.);
    Vec3d lo=host.points.front(),hi=lo;
    for(auto p:host.points) {lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
    // Fluent's conservative length scale: min(domain volume^(1/3), largest extent).
    length=std::min(std::cbrt(volume),std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z}));
    update_time(bc_speed);
    dev.fill(momentum[0]->diagonal,1);dev.fill(momentum[0]->consistent_inverse,0);
    update_boundaries();update_pressure();update_flux();
  }

  void apply_backflow_pressure() {
    kern::fv_backflow_pressure(dev,mesh.n_boundary_faces(),mesh.owner,total_pressure,mass,
      velocity.x,velocity.y,velocity.z,boundary_pressure,mesh.n_internal_faces,float(settings.density));
  }
  void update_momentum_diagonal() {
    kern::fv_momentum_diagonal(dev,mesh.n_cells,momentum[0]->diagonal,momentum[1]->diagonal,momentum[2]->diagonal,
      momentum[0]->consistent_inverse,momentum[1]->consistent_inverse,momentum[2]->consistent_inverse,
      physical_diagonal,consistent_inverse);
  }
  void turbulence_boundaries() {
    kern::fv_turbulence_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,turbulence_spec,
      turbulence_a,turbulence_b,mass,boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,
      kinetic,omega,sst->fields().wall_face_omega,boundary_k,boundary_omega,turbulence_kind,
      mesh.n_internal_faces,float(settings.density),float(settings.viscosity));
  }
  void evaluate_turbulence() {
    sst->update({velocity,{strain_gradient[0],strain_gradient[1],strain_gradient[2]},rho_cell,mu_cell,
      kinetic,omega,k_equation->grad,omega_equation->grad,boundary_velocity});
  }
  void update_turbulence() {
    if(!sst) return;
    vector_gradient(dev,mesh,stencil,strain_lsq,velocity,boundary_velocity,
      strain_base[0],strain_base[1],strain_base[2]);
    const Buffer<float>* q[]{&velocity.x,&velocity.y,&velocity.z};
    const Buffer<float>* bc[]{&boundary_velocity.x,&boundary_velocity.y,&boundary_velocity.z};
    for(int k=0;k<3;++k)
      jump_corrected_gradient(dev,mesh,stencil,jump_fit,*q[k],*bc[k],strain_base[k],strain_gradient[k]);
    evaluate_turbulence(); // bootstrap gradient-independent wall-face omega
    turbulence_boundaries();
    gradient(dev,mesh,stencil,k_lsq,kinetic,boundary_k,k_equation->grad);
    gradient(dev,mesh,stencil,omega_lsq,omega,boundary_omega,omega_equation->grad);
    evaluate_turbulence();
    kern::fv_effective_viscosity(dev,mesh.n_faces,mesh.owner,mesh.neighbour,boundary_kind,
      mesh.weight,sst->fields().mu_t,viscosity,mesh.n_internal_faces,float(settings.viscosity));
  }
  void add_transpose_stress() {
    if(!sst) return;
    reconstruct_velocity();
    const auto& u=momentum[0]->grad;const auto& v=momentum[1]->grad;const auto& w=momentum[2]->grad;
    kern::fv_transpose_stress(dev,mesh.n_faces,mesh.owner,mesh.neighbour,boundary_kind,mesh.weight,viscosity,
      mesh.area.x,mesh.area.y,mesh.area.z,u.x,u.y,u.z,v.x,v.y,v.z,w.x,w.y,w.z,
      transpose_stress.x,transpose_stress.y,transpose_stress.z,mesh.n_internal_faces);
    kern::fv_stress_source(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      mesh.volume,transpose_stress.x,transpose_stress.y,transpose_stress.z,source.x,source.y,source.z);
  }
  TransportTerms wall_terms(int component) {
    kern::fv_boundary_momentum(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,boundary_kind,
      mesh.area.x,mesh.area.y,mesh.area.z,sst?sst->fields().wall_shear_coefficient:prescribed_pressure,
      viscosity,stencil.diffusion,
      velocity.x,velocity.y,velocity.z,boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,
      wall_diagonal,wall_rhs,mesh.n_internal_faces,component,sst?1:0);
    return {nullptr,&wall_diagonal,&wall_rhs};
  }
  void assemble_turbulence(bool is_omega) {
    const auto& opt=*settings.sst;const auto& f=sst->fields();
    auto& eq=is_omega?*omega_equation:*k_equation;
    TransportSettings t;t.density=float(settings.density);t.pseudo_time_step=time_step;t.limit=settings.limit;
    t.convection=is_omega?opt.omega_convection:opt.k_convection;
    t.relaxation=float(is_omega?opt.omega_relaxation:opt.k_relaxation);
    eq.assemble_reconstructed(is_omega?omega:kinetic,is_omega?boundary_omega:boundary_k,
      turbulence_kind,mass,is_omega?f.face_diffusivity_omega:f.face_diffusivity_k,
      is_omega?f.omega_source:f.k_source,t,{is_omega?&f.omega_sink:&f.k_sink});
    if(is_omega) kern::fv_fixed_rows(dev,mesh.n_cells,eq.matrix.row_ptr,eq.matrix.diag,f.fixed_omega_mask,
      f.fixed_omega,eq.matrix.val,eq.rhs);
  }
  void solve_turbulence() {
    if(!sst) return;
    update_boundaries();update_turbulence();assemble_turbulence(false);
    report.turbulence[0]=k_equation->solve(kinetic,settings.sst->linear);
    kern::fv_floor(dev,mesh.n_cells,kinetic,settings.sst->model.k_floor);
    update_turbulence();assemble_turbulence(true);
    report.turbulence[1]=omega_equation->solve(omega,settings.sst->linear);
    kern::fv_floor(dev,mesh.n_cells,omega,settings.sst->model.omega_floor);
  }
  void monitor_turbulence() {
    if(!sst) return;
    const auto measure=[&](bool is_omega) {
      assemble_turbulence(is_omega);auto& eq=is_omega?*omega_equation:*k_equation;
      const auto& q=is_omega?omega:kinetic;
      residual(dev,eq.matrix,q,eq.rhs,residual_u);
      kern::fv_scalar_monitor(dev,mesh.n_cells,q,eq.diagonal,residual_u,
        is_omega?sst->fields().fixed_omega_mask:zero_mask,mesh.volume,absolute_residual.x,scale,
        fluent_scale.x,float(settings.density),time_step);
      const auto raw=dev.sum(absolute_residual.x,mesh.n_cells),den=dev.sum(scale,mesh.n_cells);
      const auto fluent_den=dev.sum(fluent_scale.x,mesh.n_cells);
      const auto scaled=[&](double d) {return d>0?raw/d:(raw==0?0:std::numeric_limits<double>::infinity());};
      return ResidualValue{raw,scaled(den),scaled(fluent_den),fluent_den};
    };
    report.k=measure(false);report.omega=measure(true);
  }
  void update_time(double maximum_speed) {
    time_step=settings.pseudo_time_step>0?float(settings.pseudo_time_step):
      float(.3*length*settings.pseudo_time_scale/std::max({.1*bc_speed,maximum_speed,pressure_speed,1e-10}));
  }
  void update_boundaries() {
    kern::fv_velocity_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,
      mesh.area.x,mesh.area.y,mesh.area.z,velocity.x,velocity.y,velocity.z,
      prescribed_velocity.x,prescribed_velocity.y,prescribed_velocity.z,mass,
      boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,scalar_kind,mesh.n_internal_faces,float(settings.density));
  }
  void reconstruct_velocity() {
    vector_gradient(dev,mesh,stencil,velocity_lsq,velocity,boundary_velocity,
      momentum[0]->grad,momentum[1]->grad,momentum[2]->grad);
  }
  void reconstruct_pressure(const Buffer<float>& field,const Buffer<float>& boundary,DeviceVec3& grad) {
    gradient(dev,mesh,stencil,pressure_lsq,field,boundary,grad);
  }
  void update_pressure() {
    kern::fv_pressure_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,
      mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,pressure,grad_p.x,grad_p.y,grad_p.z,
      prescribed_pressure,boundary_pressure,mesh.n_internal_faces);
    apply_backflow_pressure();
    reconstruct_pressure(pressure,boundary_pressure,grad_p);
    kern::fv_pressure_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,
      mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,pressure,grad_p.x,grad_p.y,grad_p.z,
      prescribed_pressure,boundary_pressure,mesh.n_internal_faces);
    apply_backflow_pressure();
    kern::fv_pressure_faces(dev,mesh.n_faces,mesh.owner,mesh.neighbour,
      mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,
      mesh.neighbour_to_face.x,mesh.neighbour_to_face.y,mesh.neighbour_to_face.z,
      pressure,grad_p.x,grad_p.y,grad_p.z,boundary_pressure,face_pressure,pressure_jump,mesh.n_internal_faces);
  }
  void update_flux() {
    reconstruct_velocity();
    const Buffer<float>* q[]{&velocity.x,&velocity.y,&velocity.z};
    const Buffer<float>* bc[]{&boundary_velocity.x,&boundary_velocity.y,&boundary_velocity.z};
    for(int k=0;k<3;++k) {
      auto& eq=*momentum[k];auto& g=limited_gradient[k];
      if(settings.limit) limit_gradient(dev,mesh,stencil,*q[k],*bc[k],eq.grad,eq.limiter);
      else dev.fill(eq.limiter,1);
      kern::fv_scale_gradient(dev,mesh.n_cells,eq.grad.x,eq.grad.y,eq.grad.z,eq.limiter,g.x,g.y,g.z);
    }
    const auto& gu=limited_gradient[0];const auto& gv=limited_gradient[1];const auto& gw=limited_gradient[2];
    kern::fv_momentum_faces(dev,mesh.n_faces,mesh.owner,mesh.neighbour,boundary_kind,
      mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,
      mesh.neighbour_to_face.x,mesh.neighbour_to_face.y,mesh.neighbour_to_face.z,
      velocity.x,velocity.y,velocity.z,gu.x,gu.y,gu.z,gv.x,gv.y,gv.z,gw.x,gw.y,gw.z,
      physical_diagonal,boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,
      face_velocity.x,face_velocity.y,face_velocity.z,mesh.n_internal_faces);
    kern::fv_rhie_chow(dev,mesh.n_faces,mesh.owner,mesh.neighbour,boundary_kind,
      mesh.area.x,mesh.area.y,mesh.area.z,stencil.diffusion,mesh.volume,
      physical_diagonal,consistent_inverse,face_velocity.x,face_velocity.y,face_velocity.z,
      pressure_jump,
      mass,coefficient,mesh.n_internal_faces,float(settings.density));
  }

  void step() {
    update_boundaries();update_pressure();
    update_turbulence();
    kern::fv_momentum_source(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      mesh.area.x,mesh.area.y,mesh.area.z,face_pressure,mesh.volume,body.x,body.y,body.z,
      source.x,source.y,source.z);
    add_transpose_stress();
    TransportSettings transport;
    transport.convection=settings.convection;transport.limit=settings.limit;
    transport.density=float(settings.density);transport.pseudo_time_step=time_step;
    transport.relaxation=float(settings.momentum_relaxation);
    Buffer<float>* q[]{&velocity.x,&velocity.y,&velocity.z};
    const Buffer<float>* bc[]{&boundary_velocity.x,&boundary_velocity.y,&boundary_velocity.z};
    const Buffer<float>* src[]{&source.x,&source.y,&source.z};
    reconstruct_velocity();
    for(int k=0;k<3;++k) {
      momentum[k]->assemble_reconstructed(*q[k],*bc[k],scalar_kind,mass,viscosity,*src[k],transport,wall_terms(k));
      report.momentum[k]=momentum[k]->solve(*q[k],settings.momentum_linear);
    }
    update_momentum_diagonal();
    update_boundaries();update_flux();
    kern::fv_continuity(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,mass,provisional_continuity);
    kern::fv_pressure_assemble(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      stencil.adjacent,stencil.matrix_entry,pressure_matrix->row_ptr,pressure_matrix->diag,coefficient,mass,
      pressure_matrix->val,pressure_rhs,reference);
    pressure_solver->refresh(dev,*pressure_matrix);dev.fill(correction,0);
    report.pressure=pressure_solver->solve(dev,*pressure_matrix,correction,pressure_rhs,settings.pressure_linear);
    kern::fv_correct_flux(dev,mesh.n_faces,mesh.owner,mesh.neighbour,correction,coefficient,mass,mesh.n_internal_faces);
    // Pressure correction has homogeneous Neumann on fixed-velocity boundaries
    // and homogeneous Dirichlet at pressure outlets.
    dev.fill(grad_correction.x,0);dev.fill(grad_correction.y,0);dev.fill(grad_correction.z,0);
    kern::fv_pressure_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,
      mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,correction,
      grad_correction.x,grad_correction.y,grad_correction.z,correction_boundary,
      boundary_pressure,mesh.n_internal_faces);
    reconstruct_pressure(correction,boundary_pressure,grad_correction);
    kern::fv_correct_velocity(dev,mesh.n_cells,mesh.volume,consistent_inverse,physical_diagonal,
      grad_correction.x,grad_correction.y,grad_correction.z,correction,
      velocity.x,velocity.y,velocity.z,pressure,float(settings.pressure_relaxation));
    solve_turbulence();
    ++report.iteration;
    if(report.iteration<=5 || report.iteration%settings.monitor_interval==0) monitor();
  }

  IterationReport monitor() {
    update_boundaries();update_pressure();
    update_turbulence();
    kern::fv_momentum_source(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      mesh.area.x,mesh.area.y,mesh.area.z,face_pressure,mesh.volume,body.x,body.y,body.z,
      source.x,source.y,source.z);
    add_transpose_stress();
    TransportSettings transport;
    transport.convection=settings.convection;transport.limit=settings.limit;
    transport.density=float(settings.density);transport.pseudo_time_step=time_step;
    transport.relaxation=float(settings.momentum_relaxation);
    reconstruct_velocity();
    // The components share one CSR value buffer, so each residual must follow
    // its own assembly.
    momentum[0]->assemble_reconstructed(velocity.x,boundary_velocity.x,scalar_kind,mass,viscosity,source.x,transport,wall_terms(0));
    residual(dev,momentum[0]->matrix,velocity.x,momentum[0]->rhs,residual_u);
    momentum[1]->assemble_reconstructed(velocity.y,boundary_velocity.y,scalar_kind,mass,viscosity,source.y,transport,wall_terms(1));
    residual(dev,momentum[1]->matrix,velocity.y,momentum[1]->rhs,residual_v);
    momentum[2]->assemble_reconstructed(velocity.z,boundary_velocity.z,scalar_kind,mass,viscosity,source.z,transport,wall_terms(2));
    residual(dev,momentum[2]->matrix,velocity.z,momentum[2]->rhs,residual_w);
    update_momentum_diagonal();
    kern::fv_monitor_cell(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      mass,velocity.x,velocity.y,velocity.z,physical_diagonal,residual_u,residual_v,residual_w,
      continuity,scale,speed,absolute_residual.x,absolute_residual.y,absolute_residual.z);
    kern::fv_fluent_momentum_scale(dev,mesh.n_cells,momentum[0]->diagonal,momentum[1]->diagonal,momentum[2]->diagonal,
      mesh.volume,speed,fluent_scale.x,fluent_scale.y,fluent_scale.z,float(settings.density),time_step);
    const double denominator=dev.sum(scale,mesh.n_cells);
    const double raw=dev.sum(provisional_continuity,mesh.n_cells);
    report.corrected_mass_imbalance=dev.sum(continuity,mesh.n_cells);
    if(report.iteration>0 && report.iteration<=5) first_five_continuity[report.iteration-1]=raw;
    continuity_scale=*std::max_element(std::begin(first_five_continuity),std::end(first_five_continuity));
    const auto scaled=[](double r,double d) {return d>0?r/d:(r==0?0:std::numeric_limits<double>::infinity());};
    report.u.raw=dev.sum(absolute_residual.x,mesh.n_cells);report.u.scaled=scaled(report.u.raw,denominator);
    report.v.raw=dev.sum(absolute_residual.y,mesh.n_cells);report.v.scaled=scaled(report.v.raw,denominator);
    report.w.raw=dev.sum(absolute_residual.z,mesh.n_cells);report.w.scaled=scaled(report.w.raw,denominator);
    ResidualValue* values[]{&report.u,&report.v,&report.w};
    const Buffer<float>* scales[]{&fluent_scale.x,&fluent_scale.y,&fluent_scale.z};
    for(int k=0;k<3;++k) {
      values[k]->fluent_normalization=dev.sum(*scales[k],mesh.n_cells);
      values[k]->fluent_scaled=scaled(values[k]->raw,values[k]->fluent_normalization);
    }
    report.continuity={raw,scaled(raw,continuity_scale),scaled(raw,continuity_scale),continuity_scale};
    report.maximum_velocity=dev.max_abs(speed,mesh.n_cells);
    report.pseudo_time_step=time_step;
    monitor_turbulence();
    update_time(report.maximum_velocity);
    report.residual_iteration=report.iteration;
    return report;
  }
};

SteadySolver::SteadySolver(Device& dev,PolyMesh mesh,std::vector<PatchCondition> boundaries,SolverSettings settings)
    :impl_(std::make_unique<Impl>(dev,std::move(mesh),std::move(boundaries),settings)) {}
SteadySolver::~SteadySolver()=default;
void SteadySolver::initialize(const std::function<Vec3d(Vec3d)>& velocity,const std::function<double(Vec3d)>& pressure) {
  auto& s=*impl_;KUTTA_REQUIRE(!s.sst,"FV SST initialization requires initialize_fields with k and omega");
  std::vector<Vec3d> u;std::vector<float> p;
  const double ref=s.reference>=0?pressure(s.geom.cell_centroid[s.reference]):0;
  for(auto x:s.geom.cell_centroid) {u.push_back(velocity(x));const double value=pressure(x)-ref;
    KUTTA_REQUIRE(std::isfinite(value),"FV nonfinite pressure initialization");p.push_back(float(value));}
  upload_vectors(s.dev,s.velocity,u);s.dev.upload<float>(s.pressure,p);
  s.update_boundaries();s.update_pressure();
  // Initialize a consistent advective mass flux before the first assembly.
  s.dev.fill(s.physical_diagonal,1);s.dev.fill(s.consistent_inverse,0);
  s.update_flux();
}
void SteadySolver::initialize_fields(std::span<const Vec3d> velocity,std::span<const double> pressure,
    std::span<const double> k,std::span<const double> omega) {
  auto& s=*impl_;const auto n=std::size_t(s.mesh.n_cells);
  KUTTA_REQUIRE(velocity.size()==n && pressure.size()==n,"FV initialization field size mismatch");
  KUTTA_REQUIRE(s.sst?(k.size()==n && omega.size()==n):(k.empty() && omega.empty()),"FV turbulence initialization size mismatch");
  std::vector<Vec3d> u(n);std::vector<float> p(n),kinetic(n),dissipation(n);
  const double ref=s.reference>=0?pressure[s.permutation.cells_new_to_old[s.reference]]:0;
  for(std::size_t c=0;c<n;++c) {
    const auto old=s.permutation.cells_new_to_old[c];u[c]=velocity[old];p[c]=float(pressure[old]-ref);
    KUTTA_REQUIRE(std::isfinite(p[c]),"FV nonfinite pressure initialization");
    if(s.sst) {
      KUTTA_REQUIRE(std::isfinite(k[old]) && k[old]>=0 && std::isfinite(omega[old]) && omega[old]>0,"FV invalid turbulence initialization");
      kinetic[c]=float(k[old]);dissipation[c]=float(omega[old]);
      KUTTA_REQUIRE(std::isfinite(kinetic[c]) && std::isfinite(dissipation[c]),"FV turbulence initialization exceeds float storage");
    }
  }
  upload_vectors(s.dev,s.velocity,u);s.dev.upload<float>(s.pressure,p);
  if(s.sst) {s.dev.upload<float>(s.kinetic,kinetic);s.dev.upload<float>(s.omega,dissipation);}
  s.dev.fill(s.physical_diagonal,1);s.dev.fill(s.consistent_inverse,0);
  s.update_boundaries();s.update_pressure();s.update_flux();s.update_turbulence();
  // Establish the physical momentum coefficient and its constitutive flux at
  // the supplied fields without taking a nonlinear step. In particular a
  // saved-field residual must not use the placeholder unit diagonal.
  for(int pass=0;pass<3;++pass) {
    s.monitor();s.update_flux();
    kern::fv_continuity(s.dev,s.mesh.n_cells,s.mesh.cell_face_offsets,s.mesh.cell_faces,
      s.mesh.cell_face_sign,s.mass,s.provisional_continuity);
    s.update_boundaries();s.update_pressure();
  }
}
void SteadySolver::initialize_hybrid(const HybridInitialization& h) {
  auto& s=*impl_;auto& dev=s.dev;
  KUTTA_REQUIRE(h.iterations>0,"FV hybrid initialization needs positive iteration count");
  const auto n=s.mesh.n_cells,nb=s.mesh.n_boundary_faces();
  const auto volume=std::accumulate(s.geom.cell_volume.begin(),s.geom.cell_volume.end(),0.);
  std::vector<Index> kind(nb,1),pkind(nb,1);
  std::vector<float> value(nb,0),pvalue(nb,0),normal(nb,0),flux(nb,0),source(n,0),seed(n,float(h.initial_scalar[0]));
  double inlet_area=0,outlet_area=0,pressure_sum=0,pressure_area=0,speed_sum=0;
  Vec3d direction,outlet_center;
  // Turbulence recipes: k=k0+k2*|u|^2, omega=w0+w1*|u|+w2*|u|^2.
  std::array<double,5> recipe{};
  const auto add_recipe=[&](const TurbulenceCondition& spec,double area) {
    if(const auto* v=std::get_if<TurbulenceValues>(&spec)) {recipe[0]+=area*v->k;recipe[2]+=area*v->omega;}
    else if(const auto* v=std::get_if<IntensityLength>(&spec)) {
      recipe[1]+=area*1.5*v->intensity*v->intensity;
      recipe[3]+=area*std::sqrt(1.5)*v->intensity/(std::pow(.09,.25)*v->length);
    } else {
      const auto& ratio=std::get<IntensityViscosityRatio>(spec);
      recipe[1]+=area*1.5*ratio.intensity*ratio.intensity;
      recipe[4]+=area*s.settings.density*1.5*ratio.intensity*ratio.intensity/(s.settings.viscosity*ratio.ratio);
    }
  };
  for(const auto& patch:s.host.patches) {
    const auto& bc=s.boundaries.at(patch.name);
    const auto initial_pressure=std::find_if(h.inlet_pressure.begin(),h.inlet_pressure.end(),[&](const auto& p){return p.first==patch.name;});
    for(Index f=patch.start;f<patch.start+patch.size;++f) {
      const auto b=f-s.mesh.n_internal_faces,c=s.host.owner[f];
      const auto area_vector=s.geom.face_area[f];const auto area=norm(area_vector);
      if(const auto* in=std::get_if<VelocityInlet>(&bc)) {
        const auto u=evaluate(in->velocity,s.geom.face_centroid[f]);
        flux[b]=float(dot(u,area_vector));source[c]+=flux[b]/float(s.geom.cell_volume[c]);
        normal[b]=float(dot(u,area_vector)/area*dot(s.geom.face_centroid[f]-s.geom.cell_centroid[c],area_vector)/area);
        inlet_area+=area;direction=direction+u*area;speed_sum+=norm(u)*area;
        if(s.sst) add_recipe(*in->turbulence,area);
        if(initial_pressure!=h.inlet_pressure.end()) {pkind[b]=0;pvalue[b]=float(initial_pressure->second);pressure_sum+=area*initial_pressure->second;pressure_area+=area;}
      }
      if(const auto* out=std::get_if<PressureOutlet>(&bc)) {
        kind[b]=0;pkind[b]=0;pvalue[b]=float(out->pressure);outlet_area+=area;
        outlet_center=outlet_center+s.geom.face_centroid[f]*area;
        pressure_sum+=area*out->pressure;pressure_area+=area;
      }
    }
  }
  KUTTA_REQUIRE(inlet_area==0 || outlet_area>0,"FV hybrid velocity inlets require a pressure outlet");
  if(inlet_area>0) {direction=direction/inlet_area;speed_sum/=inlet_area;for(auto& r:recipe) r/=inlet_area;}
  else if(s.sst) {add_recipe(IntensityViscosityRatio{h.turbulence_intensity,h.turbulence_viscosity_ratio},1);}
  if(outlet_area>0) outlet_center=outlet_center/outlet_area;
  std::fill(source.begin(),source.end(),0);
  for(Index f=s.mesh.n_internal_faces;f<s.mesh.n_faces;++f) {
    const auto b=f-s.mesh.n_internal_faces,c=s.host.owner[f];const auto area=s.geom.face_area[f];
    const auto r=s.geom.face_centroid[f]-s.geom.cell_centroid[c];
    if(kind[b]==0) value[b]=float(-dot(direction,s.geom.face_centroid[f]-outlet_center));
    else {
      const double q=flux[b]-dot(direction,area);
      source[c]+=float(q/s.geom.cell_volume[c]);
      normal[b]=float(q*dot(r,area)/dot(area,area));
    }
  }
  if(!h.linear_potential) for(Index c=0;c<n;++c)
    seed[c]=float(h.initial_scalar[0]-dot(direction,s.geom.cell_centroid[c]-outlet_center));
  auto potential=dev.upload<float>(seed),zero=dev.alloc<float>(s.mesh.n_faces),one=dev.alloc<float>(s.mesh.n_faces);
  dev.fill(zero,0);dev.fill(one,1);
  auto boundary=dev.alloc<float>(nb),old=dev.alloc<float>(n);
  TransportEquation equation(dev,s.mesh,s.stencil,s.host,s.geom);
  auto solve_scalar=[&](Buffer<float>& field,const std::vector<Index>& kinds,const std::vector<float>& values,
      const std::vector<float>& deltas,const std::vector<float>& sources,int scalar) {
    auto bk=dev.upload<Index>(kinds);auto bv=dev.upload<float>(values),bd=dev.upload<float>(deltas),rhs=dev.upload<float>(sources);
    std::vector<BoundaryObservation> observations(nb);
    for(Index b=0;b<nb;++b) observations[b]=kinds[b]==0?BoundaryObservation::FaceValue:BoundaryObservation::NormalFootValue;
    const auto lsq=boundary_least_squares(dev,s.host,s.geom,observations);
    TransportSettings t;t.convection=ConvectionScheme::FirstOrderUpwind;t.limit=false;t.relaxation=float(h.implicit_relaxation[scalar]);
    for(int i=0;i<h.iterations && h.solve[scalar];++i) {
      dev.copy(old,field);
      kern::fv_hybrid_boundary(dev,nb,s.mesh.owner,bk,field,bv,bd,boundary,s.mesh.n_internal_faces);
      gradient(dev,s.mesh,s.stencil,lsq,field,boundary,equation.grad);
      equation.assemble_reconstructed(field,boundary,bk,zero,one,rhs,t);
      s.pressure_solver->refresh(dev,equation.matrix);
      s.pressure_solver->solve(dev,equation.matrix,field,equation.rhs,{h.tolerance[scalar],1e-12,300});
      kern::axpbyz(dev,n,field,field,old,float(h.explicit_relaxation[scalar]),float(1-h.explicit_relaxation[scalar]));
    }
    kern::fv_hybrid_boundary(dev,nb,s.mesh.owner,bk,field,bv,bd,boundary,s.mesh.n_internal_faces);
    gradient(dev,s.mesh,s.stencil,lsq,field,boundary,equation.grad);
  };
  if(outlet_area>0) solve_scalar(potential,kind,value,normal,source,0);
  else {dev.fill(equation.grad.x,0);dev.fill(equation.grad.y,0);dev.fill(equation.grad.z,0);}
  kern::fv_hybrid_velocity(dev,n,equation.grad.x,equation.grad.y,equation.grad.z,
    s.velocity.x,s.velocity.y,s.velocity.z,float(direction.x),float(direction.y),float(direction.z),float(speed_sum),h.constant_speed?1:0);
  // Preserve the potential operator's conservative face flux when setting up
  // the first momentum equation. The gradient velocity alone is not solenoidal.
  auto potential_kind=dev.upload<Index>(kind);auto prescribed_flux=dev.upload<float>(flux);
  auto potential_flux=dev.alloc<float>(s.mesh.n_faces);
  TransportSettings t;t.convection=ConvectionScheme::FirstOrderUpwind;t.limit=false;
  auto source_device=dev.upload<float>(source);auto potential_value=dev.upload<float>(value);
  equation.assemble_reconstructed(potential,boundary,potential_kind,zero,one,source_device,t);
  kern::fv_hybrid_flux(dev,s.mesh.n_faces,s.mesh.owner,s.mesh.neighbour,potential_kind,potential,s.stencil.diffusion,
    equation.face_correction,prescribed_flux,potential_value,s.mesh.area.x,s.mesh.area.y,s.mesh.area.z,
    potential_flux,s.mesh.n_internal_faces,float(s.settings.density),float(direction.x),float(direction.y),float(direction.z));
  dev.fill(s.pressure,float(h.solve[1] && pressure_area>0?pressure_sum/pressure_area:h.initial_scalar[1]));
  if(h.solve[1] && !h.inlet_pressure.empty() && outlet_area>0)
    solve_scalar(s.pressure,pkind,pvalue,std::vector<float>(nb,0),std::vector<float>(n,0),1);
  if(s.sst) {
    kern::fv_hybrid_turbulence(dev,n,s.velocity.x,s.velocity.y,s.velocity.z,s.mesh.volume,s.kinetic,s.omega,old,s.scale,
      float(recipe[0]),float(recipe[1]),float(recipe[2]),float(recipe[3]),float(recipe[4]),
      s.settings.sst->model.k_floor,s.settings.sst->model.omega_floor);
    if(h.averaged_turbulence) {const auto k=dev.sum(old,n)/volume,w=dev.sum(s.scale,n)/volume;dev.fill(s.kinetic,float(k));dev.fill(s.omega,float(w));}
  }
  s.dev.fill(s.physical_diagonal,1);s.dev.fill(s.consistent_inverse,0);
  s.update_boundaries();s.update_pressure();s.update_flux();s.update_turbulence();
  for(int pass=0;pass<3;++pass) {s.monitor();s.update_flux();s.update_boundaries();s.update_pressure();}
  if(!h.constant_speed && h.solve[0]) dev.copy(s.mass,potential_flux);
  kern::fv_continuity(dev,n,s.mesh.cell_face_offsets,s.mesh.cell_faces,s.mesh.cell_face_sign,s.mass,s.provisional_continuity);
}
void SteadySolver::set_source(const std::function<Vec3d(Vec3d)>& source) {
  std::vector<Vec3d> values;for(auto p:impl_->geom.cell_centroid) values.push_back(source(p));
  upload_vectors(impl_->dev,impl_->body,values);
}
void SteadySolver::iterate(int count) {KUTTA_REQUIRE(count>=0,"FV negative iteration count");while(count-->0) impl_->step();}
IterationReport SteadySolver::monitor() {return impl_->monitor();}
ForceReport SteadySolver::forces(std::span<const std::string> patches,double reference_pressure) {
  auto& s=*impl_;
  KUTTA_REQUIRE(std::isfinite(reference_pressure),"FV nonfinite force reference pressure");
  std::set<std::string> selected_names(patches.begin(),patches.end());
  KUTTA_REQUIRE(selected_names.size()==patches.size(),"FV duplicate force patch");
  std::vector<float> selected(s.mesh.n_boundary_faces(),0);
  Vec3d area;
  for(const auto& name:selected_names) {
    auto it=std::find_if(s.host.patches.begin(),s.host.patches.end(),[&](const Patch& p){return p.name==name;});
    KUTTA_REQUIRE(it!=s.host.patches.end() && it->kind==PatchKind::Wall,"FV force patch must name a wall: "+name);
    for(Index f=it->start;f<it->start+it->size;++f) {
      selected[f-s.mesh.n_internal_faces]=1;
      area=area+s.geom.face_area[f];
    }
  }
  s.update_boundaries();s.update_pressure();
  s.reconstruct_velocity();
  auto mask=s.dev.upload<float>(selected);
  auto pressure=vector_field(s.dev,selected.size()),viscous=vector_field(s.dev,selected.size());
  const auto& m=s.mesh;const auto& gu=s.momentum[0]->grad;const auto& gv=s.momentum[1]->grad;const auto& gw=s.momentum[2]->grad;
  kern::fv_wall_pressure_force(s.dev,selected.size(),mask,m.area.x,m.area.y,m.area.z,s.face_pressure,
    pressure.x,pressure.y,pressure.z,m.n_internal_faces,float(reference_pressure));
  if(s.sst) {
    s.update_turbulence();const auto& t=s.sst->fields();
    kern::fv_sst_force(s.dev,selected.size(),mask,m.area.x,m.area.y,m.area.z,
      t.wall_shear.x,t.wall_shear.y,t.wall_shear.z,viscous.x,viscous.y,viscous.z,m.n_internal_faces);
  } else kern::fv_wall_viscous_force(s.dev,selected.size(),m.owner,mask,m.area.x,m.area.y,m.area.z,
    m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,s.velocity.x,s.velocity.y,s.velocity.z,
    gu.x,gu.y,gu.z,gv.x,gv.y,gv.z,gw.x,gw.y,gw.z,
    s.boundary_velocity.x,s.boundary_velocity.y,s.boundary_velocity.z,viscous.x,viscous.y,viscous.z,
    m.n_internal_faces,float(s.settings.viscosity));
  return {{s.dev.sum(pressure.x,selected.size()),s.dev.sum(pressure.y,selected.size()),s.dev.sum(pressure.z,selected.size())},
    {s.dev.sum(viscous.x,selected.size()),s.dev.sum(viscous.y,selected.size()),s.dev.sum(viscous.z,selected.size())},area,reference_pressure};
}
const PolyMesh& SteadySolver::mesh() const {return impl_->host;}
const MeshGeometry& SteadySolver::geometry() const {return impl_->geom;}
const DeviceMesh& SteadySolver::device_mesh() const {return impl_->mesh;}
const DeviceVec3& SteadySolver::velocity() const {return impl_->velocity;}
const Buffer<float>& SteadySolver::pressure() const {return impl_->pressure;}
const Buffer<float>& SteadySolver::mass_flux() const {return impl_->mass;}
const Buffer<float>& SteadySolver::kinetic_energy() const {KUTTA_REQUIRE(impl_->sst,"FV SST disabled");return impl_->kinetic;}
const Buffer<float>& SteadySolver::specific_dissipation() const {KUTTA_REQUIRE(impl_->sst,"FV SST disabled");return impl_->omega;}
const SstFields& SteadySolver::turbulence() const {KUTTA_REQUIRE(impl_->sst,"FV SST disabled");return impl_->sst->fields();}
const MeshPermutation& SteadySolver::permutation() const {return impl_->permutation;}
const Buffer<float>& SteadySolver::boundary_pressure() const {return impl_->boundary_pressure;}
const IterationReport& SteadySolver::report() const {return impl_->report;}
Vec3d ForceReport::force() const {return pressure+viscous;}
double ForceReport::coefficient(Vec3d direction,const ReferenceValues& ref) const {
  KUTTA_REQUIRE(ref.density>0 && ref.velocity>0 && ref.area>0,"FV positive force references required");
  const auto adjusted=force()+area_vector*(reference_pressure-ref.pressure);
  return dot(adjusted,direction)/(.5*ref.density*ref.velocity*ref.velocity*ref.area);
}
}  // namespace kutta::fv
