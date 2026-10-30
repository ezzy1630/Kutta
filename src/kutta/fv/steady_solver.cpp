#include "kutta/fv/steady_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>

#include "kernels.hpp"
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
  DeviceVec3 velocity,body,boundary_velocity,prescribed_velocity,source,face_velocity;
  Buffer<float> pressure,mass,viscosity,boundary_pressure,prescribed_pressure;
  Buffer<Index> boundary_kind,scalar_kind;
  DeviceVec3 grad_p;
  std::unique_ptr<TransportEquation> momentum[3];
  std::shared_ptr<CsrMatrix> pressure_matrix;
  std::unique_ptr<LinearSolver> pressure_solver;
  Buffer<float> pressure_rhs,correction,correction_boundary,coefficient,face_pressure,pressure_jump;
  DeviceVec3 grad_correction;
  Buffer<float> residual_u,residual_v,residual_w,continuity,provisional_continuity,scale,speed;
  DeviceVec3 absolute_residual;
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
    renumber_cells_rcm(host);geom=compute_geometry(host);
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
    grad_p=vector_field(dev,n);grad_correction=vector_field(dev,n);absolute_residual=vector_field(dev,n);
    pressure=dev.alloc<float>(n);dev.fill(pressure,0);
    mass=dev.alloc<float>(nf);dev.fill(mass,0);viscosity=dev.alloc<float>(nf);dev.fill(viscosity,float(settings.viscosity));
    boundary_pressure=dev.alloc<float>(nb);dev.fill(boundary_pressure,0);
    prescribed_pressure=dev.alloc<float>(nb);boundary_kind=dev.alloc<Index>(nb);scalar_kind=dev.alloc<Index>(nb);
    pressure_rhs=dev.alloc<float>(n);correction=dev.alloc<float>(n);dev.fill(correction,0);
    correction_boundary=dev.alloc<float>(nb);dev.fill(correction_boundary,0);
    coefficient=dev.alloc<float>(nf);face_pressure=dev.alloc<float>(nf);pressure_jump=dev.alloc<float>(nf);
    residual_u=dev.alloc<float>(n);residual_v=dev.alloc<float>(n);residual_w=dev.alloc<float>(n);
    continuity=dev.alloc<float>(n);provisional_continuity=dev.alloc<float>(n);dev.fill(provisional_continuity,0);scale=dev.alloc<float>(n);speed=dev.alloc<float>(n);
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
        pressures[b]=float(out->pressure);pmin=std::min(pmin,out->pressure);pmax=std::max(pmax,out->pressure);
      }
      if(std::holds_alternative<Symmetry>(bc)) kinds[b]=3;
      bc_speed=std::max(bc_speed,norm(velocities[b]));
    }
    pressure_speed=outlet?.1*std::sqrt((pmax-pmin)/settings.density):0;
    dev.upload<Index>(boundary_kind,kinds);dev.upload<float>(prescribed_pressure,pressures);
    upload_vectors(dev,prescribed_velocity,velocities);
    pressure_matrix=std::make_shared<CsrMatrix>(upload(dev,stencil.pattern));
    for(auto& eq:momentum) eq=std::make_unique<TransportEquation>(dev,mesh,stencil,host,geom,pressure_matrix);
    CsrMatrixHost initial{stencil.pattern,std::vector<double>(stencil.pattern.nnz(),0)};
    const auto cf=build_cell_faces(host);
    for(Index c=0;c<n;++c) {
      double ap=0;
      for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
        auto f=cf.faces[j];auto other=f<mesh.n_internal_faces?(host.owner[f]==c?host.neighbour[f]:host.owner[f]):-1;
        if(other<0 && kinds[f-mesh.n_internal_faces]!=2 && kinds[f-mesh.n_internal_faces]!=4) continue;
        auto d=(other>=0?geom.cell_centroid[other]:geom.face_centroid[f])-geom.cell_centroid[c];
        double a=dot(geom.face_area[f],geom.face_area[f])/std::abs(dot(geom.face_area[f],d));
        ap+=a;
        if(other>=0 && other!=reference && c!=reference) {
          auto index=std::lower_bound(initial.pattern.col.begin()+initial.pattern.row_ptr[c],
            initial.pattern.col.begin()+initial.pattern.row_ptr[c+1],other)-initial.pattern.col.begin();
          initial.val[index]-=a;
        }
      }
      initial.val[initial.pattern.diag[c]]=c==reference?1:ap;
    }
    pressure_solver=make_solver(pressure_solver_kind,dev,initial);
    const double volume=std::accumulate(geom.cell_volume.begin(),geom.cell_volume.end(),0.);
    Vec3d lo=host.points.front(),hi=lo;
    for(auto p:host.points) {lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
    // Fluent's conservative length scale: min(domain volume^(1/3), largest extent).
    length=std::min(std::cbrt(volume),std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z}));
    update_time(bc_speed);
    dev.fill(momentum[0]->diagonal,1);dev.fill(momentum[0]->consistent_inverse,0);
    update_boundaries();update_pressure();update_flux();
  }

  void update_time(double maximum_speed) {
    time_step=settings.pseudo_time_step>0?float(settings.pseudo_time_step):
      float(.3*length*settings.pseudo_time_scale/std::max({.1*bc_speed,maximum_speed,pressure_speed,1e-10}));
  }
  void update_boundaries() {
    auto& gu=momentum[0]->grad;auto& gv=momentum[1]->grad;auto& gw=momentum[2]->grad;
    kern::fv_velocity_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,
      mesh.area.x,mesh.area.y,mesh.area.z,mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,
      velocity.x,velocity.y,velocity.z,gu.x,gu.y,gu.z,gv.x,gv.y,gv.z,gw.x,gw.y,gw.z,
      prescribed_velocity.x,prescribed_velocity.y,prescribed_velocity.z,mass,
      boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,scalar_kind,mesh.n_internal_faces,float(settings.density));
  }
  void update_pressure() {
    kern::fv_pressure_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,
      mesh.area.x,mesh.area.y,mesh.area.z,mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,pressure,grad_p.x,grad_p.y,grad_p.z,
      prescribed_pressure,boundary_pressure,mesh.n_internal_faces);
    gradient(dev,mesh,stencil,pressure,boundary_pressure,grad_p);
    kern::fv_pressure_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,
      mesh.area.x,mesh.area.y,mesh.area.z,mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,pressure,grad_p.x,grad_p.y,grad_p.z,
      prescribed_pressure,boundary_pressure,mesh.n_internal_faces);
    kern::fv_pressure_faces(dev,mesh.n_faces,mesh.owner,mesh.neighbour,
      mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,
      mesh.neighbour_to_face.x,mesh.neighbour_to_face.y,mesh.neighbour_to_face.z,
      pressure,grad_p.x,grad_p.y,grad_p.z,boundary_pressure,face_pressure,pressure_jump,mesh.n_internal_faces);
  }
  void update_flux() {
    vector_gradient(dev,mesh,stencil,velocity,boundary_velocity,momentum[0]->grad,momentum[1]->grad,momentum[2]->grad);
    const auto& gu=momentum[0]->grad;const auto& gv=momentum[1]->grad;const auto& gw=momentum[2]->grad;
    kern::fv_momentum_faces(dev,mesh.n_faces,mesh.owner,mesh.neighbour,boundary_kind,
      mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,
      mesh.neighbour_to_face.x,mesh.neighbour_to_face.y,mesh.neighbour_to_face.z,
      velocity.x,velocity.y,velocity.z,gu.x,gu.y,gu.z,gv.x,gv.y,gv.z,gw.x,gw.y,gw.z,
      momentum[0]->diagonal,boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,
      face_velocity.x,face_velocity.y,face_velocity.z,mesh.n_internal_faces);
    kern::fv_rhie_chow(dev,mesh.n_faces,mesh.owner,mesh.neighbour,boundary_kind,
      mesh.area.x,mesh.area.y,mesh.area.z,stencil.diffusion,mesh.volume,
      momentum[0]->diagonal,momentum[0]->consistent_inverse,face_velocity.x,face_velocity.y,face_velocity.z,
      pressure_jump,
      mass,coefficient,mesh.n_internal_faces,float(settings.density));
  }

  void step() {
    update_boundaries();update_pressure();
    kern::fv_momentum_source(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      mesh.area.x,mesh.area.y,mesh.area.z,face_pressure,mesh.volume,body.x,body.y,body.z,
      source.x,source.y,source.z);
    TransportSettings transport;
    transport.convection=settings.convection;transport.limit=settings.limit;
    transport.density=float(settings.density);transport.pseudo_time_step=time_step;
    transport.relaxation=float(settings.momentum_relaxation);
    Buffer<float>* q[]{&velocity.x,&velocity.y,&velocity.z};
    const Buffer<float>* bc[]{&boundary_velocity.x,&boundary_velocity.y,&boundary_velocity.z};
    const Buffer<float>* src[]{&source.x,&source.y,&source.z};
    vector_gradient(dev,mesh,stencil,velocity,boundary_velocity,momentum[0]->grad,momentum[1]->grad,momentum[2]->grad);
    for(int k=0;k<3;++k) {
      momentum[k]->assemble_reconstructed(*q[k],*bc[k],scalar_kind,mass,viscosity,*src[k],transport);
      report.momentum[k]=momentum[k]->solve(*q[k],settings.momentum_linear);
    }
    update_boundaries();update_flux();
    kern::fv_continuity(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,mass,provisional_continuity);
    kern::fv_pressure_assemble(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      stencil.adjacent,stencil.matrix_entry,pressure_matrix->diag,coefficient,mass,
      pressure_matrix->val,pressure_rhs,reference);
    pressure_solver->refresh(dev,*pressure_matrix);dev.fill(correction,0);
    report.pressure=pressure_solver->solve(dev,*pressure_matrix,correction,pressure_rhs,settings.pressure_linear);
    kern::fv_correct_flux(dev,mesh.n_faces,mesh.owner,mesh.neighbour,correction,coefficient,mass,mesh.n_internal_faces);
    // Pressure correction has homogeneous Neumann on fixed-velocity boundaries
    // and homogeneous Dirichlet at pressure outlets.
    dev.fill(grad_correction.x,0);dev.fill(grad_correction.y,0);dev.fill(grad_correction.z,0);
    kern::fv_pressure_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,
      mesh.area.x,mesh.area.y,mesh.area.z,mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,correction,
      grad_correction.x,grad_correction.y,grad_correction.z,correction_boundary,
      boundary_pressure,mesh.n_internal_faces);
    gradient(dev,mesh,stencil,correction,boundary_pressure,grad_correction);
    kern::fv_correct_velocity(dev,mesh.n_cells,mesh.volume,momentum[0]->consistent_inverse,momentum[0]->diagonal,
      grad_correction.x,grad_correction.y,grad_correction.z,correction,
      velocity.x,velocity.y,velocity.z,pressure,float(settings.pressure_relaxation));
    ++report.iteration;
    if(report.iteration<=5 || report.iteration%settings.monitor_interval==0) monitor();
  }

  IterationReport monitor() {
    update_boundaries();update_pressure();
    kern::fv_momentum_source(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      mesh.area.x,mesh.area.y,mesh.area.z,face_pressure,mesh.volume,body.x,body.y,body.z,
      source.x,source.y,source.z);
    TransportSettings transport;
    transport.convection=settings.convection;transport.limit=settings.limit;
    transport.density=float(settings.density);transport.pseudo_time_step=time_step;
    transport.relaxation=float(settings.momentum_relaxation);
    vector_gradient(dev,mesh,stencil,velocity,boundary_velocity,momentum[0]->grad,momentum[1]->grad,momentum[2]->grad);
    // The components share one CSR value buffer, so each residual must follow
    // its own assembly.
    momentum[0]->assemble_reconstructed(velocity.x,boundary_velocity.x,scalar_kind,mass,viscosity,source.x,transport);
    residual(dev,momentum[0]->matrix,velocity.x,momentum[0]->rhs,residual_u);
    momentum[1]->assemble_reconstructed(velocity.y,boundary_velocity.y,scalar_kind,mass,viscosity,source.y,transport);
    residual(dev,momentum[1]->matrix,velocity.y,momentum[1]->rhs,residual_v);
    momentum[2]->assemble_reconstructed(velocity.z,boundary_velocity.z,scalar_kind,mass,viscosity,source.z,transport);
    residual(dev,momentum[2]->matrix,velocity.z,momentum[2]->rhs,residual_w);
    kern::fv_monitor_cell(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      mass,velocity.x,velocity.y,velocity.z,momentum[0]->diagonal,residual_u,residual_v,residual_w,
      continuity,scale,speed,absolute_residual.x,absolute_residual.y,absolute_residual.z);
    const double denominator=dev.sum(scale,mesh.n_cells);
    const double raw=dev.sum(provisional_continuity,mesh.n_cells);
    report.corrected_mass_imbalance=dev.sum(continuity,mesh.n_cells);
    if(report.iteration>0 && report.iteration<=5) first_five_continuity[report.iteration-1]=raw;
    continuity_scale=*std::max_element(std::begin(first_five_continuity),std::end(first_five_continuity));
    const auto scaled=[](double r,double d) {return d>0?r/d:(r==0?0:std::numeric_limits<double>::infinity());};
    report.u.raw=dev.sum(absolute_residual.x,mesh.n_cells);report.u.scaled=scaled(report.u.raw,denominator);
    report.v.raw=dev.sum(absolute_residual.y,mesh.n_cells);report.v.scaled=scaled(report.v.raw,denominator);
    report.w.raw=dev.sum(absolute_residual.z,mesh.n_cells);report.w.scaled=scaled(report.w.raw,denominator);
    report.continuity={raw,scaled(raw,continuity_scale)};
    report.maximum_velocity=dev.max_abs(speed,mesh.n_cells);update_time(report.maximum_velocity);
    report.pseudo_time_step=time_step;
    report.residual_iteration=report.iteration;
    return report;
  }
};

SteadySolver::SteadySolver(Device& dev,PolyMesh mesh,std::vector<PatchCondition> boundaries,SolverSettings settings)
    :impl_(std::make_unique<Impl>(dev,std::move(mesh),std::move(boundaries),settings)) {}
SteadySolver::~SteadySolver()=default;
void SteadySolver::initialize(const std::function<Vec3d(Vec3d)>& velocity,const std::function<double(Vec3d)>& pressure) {
  auto& s=*impl_;std::vector<Vec3d> u;std::vector<float> p;
  const double ref=s.reference>=0?pressure(s.geom.cell_centroid[s.reference]):0;
  for(auto x:s.geom.cell_centroid) {u.push_back(velocity(x));const double value=pressure(x)-ref;
    KUTTA_REQUIRE(std::isfinite(value),"FV nonfinite pressure initialization");p.push_back(float(value));}
  upload_vectors(s.dev,s.velocity,u);s.dev.upload<float>(s.pressure,p);
  s.update_boundaries();s.update_pressure();
  // Initialize a consistent advective mass flux before the first assembly.
  s.dev.fill(s.momentum[0]->diagonal,1);s.dev.fill(s.momentum[0]->consistent_inverse,0);
  s.update_flux();
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
  vector_gradient(s.dev,s.mesh,s.stencil,s.velocity,s.boundary_velocity,s.momentum[0]->grad,s.momentum[1]->grad,s.momentum[2]->grad);
  auto mask=s.dev.upload<float>(selected);
  auto pressure=vector_field(s.dev,selected.size()),viscous=vector_field(s.dev,selected.size());
  const auto& m=s.mesh;const auto& gu=s.momentum[0]->grad;const auto& gv=s.momentum[1]->grad;const auto& gw=s.momentum[2]->grad;
  kern::fv_wall_pressure_force(s.dev,selected.size(),mask,m.area.x,m.area.y,m.area.z,s.face_pressure,
    pressure.x,pressure.y,pressure.z,m.n_internal_faces,float(reference_pressure));
  kern::fv_wall_viscous_force(s.dev,selected.size(),m.owner,mask,m.area.x,m.area.y,m.area.z,
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
const IterationReport& SteadySolver::report() const {return impl_->report;}
Vec3d ForceReport::force() const {return pressure+viscous;}
double ForceReport::coefficient(Vec3d direction,const ReferenceValues& ref) const {
  KUTTA_REQUIRE(ref.density>0 && ref.velocity>0 && ref.area>0,"FV positive force references required");
  const auto adjusted=force()+area_vector*(reference_pressure-ref.pressure);
  return dot(adjusted,direction)/(.5*ref.density*ref.velocity*ref.velocity*ref.area);
}
}  // namespace kutta::fv
