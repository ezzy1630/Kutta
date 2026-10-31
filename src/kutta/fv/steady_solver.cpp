#include "kutta/fv/steady_solver.hpp"
#include "kutta/fv/turbulence/sa_transport.hpp"
#include "kutta/fv/turbulence/realizable_ke_transport.hpp"
#include "kutta/fv/turbulence/intermittency.hpp"
#include "kutta/fv/coupled_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>

#include "kernels.hpp"
#include "kutta/linalg/amg.hpp"
#include "kutta/kernels/prelude.h"
#include "kutta/fv/boundary_physics.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
#include "kutta/core/memory.hpp"

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

namespace {
// A quadratic cloud assumes a continuous trace. A moving lid meeting a fixed
// wall violates that assumption; retain the existing scheme in its support.
std::vector<Index> discontinuous_wall_cells(const PolyMesh& mesh,
    const std::map<std::string,BoundaryCondition>& conditions,std::span<const RigidFrame> frames) {
  struct Trace {Index point,cell;Vec3d velocity;};std::vector<Trace> traces;
  for(const auto& patch:mesh.patches) if(const auto* wall=std::get_if<Wall>(&conditions.at(patch.name)))
    for(Index f=patch.start;f<patch.start+patch.size;++f) {
      const auto frame=frames.empty()?RigidFrame{}:frames[mesh.owner[f]];
      for(auto point:detail::nodes(mesh,f)) traces.push_back({point,mesh.owner[f],wall_velocity(*wall,mesh.points[point],frame)});
    }
  std::sort(traces.begin(),traces.end(),[](const auto& a,const auto& b){return a.point<b.point;});
  std::vector<Index> cells(mesh.n_cells,0);
  for(std::size_t first=0;first<traces.size();) {
    std::size_t end=first+1;while(end<traces.size() && traces[end].point==traces[first].point) ++end;
    bool discontinuous=false;
    for(auto j=first+1;j<end;++j) {
      const auto a=traces[first].velocity,b=traces[j].velocity;
      discontinuous|=norm(a-b)>64*std::numeric_limits<double>::epsilon()*std::max({1.,norm(a),norm(b)});
    }
    if(discontinuous) for(auto j=first;j<end;++j) cells[traces[j].cell]=1;
    first=end;
  }
  return cells;
}
}

struct SteadySolver::Impl {
  Device& dev;
  PolyMesh host;
  MeshGeometry geom;
  SolverSettings settings;
  DeviceMesh mesh;
  Stencil stencil;
  BoundaryLeastSquares velocity_lsq,pressure_lsq,k_lsq,omega_lsq,strain_lsq;FaceJumpFit jump_fit;
  DeviceVec3 strain_base[3],strain_gradient[3];
  MeshPermutation permutation;
  std::unique_ptr<SstModel> sst;
  std::unique_ptr<SaTransport> sa;
  std::unique_ptr<RealizableKeTransport> rke;
  std::unique_ptr<IntermittencyTransport> transition;
  std::unique_ptr<TransportEquation> k_equation,omega_equation;
  Buffer<float> kinetic,omega,rho_cell,mu_cell,boundary_k,boundary_omega,turbulence_a,turbulence_b;
  Buffer<Index> turbulence_spec,turbulence_kind,zero_mask,total_pressure;
  Buffer<float> wall_diagonal,wall_rhs,physical_diagonal,consistent_inverse;
  // Distance-based Rhie-Chow (TG 23-48) is the momentum-weighted face
  // interpolation with equal weights.
  Buffer<float> unit_diagonal;
  Buffer<float> wall_mask;DeviceVec3 laminar_wall_shear;
  DeviceVec3 transpose_stress;
  DeviceVec3 velocity,body,boundary_velocity,prescribed_velocity,source,face_velocity;
  Buffer<float> pressure,mass,viscosity,boundary_pressure,prescribed_pressure;
  Buffer<Index> boundary_kind,scalar_kind;
  DeviceVec3 grad_p,limited_gradient[3];
  std::unique_ptr<TransportEquation> momentum[3];
  std::shared_ptr<const ProjectedFaceReconstruction> diffusion_reconstruction;
  std::shared_ptr<const WalkFaceReconstruction> turbulent_faces;
  std::shared_ptr<const FaceReconstruction> velocity_faces;
  DeviceVec3 projected_gradient[3];
  std::shared_ptr<CsrMatrix> pressure_matrix;
  // Momentum, turbulence and pressure solves never interleave, so their
  // Krylov vectors are one set instead of one per equation.
  std::shared_ptr<KrylovWorkspace> krylov_workspace;
  std::unique_ptr<LinearSolver> pressure_solver;
  std::unique_ptr<CoupledSolver> coupled;
  Buffer<float> pressure_rhs,correction,correction_boundary,coefficient,face_pressure,pressure_jump;
  DeviceVec3 grad_correction;
  Buffer<float> residual_u,residual_v,residual_w,continuity,provisional_continuity,scale,speed;
  DeviceVec3 absolute_residual,fluent_scale;
  IterationReport report;
  std::map<std::string,BoundaryCondition> boundaries;
  // BEGIN FSAE physics: isolated volume sources and open-boundary overrides.
  std::unique_ptr<ZonePhysics> zone_physics;
  std::unique_ptr<InteriorJumps> interior_jumps;
  Buffer<Index> open_boundary_kind;
  DeviceVec3 inlet_direction;
  bool pressure_inlets=false;
  // END FSAE physics
  Index reference=-1;
  double length=0,bc_speed=0,pressure_speed=0,continuity_scale=0;
  double first_five_continuity[5]{};
  float time_step=0;
  double velocity_scale=1; // the automatic time step's velocity scale

  Impl(Device& device,PolyMesh input,std::vector<PatchCondition> conditions,SolverSettings options)
      :dev(device),host(std::move(input)),settings(options) {
    KUTTA_REQUIRE(settings.pseudo_transient || settings.coupling!=PressureVelocityCoupling::Coupled,
      "FV: the Coupled solver requires pseudo-transient stepping");
    KUTTA_REQUIRE(int(settings.sa.has_value())+int(settings.sst.has_value())+int(settings.realizable_ke.has_value())<=1,"FV: select exactly one turbulence model");
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
      if(std::holds_alternative<MassFlowInlet>(it->second)) patch.kind=PatchKind::VelocityInlet;
      if(std::holds_alternative<PressureInlet>(it->second)) {patch.kind=PatchKind::PressureInlet;outlet=true;}
      if(std::holds_alternative<PressureOutlet>(it->second)) {patch.kind=PatchKind::PressureOutlet;outlet=true;}
      if(std::holds_alternative<Symmetry>(it->second)) patch.kind=PatchKind::Symmetry;
    }
    permutation=renumber_cells_rcm(host);memory::phase("renumber");geom=compute_geometry(host);memory::phase("geometry");
    KUTTA_REQUIRE(host.n_cells>0,"FV steady solver needs a nonempty mesh");
    if(!outlet) reference=locate_cell(host,geom,settings.reference_point);
    // A single pin cannot constrain disconnected Neumann components.
    const auto neighbours=build_cell_neighbours(host);
    std::vector<bool> seen(host.n_cells);std::vector<Index> stack{0};seen[0]=true;
    for(std::size_t k=0;k<stack.size();++k) for(Index j=neighbours.offsets[stack[k]];j<neighbours.offsets[stack[k]+1];++j)
      if(!seen[neighbours.cells[j]]) {seen[neighbours.cells[j]]=true;stack.push_back(neighbours.cells[j]);}
    KUTTA_REQUIRE(stack.size()==std::size_t(host.n_cells),"FV disconnected domains require separate solvers");
    memory::phase("connectivity");
    mesh=upload_mesh(dev,host,geom);memory::phase("upload-mesh");stencil=make_stencil(dev,host,geom);memory::phase("stencil");
    // BEGIN FSAE physics
    zone_physics=std::make_unique<ZonePhysics>(dev,host,geom,settings.cell_zone_physics,settings.density,settings.viscosity);
    interior_jumps=std::make_unique<InteriorJumps>(dev,host,geom,settings.interior_conditions,settings.density,settings.viscosity);
    const auto frames=settings.cell_zone_physics.empty()?std::vector<RigidFrame>{}:cell_frames(host,settings.cell_zone_physics);
    const auto adjacent_frame=[&](Index face) {return frames.empty()?RigidFrame{}:frames[host.owner[face]];};
    // END FSAE physics
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
    if(settings.face_velocity==FaceVelocityWeights::Distance) {unit_diagonal=dev.alloc<float>(n);dev.fill(unit_diagonal,1);}
    std::vector<Index> kinds(nb);std::vector<float> pressures(nb,0);std::vector<Vec3d> velocities(nb);
    std::vector<Index> open_kinds(nb,0);std::vector<Vec3d> directions(nb);
    std::map<std::string,std::pair<double,double>> inlet_areas;
    for(const auto& patch:host.patches) if(const auto* inlet=std::get_if<MassFlowInlet>(&boundaries.at(patch.name))) {
      auto& [area,projected]=inlet_areas[patch.name];
      for(Index f=patch.start;f<patch.start+patch.size;++f) {
        area+=norm(geom.face_area[f]);
        if(inlet->direction && norm(*inlet->direction)>0) projected-=dot(*inlet->direction/norm(*inlet->direction),geom.face_area[f]);
      }
    }
    double pmin=std::numeric_limits<double>::infinity(),pmax=-pmin;
    for(const auto& patch:host.patches) for(Index f=patch.start;f<patch.start+patch.size;++f) {
      const auto b=f-mesh.n_internal_faces;const auto& bc=boundaries.at(patch.name);
      if(const auto* wall=std::get_if<Wall>(&bc)) {kinds[b]=0;velocities[b]=wall_velocity(*wall,geom.face_centroid[f],adjacent_frame(f));}
      if(const auto* inlet=std::get_if<VelocityInlet>(&bc)) {kinds[b]=1;velocities[b]=evaluate(inlet->velocity,geom.face_centroid[f]);}
      // BEGIN FSAE physics: mass flow is integrated over the whole patch.
      if(const auto* inlet=std::get_if<MassFlowInlet>(&bc)) {
        const auto [area,projected]=inlet_areas.at(patch.name);
        kinds[b]=1;velocities[b]=mass_flow_velocity(*inlet,geom.face_area[f],area,projected,settings.density);
      }
      if(const auto* inlet=std::get_if<PressureInlet>(&bc)) {
        KUTTA_REQUIRE(std::isfinite(inlet->total_pressure),"pressure-inlet: nonfinite total pressure");
        kinds[b]=4;totals[b]=1;pressures[b]=float(inlet->total_pressure);open_kinds[b]=1;pressure_inlets=true;
        directions[b]=inlet->direction.value_or(geom.face_area[f]*(-1/norm(geom.face_area[f])));
        KUTTA_REQUIRE(finite(directions[b]) && dot(directions[b],geom.face_area[f])<0,"pressure-inlet: direction must enter every face");
        pmin=std::min(pmin,inlet->total_pressure);pmax=std::max(pmax,inlet->total_pressure);
        KUTTA_REQUIRE(norm(adjacent_frame(f).velocity(geom.face_centroid[f]))==0,"pressure-inlet: moving adjacent frame is not yet verified");
      }
      // END FSAE physics
      if(const auto* out=std::get_if<PressureOutlet>(&bc)) {
        KUTTA_REQUIRE(norm(adjacent_frame(f).velocity(geom.face_centroid[f]))==0,
          "pressure-outlet: moving adjacent frame/backflow conventions are not yet verified");
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
    open_boundary_kind=dev.upload<Index>(open_kinds);inlet_direction=vector_field(dev,nb);upload_vectors(dev,inlet_direction,directions);
    using enum BoundaryObservation;
    std::vector<BoundaryObservation> independent_velocity(nb),independent_pressure(nb),independent_k(nb),independent_omega(nb);
    for(Index b=0;b<nb;++b) {
      independent_velocity[b]=kinds[b]==3?NormalFootValue:((kinds[b]==0 || kinds[b]==1)?FaceValue:Omitted);
      // The standard scheme's unprescribed boundary rows hold the cell
      // pressure, the value Fluent stores there, at the normal foot: its zero
      // normal derivative. Second-order pressure extrapolates them, so they
      // carry no information and are omitted, except where that
      // extrapolation amplifies (below).
      const bool zero_normal=settings.pressure_interpolation==PressureInterpolation::Standard;
      independent_pressure[b]=(kinds[b]==2 || kinds[b]==4)?FaceValue:(zero_normal?NormalFootValue:Omitted);
      independent_k[b]=kinds[b]==1?FaceValue:Omitted;independent_omega[b]=(kinds[b]==0 || kinds[b]==1)?FaceValue:Omitted;
      // BEGIN FSAE physics: pressure inlets prescribe transported inflow data.
      if(open_kinds[b]==1) {independent_velocity[b]=FaceValue;independent_k[b]=FaceValue;independent_omega[b]=FaceValue;}
      // END FSAE physics
    }
    dev.upload<Index>(boundary_kind,kinds);dev.upload<float>(prescribed_pressure,pressures);
    upload_vectors(dev,prescribed_velocity,velocities);
    memory::phase("fields");
    pressure_matrix=std::make_shared<CsrMatrix>(upload(dev,stencil.pattern));
    krylov_workspace=std::make_shared<KrylovWorkspace>(n);
    const bool green_gauss=settings.gradient==GradientScheme::GreenGaussCell;
    // SST's strain and SA's gradients have only been verified with LSQ fits.
    KUTTA_REQUIRE(!green_gauss || (!settings.sst && !settings.sa),"FV: Green-Gauss gradients serve laminar and k-epsilon flows only");
    const auto fit=[&](std::span<const BoundaryObservation> observations) {
      return green_gauss?green_gauss_cell(dev,host,geom):boundary_least_squares(dev,host,geom,observations,&stencil);
    };
    // The projected reconstruction is a quadratic least-squares extension of
    // the linear primary/secondary split; a Green-Gauss case keeps Fluent's
    // documented split with its own cell gradients.
    if(!settings.sst && !settings.sa && !settings.realizable_ke && !green_gauss) {
      const auto nonsmooth=discontinuous_wall_cells(host,boundaries,frames);
      diffusion_reconstruction=std::make_shared<ProjectedFaceReconstruction>(dev,host,geom,independent_velocity,nonsmooth);
      if(diffusion_reconstruction->corrected_faces()==0) diffusion_reconstruction.reset();
      velocity_faces=diffusion_reconstruction;
    } else if(!green_gauss) {
      // Turbulent models use the compact walk-weighted fit. A wall function
      // owns the wall-to-cell profile, so SST and k-epsilon wall values stay
      // out of the fit, as they do from SST's strain LSQ. Resolved SA walls
      // are Dirichlet data, like laminar walls.
      auto observed=independent_velocity;
      if(!settings.sa) for(Index b=0;b<nb;++b) if(kinds[b]==0) observed[b]=Omitted;
      const auto nonsmooth=settings.sa?discontinuous_wall_cells(host,boundaries,frames):std::vector<Index>{};
      turbulent_faces=std::make_shared<WalkFaceReconstruction>(dev,host,geom,mesh,stencil,observed,nonsmooth);
      if(turbulent_faces->corrected_faces()==0) turbulent_faces.reset();
      velocity_faces=turbulent_faces;
    }
    memory::phase("face-reconstruction");
    if(diffusion_reconstruction) {
      for(auto& g:projected_gradient) g=vector_field(dev,nf);
    }
    // Under wall functions only the momentum predictor is projected: on V6
    // the projected mu+mu_t diffusion moved no field by more than 0.03 points
    // (validation/results/turbulent_poly) and cost three fits per iteration.
    // Resolved SA walls keep it, because their wall traction is the projected
    // wall flux.
    const bool momentum_diffusion=velocity_faces && !settings.sst && !settings.realizable_ke;
    for(auto& eq:momentum) eq=std::make_unique<TransportEquation>(dev,mesh,stencil,host,geom,pressure_matrix,
      momentum_diffusion?velocity_faces:nullptr,krylov_workspace);
    memory::phase("transport");
    velocity_lsq=fit(independent_velocity);
    if(settings.pressure_interpolation!=PressureInterpolation::Standard && !green_gauss) {
      // Where the fit's extrapolation to an unprescribed face can amplify
      // neighbour differences (gain > 1), as at tet wall cells whose
      // neighbours lie in one plane, it feeds back through momentum's
      // boundary pressure and diverges (the mixing-elbow tets). There the
      // face holds the cell pressure, Fluent's stored value, at its centroid.
      // Elsewhere a zero normal derivative would bias the gradient wherever
      // the face carries a pressure gradient: an inlet under a Darcy
      // gradient (FSAE V8), a rotating wall (FSAE V7.2)
      // (validation/results/tet_robustness/REPORT.md).
      const auto gain=boundary_extrapolation_gain(host,geom,independent_pressure);
      for(Index b=0;b<nb;++b) if(independent_pressure[b]==Omitted && gain[b]>1) independent_pressure[b]=FaceValue;
    }
    pressure_lsq=fit(independent_pressure);
    memory::phase("lsq");
    if(settings.sst) {
      std::vector<BoundaryObservation> strain_rows(nb);
      for(Index b=0;b<nb;++b) strain_rows[b]=kinds[b]!=0?FaceValue:Omitted;
      strain_lsq=boundary_least_squares(dev,host,geom,strain_rows,&stencil);
      jump_fit=face_jump_least_squares(dev,host,geom);
      for(int k=0;k<3;++k) {strain_base[k]=vector_field(dev,n);strain_gradient[k]=vector_field(dev,n);}
      const auto& t=*settings.sst;
      KUTTA_REQUIRE(t.k_relaxation>0 && t.k_relaxation<=1 && t.omega_relaxation>0 && t.omega_relaxation<=1,"FV invalid SST relaxation");
      k_lsq=boundary_least_squares(dev,host,geom,independent_k,&stencil);
      omega_lsq=boundary_least_squares(dev,host,geom,independent_omega,&stencil);
      memory::phase("turbulence-lsq");
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
        if(const auto* in=std::get_if<PressureInlet>(&bc)) spec=&in->turbulence;
        if(const auto* in=std::get_if<MassFlowInlet>(&bc)) spec=&in->turbulence;
        if(!spec) continue;
        KUTTA_REQUIRE(spec->has_value(),"FV missing turbulence specification: "+patch.name);
        KUTTA_REQUIRE(std::holds_alternative<TurbulenceValues>(**spec) || std::holds_alternative<IntensityLength>(**spec) || std::holds_alternative<IntensityViscosityRatio>(**spec),"FV: wrong boundary turbulence model supplied to SST: "+patch.name);
        for(Index f=patch.start;f<patch.start+patch.size;++f) {
          const auto bf=f-mesh.n_internal_faces;
          if(const auto* v=std::get_if<TurbulenceValues>(&**spec)) {a[bf]=float(v->k);b[bf]=float(v->omega);}
          if(const auto* v=std::get_if<IntensityLength>(&**spec)) {specs[bf]=1;a[bf]=float(v->intensity);b[bf]=float(v->length);}
          if(const auto* v=std::get_if<IntensityViscosityRatio>(&**spec)) {specs[bf]=2;a[bf]=float(v->intensity);b[bf]=float(v->ratio);}
          KUTTA_REQUIRE(a[bf]>=0 && b[bf]>0 && std::isfinite(a[bf]) && std::isfinite(b[bf]),"FV invalid turbulence specification: "+patch.name);
        }
      }
      turbulence_spec=dev.upload<Index>(specs);turbulence_a=dev.upload<float>(a);turbulence_b=dev.upload<float>(b);
      // Wall omega is the wall treatment's value, not a smooth trace.
      auto omega_observed=independent_omega;
      for(Index b=0;b<nb;++b) if(kinds[b]==0) omega_observed[b]=Omitted;
      // Projected k diffusion had no measurable effect on V6; omega diffusion
      // on skewed faces was V6's largest turbulent-operator error.
      k_equation=std::make_unique<TransportEquation>(dev,mesh,stencil,host,geom,pressure_matrix,nullptr,krylov_workspace);
      omega_equation=std::make_unique<TransportEquation>(dev,mesh,stencil,host,geom,pressure_matrix,
        walk_reconstruction(turbulent_faces,dev,host,geom,mesh,stencil,omega_observed),krylov_workspace);
      memory::phase("turbulence-transport");
      KUTTA_REQUIRE(!t.model.intermittency,"FV: enable transition with SstSettings::transition, not SstOptions::intermittency");
      auto sst_options=t.model;
      if(t.transition) {
        sst_options.intermittency=SstIntermittency{t.transition->model.ck,t.transition->model.csep};
        transition=std::make_unique<IntermittencyTransport>(dev,mesh,stencil,host,geom,boundaries,settings,
          pressure_matrix,k_lsq,krylov_workspace);
        memory::phase("transition");
      }
      sst=std::make_unique<SstModel>(dev,mesh,sst_options);
    }
    if(settings.sa) {
      sa=std::make_unique<SaTransport>(dev,mesh,stencil,host,geom,boundaries,settings,pressure_matrix,turbulent_faces);
      transpose_stress=vector_field(dev,nf);
    }
    if(settings.realizable_ke) {
      rke=std::make_unique<RealizableKeTransport>(dev,mesh,stencil,host,geom,boundaries,settings,pressure_matrix,turbulent_faces);
      transpose_stress=vector_field(dev,nf);
    }
    memory::phase("turbulence");
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
    if(settings.sst) {
      AmgOptions options;options.pre_sweeps=options.post_sweeps=2;
      options.structure=OperatorStructure::SymmetricMMatrix;options.workspace=krylov_workspace;
      pressure_solver=make_amg_solver(pressure_solver_kind,dev,initial,options);
    } else pressure_solver=make_solver(pressure_solver_kind,dev,initial,
      {.structure=OperatorStructure::SymmetricMMatrix,.workspace=krylov_workspace});
    memory::phase("pressure-solver");
    const double volume=std::accumulate(geom.cell_volume.begin(),geom.cell_volume.end(),0.);
    Vec3d lo=host.points.front(),hi=lo;
    for(auto p:host.points) {lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
    // Fluent's conservative length scale: min(domain volume^(1/3), largest extent).
    length=std::min(std::cbrt(volume),std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z}));
    if(settings.coupling==PressureVelocityCoupling::Coupled)
      coupled=std::make_unique<CoupledSolver>(dev,mesh,stencil,host,geom,boundary_kind,pressure_lsq,reference);
    memory::phase("coupled-solver");
    update_time(bc_speed,bc_speed);
    dev.fill(momentum[0]->diagonal,1);dev.fill(momentum[0]->consistent_inverse,0);
    update_boundaries();update_pressure();update_flux();
    memory::phase("solver-ready");
  }

  void apply_backflow_pressure() {
    kern::fv_backflow_pressure(dev,mesh.n_boundary_faces(),mesh.owner,total_pressure,mass,
      velocity.x,velocity.y,velocity.z,boundary_pressure,mesh.n_internal_faces,float(settings.density));
    // BEGIN FSAE physics: inlet Bernoulli uses the inlet-plane speed.
    if(pressure_inlets) kern::fv_pressure_inlet_pressure(dev,mesh.n_boundary_faces(),open_boundary_kind,mass,
      prescribed_pressure,boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,boundary_pressure,mesh.n_internal_faces,float(settings.density));
    // END FSAE physics
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
      kinetic,omega,k_equation->grad,omega_equation->grad,boundary_velocity,
      transition?&transition->gamma():nullptr});
  }
  void update_turbulence(bool velocity_changed=true) {
    if(rke) {
      reconstruct_velocity();
      rke->update(velocity,boundary_velocity,{momentum[0]->grad,momentum[1]->grad,momentum[2]->grad},boundary_kind,mass,viscosity);
      return;
    }
    if(sa) {
      reconstruct_velocity();
      sa->update(velocity,boundary_velocity,momentum[0]->grad,momentum[1]->grad,momentum[2]->grad,boundary_kind,mass,viscosity);
      return;
    }
    if(!sst) return;
    if(velocity_changed) {
      vector_gradient(dev,mesh,stencil,strain_lsq,velocity,boundary_velocity,
        strain_base[0],strain_base[1],strain_base[2]);
      const Buffer<float>* q[]{&velocity.x,&velocity.y,&velocity.z};
      const Buffer<float>* bc[]{&boundary_velocity.x,&boundary_velocity.y,&boundary_velocity.z};
      for(int k=0;k<3;++k)
        jump_corrected_gradient(dev,mesh,stencil,jump_fit,*q[k],*bc[k],strain_base[k],strain_gradient[k]);
    }
    evaluate_turbulence(); // bootstrap gradient-independent wall-face omega
    turbulence_boundaries();
    gradient(dev,mesh,stencil,k_lsq,kinetic,boundary_k,k_equation->grad);
    gradient(dev,mesh,stencil,omega_lsq,omega,boundary_omega,omega_equation->grad);
    evaluate_turbulence();
    if(transition) transition->update(sst->fields(),rho_cell,mu_cell,kinetic,omega,strain_gradient,boundary_kind,mass);
    kern::fv_effective_viscosity(dev,mesh.n_faces,mesh.owner,mesh.neighbour,boundary_kind,
      mesh.weight,sst->fields().mu_t,viscosity,mesh.n_internal_faces,float(settings.viscosity));
  }
  void add_transpose_stress() {
    // BEGIN FSAE physics: each side sees its own discontinuous face pressure.
    interior_jumps->momentum(mesh,source);
    // END FSAE physics
    // Callers reconstruct the velocity gradients once before this.
    if(!sst && !sa && !rke) return;
    const auto& u=momentum[0]->grad;const auto& v=momentum[1]->grad;const auto& w=momentum[2]->grad;
    kern::fv_transpose_stress(dev,mesh.n_faces,mesh.owner,mesh.neighbour,boundary_kind,mesh.weight,viscosity,
      mesh.area.x,mesh.area.y,mesh.area.z,u.x,u.y,u.z,v.x,v.y,v.z,w.x,w.y,w.z,
      transpose_stress.x,transpose_stress.y,transpose_stress.z,mesh.n_internal_faces);
    kern::fv_stress_source(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      mesh.volume,transpose_stress.x,transpose_stress.y,transpose_stress.z,source.x,source.y,source.z);
  }
  TransportTerms wall_terms(int component) {
    kern::fv_boundary_momentum(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,boundary_kind,
      mesh.area.x,mesh.area.y,mesh.area.z,sst?sst->fields().wall_shear_coefficient:(rke?rke->fields().wall_coefficient:prescribed_pressure),
      viscosity,stencil.diffusion,
      velocity.x,velocity.y,velocity.z,boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,
      wall_diagonal,wall_rhs,mesh.n_internal_faces,component,(sst || rke)?1:0);
    // BEGIN FSAE physics: same implicit terms in iteration and residual assembly.
    zone_physics->momentum(component,velocity,mesh.volume,wall_diagonal,wall_rhs);
    // END FSAE physics
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
    if(rke) {
      update_boundaries();update_turbulence();report.turbulence[0]=rke->solve(false,time_step);
      update_turbulence();report.turbulence[1]=rke->solve(true,time_step);return;
    }
    if(sa) {update_boundaries();update_turbulence();report.turbulence[0]=sa->solve(time_step);return;}
    if(!sst) return;
    update_boundaries();update_turbulence();assemble_turbulence(false);
    report.turbulence[0]=k_equation->solve(kinetic,settings.sst->linear);
    kern::fv_floor(dev,mesh.n_cells,kinetic,settings.sst->model.k_floor);
    // k changed; velocity and its boundary values did not. Keep the derivative
    // gradients used for strain/production, while refreshing both scalar fits.
    update_turbulence(false);assemble_turbulence(true);
    report.turbulence[1]=omega_equation->solve(omega,settings.sst->linear);
    kern::fv_floor(dev,mesh.n_cells,omega,settings.sst->model.omega_floor);
    // Fluent order: k, omega, then intermittency, each from the latest fields.
    if(transition) {update_turbulence(false);report.transition=transition->solve(time_step);}
  }
  void monitor_turbulence() {
    if(rke) {auto r=rke->monitor(time_step);report.k=r[0];report.epsilon=r[1];return;}
    if(sa) {report.nu_tilde=sa->monitor(time_step);return;}
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
    if(transition) report.intermittency=transition->monitor(time_step);
  }
  // Fluent's automatic global time step (Theory Guide 23.6.2, eqs. 23-102 to
  // 23-105): 0.3 L / max(0.1 U_bc, U_domain, U_dp). Fluent defines it for the
  // Coupled solver only, where U_domain is the arithmetic mean cell speed.
  // SIMPLE and SIMPLEC keep Kutta's older, more conservative maximum cell
  // speed. The Coupled block scaling wc keeps the maximum speed either way.
  void update_time(double maximum_speed,double mean_speed) {
    velocity_scale=std::max({.1*bc_speed,maximum_speed,pressure_speed,1e-10});
    if(!settings.pseudo_transient) {time_step=0;return;}
    const double time_speed=coupled?std::max({.1*bc_speed,mean_speed,pressure_speed,1e-10}):velocity_scale;
    time_step=settings.pseudo_time_step>0?float(settings.pseudo_time_step):
      float(.3*length*settings.pseudo_time_scale/time_speed);
  }
  // The automatic pseudo time step follows the current flow. step() and
  // monitor() both refresh it from the same state, so monitoring never changes
  // the step a run takes; it used to be refreshed only in monitor(), which
  // made the solution depend on the monitor interval.
  void refresh_time_step() {
    if(settings.pseudo_time_step>0 || !settings.pseudo_transient) return;
    kern::fv_speed(dev,mesh.n_cells,velocity.x,velocity.y,velocity.z,speed);
    if(!coupled) {const auto maximum=dev.max_abs(speed,mesh.n_cells);update_time(maximum,maximum);return;}
    const auto speeds=dev.reduce_many({Device::HostReduction::max_abs(speed,mesh.n_cells),
      Device::HostReduction::sum(speed,mesh.n_cells)});
    update_time(speeds[0],speeds[1]/double(mesh.n_cells));
  }
  void update_boundaries() {
    kern::fv_velocity_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,
      mesh.area.x,mesh.area.y,mesh.area.z,velocity.x,velocity.y,velocity.z,
      prescribed_velocity.x,prescribed_velocity.y,prescribed_velocity.z,mass,
      boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,scalar_kind,mesh.n_internal_faces,float(settings.density));
    // BEGIN FSAE physics
    if(pressure_inlets) pressure_inlet_velocity(dev,mesh,open_boundary_kind,inlet_direction,mass,boundary_velocity,settings.density);
    if(zone_physics->has_frame_motion()) kern::fv_frame_wall_velocity(dev,mesh.n_boundary_faces(),boundary_kind,
      zone_physics->frame_normal_flux(),mesh.area.x,mesh.area.y,mesh.area.z,
      boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,mesh.n_internal_faces,float(settings.density));
    // END FSAE physics
  }
  void reconstruct_velocity() {
    vector_gradient(dev,mesh,stencil,velocity_lsq,velocity,boundary_velocity,
      momentum[0]->grad,momentum[1]->grad,momentum[2]->grad);
  }
  void reconstruct_pressure(const Buffer<float>& field,const Buffer<float>& boundary,DeviceVec3& grad) {
    gradient(dev,mesh,stencil,pressure_lsq,field,boundary,grad);
    // BEGIN FSAE physics: the frozen physical jump is absent from p-correction.
    if(&field==&pressure) interior_jumps->gradient(mesh,pressure_lsq,grad);
    // END FSAE physics
  }
  void update_pressure() {
    // BEGIN FSAE physics
    interior_jumps->update(mass);
    // END FSAE physics
    kern::fv_pressure_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,
      mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,pressure,grad_p.x,grad_p.y,grad_p.z,
      prescribed_pressure,boundary_pressure,mesh.n_internal_faces);
    kern::fv_cell_pressure_boundary(dev,mesh.n_boundary_faces(),mesh.owner,boundary_kind,pressure,
      boundary_pressure,mesh.n_internal_faces);
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
    if(settings.pressure_interpolation==PressureInterpolation::Standard)
      kern::fv_standard_pressure_faces(dev,mesh.n_faces,mesh.owner,mesh.neighbour,boundary_kind,
        pressure,physical_diagonal,face_pressure,boundary_pressure,mesh.n_internal_faces);
  }
  const Buffer<float>& face_velocity_weights() const {
    return settings.face_velocity==FaceVelocityWeights::Distance?unit_diagonal:physical_diagonal;
  }
  // SIMPLE replaces SIMPLEC's consistent inverse with the relaxed diagonal's.
  void momentum_correction_coefficient(int component) {
    if(settings.coupling!=PressureVelocityCoupling::Simple) return;
    auto& eq=*momentum[component];
    kern::fv_simple_inverse(dev,mesh.n_cells,eq.matrix.diag,eq.matrix.val,eq.consistent_inverse);
  }
  void update_flux() {
    reconstruct_velocity();
    const Buffer<float>* q[]{&velocity.x,&velocity.y,&velocity.z};
    const Buffer<float>* bc[]{&boundary_velocity.x,&boundary_velocity.y,&boundary_velocity.z};
    for(int k=0;k<3;++k) {
      auto& eq=*momentum[k];auto& g=limited_gradient[k];
      if(settings.limit) kern::fv_limiter_scale(dev,mesh.n_cells,mesh.cell_face_offsets,stencil.adjacent,
        stencil.face_delta.x,stencil.face_delta.y,stencil.face_delta.z,*q[k],*bc[k],
        eq.grad.x,eq.grad.y,eq.grad.z,eq.limiter,g.x,g.y,g.z);
      else {
        dev.fill(eq.limiter,1);
        kern::fv_scale_gradient(dev,mesh.n_cells,eq.grad.x,eq.grad.y,eq.grad.z,eq.limiter,g.x,g.y,g.z);
      }
    }
    const auto& gu=limited_gradient[0];const auto& gv=limited_gradient[1];const auto& gw=limited_gradient[2];
    kern::fv_momentum_faces(dev,mesh.n_faces,mesh.owner,mesh.neighbour,boundary_kind,
      mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,
      mesh.neighbour_to_face.x,mesh.neighbour_to_face.y,mesh.neighbour_to_face.z,
      velocity.x,velocity.y,velocity.z,gu.x,gu.y,gu.z,gv.x,gv.y,gv.z,gw.x,gw.y,gw.z,
      face_velocity_weights(),boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,
      face_velocity.x,face_velocity.y,face_velocity.z,mesh.n_internal_faces);
    if(velocity_faces) {
      velocity_faces->momentum_values(velocity.x,boundary_velocity.x,face_velocity_weights(),face_velocity.x,settings.limit?&momentum[0]->limiter:nullptr);
      velocity_faces->momentum_values(velocity.y,boundary_velocity.y,face_velocity_weights(),face_velocity.y,settings.limit?&momentum[1]->limiter:nullptr);
      velocity_faces->momentum_values(velocity.z,boundary_velocity.z,face_velocity_weights(),face_velocity.z,settings.limit?&momentum[2]->limiter:nullptr);
    }
    kern::fv_rhie_chow(dev,mesh.n_faces,mesh.owner,mesh.neighbour,boundary_kind,
      mesh.area.x,mesh.area.y,mesh.area.z,stencil.diffusion,mesh.volume,
      physical_diagonal,consistent_inverse,face_velocity.x,face_velocity.y,face_velocity.z,
      pressure_jump,
      mass,coefficient,mesh.n_internal_faces,float(settings.density));
    // BEGIN FSAE physics
    zone_physics->relative_flux(boundary_kind,mesh.n_internal_faces, mass);
    interior_jumps->flux(mesh,stencil,physical_diagonal,mass,coefficient);
    // END FSAE physics
  }

  // Coupled has no implicit momentum relaxation: its relaxation is explicit.
  TransportSettings momentum_transport() const {
    TransportSettings transport;
    transport.convection=settings.convection;transport.limit=settings.limit;
    transport.density=float(settings.density);transport.pseudo_time_step=time_step;
    transport.relaxation=coupled?1.f:float(settings.momentum_relaxation);
    return transport;
  }
  void momentum_sources() {
    kern::fv_momentum_source(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      mesh.area.x,mesh.area.y,mesh.area.z,face_pressure,mesh.volume,body.x,body.y,body.z,
      source.x,source.y,source.z);
    reconstruct_velocity();add_transpose_stress();
  }
  // One pseudo-time Coupled iteration (Fluent Theory Guide 23.4.3.2 and the
  // pseudo-time form 23-101). The block RHS is the steady residual at the
  // current state: momentum b - A u with the explicit reconstructed face
  // pressure, and the Rhie-Chow mass imbalance. Turbulence follows segregated.
  void coupled_step() {
    refresh_time_step();
    update_boundaries();update_pressure();update_turbulence();momentum_sources();
    const auto transport=momentum_transport();
    Buffer<float>* q[]{&velocity.x,&velocity.y,&velocity.z};
    const Buffer<float>* bc[]{&boundary_velocity.x,&boundary_velocity.y,&boundary_velocity.z};
    const Buffer<float>* src[]{&source.x,&source.y,&source.z};
    Buffer<float>* res[]{&residual_u,&residual_v,&residual_w};
    for(int k=0;k<3;++k) {
      momentum[k]->assemble_reconstructed(*q[k],*bc[k],scalar_kind,mass,viscosity,*src[k],transport,wall_terms(k));
      residual(dev,momentum[k]->matrix,*q[k],momentum[k]->rhs,*res[k]);
      coupled->capture_diagonal(k,momentum[k]->matrix,momentum[k]->diagonal);
      report.momentum[k]={};
    }
    update_momentum_diagonal();
    report.pressure=coupled->solve({residual_u,residual_v,residual_w,mass,viscosity,physical_diagonal,
      face_velocity_weights(),settings.density,velocity_scale},settings.coupled_linear);
    coupled->update(velocity,pressure,settings.momentum_relaxation,settings.pressure_relaxation);
    update_boundaries();update_pressure();update_flux();
    kern::fv_continuity(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,mass,provisional_continuity);
    solve_turbulence();
    ++report.iteration;
    if(report.iteration<=5 || report.iteration%settings.monitor_interval==0) monitor();
  }

  void step() {
    if(coupled) {coupled_step();return;}
    refresh_time_step();
    update_boundaries();update_pressure();
    update_turbulence();
    kern::fv_momentum_source(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      mesh.area.x,mesh.area.y,mesh.area.z,face_pressure,mesh.volume,body.x,body.y,body.z,
      source.x,source.y,source.z);
    reconstruct_velocity();add_transpose_stress();
    TransportSettings transport;
    transport.convection=settings.convection;transport.limit=settings.limit;
    transport.density=float(settings.density);transport.pseudo_time_step=time_step;
    transport.relaxation=float(settings.momentum_relaxation);
    Buffer<float>* q[]{&velocity.x,&velocity.y,&velocity.z};
    const Buffer<float>* bc[]{&boundary_velocity.x,&boundary_velocity.y,&boundary_velocity.z};
    const Buffer<float>* src[]{&source.x,&source.y,&source.z};
    for(int k=0;k<3;++k) {
      momentum[k]->assemble_reconstructed(*q[k],*bc[k],scalar_kind,mass,viscosity,*src[k],transport,wall_terms(k));
      momentum_correction_coefficient(k);
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

  // Viscous force on each masked wall face for the active closure. Callers
  // refresh boundaries, pressure and the velocity gradients first.
  void viscous_face_forces(const Buffer<float>& mask,DeviceVec3& viscous) {
    const auto n=mask.size();const auto& m=mesh;const auto& gu=momentum[0]->grad;const auto& gv=momentum[1]->grad;const auto& gw=momentum[2]->grad;
    if(sst || rke) {
      update_turbulence();const auto& shear=sst?sst->fields().wall_shear:rke->fields().wall_shear;
      kern::fv_sst_force(dev,n,mask,m.area.x,m.area.y,m.area.z,
        shear.x,shear.y,shear.z,viscous.x,viscous.y,viscous.z,m.n_internal_faces);
    } else if(diffusion_reconstruction) {
      kern::fv_wall_viscous_force(dev,n,m.owner,mask,m.area.x,m.area.y,m.area.z,
        m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,velocity.x,velocity.y,velocity.z,
        gu.x,gu.y,gu.z,gv.x,gv.y,gv.z,gw.x,gw.y,gw.z,
        boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,viscous.x,viscous.y,viscous.z,
        m.n_internal_faces,float(settings.viscosity));
      diffusion_reconstruction->flux_gradient(velocity.x,boundary_velocity.x,m,stencil,projected_gradient[0]);
      diffusion_reconstruction->flux_gradient(velocity.y,boundary_velocity.y,m,stencil,projected_gradient[1]);
      diffusion_reconstruction->flux_gradient(velocity.z,boundary_velocity.z,m,stencil,projected_gradient[2]);
      const auto& u=projected_gradient[0];const auto& v=projected_gradient[1];const auto& w=projected_gradient[2];
      kern::fv_wall_projected_viscous_force(dev,n,mask,diffusion_reconstruction->mask(),m.area.x,m.area.y,m.area.z,
        u.x,u.y,u.z,v.x,v.y,v.z,w.x,w.y,w.z,viscous.x,viscous.y,viscous.z,m.n_internal_faces,float(settings.viscosity));
    } else {
      kern::fv_wall_viscous_force(dev,n,m.owner,mask,m.area.x,m.area.y,m.area.z,
        m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,velocity.x,velocity.y,velocity.z,
        gu.x,gu.y,gu.z,gv.x,gv.y,gv.z,gw.x,gw.y,gw.z,
        boundary_velocity.x,boundary_velocity.y,boundary_velocity.z,viscous.x,viscous.y,viscous.z,
        m.n_internal_faces,float(settings.viscosity));
      // Resolved SA walls carry the projected momentum flux; keep the
      // traction's normal derivative equal to it.
      if(sa && turbulent_faces) {
        auto correction=dev.alloc<float>(m.n_faces);auto shift=vector_field(dev,n);
        const Buffer<float>* q[]{&velocity.x,&velocity.y,&velocity.z};
        const Buffer<float>* bc[]{&boundary_velocity.x,&boundary_velocity.y,&boundary_velocity.z};
        const DeviceVec3* g[]{&gu,&gv,&gw};Buffer<float>* out[]{&shift.x,&shift.y,&shift.z};
        for(int k=0;k<3;++k) {
          turbulent_faces->diffusion_correction(*q[k],*bc[k],correction);
          kern::fv_wall_projected_shift(dev,n,m.owner,turbulent_faces->mask(),m.area.x,m.area.y,m.area.z,
            stencil.nonorthogonal.x,stencil.nonorthogonal.y,stencil.nonorthogonal.z,g[k]->x,g[k]->y,g[k]->z,correction,*out[k],m.n_internal_faces);
        }
        kern::fv_wall_normal_force_shift(dev,n,mask,m.area.x,m.area.y,m.area.z,shift.x,shift.y,shift.z,
          viscous.x,viscous.y,viscous.z,m.n_internal_faces,float(settings.viscosity));
      }
    }
  }
  IterationReport monitor() {
    refresh_time_step();
    update_boundaries();update_pressure();
    update_turbulence();
    kern::fv_momentum_source(dev,mesh.n_cells,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,
      mesh.area.x,mesh.area.y,mesh.area.z,face_pressure,mesh.volume,body.x,body.y,body.z,
      source.x,source.y,source.z);
    reconstruct_velocity();add_transpose_stress();
    const auto transport=momentum_transport();
    // The components share one CSR value buffer, so each residual must follow
    // its own assembly.
    momentum[0]->assemble_reconstructed(velocity.x,boundary_velocity.x,scalar_kind,mass,viscosity,source.x,transport,wall_terms(0));
    momentum_correction_coefficient(0);
    residual(dev,momentum[0]->matrix,velocity.x,momentum[0]->rhs,residual_u);
    momentum[1]->assemble_reconstructed(velocity.y,boundary_velocity.y,scalar_kind,mass,viscosity,source.y,transport,wall_terms(1));
    momentum_correction_coefficient(1);
    residual(dev,momentum[1]->matrix,velocity.y,momentum[1]->rhs,residual_v);
    momentum[2]->assemble_reconstructed(velocity.z,boundary_velocity.z,scalar_kind,mass,viscosity,source.z,transport,wall_terms(2));
    momentum_correction_coefficient(2);
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
    std::span<const double> k,std::span<const double> omega,std::span<const double> intermittency) {
  auto& s=*impl_;const auto n=std::size_t(s.mesh.n_cells);
  KUTTA_REQUIRE(intermittency.empty() || s.transition,"FV intermittency initialization requires the transition model");
  if(!intermittency.empty()) s.transition->initialize(intermittency,s.permutation);
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
  KUTTA_REQUIRE(!impl_->sa,"SA: hybrid turbulence initialization is not implemented");
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
    // k-epsilon specifications enter through omega = epsilon/(Cmu k), Cmu=.09.
    if(const auto* v=std::get_if<TurbulenceValues>(&spec)) {recipe[0]+=area*v->k;recipe[2]+=area*v->omega;}
    else if(const auto* v=std::get_if<KeValues>(&spec)) {recipe[0]+=area*v->k;recipe[2]+=area*v->epsilon/(.09*v->k);}
    else if(const auto* v=std::get_if<IntensityDiameter>(&spec)) {
      recipe[1]+=area*1.5*v->intensity*v->intensity;
      recipe[3]+=area*std::sqrt(1.5)*v->intensity/(std::pow(.09,.25)*.07*v->diameter);
    } else if(const auto* v=std::get_if<IntensityLength>(&spec)) {
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
    double mass_patch_area=0,mass_projected=0;
    if(const auto* in=std::get_if<MassFlowInlet>(&bc)) for(Index j=patch.start;j<patch.start+patch.size;++j) {
      mass_patch_area+=norm(s.geom.face_area[j]);
      if(in->direction && norm(*in->direction)>0) mass_projected-=dot(*in->direction/norm(*in->direction),s.geom.face_area[j]);
    }
    const auto initial_pressure=std::find_if(h.inlet_pressure.begin(),h.inlet_pressure.end(),[&](const auto& p){return p.first==patch.name;});
    for(Index f=patch.start;f<patch.start+patch.size;++f) {
      const auto b=f-s.mesh.n_internal_faces,c=s.host.owner[f];
      const auto area_vector=s.geom.face_area[f];const auto area=norm(area_vector);
      if(const auto* in=std::get_if<VelocityInlet>(&bc)) {
        const auto u=evaluate(in->velocity,s.geom.face_centroid[f]);
        flux[b]=float(dot(u,area_vector));source[c]+=flux[b]/float(s.geom.cell_volume[c]);
        normal[b]=float(dot(u,area_vector)/area*dot(s.geom.face_centroid[f]-s.geom.cell_centroid[c],area_vector)/area);
        inlet_area+=area;direction=direction+u*area;speed_sum+=norm(u)*area;
        if(s.sst || s.rke) add_recipe(*in->turbulence,area);
        if(initial_pressure!=h.inlet_pressure.end()) {pkind[b]=0;pvalue[b]=float(initial_pressure->second);pressure_sum+=area*initial_pressure->second;pressure_area+=area;}
      }
      // BEGIN FSAE physics: constant-density mass flow supplies potential flux;
      // pressure inlets supply a pressure anchor, with speed solved afterward.
      if(const auto* in=std::get_if<MassFlowInlet>(&bc)) {
        const auto u=mass_flow_velocity(*in,area_vector,mass_patch_area,mass_projected,s.settings.density);
        flux[b]=float(dot(u,area_vector));inlet_area+=area;direction=direction+u*area;speed_sum+=norm(u)*area;
        if(s.sst || s.rke) add_recipe(*in->turbulence,area);
      }
      if(const auto* in=std::get_if<PressureInlet>(&bc)) {
        kind[b]=0;pkind[b]=0;pvalue[b]=float(in->total_pressure);
        pressure_sum+=area*in->total_pressure;pressure_area+=area;
      }
      // END FSAE physics
      if(const auto* out=std::get_if<PressureOutlet>(&bc)) {
        kind[b]=0;pkind[b]=0;pvalue[b]=float(out->pressure);outlet_area+=area;
        outlet_center=outlet_center+s.geom.face_centroid[f]*area;
        pressure_sum+=area*out->pressure;pressure_area+=area;
      }
    }
  }
  KUTTA_REQUIRE(inlet_area==0 || outlet_area>0,"FV hybrid velocity inlets require a pressure outlet");
  if(inlet_area>0) {direction=direction/inlet_area;speed_sum/=inlet_area;for(auto& r:recipe) r/=inlet_area;}
  else if(s.sst || s.rke) {add_recipe(IntensityViscosityRatio{h.turbulence_intensity,h.turbulence_viscosity_ratio},1);}
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
  if(h.solve[1] && (!h.inlet_pressure.empty() || s.pressure_inlets) && outlet_area>0)
    solve_scalar(s.pressure,pkind,pvalue,std::vector<float>(nb,0),std::vector<float>(n,0),1);
  if(s.sst) {
    kern::fv_hybrid_turbulence(dev,n,s.velocity.x,s.velocity.y,s.velocity.z,s.mesh.volume,s.kinetic,s.omega,old,s.scale,
      float(recipe[0]),float(recipe[1]),float(recipe[2]),float(recipe[3]),float(recipe[4]),
      s.settings.sst->model.k_floor,s.settings.sst->model.omega_floor);
    if(h.averaged_turbulence) {const auto k=dev.sum(old,n)/volume,w=dev.sum(s.scale,n)/volume;dev.fill(s.kinetic,float(k));dev.fill(s.omega,float(w));}
  }
  if(s.rke) {
    // The SST recipe gives domain-averaged k and omega; epsilon = Cmu k omega.
    // Fluent's k-epsilon averaging is undocumented; this is startup only.
    KUTTA_REQUIRE(h.averaged_turbulence,"k-epsilon hybrid initialization requires domain-averaged turbulence");
    const auto& o=s.settings.realizable_ke->model;
    auto kinetic=dev.alloc<float>(n),omega=dev.alloc<float>(n);
    kern::fv_hybrid_turbulence(dev,n,s.velocity.x,s.velocity.y,s.velocity.z,s.mesh.volume,kinetic,omega,old,s.scale,
      float(recipe[0]),float(recipe[1]),float(recipe[2]),float(recipe[3]),float(recipe[4]),o.k_floor,o.epsilon_floor);
    const auto k=std::max(dev.sum(old,n)/volume,double(o.k_floor)),w=dev.sum(s.scale,n)/volume;
    const std::vector<double> ks(n,k),es(n,std::max(.09*k*w,double(o.epsilon_floor)));
    s.rke->initialize(ks,es,s.permutation);
  }
  s.dev.fill(s.physical_diagonal,1);s.dev.fill(s.consistent_inverse,0);
  s.update_boundaries();s.update_pressure();s.update_flux();s.update_turbulence();
  for(int pass=0;pass<3;++pass) {s.monitor();s.update_flux();s.update_boundaries();s.update_pressure();}
  if(!h.constant_speed && h.solve[0]) {
    dev.copy(s.mass,potential_flux);
    // BEGIN FSAE physics: the potential field is expressed in absolute velocity.
    s.zone_physics->relative_flux(s.boundary_kind,s.mesh.n_internal_faces,s.mass);
    // END FSAE physics
  }
  kern::fv_continuity(dev,n,s.mesh.cell_face_offsets,s.mesh.cell_faces,s.mesh.cell_face_sign,s.mass,s.provisional_continuity);
}
void SteadySolver::set_source(const std::function<Vec3d(Vec3d)>& source) {
  std::vector<Vec3d> values;for(auto p:impl_->geom.cell_centroid) values.push_back(source(p));
  upload_vectors(impl_->dev,impl_->body,values);
}
void SteadySolver::initialize_sa_fields(std::span<const Vec3d> velocity,std::span<const double> pressure,
    std::span<const double> modified_viscosity) {
  auto& s=*impl_;KUTTA_REQUIRE(s.sa,"FV: SA disabled");
  s.sa->initialize(modified_viscosity,s.permutation);
  initialize_fields(velocity,pressure);
}
void SteadySolver::initialize_ke_fields(std::span<const Vec3d> velocity,std::span<const double> pressure,
    std::span<const double> k,std::span<const double> epsilon) {
  auto& s=*impl_;KUTTA_REQUIRE(s.rke,"FV: realizable k-epsilon disabled");
  s.rke->initialize(k,epsilon,s.permutation);initialize_fields(velocity,pressure);
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
  const auto& m=s.mesh;
  kern::fv_wall_pressure_force(s.dev,selected.size(),mask,m.area.x,m.area.y,m.area.z,s.face_pressure,
    pressure.x,pressure.y,pressure.z,m.n_internal_faces,float(reference_pressure));
  s.viscous_face_forces(mask,viscous);
  return {{s.dev.sum(pressure.x,selected.size()),s.dev.sum(pressure.y,selected.size()),s.dev.sum(pressure.z,selected.size())},
    {s.dev.sum(viscous.x,selected.size()),s.dev.sum(viscous.y,selected.size()),s.dev.sum(viscous.z,selected.size())},area,reference_pressure};
}
const PolyMesh& SteadySolver::mesh() const {return impl_->host;}
const MeshGeometry& SteadySolver::geometry() const {return impl_->geom;}
const DeviceMesh& SteadySolver::device_mesh() const {return impl_->mesh;}
const DeviceVec3& SteadySolver::velocity() const {return impl_->velocity;}
const Buffer<float>& SteadySolver::pressure() const {return impl_->pressure;}
const Buffer<float>& SteadySolver::mass_flux() const {return impl_->mass;}
const Buffer<float>& SteadySolver::kinetic_energy() const {KUTTA_REQUIRE(impl_->sst || impl_->rke,"FV two-equation turbulence disabled");return impl_->rke?impl_->rke->kinetic():impl_->kinetic;}
const Buffer<float>& SteadySolver::specific_dissipation() const {KUTTA_REQUIRE(impl_->sst,"FV SST disabled");return impl_->omega;}
const DeviceVec3& SteadySolver::wall_shear() const {
  auto& s=*impl_;
  if(s.sst) return s.sst->fields().wall_shear;
  if(s.rke) return s.rke->fields().wall_shear;
  // Laminar: the same face forces as forces(), on every wall face.
  const auto nb=s.mesh.n_boundary_faces();
  if(s.laminar_wall_shear.x.size()!=std::size_t(nb)) {
    std::vector<float> walls(nb,0);
    for(const auto& p:s.host.patches) if(p.kind==PatchKind::Wall) for(Index f=p.start;f<p.start+p.size;++f) walls[f-s.mesh.n_internal_faces]=1;
    s.wall_mask=s.dev.upload<float>(walls);s.laminar_wall_shear=vector_field(s.dev,nb);
  }
  s.update_boundaries();s.update_pressure();s.reconstruct_velocity();
  auto force=vector_field(s.dev,nb);s.viscous_face_forces(s.wall_mask,force);
  kern::fv_force_to_shear(s.dev,nb,s.wall_mask,s.mesh.area.x,s.mesh.area.y,s.mesh.area.z,force.x,force.y,force.z,
    s.laminar_wall_shear.x,s.laminar_wall_shear.y,s.laminar_wall_shear.z,s.mesh.n_internal_faces);
  return s.laminar_wall_shear;
}
const SstFields& SteadySolver::turbulence() const {KUTTA_REQUIRE(impl_->sst,"FV SST disabled");return impl_->sst->fields();}
const Buffer<float>& SteadySolver::intermittency() const {KUTTA_REQUIRE(impl_->transition,"FV transition model disabled");return impl_->transition->gamma();}
const IntermittencyFields& SteadySolver::transition() const {KUTTA_REQUIRE(impl_->transition,"FV transition model disabled");return impl_->transition->fields();}
const SaFields& SteadySolver::sa_turbulence() const {KUTTA_REQUIRE(impl_->sa,"FV SA disabled");return impl_->sa->fields();}
const SaWallFields& SteadySolver::sa_wall() const {KUTTA_REQUIRE(impl_->sa,"FV SA disabled");return impl_->sa->wall_fields();}
const RealizableKeFields& SteadySolver::realizable_turbulence() const {KUTTA_REQUIRE(impl_->rke,"FV realizable k-epsilon disabled");return impl_->rke->fields();}
const Buffer<float>& SteadySolver::dissipation() const {KUTTA_REQUIRE(impl_->rke,"FV realizable k-epsilon disabled");return impl_->rke->epsilon();}
const Buffer<float>& SteadySolver::modified_viscosity() const {KUTTA_REQUIRE(impl_->sa,"FV SA disabled");return impl_->sa->modified_viscosity();}
const MeshPermutation& SteadySolver::permutation() const {return impl_->permutation;}
const Buffer<float>& SteadySolver::boundary_pressure() const {return impl_->boundary_pressure;}
const IterationReport& SteadySolver::report() const {return impl_->report;}
std::size_t SteadySolver::coupled_solver_bytes() const {return impl_->coupled?impl_->coupled->bytes():0;}
Vec3d ForceReport::force() const {return pressure+viscous;}
double ForceReport::coefficient(Vec3d direction,const ReferenceValues& ref) const {
  KUTTA_REQUIRE(ref.density>0 && ref.velocity>0 && ref.area>0,"FV positive force references required");
  const auto adjusted=force()+area_vector*(reference_pressure-ref.pressure);
  return dot(adjusted,direction)/(.5*ref.density*ref.velocity*ref.velocity*ref.area);
}
}  // namespace kutta::fv
