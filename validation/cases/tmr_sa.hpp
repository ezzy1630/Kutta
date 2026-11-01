#pragma once
#include "tmr.hpp"

namespace kutta::validation::tmr::sa {
inline Sample sample(Device& dev,SteadySolver& solver,Case c) {
  const auto& t=solver.sa_wall();const auto tx=dev.download(t.shear.x),ty=dev.download(t.shear.y),yp=dev.download(t.y_plus);
  const auto pressure=dev.download(solver.boundary_pressure());
  std::vector<std::array<double,2>> cf,cp;std::vector<std::string> walls;Sample result;
  for(const auto& patch:solver.mesh().patches) if(patch.kind==PatchKind::Wall) {
    walls.push_back(patch.name);
    for(Index f=patch.start;f<patch.start+patch.size;++f) {
      const auto b=f-solver.mesh().n_internal_faces();const auto a=solver.geometry().face_area[f];
      const double x=solver.geometry().face_centroid[f].x;
      cf.push_back({x,2*(a.y*tx[b]-a.x*ty[b])/std::hypot(a.x,a.y)});cp.push_back({x,2*pressure[b]});
      result.max_yplus=std::max(result.max_yplus,double(yp[b]));
    }
  }
  std::sort(cf.begin(),cf.end());std::sort(cp.begin(),cp.end());
  for(int i=0;i<(c==Case::FlatPlate?1:3);++i) {
    const double x=c==Case::FlatPlate?plate_station:bump_stations[i];
    result.cf[i]=interpolate(cf,x);result.cp[i]=interpolate(cp,x);
  }
  result.drag=solver.forces(walls).coefficient({1,0,0},{1,1,(c==Case::FlatPlate?2.:1.5)*span,0});return result;
}
// Stops on convergence or after max_iterations, never on wall time.
inline Result run(Device& dev,const std::filesystem::path& root,Case c,int level,int max_iterations=3000) {
  using Clock=std::chrono::steady_clock;auto start=Clock::now();auto path=grid(root,c,level);
  KUTTA_REQUIRE(std::filesystem::exists(path),"SA TMR grid missing; run tools/prepare_tmr.py");
  auto mesh=read_plot3d(path,span);const double nu=c==Case::FlatPlate?1./5e6:1./3e6;
  std::vector<PatchCondition> bc;
  for(const auto& p:mesh.patches) {
    if(p.kind==PatchKind::Wall) bc.push_back({p.name,Wall{}});
    else if(p.source_type=="subsonic_inflow_pt") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0},ModifiedViscosity{3*nu}}});
    else if(p.source_type=="back_pressure" || p.source_type=="farfield_riem")
      bc.push_back({p.name,PressureOutlet{0,Vec3d{1,0,0},BackflowPressure::Static,ModifiedViscosity{3*nu}}});
    else if(p.kind==PatchKind::Symmetry) bc.push_back({p.name,Symmetry{}});
    else throw Error("SA TMR unmapped patch: "+p.name);
  }
  SolverSettings settings;settings.viscosity=nu;settings.sa=SaSettings{};
  settings.sa->model.formulation=SaFormulation::TmrStandard;
  settings.sa->convection=ConvectionScheme::FirstOrderUpwind;settings.monitor_interval=100;
  if(c==Case::Bump) {settings.pseudo_time_step=.01;settings.momentum_relaxation=.5;settings.pressure_relaxation=.3;}
  SteadySolver solver(dev,std::move(mesh),std::move(bc),settings);const auto n=solver.mesh().n_cells;
  solver.initialize_sa_fields(std::vector<Vec3d>(n,{1,0,0}),std::vector<double>(n,0),std::vector<double>(n,3*nu));
  dev.synchronize();Result result;result.cells=n;result.setup_seconds=std::chrono::duration<double>(Clock::now()-start).count();
  start=Clock::now();Sample previous;bool has_previous=false;
  while(solver.report().iteration<max_iterations) {
    try {solver.iterate(std::min(100,max_iterations-solver.report().iteration));result.residual=solver.monitor();result.values=sa::sample(dev,solver,c);}
    catch(const Error& e) {result.failure=e.what();result.iterations=solver.report().iteration;break;}
    result.iterations=result.residual.iteration;result.load_change=has_previous?change(previous,result.values,c):std::numeric_limits<double>::infinity();
    const auto& r=result.residual;const double defect=std::max({r.u.scaled,r.v.scaled,r.nu_tilde.scaled,r.continuity.scaled});
    result.converged=has_previous && defect<=1e-5 && result.load_change<=.001;
    std::cout<<std::setprecision(12)<<"SA TMR progress "<<backend_name(dev.backend())<<' '<<(c==Case::FlatPlate?"plate":"bump")<<' '<<level<<" it="<<r.iteration<<" defect="<<defect<<" load_change="<<result.load_change<<" Cf="<<result.values.cf[0]<<" Cd="<<result.values.drag<<std::endl;
    if(result.converged || !std::isfinite(defect) || !std::isfinite(result.values.drag)) break;
    previous=result.values;has_previous=true;
  }
  dev.synchronize();result.solve_seconds=std::chrono::duration<double>(Clock::now()-start).count();return result;
}
} // namespace kutta::validation::tmr::sa
