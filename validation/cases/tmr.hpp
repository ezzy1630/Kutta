#pragma once

// TMR setup and host-only diagnostics. The solve remains in SteadySolver/Device.
#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string_view>

#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"

namespace kutta::validation::tmr {
using namespace fv;
enum class Case { FlatPlate, Bump };
inline constexpr double span=.01;
inline constexpr double plate_station=.970084071;
inline constexpr std::array<double,3> bump_stations{.6321975,.75,.8678025};
struct Zone {std::string name;std::vector<std::vector<double>> rows;};
inline std::vector<Zone> read_reference(const std::filesystem::path& path) {
  std::ifstream in(path);KUTTA_REQUIRE(in,"TMR cannot open reference: "+path.string());
  std::vector<Zone> zones;std::string line;
  while(std::getline(in,line)) {
    const auto first=line.find_first_not_of(" \t\r");if(first==std::string::npos) continue;
    line.erase(0,first);auto lower=line;
    std::transform(lower.begin(),lower.end(),lower.begin(),[](unsigned char ch){return static_cast<char>(std::tolower(ch));});
    if(lower.starts_with("zone")) {zones.push_back({line,{}});continue;}
    if(line[0]=='#' || lower.starts_with("variables")) continue;
    std::istringstream row(line);std::vector<double> values;double value;
    while(row>>value) {KUTTA_REQUIRE(std::isfinite(value),"TMR nonfinite reference");values.push_back(value);}
    KUTTA_REQUIRE(!values.empty() && row.eof() && !zones.empty(),"TMR invalid reference row: "+line);
    KUTTA_REQUIRE(zones.back().rows.empty() || values.size()==zones.back().rows.front().size(),"TMR inconsistent reference columns");
    zones.back().rows.push_back(std::move(values));
  }
  KUTTA_REQUIRE(!zones.empty(),"TMR missing reference zones");return zones;
}
inline double interpolate(const std::vector<std::array<double,2>>& points,double x) {
  KUTTA_REQUIRE(points.size()>=2 && x>=points.front()[0] && x<=points.back()[0],"TMR station outside sampled range");
  const auto high=std::upper_bound(points.begin(),points.end(),x,[](double a,const auto& b){return a<b[0];});
  if(high==points.end()) return points.back()[1];
  const auto& lo=*(high-1);return lo[1]+(high->at(1)-lo[1])*(x-lo[0])/(high->at(0)-lo[0]);
}
inline double reference_station(const Zone& zone,double x) {
  std::vector<std::array<double,2>> points;
  for(const auto& row:zone.rows) {KUTTA_REQUIRE(row.size()==2,"TMR expected two-column distribution");points.push_back({row[0],row[1]});}
  std::sort(points.begin(),points.end());return interpolate(points,x);
}
inline std::filesystem::path grid(const std::filesystem::path& root,Case c,int level) {
  // Plate level 3 is the 273x193 grid (tools/prepare_tmr.py expands it); used
  // for the transition model's leading-edge mesh check.
  if(c==Case::FlatPlate && level==3) return root/"2DZP_flatplate/grids/flatplate_clust2_1leveldown_273x193.p2dfmt";
  KUTTA_REQUIRE(level>=0 && level<=2,"TMR level must be 0..2 (plate: 0..3)");
  const std::array<std::string,3> plate{"flatplate_clust2_4levelsdown_35x25","flatplate_clust2_3levelsdown_69x49","flatplate_clust2_2levelsdown_137x97"};
  const std::array<std::string,3> bump{"bump_4levelsdown_89x41","bump_3levelsdown_177x81","bump_2levelsdown_353x161"};
  return root/(c==Case::FlatPlate?"2DZP_flatplate":"2DBC_bump")/"grids"/((c==Case::FlatPlate?plate:bump)[level]+".p2dfmt");
}
struct Options {
  int max_iterations=3000;
  double max_seconds=600;
  bool symmetry_top=false;
  SstOptions model;
  // Inflow k and omega instead of the TMR freestream (the fully turbulent
  // transition regression raises the freestream turbulence).
  std::optional<TurbulenceValues> inflow;
  std::optional<IntermittencySettings> transition;
  std::filesystem::path output;
};
struct Sample {
  std::array<double,3> cf{},cp{};
  double drag=0,max_yplus=0;
};
inline Sample sample(Device& dev,SteadySolver& solver,Case c) {
  const auto tau_x=dev.download(solver.turbulence().wall_shear.x);
  const auto tau_y=dev.download(solver.turbulence().wall_shear.y);
  const auto yp=dev.download(solver.turbulence().wall_y_plus_utau);
  const auto pressure=dev.download(solver.boundary_pressure());
  std::vector<std::array<double,2>> cf,cp;std::vector<std::string> walls;
  Sample result;
  for(const auto& patch:solver.mesh().patches) if(patch.kind==PatchKind::Wall) {
    walls.push_back(patch.name);
    for(Index f=patch.start;f<patch.start+patch.size;++f) {
      const auto b=f-solver.mesh().n_internal_faces();const auto a=solver.geometry().face_area[f];
      const double length=std::hypot(a.x,a.y),x=solver.geometry().face_centroid[f].x;
      // Positive tangent follows increasing x on the lower wall (outward normal down).
      cf.push_back({x,2*(a.y*tau_x[b]-a.x*tau_y[b])/length});
      cp.push_back({x,2*pressure[b]});result.max_yplus=std::max(result.max_yplus,double(yp[b]));
    }
  }
  std::sort(cf.begin(),cf.end());std::sort(cp.begin(),cp.end());
  const int count=c==Case::FlatPlate?1:3;
  for(int i=0;i<count;++i) {const double x=c==Case::FlatPlate?plate_station:bump_stations[i];result.cf[i]=interpolate(cf,x);result.cp[i]=interpolate(cp,x);}
  result.drag=solver.forces(walls).coefficient({1,0,0},{1,1,(c==Case::FlatPlate?2.:1.5)*span,0});
  return result;
}
inline double change(const Sample& a,const Sample& b,Case c) {
  double result=std::abs(a.drag-b.drag)/std::max(std::abs(b.drag),1e-12);
  for(int i=0;i<(c==Case::FlatPlate?1:3);++i) {
    result=std::max(result,std::abs(a.cf[i]-b.cf[i])/std::max(std::abs(b.cf[i]),1e-12));
    if(c==Case::Bump) result=std::max(result,std::abs(a.cp[i]-b.cp[i])/std::max(std::abs(b.cp[i]),.03));
  }
  return result;
}
inline double defect(const IterationReport& r) {
  return std::max({r.u.scaled,r.v.scaled,r.k.scaled,r.omega.scaled,r.continuity.scaled,r.intermittency.scaled});
}
struct Result {
  Index cells=0;int iterations=0;double setup_seconds=0,solve_seconds=0,load_change=0;
  bool converged=false,stopped_for_budget=false;std::string failure;Sample values;IterationReport residual;
};
inline Result run(Device& dev,const std::filesystem::path& root,Case c,int level,const Options& options={}) {
  using Clock=std::chrono::steady_clock;auto start=Clock::now();
  auto path=grid(root,c,level);KUTTA_REQUIRE(std::filesystem::exists(path),"TMR grid missing; run python3 tools/prepare_tmr.py: "+path.string());auto mesh=read_plot3d(path,span);
  std::vector<PatchCondition> bc;const double nu=c==Case::FlatPlate?1./5e6:1./3e6;
  const auto turbulence=options.inflow.value_or(TurbulenceValues{2.25e-7,25e-6/nu});
  for(const auto& p:mesh.patches) {
    if(p.kind==PatchKind::Wall) bc.push_back({p.name,Wall{}});
    else if(p.source_type=="subsonic_inflow_pt") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0},turbulence}});
    else if(p.source_type=="back_pressure" || (p.source_type=="farfield_riem" && !options.symmetry_top))
      bc.push_back({p.name,PressureOutlet{0,Vec3d{1,0,0},BackflowPressure::Static,turbulence}});
    else if(p.kind==PatchKind::Symmetry || (p.source_type=="farfield_riem" && options.symmetry_top)) bc.push_back({p.name,Symmetry{}});
    else throw Error("TMR unmapped boundary: "+p.name+" ("+p.source_type+")");
  }
  SolverSettings settings;settings.viscosity=nu;settings.sst=SstSettings{};settings.sst->model=options.model;settings.sst->transition=options.transition;
  settings.sst->k_convection=ConvectionScheme::FirstOrderUpwind;settings.sst->omega_convection=ConvectionScheme::FirstOrderUpwind;
  settings.monitor_interval=100;
  // The bump domain extends 25 lengths upstream/downstream. Its global
  // automatic time scale is too large for uniform no-slip-wall startup.
  if(c==Case::Bump) {settings.pseudo_time_step=.01;settings.momentum_relaxation=.5;settings.pressure_relaxation=.3;}
  SteadySolver solver(dev,std::move(mesh),std::move(bc),settings);
  const auto n=solver.mesh().n_cells;std::vector<Vec3d> u(n,Vec3d{1,0,0});
  std::vector<double> p(n,0),k(n,turbulence.k),w(n,turbulence.omega);solver.initialize_fields(u,p,k,w);dev.synchronize();
  Result result;result.cells=n;result.setup_seconds=std::chrono::duration<double>(Clock::now()-start).count();start=Clock::now();
  std::ofstream history;
  if(!options.output.empty()) {std::filesystem::create_directories(options.output);history.open(options.output/"history.csv");history<<"iteration,u,v,k,omega,continuity,load_change,drag,cf0,cf1,cf2,cp0,cp1,cp2,pressure_linear_ratio,continuity_raw,continuity_scale,corrected_mass,maximum_velocity,pressure_iterations\n"<<std::setprecision(12);}
  Sample previous;bool has_previous=false;
  while(solver.report().iteration<options.max_iterations) {
    try {
      solver.iterate(std::min(100,options.max_iterations-solver.report().iteration));
      result.residual=solver.monitor();result.values=sample(dev,solver,c);
    } catch(const Error& error) {
      result.failure=error.what();result.iterations=solver.report().iteration;result.load_change=std::numeric_limits<double>::quiet_NaN();
      result.residual=solver.report();result.values.drag=std::numeric_limits<double>::quiet_NaN();
      result.values.cf.fill(std::numeric_limits<double>::quiet_NaN());result.values.cp.fill(std::numeric_limits<double>::quiet_NaN());
      break;
    }
    result.iterations=result.residual.iteration;
    result.load_change=has_previous?change(previous,result.values,c):std::numeric_limits<double>::infinity();
    result.converged=has_previous && defect(result.residual)<=1e-5 && result.load_change<=.001;
    if(history) {
      const auto& r=result.residual;history<<r.iteration<<','<<r.u.scaled<<','<<r.v.scaled<<','<<r.k.scaled<<','<<r.omega.scaled<<','<<r.continuity.scaled<<','<<result.load_change<<','<<result.values.drag;
      for(auto q:result.values.cf) history<<','<<q;
      for(auto q:result.values.cp) history<<','<<q;
      history<<','<<r.pressure.final_residual/std::max(r.pressure.initial_residual,1e-30)<<','<<r.continuity.raw<<','<<r.continuity.fluent_normalization<<','<<r.corrected_mass_imbalance<<','<<r.maximum_velocity<<','<<r.pressure.iterations<<'\n';history.flush();
    }
    std::cout<<std::setprecision(12)<<"TMR progress "<<backend_name(dev.backend())<<' '<<(c==Case::FlatPlate?"plate":"bump")<<' '<<level<<" it="<<result.iterations<<" defect="<<defect(result.residual)<<" change="<<result.load_change<<" cf="<<result.values.cf[0]<<" drag="<<result.values.drag<<std::endl;
    if(!std::isfinite(defect(result.residual)) || !std::isfinite(result.values.drag)) break;
    if(result.converged) break;
    if(std::chrono::duration<double>(Clock::now()-start).count()>=options.max_seconds) {result.stopped_for_budget=true;break;}
    previous=result.values;has_previous=true;
  }
  dev.synchronize();result.solve_seconds=std::chrono::duration<double>(Clock::now()-start).count();
  if(!options.output.empty() && result.failure.empty()) {
    auto x=dev.download(solver.velocity().x),y=dev.download(solver.velocity().y),mt=dev.download(solver.turbulence().mu_t);
    auto ke=dev.download(solver.kinetic_energy()),om=dev.download(solver.specific_dissipation());
    const auto gamma=options.transition?dev.download(solver.intermittency()):std::vector<float>{};
    std::ofstream fields(options.output/"cells.csv");fields<<"original_cell,x,y,u,v,mut_over_mu,k,omega"<<(options.transition?",intermittency":"")<<'\n'<<std::setprecision(12);
    for(Index i=0;i<n;++i) {
      const auto point=solver.geometry().cell_centroid[i];fields<<solver.permutation().cells_new_to_old[i]<<','<<point.x<<','<<point.y<<','<<x[i]<<','<<y[i]<<','<<mt[i]/nu<<','<<ke[i]<<','<<om[i];
      if(options.transition) fields<<','<<gamma[i];
      fields<<'\n';
    }
    const auto tau_x=dev.download(solver.turbulence().wall_shear.x),tau_y=dev.download(solver.turbulence().wall_shear.y);
    std::ofstream wall(options.output/"wall_cf.csv");wall<<"patch,x,cf\n"<<std::setprecision(12);
    for(const auto& patch:solver.mesh().patches) if(patch.kind==PatchKind::Wall)
      for(Index f=patch.start;f<patch.start+patch.size;++f) {
        const auto b=f-solver.mesh().n_internal_faces();const auto a=solver.geometry().face_area[f];
        wall<<patch.name<<','<<solver.geometry().face_centroid[f].x<<','<<2*(a.y*tau_x[b]-a.x*tau_y[b])/std::hypot(a.x,a.y)<<'\n';
      }
  }
  return result;
}
} // namespace kutta::validation::tmr
