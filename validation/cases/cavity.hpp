#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"

namespace kutta::validation {
using Profile=std::vector<std::pair<double,double>>;
struct CavityResult {
  Profile u,v;
  fv::IterationReport report;
  double update=0;
  bool converged=false;
};

inline double interpolate(const Profile& profile,double coordinate) {
  auto next=std::lower_bound(profile.begin(),profile.end(),coordinate,
    [](const auto& p,double x){return p.first<x;});
  if(next==profile.begin()) return next->second;
  if(next==profile.end()) return profile.back().second;
  const auto& previous=*(next-1);
  const double weight=(coordinate-previous.first)/(next->first-previous.first);
  return previous.second+weight*(next->second-previous.second);
}

inline Profile centerline(const MeshGeometry& g,const std::vector<float>& value,bool vertical) {
  std::map<long long,Profile> lines;
  for(std::size_t c=0;c<value.size();++c) {
    auto p=g.cell_centroid[c];
    lines[std::llround((vertical?p.y:p.x)*1e10)].push_back({vertical?p.x:p.y,value[c]});
  }
  Profile result;
  for(auto& [key,line]:lines) {
    std::sort(line.begin(),line.end());
    result.push_back({key/1e10,interpolate(line,.5)});
  }
  result.insert(result.begin(),{0,0});result.push_back({1,vertical?1.:0.});
  return result;
}

inline std::vector<std::vector<double>> read_csv(const std::filesystem::path& path) {
  std::ifstream input(path);KUTTA_REQUIRE(input.good(),"Missing cavity data: "+path.string());
  std::string line;std::getline(input,line);
  std::vector<std::vector<double>> rows;
  while(std::getline(input,line)) {
    std::stringstream stream(line);std::string token;std::vector<double> row;
    while(std::getline(stream,token,',')) row.push_back(std::stod(token));
    if(!row.empty()) rows.push_back(std::move(row));
  }
  return rows;
}

inline double reference_error(const Profile& profile,const std::filesystem::path& file,
                              int coordinate,int value,bool omit_suspect=false) {
  double error=0;
  for(const auto& row:read_csv(file)) {
    if(omit_suspect && std::abs(row[coordinate]-.90625)<1e-9) continue;
    error=std::max(error,std::abs(interpolate(profile,row[coordinate])-row[value]));
  }
  return error;
}

inline CavityResult cavity(Device& dev,PolyMesh mesh,int reynolds,int maximum_iterations=6000) {
  std::vector<fv::PatchCondition> boundaries;
  for(const auto& patch:mesh.patches) {
    if(patch.kind==PatchKind::Symmetry || patch.name=="zmin" || patch.name=="zmax")
      boundaries.push_back({patch.name,fv::Symmetry{}});
    else boundaries.push_back({patch.name,fv::Wall{{patch.name=="ymax" || patch.name.find("moving")!=std::string::npos?1.:0.,0,0}}});
  }
  fv::SolverSettings settings;settings.viscosity=1./reynolds;
  settings.pseudo_time_scale=5;settings.momentum_relaxation=1;settings.pressure_relaxation=1;settings.monitor_interval=100;
  fv::SteadySolver solver(dev,std::move(mesh),boundaries,settings);
  CavityResult result;std::vector<float> old_u,old_v;
  for(int iteration=0;iteration<maximum_iterations;iteration+=100) {
    solver.iterate(std::min(100,maximum_iterations-iteration));
    auto u=dev.download(solver.velocity().x),v=dev.download(solver.velocity().y);
    result.report=solver.report();result.update=1;
    if(!old_u.empty()) {
      result.update=0;
      for(std::size_t c=0;c<u.size();++c)
        result.update=std::max({result.update,std::abs(double(u[c]-old_u[c])),std::abs(double(v[c]-old_v[c]))});
    }
    std::cout<<"FV_CAVITY_PROGRESS backend="<<backend_name(dev.backend())<<" Re="<<reynolds
      <<" cells="<<solver.mesh().n_cells<<" iterations="<<result.report.iteration
      <<" ru="<<result.report.u.scaled<<" rv="<<result.report.v.scaled
      <<" rc="<<result.report.continuity.scaled<<" mass="<<result.report.corrected_mass_imbalance
      <<" update="<<result.update<<std::endl;
    const bool steady=std::max(result.report.u.scaled,result.report.v.scaled)<5e-7 &&
      result.update<2e-5 && result.report.corrected_mass_imbalance<1e-5;
    if(steady || iteration+100>=maximum_iterations) {
      result.converged=steady;
      result.u=centerline(solver.geometry(),u,true);result.v=centerline(solver.geometry(),v,false);
      break;
    }
    old_u=std::move(u);old_v=std::move(v);
  }
  return result;
}

inline void compare_cavity(const CavityResult& result,const std::filesystem::path& data,int reynolds,
                           Backend backend,const std::string& label) {
  const int column=reynolds==100?3:(reynolds==400?4:5);
  const double eu=reference_error(result.u,data/"ghia1982_cavity/ghia1982_u_vertical_centerline.csv",2,column);
  const double ev=reference_error(result.v,data/"ghia1982_cavity/ghia1982_v_horizontal_centerline.csv",2,column,reynolds==400);
  std::cout<<"FV_CAVITY_RESULT backend="<<backend_name(backend)<<" case="<<label<<" Re="<<reynolds
    <<" iterations="<<result.report.iteration<<" converged="<<result.converged
    <<" Ghia_u_max="<<eu<<" Ghia_v_max="<<ev;
  if(reynolds==1000) {
    std::cout<<" Erturk_u_max="<<reference_error(result.u,data/"erturk2005_cavity/erturk2005_u_vertical_centerline_Re1000.csv",0,1)
      <<" Erturk_v_max="<<reference_error(result.v,data/"erturk2005_cavity/erturk2005_v_horizontal_centerline_Re1000.csv",0,1);
  }
  std::cout<<std::endl;
}
}  // namespace kutta::validation
