#pragma once

#include <cmath>
#include <numbers>
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"

namespace kutta::bench {
// Freudenthal tetrahedra in smoothly graded wall layers. This exercises the
// unstructured gather path, not the poly-hexcore/prism target topology.
inline PolyMesh sst_channel(int nx,int ny,int nz) {
  auto mesh=make_box_tet(nx,ny,nz,{0,0,0},{4,1,.5});
  for(auto& p:mesh.points) {
    p.y=.5*(1+std::tanh(1.2*(2*p.y-1))/std::tanh(1.2));
    // Interior skew leaves all six planar boundary patches intact.
    p.x+=.015*std::sin(std::numbers::pi*p.x/4)*std::sin(2*std::numbers::pi*p.y)*
      std::sin(2*std::numbers::pi*p.z);
  }
  return mesh;
}
inline std::vector<fv::PatchCondition> channel_boundaries(const PolyMesh& mesh) {
  std::vector<fv::PatchCondition> bc;
  for(const auto& p:mesh.patches) {
    if(p.name=="xmin") bc.push_back({p.name,fv::VelocityInlet{Vec3d{1,0,0},fv::TurbulenceValues{.0015,10}}});
    else if(p.name=="xmax") bc.push_back({p.name,fv::PressureOutlet{0,fv::NormalBackflow{},fv::BackflowPressure::Static,
      fv::TurbulenceValues{.0015,10}}});
    else if(p.name=="ymin" || p.name=="ymax") bc.push_back({p.name,fv::Wall{}});
    else bc.push_back({p.name,fv::Symmetry{}});
  }
  return bc;
}
inline fv::SolverSettings channel_settings() {
  fv::SolverSettings s;s.viscosity=.001;s.sst=fv::SstSettings{};
  s.pseudo_time_step=.01;s.monitor_interval=20;
  s.momentum_relaxation=.3;s.pressure_relaxation=.05;
  s.sst->k_convection=fv::ConvectionScheme::FirstOrderUpwind;
  s.sst->omega_convection=fv::ConvectionScheme::FirstOrderUpwind;
  return s;
}
inline void initialize_channel(fv::SteadySolver& solver) {
  const auto n=solver.mesh().n_cells;
  solver.initialize_fields(std::vector<Vec3d>(n,{1,0,0}),std::vector<double>(n,0),
    std::vector<double>(n,.0015),std::vector<double>(n,10));
}
} // namespace kutta::bench
