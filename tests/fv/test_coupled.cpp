#include <cmath>

#include "../backend_helpers.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"

using namespace kutta;
using namespace kutta::fv;

namespace {
struct Fields {std::vector<float> u,v,p;int iterations=0;};
// Iterates until the steady momentum and raw continuity defects stop
// changing; both schemes share this exit, so neither is favoured.
Fields converge(Device& dev,const PolyMesh& mesh,const std::vector<PatchCondition>& bc,SolverSettings s,
                Vec3d initial) {
  s.monitor_interval=25;
  SteadySolver solver(dev,mesh,bc,s);
  solver.initialize([initial](Vec3d){return initial;},[](Vec3d){return 0.;});
  for(int i=0;i<4000;i+=25) {
    solver.iterate(25);const auto r=solver.monitor();
    REQUIRE(std::isfinite(r.u.scaled));
    if(std::max(r.u.scaled,r.v.scaled)<2e-7 && r.continuity.raw<1e-7) break;
  }
  return {dev.download(solver.velocity().x),dev.download(solver.velocity().y),dev.download(solver.pressure()),
          solver.report().iteration};
}
double relative(const std::vector<float>& a,const std::vector<float>& b) {
  double d=0,r=0;for(std::size_t i=0;i<a.size();++i) {d+=std::pow(double(a[i])-b[i],2);r+=double(b[i])*b[i];}
  return std::sqrt(d/std::max(r,1e-300));
}
}

TEST_CASE("FV Coupled reaches SIMPLEC's fixed point in an open channel and a closed cavity") {
  for(bool cavity:{false,true}) {
    CAPTURE(cavity);
    auto mesh=make_box_hex(12,cavity?12:6,1,{0,0,0},{1,cavity?1.:.5,.1});
    std::vector<PatchCondition> bc;
    for(const auto& p:mesh.patches) {
      if(p.name=="zmin" || p.name=="zmax") bc.push_back({p.name,Symmetry{}});
      else if(cavity) bc.push_back({p.name,p.name=="ymax"?Wall{Vec3d{1,0,0}}:Wall{}});
      else if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0}}});
      else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0}});
      else bc.push_back({p.name,p.name=="ymin"?BoundaryCondition{Wall{}}:BoundaryCondition{Symmetry{}}});
    }
    SolverSettings simplec;simplec.viscosity=.01;simplec.reference_point={.5,.5,.05};
    auto coupled=simplec;coupled.coupling=PressureVelocityCoupling::Coupled;
    coupled.momentum_relaxation=coupled.pressure_relaxation=.75;
    const Vec3d start=cavity?Vec3d{}:Vec3d{1,0,0};
    std::vector<float> first;
    KUTTA_FOR_EACH_BACKEND(dev) {
      const auto a=converge(*dev,mesh,bc,simplec,start),b=converge(*dev,mesh,bc,coupled,start);
      MESSAGE("iterations simplec="<<a.iterations<<" coupled="<<b.iterations);
      CHECK(b.iterations<4000);
      // One discretization: only the float convergence floors differ.
      CHECK(relative(b.u,a.u)<2e-5);CHECK(relative(b.v,a.v)<2e-4);
      double scale=0;for(auto x:a.p) scale=std::max(scale,double(std::abs(x)));
      for(std::size_t c=0;c<a.p.size();++c) CHECK(std::abs(b.p[c]-a.p[c])<1e-4*scale);
      if(first.empty()) first=b.u;
      else CHECK(relative(b.u,first)<1e-5);  // backends agree
    }
  }
}

// Tetrahedral boundary cells feel their own pressure through the least-squares
// face extrapolation; without that term in the Jacobian this case stalls at
// raw continuity ~1e-4 and larger tetrahedral channels diverge.
TEST_CASE("FV Coupled converges on a tetrahedral channel") {
  const auto mesh=make_box_tet(12,4,2,{0,0,0},{4,1,.5});
  std::vector<PatchCondition> bc;
  for(const auto& p:mesh.patches) {
    if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0}}});
    else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0}});
    else bc.push_back({p.name,p.name=="ymin"?BoundaryCondition{Wall{}}:BoundaryCondition{Symmetry{}}});
  }
  SolverSettings s;s.viscosity=.01;s.coupling=PressureVelocityCoupling::Coupled;
  s.momentum_relaxation=s.pressure_relaxation=.5;s.monitor_interval=50;
  std::vector<float> first;
  KUTTA_FOR_EACH_BACKEND(dev) {
    SteadySolver solver(*dev,mesh,bc,s);
    solver.initialize([](Vec3d){return Vec3d{1,0,0};},[](Vec3d){return 0.;});
    solver.iterate(400);const auto r=solver.monitor();
    CHECK(std::max(r.u.scaled,r.v.scaled)<2e-7);CHECK(r.continuity.raw<3e-6);
    const auto u=dev->download(solver.velocity().x);
    if(first.empty()) first=u;
    else CHECK(relative(u,first)<1e-5);
  }
}

// The flow Courant number (SolverSettings::coupled_courant) relaxes the
// momentum Jacobian only: it changes the path, never the fixed point.
TEST_CASE("FV Coupled flow Courant number changes the path, not the fixed point") {
  auto mesh=make_box_hex(12,6,1,{0,0,0},{1,.5,.1});
  std::vector<PatchCondition> bc;
  for(const auto& p:mesh.patches) {
    if(p.name=="zmin" || p.name=="zmax" || p.name=="ymax") bc.push_back({p.name,Symmetry{}});
    else if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0}}});
    else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0}});
    else bc.push_back({p.name,Wall{}});
  }
  SolverSettings s;s.viscosity=.01;s.coupling=PressureVelocityCoupling::Coupled;
  s.momentum_relaxation=s.pressure_relaxation=.75;
  auto low=s;low.coupled_courant=5;
  KUTTA_FOR_EACH_BACKEND(dev) {
    const auto a=converge(*dev,mesh,bc,s,{1,0,0}),b=converge(*dev,mesh,bc,low,{1,0,0});
    CHECK(relative(b.u,a.u)<2e-5);CHECK(relative(b.v,a.v)<2e-4);
    std::vector<float> early[2];
    for(int k=0;k<2;++k) {
      SteadySolver solver(*dev,mesh,bc,k?low:s);
      solver.initialize([](Vec3d){return Vec3d{1,0,0};},[](Vec3d){return 0.;});
      solver.iterate(5);early[k]=dev->download(solver.velocity().x);
    }
    CHECK(relative(early[1],early[0])>1e-3);
  }
}
