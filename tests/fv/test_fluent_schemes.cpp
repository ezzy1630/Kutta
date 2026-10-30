#include <algorithm>
#include <cmath>
#include "../backend_helpers.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"

using namespace kutta;
using namespace kutta::fv;

namespace {
struct Channel {
  std::vector<float> u,v;
  IterationReport report;
};
// Developing channel flow, Re 20, on tets so wall cells have uneven stencils.
Channel channel(Device& dev,const SolverSettings& settings,int iterations) {
  auto mesh=make_box_tet(8,4,1,{0,0,0},{2,1,.2});
  std::vector<PatchCondition> bc;
  for(const auto& p:mesh.patches) {
    if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0}}});
    else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0}});
    else if(p.name=="zmin" || p.name=="zmax") bc.push_back({p.name,Symmetry{}});
    else bc.push_back({p.name,Wall{}});
  }
  SteadySolver solver(dev,std::move(mesh),bc,settings);
  solver.initialize([](Vec3d){return Vec3d{};},[](Vec3d){return 0.;});
  solver.iterate(iterations);
  return {dev.download(solver.velocity().x),dev.download(solver.velocity().y),solver.monitor()};
}
double max_difference(const Channel& a,const Channel& b) {
  double out=0;
  for(std::size_t i=0;i<a.u.size();++i) out=std::max({out,double(std::abs(a.u[i]-b.u[i])),double(std::abs(a.v[i]-b.v[i]))});
  return out;
}
}

TEST_CASE("FV SIMPLE and classic steady iteration reach the pseudo-time SIMPLEC fixed point") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    SolverSettings base;base.viscosity=.05;
    const auto reference=channel(*dev,base,400);
    REQUIRE(reference.report.continuity.raw<1e-6);
    // Fluent's classic steady SIMPLE: no pseudo-time term, implicit momentum
    // relaxation and explicit pressure relaxation only. Rhie-Chow uses the
    // unrelaxed diagonal, so the fixed point is the same.
    auto classic=base;classic.pseudo_transient=false;classic.pressure_correction=PressureCorrection::Simple;
    classic.momentum_relaxation=.7;classic.pressure_relaxation=.3;
    const auto simple=channel(*dev,classic,1200);
    CHECK(simple.report.continuity.raw<1e-6);
    CHECK(simple.report.pseudo_time_step==0);
    CHECK(max_difference(simple,reference)<1e-5);
    // Standard pressure (TG 23-42) and distance-weighted Rhie-Chow (TG 23-48)
    // are different discretizations: they must converge and move the result.
    for(const bool standard:{true,false}) {
      INFO("standard pressure: "<<standard);
      auto fluent=classic;
      if(standard) fluent.pressure_interpolation=PressureInterpolation::Standard;
      else fluent.face_velocity=FaceVelocityWeights::Distance;
      const auto alternative=channel(*dev,fluent,1200);
      CHECK(alternative.report.continuity.raw<1e-6);
      CHECK(max_difference(alternative,reference)>1e-4);
    }
  }
}
