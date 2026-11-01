#include <iostream>
#include "../backend_helpers.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"

using namespace kutta;
using namespace kutta::fv;

TEST_CASE("FV symmetry constrains normal velocity in a closed box") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=make_box_hex(1,1,1,{0,0,0},{1,1,1});
    std::vector<PatchCondition> bc;
    for(const auto& patch:mesh.patches) bc.push_back({patch.name,Symmetry{}});
    SolverSettings settings;settings.viscosity=.1;
    SteadySolver solver(*dev,std::move(mesh),bc,settings);
    solver.initialize([](Vec3d) {return Vec3d{1,.4,-.2};},[](Vec3d) {return 0.;});
    solver.iterate(100);
    CHECK(solver.report().maximum_velocity<1e-5);
  }
}

TEST_CASE("FV SIMPLEC preserves rest and uniform inlet outlet flow") {
  for(int through:{0,1,-1}) {
    auto mesh=make_box_hex(8,4,1,{0,0,0},{1,1,.1});
    std::vector<PatchCondition> conditions;
    for(const auto& p:mesh.patches) {
      if(through && p.name=="xmin") conditions.push_back({p.name,VelocityInlet{Vec3d{double(through),0,0}}});
      else if(through && p.name=="xmax") conditions.push_back({p.name,PressureOutlet{0,Vec3d{double(through),0,0}}});
      else if(through || p.name=="zmin" || p.name=="zmax") conditions.push_back({p.name,Symmetry{}});
      else conditions.push_back({p.name,Wall{}});
    }
    std::vector<float> baseline;
    KUTTA_FOR_EACH_BACKEND(dev) {
      SteadySolver solver(*dev,mesh,conditions);
      solver.initialize([through](Vec3d){return Vec3d{double(through),0,0};},[](Vec3d){return 0.;});
      solver.iterate(5);auto report=solver.monitor();
      CHECK(report.residual_iteration==5);
      solver.iterate(1);
      CHECK(solver.report().iteration==6);CHECK(solver.report().residual_iteration==5);
      CHECK(solver.monitor().residual_iteration==6);
      auto u=dev->download(solver.velocity().x),v=dev->download(solver.velocity().y),p=dev->download(solver.pressure());
      for(auto x:u) CHECK(x==doctest::Approx(through).epsilon(1e-5));
      for(auto x:v) CHECK(std::abs(x)<1e-6);
      for(auto x:p) CHECK(std::abs(x)<1e-6);
      CHECK(report.continuity.raw<1e-6);
      if(baseline.empty()) baseline=u;else CHECK(u==baseline);
    }
  }
}

TEST_CASE("FV SIMPLEC develops planar Couette flow with moving and slip walls" * doctest::test_suite("slow")) {
  auto mesh=make_box_hex(4,12,1,{0,0,0},{1,1,.1});
  std::vector<PatchCondition> boundaries;
  for(auto& p:mesh.patches) {
    if(p.name=="ymax") boundaries.push_back({p.name,Wall{{1,4,0}}});
    else if(p.name=="ymin") boundaries.push_back({p.name,Wall{}});
    else if(p.name=="xmin" || p.name=="xmax") boundaries.push_back({p.name,VelocityInlet{
      std::function<Vec3d(Vec3d)>([](Vec3d p){return Vec3d{p.y,0,0};})}});
    else boundaries.push_back({p.name,Symmetry{}});
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    SolverSettings settings;settings.viscosity=1;
    SteadySolver solver(*dev,mesh,boundaries,settings);
    for(int k=0;k<20;++k) solver.iterate(20);
    auto u=dev->download(solver.velocity().x),v=dev->download(solver.velocity().y);
    double error=0;
    for(Index c=0;c<mesh.n_cells;++c) error=std::max(error,std::abs(u[c]-solver.geometry().cell_centroid[c].y));
    std::cout<<"FV_COUETTE backend="<<backend_name(dev->backend())<<" max_error="<<error
      <<" residual="<<solver.report().u.scaled<<std::endl;
    CHECK(error<2e-4);
    CHECK(dev->max_abs(solver.velocity().y,v.size())<2e-4);
    const std::vector<std::string> top{"ymax"};
    auto force=solver.forces(top);
    CHECK(force.viscous.x==doctest::Approx(-.1).epsilon(2e-4));
    CHECK(std::abs(force.pressure.x)<1e-6);
    CHECK(force.coefficient({1,0,0},{1,1,.1,0})==doctest::Approx(-2).epsilon(2e-4));
    auto shifted=solver.forces(top,2);
    CHECK(shifted.pressure.y-force.pressure.y==doctest::Approx(-.2).epsilon(1e-6));
    CHECK(std::abs(shifted.coefficient({0,1,0},{1,1,.1,0})-
                   force.coefficient({0,1,0},{1,1,.1,0}))<1e-6);
  }
}

TEST_CASE("FV pressure outlet defaults to flux-consistent normal backflow") {
  auto mesh=make_box_hex(4,3,1,{0,0,0},{1,1,.1});
  std::vector<PatchCondition> conditions;
  for(const auto& p:mesh.patches) {
    if(p.name=="xmin") conditions.push_back({p.name,VelocityInlet{Vec3d{-1,0,0}}});
    else if(p.name=="xmax") conditions.push_back({p.name,PressureOutlet{}});
    else conditions.push_back({p.name,Symmetry{}});
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    SolverSettings settings;settings.density=2;
    SteadySolver solver(*dev,mesh,conditions,settings);
    solver.initialize([](Vec3d){return Vec3d{-1,0,0};},[](Vec3d){return 0.;});
    solver.iterate(10);
    for(auto u:dev->download(solver.velocity().x)) CHECK(u==doctest::Approx(-1).epsilon(1e-6));
    CHECK(dev->max_abs(solver.velocity().y,mesh.n_cells)<1e-6);
    CHECK(dev->max_abs(solver.pressure(),mesh.n_cells)<1e-6);
  }
}

TEST_CASE("FV closed-domain gauge pressure is zero at the reference point's cell") {
  auto mesh=make_box_hex(6,6,1,{0,0,0},{1,1,.1});
  std::vector<PatchCondition> conditions;
  for(const auto& p:mesh.patches) {
    if(p.name=="ymax") conditions.push_back({p.name,Wall{{1,0,0}}});
    else if(p.name=="zmin" || p.name=="zmax") conditions.push_back({p.name,Symmetry{}});
    else conditions.push_back({p.name,Wall{}});
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    SolverSettings settings;settings.reference_point={.75,.25,.05};  // cell centre (5/12, 1/12 grid)
    SteadySolver solver(*dev,mesh,conditions,settings);
    solver.initialize([](Vec3d){return Vec3d{};},[](Vec3d x){return x.x+2*x.y;});
    solver.iterate(10);
    const auto p=dev->download(solver.pressure());
    const auto& centres=solver.geometry().cell_centroid;
    Index pinned=-1;
    for(Index c=0;c<mesh.n_cells;++c)
      if(std::abs(centres[c].x-.75)<1e-9 && std::abs(centres[c].y-.25)<1e-9) pinned=c;
    REQUIRE(pinned>=0);
    CHECK(p[pinned]==0.f);
    CHECK(dev->max_abs(solver.pressure(),mesh.n_cells)>1e-3);
  }
}

// The float32 state's continuity plateaus at its velocity-rounding floor
// (validation/results/residual_floor). The float-pair state continues to a
// floor about 100x lower on this channel (validation/results/paired_state),
// with the same fields to float precision.
TEST_CASE("FV SIMPLEC float-pair state falls below the float32 continuity floor") {
  auto mesh=make_box_hex(24,8,1,{0,0,0},{3,1,.1});
  std::vector<PatchCondition> bc;
  for(const auto& p:mesh.patches) {
    if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0}}});
    else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0,Vec3d{1,0,0}}});
    else if(p.name=="zmin" || p.name=="zmax") bc.push_back({p.name,Symmetry{}});
    else bc.push_back({p.name,Wall{}});
  }
  std::vector<float> baseline;
  KUTTA_FOR_EACH_BACKEND(dev) {
    std::vector<float> u[2];double continuity[2];
    for(bool paired:{false,true}) {
      SolverSettings settings;settings.viscosity=.02;settings.paired_state=paired;
      SteadySolver solver(*dev,mesh,bc,settings);
      solver.initialize([](Vec3d){return Vec3d{1,0,0};},[](Vec3d){return 0.;});
      solver.iterate(200);
      continuity[paired]=solver.monitor().continuity.scaled;u[paired]=dev->download(solver.velocity().x);
    }
    CHECK(continuity[0]>1e-6);CHECK(continuity[1]<1e-7);
    for(std::size_t i=0;i<u[0].size();++i) CHECK(u[1][i]==doctest::Approx(u[0][i]).epsilon(1e-5));
    if(baseline.empty()) baseline=u[1];else CHECK(u[1]==baseline);
  }
}

// Coupled carries the same pairs: its block RHS is the pair's residual.
TEST_CASE("FV Coupled float-pair state falls below the float32 continuity floor") {
  auto mesh=make_box_hex(24,8,1,{0,0,0},{3,1,.1});
  std::vector<PatchCondition> bc;
  for(const auto& p:mesh.patches) {
    if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0}}});
    else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0,Vec3d{1,0,0}}});
    else if(p.name=="zmin" || p.name=="zmax") bc.push_back({p.name,Symmetry{}});
    else bc.push_back({p.name,Wall{}});
  }
  std::vector<float> baseline;
  KUTTA_FOR_EACH_BACKEND(dev) {
    std::vector<float> u[2];double continuity[2];
    for(bool paired:{false,true}) {
      SolverSettings settings;settings.viscosity=.02;settings.paired_state=paired;
      settings.coupling=PressureVelocityCoupling::Coupled;settings.momentum_relaxation=settings.pressure_relaxation=.5;
      SteadySolver solver(*dev,mesh,bc,settings);
      solver.initialize([](Vec3d){return Vec3d{1,0,0};},[](Vec3d){return 0.;});
      solver.iterate(1400);
      continuity[paired]=solver.monitor().continuity.scaled;u[paired]=dev->download(solver.velocity().x);
    }
    // Float32 stalls near 1e-5 from iteration 800; the pair reaches 1.9e-8.
    CHECK(continuity[0]>5e-6);CHECK(continuity[1]<1e-7);
    for(std::size_t i=0;i<u[0].size();++i) CHECK(u[1][i]==doctest::Approx(u[0][i]).epsilon(1e-5));
    if(baseline.empty()) baseline=u[1];else CHECK(u[1]==baseline);
  }
}
