#include <cmath>
#include "kernels.hpp"
#include "../backend_helpers.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"

using namespace kutta;
using namespace kutta::fv;

TEST_CASE("FV SST coupled transport and wall force agree with CPU oracle") {
  auto mesh=make_box_hex(4,3,2,{0,0,0},{1,.4,.3});
  std::vector<PatchCondition> bc;
  for(const auto& p:mesh.patches) {
    if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0},TurbulenceValues{.01,4}}});
    else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0,NormalBackflow{},BackflowPressure::Total,IntensityViscosityRatio{.05,4}}});
    else if(p.name=="zmin") bc.push_back({p.name,Symmetry{}});
    else bc.push_back({p.name,Wall{}});
  }
  std::vector<float> oracle;
  KUTTA_FOR_EACH_BACKEND(dev) {
    SolverSettings settings;settings.viscosity=.001;settings.sst=SstSettings{};settings.monitor_interval=10;
    SteadySolver solver(*dev,mesh,bc,settings);
    std::vector<Vec3d> u(mesh.n_cells,Vec3d{.8,.01,.02});
    std::vector<double> p(mesh.n_cells,.01),k(mesh.n_cells,.01),w(mesh.n_cells,4);
    auto overflow=k;overflow[0]=1e100;
    CHECK_THROWS_WITH(solver.initialize_fields(u,p,overflow,w),doctest::Contains("exceeds float storage"));
    solver.initialize_fields(u,p,k,w);
    const auto initial_flux=dev->download(solver.mass_flux());
    auto changed=settings;changed.momentum_relaxation=.9;changed.pressure_relaxation=.2;changed.pseudo_time_scale=13;
    SteadySolver different_path(*dev,mesh,bc,changed);different_path.initialize_fields(u,p,k,w);
    CHECK(dev->download(different_path.mass_flux())==initial_flux);
    // Monitoring must not impose fixed omega rows on the caller's saved state.
    const auto before=dev->download(solver.specific_dissipation());solver.monitor();
    CHECK(dev->download(solver.specific_dissipation())==before);
    solver.iterate(3);const auto r=solver.monitor();
    CHECK(std::isfinite(r.k.scaled));CHECK(std::isfinite(r.omega.scaled));
    CHECK(r.k.raw>0);CHECK(r.omega.raw>0);
    const std::vector<std::string> walls{"ymin","ymax","zmax"};
    const auto f=solver.forces(walls);
    CHECK(f.viscous.x>0);
    std::vector<float> result;
    for(const auto* q:{&solver.velocity().x,&solver.velocity().y,&solver.velocity().z,
                      &solver.pressure(),&solver.kinetic_energy(),&solver.specific_dissipation()}) {
      auto v=dev->download(*q);result.insert(result.end(),v.begin(),v.end());
      for(auto x:v) REQUIRE(std::isfinite(x));
    }
    if(oracle.empty()) oracle=result;
    else for(std::size_t i=0;i<result.size();++i)
      CHECK(result[i]==doctest::Approx(oracle[i]).epsilon(3e-4).scale(.01));
  }
}

TEST_CASE("FV total-pressure outlet backflow uses adjacent-cell dynamic pressure") {
  auto mesh=make_box_hex(3,2,2,{0,0,0},{1,1,1});std::vector<PatchCondition> bc;
  for(const auto& p:mesh.patches) {
    if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{-2,0,0}}});
    else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{7,NormalBackflow{},BackflowPressure::Total}});
    else bc.push_back({p.name,Symmetry{}});
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    SolverSettings settings;settings.density=3;
    SteadySolver solver(*dev,mesh,bc,settings);
    solver.initialize([](Vec3d){return Vec3d{-2,0,0};},[](Vec3d){return 1.;});
    solver.monitor();const auto pressure=dev->download(solver.boundary_pressure());
    for(const auto& patch:solver.mesh().patches) if(patch.name=="xmax")
      for(Index f=patch.start;f<patch.start+patch.size;++f)
        CHECK(pressure[f-solver.mesh().n_internal_faces()]==doctest::Approx(1));
  }
}

TEST_CASE("FV transport integrates implicit sinks and component wall terms before relaxation") {
  const auto host=make_box_hex(1,1,1,{0,0,0},{1,2,3});const auto geometry=compute_geometry(host);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,host,geometry);auto s=make_stencil(*dev,host,geometry);
    TransportEquation eq(*dev,m,s,host,geometry);
    auto field=dev->alloc<float>(1),source=dev->alloc<float>(1),sink=dev->alloc<float>(1);
    auto extra=dev->alloc<float>(1),rhs=dev->alloc<float>(1);
    dev->fill(field,2);dev->fill(source,7);dev->fill(sink,3);dev->fill(extra,2);dev->fill(rhs,5);
    auto boundary=dev->alloc<float>(6),mass=dev->alloc<float>(6),mu=dev->alloc<float>(6);
    auto kind=dev->alloc<Index>(6);dev->fill(boundary,0);dev->fill(mass,0);dev->fill(mu,0);dev->fill(kind,1);
    TransportSettings settings;settings.relaxation=.5;
    eq.assemble(field,boundary,kind,mass,mu,source,settings,{&sink,&extra,&rhs});
    CHECK(dev->download(eq.diagonal)[0]==doctest::Approx(20));
    // (physical A=3*6+2), b=7*6+5; URF contributes 20*old.
    eq.solve(field,{1e-7,1e-12,10});CHECK(dev->download(field)[0]==doctest::Approx(87./40));
    auto mask=dev->alloc<Index>(1);auto fixed=dev->alloc<float>(1);dev->fill(mask,1);dev->fill(fixed,4.2f);
    kern::fv_fixed_rows(*dev,1,eq.matrix.row_ptr,eq.matrix.diag,mask,fixed,eq.matrix.val,eq.rhs);
    CHECK(dev->download(eq.matrix.val)[0]==doctest::Approx(40));
    CHECK(dev->download(field)[0]==doctest::Approx(87./40));
    eq.solve(field,{1e-7,1e-12,10});CHECK(dev->download(field)[0]==doctest::Approx(4.2));
  }
}
