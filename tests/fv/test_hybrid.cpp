#include "../backend_helpers.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/setup/fluent.hpp"

using namespace kutta;
using namespace kutta::fv;

TEST_CASE("Hybrid potential obeys inlet flux and outlet potential, with constant and local turbulence recipes") {
  const auto mesh=make_box_hex(6,3,2,{0,0,0},{1,.3,.2});
  std::vector<PatchCondition> bc;
  for(const auto& p:mesh.patches) {
    if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{2,0,0},IntensityViscosityRatio{.1,5}}});
    else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{3,NormalBackflow{},BackflowPressure::Static,TurbulenceValues{.06,12}}});
    else bc.push_back({p.name,Symmetry{}});
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    SolverSettings settings;settings.viscosity=.001;settings.sst=SstSettings{};
    SteadySolver solver(*dev,mesh,bc,settings);
    HybridInitialization h;h.iterations=3;solver.initialize_hybrid(h);
    for(auto x:dev->download(solver.velocity().x)) CHECK(x==doctest::Approx(2).epsilon(2e-4));
    for(auto y:dev->download(solver.velocity().y)) CHECK(std::abs(y)<1e-4);
    for(auto p:dev->download(solver.pressure())) CHECK(p==3);
    for(auto k:dev->download(solver.kinetic_energy())) CHECK(k==doctest::Approx(.06).epsilon(3e-4));
    for(auto w:dev->download(solver.specific_dissipation())) CHECK(w==doctest::Approx(12).epsilon(3e-4));
    CHECK(solver.monitor().continuity.raw<2e-6);
    h.averaged_turbulence=false;h.constant_speed=true;h.linear_potential=true;
    h.inlet_pressure={{"xmin",5}};solver.initialize_hybrid(h);
    const auto p=dev->download(solver.pressure());
    for(Index c=0;c<mesh.n_cells;++c)
      CHECK(p[c]==doctest::Approx(5-2*solver.geometry().cell_centroid[c].x).epsilon(2e-4));
  }
}

TEST_CASE("Fluent hybrid controls are parsed strictly and do not require stale standard defaults") {
  fluent::FluentCase c;
  c.variables.alist=fluent::parse_sexpr(R"(((hyb-init? #t) (fmg-init/enable? #f)
    (hyb-init/niter 7) (hyb-init/explicit-relax (.6 .8)) (hyb-init/relax (.9 1))
    (hyb-init/residual-tolerance (1e-5 2e-5)) (hyb-init/solve? (#t #f)) (hyb-init/patch (2 3))
    (hyb-init/turb-domain-averaged? #f) (hyb-init/patch-ads? #t) (hyb-init/constant-vel-vector? #t)
    (hyb-init/turb-intens .03) (hyb-init/turb-vr 4) (hyb-init/pressure-from-bc? #f)))");
  const auto mesh=make_box_hex(1,1,1,{0,0,0},{1,1,1});
  const auto init=setup::initial_fields(c,mesh);
  REQUIRE(std::holds_alternative<HybridInitialization>(init));
  const auto& h=std::get<HybridInitialization>(init);
  CHECK(h.iterations==7);CHECK(h.explicit_relaxation[0]==.6);CHECK(h.implicit_relaxation[0]==.9);
  CHECK(h.tolerance[1]==2e-5);CHECK_FALSE(h.solve[1]);CHECK(h.initial_scalar[0]==2);
  CHECK_FALSE(h.averaged_turbulence);CHECK(h.linear_potential);CHECK(h.constant_speed);
  CHECK(h.turbulence_intensity==.03);CHECK(h.turbulence_viscosity_ratio==4);
  for(auto& e:c.variables.alist.elements) if(e.elements.front().atom_text()=="hyb-init/niter")
    e.elements[1]=fluent::parse_sexpr("1.5");
  CHECK_THROWS_WITH(setup::initial_fields(c,mesh),doctest::Contains("hyb-init/niter"));
}
