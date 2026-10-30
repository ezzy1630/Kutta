#include "../backend_helpers.hpp"
#include <algorithm>
#include "../fluent_examples.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/setup/fluent.hpp"

using namespace kutta;
namespace {
void set(fluent::SExpr& alist,const std::string& settings) {
  for(auto replacement:fluent::parse_sexprs(settings)) {
    const auto key=replacement.elements.front().atom_text();
    auto it=std::find_if(alist.elements.begin(),alist.elements.end(),[&](const auto& e){return !e.elements.empty() && e.elements.front().atom_text()==key;});
    if(it==alist.elements.end()) alist.elements.push_back(std::move(replacement));else *it=std::move(replacement);
  }
}
fluent::FluentCase change(const fluent::FluentCase& original,std::int64_t id,const std::string& type,
    const std::string& settings,const std::string& variables="") {
  auto sections=original.raw_sections;
  for(auto& thread:sections.at("Thread Variables")) if(thread.elements[1].elements[0].integer==id) {
    thread.elements[1].elements[1]=fluent::SExpr::symbol(type);set(thread.elements[2],settings);
  }
  set(sections.at("Rampant Variables")[0].elements[1],variables);
  auto c=fluent::extract_case(std::move(sections),original.solver_version,original.origin);c.mesh=original.mesh;return c;
}
}
TEST_CASE("FSAE strict Fluent adapter accepts verified settings and holds unmet physics gates") {
#if KUTTA_HAS_HDF5
  const auto path=fluent_examples()/"mixing_elbow/mixing_elbow.cas.h5";
  if(!std::filesystem::exists(path)) {MESSAGE("SKIP external elbow fixture absent; parser conventions have independent synthetic coverage");return;}
  const auto c=fluent::read_case(path);const auto mesh=read_fluent_cff(path);const auto allowed=setup::parse_approximations("wfgc");
  const auto wall=std::find_if(c.boundaries.begin(),c.boundaries.end(),[](const auto& z){return z.type==fluent::ZoneType::wall;});REQUIRE(wall!=c.boundaries.end());
  const auto inlet=std::find_if(c.boundaries.begin(),c.boundaries.end(),[](const auto& z){return z.type==fluent::ZoneType::velocity_inlet;});REQUIRE(inlet!=c.boundaries.end());
  auto moving=change(c,wall->id,"wall","(motion-bc . 1) (rotating? . #t) (relative? . #f) (omega . 2) (ai . 0) (aj . 0) (ak . 2) (x-origin . 0) (y-origin . 0) (z-origin . 0)");
  const auto rotating=setup::fluent_setup(moving,mesh,allowed);
  const auto adapted=std::find_if(rotating.boundaries.begin(),rotating.boundaries.end(),[&](const auto& p){return p.patch==wall->name;});
  REQUIRE(adapted!=rotating.boundaries.end());
  const auto& native_wall=std::get<fv::Wall>(adapted->condition);
  REQUIRE(native_wall.rotation.has_value());
  CHECK(native_wall.rotation->angular_velocity.z==2);
  moving=change(c,wall->id,"wall","(motion-bc . 1) (rotating? . #f) (components? . #t) (relative? . #f) (u . 12) (v . 0) (w . 0)");
  CHECK_NOTHROW(setup::fluent_setup(moving,mesh,allowed));
  auto pressure=change(c,inlet->id,"pressure-inlet","(p0 . 1) (direction-spec . 1) (flow-spec . 0)");
  const auto pressure_setup=setup::fluent_setup(pressure,mesh,allowed);
  CHECK(std::any_of(pressure_setup.boundaries.begin(),pressure_setup.boundaries.end(),[](const auto& p){return std::holds_alternative<fv::PressureInlet>(p.condition);}));
  auto mass=change(c,inlet->id,"mass-flow-inlet","(mass-flow . .2) (direction-spec . 1) (flow-spec . 0)");
  CHECK_NOTHROW(setup::fluent_setup(mass,mesh,allowed));
  mass=change(c,inlet->id,"mass-flow-inlet","(mass-flux . 2) (direction-spec . 1) (flow-spec . 1)");
  CHECK_NOTHROW(setup::fluent_setup(mass,mesh,allowed));
  const auto fluid=c.cell_zones.front().id;
  auto porous=change(c,fluid,"fluid",R"((porous? . #t) (rel-vel-resistance? . #t)
    (porous-r/1 . 30) (porous-r/2 . 200) (porous-r/3 . 400) (porous-c/1 . 4) (porous-c/2 . 20) (porous-c/3 . 40)
    (direction-1/x . 1) (direction-1/y . 0) (direction-1/z . 0) (direction-2/x . 0) (direction-2/y . 1) (direction-2/z . 0)
    (porosity . .35) (alt-inertial-form? . #f) (c0 . 0) (c1 . 0))","(porous/physical-velocity-form? . #f)");
  const auto porous_setup=setup::fluent_setup(porous,mesh,allowed);
  REQUIRE(porous_setup.settings.cell_zone_physics.size()==1);CHECK(porous_setup.settings.cell_zone_physics[0].porous->porosity==.35);
  auto frame=change(c,fluid,"fluid","(mrf-motion? . #t) (mrf-omega . .4)","(solve-absolute-velocities? . #t)");
  CHECK_NOTHROW(setup::fluent_setup(frame,mesh,allowed));
  auto fan=change(c,wall->id,"fan","(pressure-jump (polynomial velocity 5 -1 -1))");
  CHECK_THROWS_WITH(setup::fluent_setup(fan,mesh,allowed),doctest::Contains("native Fluent fan polynomial input conventions"));
  auto farfield=change(c,inlet->id,"pressure-far-field","");
  CHECK_THROWS_WITH(setup::fluent_setup(farfield,mesh,allowed),doctest::Contains("incompressible constant-density"));
#endif
}
