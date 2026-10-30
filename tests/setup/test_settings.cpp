#include "kutta/setup/settings.hpp"
#include "kutta/core/json.hpp"
#include "kutta/io/mesh_io.hpp"
#include "../fluent_examples.hpp"
#include "../io/test_support.hpp"
#include <doctest/doctest.h>
#include <algorithm>
#include <map>
#include <set>

TEST_CASE("Mixing elbow inventories every retained case and zone setting with Fluent groups") {
#if KUTTA_HAS_HDF5
  const auto path=fluent_examples()/"mixing_elbow/mixing_elbow.cas.h5";
  if(!std::filesystem::exists(path)){MESSAGE("SKIP: optional Fluent mixing elbow absent");return;}
  auto c=kutta::fluent::read_case(path);kutta::setup::RunControls run;
  const auto settings=kutta::setup::settings_inventory(c,run);
  std::set<std::string> paths,groups;
  for(const auto& s:settings) {
    CHECK(paths.insert(s.path).second);groups.insert(s.group);
    CHECK(!s.label.empty());CHECK(s.path.find(' ')==std::string::npos);
    CHECK(std::none_of(s.path.begin(),s.path.end(),[](char c){return c>='A' && c<='Z';}));
    CHECK(s.editable==bool(s.apply));if(!s.editable)CHECK(!s.reason.empty());
  }
  CHECK(settings.size()>=c.variables.alist.elements.size());
  for(auto g:{"General","Models","Materials","Cell Zone Conditions","Boundary Conditions","Reference Values","Methods","Controls","Monitors","Initialization","Run Calculation"})CHECK(groups.contains(g));
  for(const auto& z:c.boundaries) {
    // Fixture names are lowercase and path-safe; canonical identities are
    // asserted from the actual case rather than a guessed numeric zone ID.
    CHECK(paths.contains("boundary/"+z.name+"/type"));
    for(const auto& entry:z.raw.alist.elements)if(!entry.elements.empty() && entry.elements[0].atom_text()) {
      const std::string key(*entry.elements[0].atom_text());
      if(key=="vmag")CHECK(paths.contains("boundary/"+z.name+"/velocity-magnitude"));
    }
  }
  MeshTemp t;
  mesh_text(t.path("bad.json"),"{\"v\":2,\"overrides\":{\"missing/path\":1}}");
  CHECK_THROWS_AS(kutta::setup::apply_setup(t.path("bad.json"),settings),kutta::setup::RejectedSettings);
#else
  MESSAGE("SKIP: HDF5 disabled");
#endif
}
TEST_CASE("Mixing elbow settings carry Fluent labels, activity and a guarded turbulence-method switch") {
#if KUTTA_HAS_HDF5
  using kutta::Json;
  const auto path=fluent_examples()/"mixing_elbow/mixing_elbow.cas.h5";
  if(!std::filesystem::exists(path)){MESSAGE("SKIP: optional Fluent mixing elbow absent");return;}
  auto c=kutta::fluent::read_case(path);kutta::setup::RunControls run;
  const auto settings=kutta::setup::settings_inventory(c,run);
  std::map<std::string,const kutta::setup::Setting*> by;for(const auto& s:settings)by[s.path]=&s;
  for(const auto& s:settings) {
    CAPTURE(s.path);
    if(!s.option_labels.empty())CHECK(s.option_labels.size()==s.options.size());
    if(s.depends_on) {
      // Activity is exactly the controlling value's membership in active_when.
      REQUIRE(by.contains(*s.depends_on));const auto& v=by[*s.depends_on]->value;
      CHECK(s.active==(std::find(s.active_when.begin(),s.active_when.end(),v)!=s.active_when.end()));
    } else CHECK(s.active_when.empty());
  }
  const auto& spec=*by.at("boundary/cold-inlet/turbulence-specification");
  CHECK(spec.value==Json(3));CHECK(spec.value_label=="Intensity and Hydraulic Diameter");CHECK(spec.editable);CHECK(spec.active);
  CHECK(spec.options==Json::Array{0,1,2,3});CHECK(spec.option_labels.size()==spec.options.size());CHECK(spec.option_labels[0]==Json("K and Omega"));
  for(auto [key,active]:std::vector<std::pair<std::string,bool>>{{"turb-length-scale",false},{"k",false},{"omega",false},{"turb-viscosity-ratio",false},{"turb-hydraulic-diam",true},{"turb-intensity",true}}) {
    const auto& s=*by.at("boundary/cold-inlet/"+key);CAPTURE(key);
    CHECK(s.active==active);CHECK(s.depends_on==std::optional<std::string>("boundary/cold-inlet/turbulence-specification"));CHECK(s.editable);
  }
  CHECK(by.at("boundary/cold-inlet/turb-length-scale")->active_when==Json::Array{1});
  CHECK(by.at("boundary/cold-inlet/turb-intensity")->active_when==Json::Array{1,2,3});
  CHECK(by.at("boundary/cold-inlet/k")->active_when==Json::Array{0});CHECK(by.at("boundary/cold-inlet/omega")->active_when==Json::Array{0});
  CHECK(by.at("boundary/cold-inlet/epsilon")->active_when.empty());CHECK_FALSE(by.at("boundary/cold-inlet/epsilon")->editable);
  CHECK(by.at("boundary/cold-inlet/turb-hydraulic-diam")->label=="Hydraulic diameter");
  const auto& velocity=*by.at("boundary/cold-inlet/velocity-specification");
  CHECK(velocity.value==Json(2));CHECK(velocity.value_label=="Magnitude, Normal to Boundary");CHECK_FALSE(velocity.editable);CHECK(velocity.options.size()==3);
  CHECK(by.at("boundary/cold-inlet/velocity-magnitude")->active);CHECK_FALSE(by.at("boundary/cold-inlet/direction")->active);
  CHECK_FALSE(by.at("boundary/wall-elbow/turbulence-specification")->active);CHECK_FALSE(by.at("boundary/wall-elbow/gauge-pressure")->active);
  CHECK_FALSE(by.at("boundary/wall-elbow/velocity-magnitude")->depends_on);
  CHECK(by.at("boundary/outlet/gauge-pressure")->active);CHECK_FALSE(by.contains("boundary/outlet/velocity-specification"));
  CHECK(by.at("models/viscous")->value_label=="k-omega SST");CHECK(by.at("boundary/wall-elbow/type")->value_label=="wall");
  CHECK(by.at("methods/momentum-scheme")->value_label=="Second Order Upwind");
  CHECK(by.at("initialization/method")->value_label=="Hybrid Initialization");
  const auto json=kutta::dump_json(kutta::setup::settings_json(settings));
  for(auto key:{"\"option_labels\":","\"value_label\":","\"active\":","\"depends_on\":","\"active_when\":"})CHECK(json.find(key)!=std::string::npos);

  // A switch never reaches a spec whose values were not saved.
  MeshTemp t;
  auto missing=c;auto& absent=*std::find_if(missing.boundaries.begin(),missing.boundaries.end(),[](const auto& z){return z.name=="cold-inlet";});absent.turbulence.erase("turb-length-scale");
  const auto guarded=kutta::setup::settings_inventory(missing,run);
  const auto& restricted=*std::find_if(guarded.begin(),guarded.end(),[](const auto& s){return s.path=="boundary/cold-inlet/turbulence-specification";});
  CHECK(restricted.options==Json::Array{0,2,3});CHECK(restricted.option_labels.size()==3);
  mesh_text(t.path("absent.json"),"{\"v\":1,\"overrides\":{\"boundary/cold-inlet/turbulence-specification\":1}}");
  CHECK_THROWS_AS(kutta::setup::apply_setup(t.path("absent.json"),guarded),kutta::setup::RejectedSettings);
  CHECK(absent.turbulence_spec==kutta::fluent::TurbulenceSpec::intensity_diameter);
  // An inactive-but-editable value may be set alongside the switch; the
  // solver setup then reads the new method and value.
  mesh_text(t.path("switch.json"),"{\"v\":1,\"overrides\":{\"boundary/cold-inlet/turbulence-specification\":2,\"boundary/cold-inlet/turb-viscosity-ratio\":7}}");
  CHECK(kutta::setup::apply_setup(t.path("switch.json"),settings).size()==2);
  const auto& zone=*std::find_if(c.boundaries.begin(),c.boundaries.end(),[](const auto& z){return z.name=="cold-inlet";});
  CHECK(zone.turbulence_spec==kutta::fluent::TurbulenceSpec::intensity_viscosity_ratio);
  const kutta::setup::Approximation wfgc[]={kutta::setup::Approximation::wfgc};
  const auto setup=kutta::setup::fluent_setup(c,kutta::read_mesh(path),wfgc);
  const auto bc=std::find_if(setup.boundaries.begin(),setup.boundaries.end(),[](const auto& b){return b.patch=="cold-inlet";});REQUIRE(bc!=setup.boundaries.end());
  const auto* inlet=std::get_if<kutta::fv::VelocityInlet>(&bc->condition);REQUIRE(inlet);REQUIRE(inlet->turbulence);
  const auto* ratio=std::get_if<kutta::fv::IntensityViscosityRatio>(&*inlet->turbulence);REQUIRE(ratio);
  CHECK(ratio->intensity==doctest::Approx(.05));CHECK(ratio->ratio==doctest::Approx(7));
#else
  MESSAGE("SKIP: HDF5 disabled");
#endif
}
