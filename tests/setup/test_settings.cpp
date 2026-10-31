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
TEST_CASE("Mixing elbow convergence controls read Fluent residual criteria and stay run controls") {
#if KUTTA_HAS_HDF5
  using kutta::Json;
  const auto path=fluent_examples()/"mixing_elbow/mixing_elbow.cas.h5";
  if(!std::filesystem::exists(path)){MESSAGE("SKIP: optional Fluent mixing elbow absent");return;}
  auto c=kutta::fluent::read_case(path);kutta::setup::RunControls run;
  auto settings=kutta::setup::settings_inventory(c,run);
  std::map<std::string,const kutta::setup::Setting*> by;for(const auto& s:settings)by[s.path]=&s;
  CHECK_FALSE(by.contains("monitors/convergence-criteria"));
  const auto& check=*by.at("monitors/check-convergence");
  CHECK(check.value==Json(false));CHECK(check.kind=="bool");CHECK(check.editable);CHECK(check.group=="Monitors");
  const auto& condition=*by.at("monitors/condition");
  CHECK(condition.options==Json::Array{"all","any"});CHECK(condition.option_labels==Json::Array{"All Conditions are Met","Any Condition is Met"});
  CHECK(condition.value_label=="All Conditions are Met");CHECK(condition.editable);
  // The case saves every residual checked at 1e-3 (energy at 1e-6, which
  // Kutta does not solve, so the raw row is retained read-only).
  for(auto eq:kutta::setup::residual_equations) {
    const std::string prefix="monitors/residual/"+std::string(eq)+"/";CAPTURE(prefix);
    CHECK(by.at(prefix+"check")->value==Json(true));CHECK(by.at(prefix+"check")->editable);
    const auto& criterion=*by.at(prefix+"absolute-criterion");
    CHECK(criterion.value==Json(1e-3));CHECK(criterion.editable);CHECK(criterion.min==0.0);CHECK(criterion.label=="Absolute Criteria ("+std::string(eq)+")");
  }
  CHECK(run.convergence.residuals.size()==6);CHECK_FALSE(run.convergence.check);
  CHECK(std::all_of(run.convergence.residuals.begin(),run.convergence.residuals.end(),[](const auto& r){return r.check && r.absolute==1e-3;}));
  CHECK(by.contains("case/residuals/settings"));CHECK_FALSE(by.at("case/residuals/settings")->editable);
  CHECK(run.convergence.reports.empty()); // the elbow has no report definitions
  CHECK(by.at("reference/length")->editable);CHECK(by.at("reference/length")->value==Json(1));
  for(const auto& s:settings)if(s.path.starts_with("monitors/"))CHECK(kutta::setup::run_control(s.path));
  CHECK(kutta::setup::run_control("run/iterations"));CHECK_FALSE(kutta::setup::run_control("reference/length"));

  MeshTemp t;
  mesh_text(t.path("convergence.json"),R"({"v":1,"overrides":{"monitors/check-convergence":true,"monitors/condition":"any",
    "monitors/residual/k/check":false,"monitors/residual/omega/absolute-criterion":0.0001,"reference/length":2.5}})");
  CHECK(kutta::setup::apply_setup(t.path("convergence.json"),settings).size()==5);
  CHECK(run.convergence.check);CHECK(run.convergence.combine==kutta::setup::ConvergenceControls::Combine::any);
  CHECK_FALSE(run.convergence.residuals[4].check);CHECK(run.convergence.residuals[5].absolute==1e-4);CHECK(c.reference.length==2.5);
  for(auto [key,value]:std::vector<std::pair<std::string,std::string>>{{"monitors/condition","\"some\""},{"monitors/residual/k/absolute-criterion","-1"},
      {"monitors/check-convergence","1"},{"monitors/residual/energy/check","true"},{"reference/length","0"}}) {
    CAPTURE(key);mesh_text(t.path("bad.json"),"{\"v\":1,\"overrides\":{\""+key+"\":"+value+"}}");
    try {kutta::setup::apply_setup(t.path("bad.json"),settings);FAIL("accepted");}
    catch(const kutta::setup::RejectedSettings& e){REQUIRE(e.issues().size()==1);CHECK(e.issues()[0].path==key);}
  }

  // Residual rows only for Kutta's equations are fully represented, so the raw
  // value is consumed. Any other row layout leaves checks off and criteria unknown.
  const auto with_residuals=[&](const std::string& rows) {
    auto modified=kutta::fluent::read_case(path);
    for(auto& e:modified.variables.alist.elements)if(!e.elements.empty() && e.elements[0].atom_text()=="residuals/settings")e=kutta::fluent::parse_sexpr("(residuals/settings "+rows+")");
    kutta::setup::RunControls controls;const auto inventory=kutta::setup::settings_inventory(modified,controls);
    const bool raw=std::any_of(inventory.begin(),inventory.end(),[](const auto& s){return s.path=="case/residuals/settings";});
    return std::pair{controls.convergence,raw};
  };
  auto [parsed,raw]=with_residuals("((continuity #t 0.0 #f 0.01 0.0 0.05) (k #t 0.0 #t 0.002 0.0 0.05))");
  CHECK_FALSE(raw);
  CHECK_FALSE(parsed.residuals[0].check);CHECK(parsed.residuals[0].absolute==0.01);
  CHECK(parsed.residuals[4].check);CHECK(parsed.residuals[4].absolute==0.002);
  CHECK_FALSE(parsed.residuals[1].check);CHECK_FALSE(parsed.residuals[1].absolute); // absent row
  std::tie(parsed,raw)=with_residuals("((continuity #t 0 #f 0.001 0.001 1e-10 0.0001 0.0007 0.05))");
  CHECK(raw);CHECK(std::none_of(parsed.residuals.begin(),parsed.residuals.end(),[](const auto& r){return r.check || r.absolute;}));
#else
  MESSAGE("SKIP: HDF5 disabled");
#endif
}
TEST_CASE("Mixing elbow runs Fluent's Coupled scheme natively and can switch to SIMPLEC") {
#if KUTTA_HAS_HDF5
  using kutta::Json;
  const auto path=fluent_examples()/"mixing_elbow/mixing_elbow.cas.h5";
  if(!std::filesystem::exists(path)){MESSAGE("SKIP: optional Fluent mixing elbow absent");return;}
  auto c=kutta::fluent::read_case(path);kutta::setup::RunControls run;
  const auto mesh=kutta::read_mesh(path);
  const kutta::setup::Approximation wfgc[]={kutta::setup::Approximation::wfgc};
  const auto notes=[](const kutta::setup::FluentSetup& s) {
    std::string all;for(const auto& note:s.convergence_notes)all+=note+'\n';return all;
  };
  const auto saved=kutta::setup::fluent_setup(c,mesh,wfgc);
  CHECK(saved.settings.coupling==kutta::fv::PressureVelocityCoupling::Coupled);
  CHECK(notes(saved).find("flow/scheme")==std::string::npos);
  // Fluent's pseudo-time explicit relaxation factors of the Coupled scheme.
  CHECK(saved.settings.momentum_relaxation==.5);CHECK(saved.settings.pressure_relaxation==.5);
  const auto settings=kutta::setup::settings_inventory(c,run);
  const auto it=std::find_if(settings.begin(),settings.end(),[](const auto& s){return s.path=="methods/pressure-velocity-coupling";});
  REQUIRE(it!=settings.end());
  CHECK(it->value==Json("coupled"));CHECK(it->value_label=="Coupled");CHECK(it->editable);
  CHECK(it->options==Json::Array{"simple","simplec","coupled"});
  std::map<std::string,const kutta::setup::Setting*> by;for(const auto& s:settings)by[s.path]=&s;
  CHECK(by.at("methods/pseudo-transient")->value==Json(true));CHECK_FALSE(by.at("methods/pseudo-transient")->editable);
  CHECK(by.at("methods/gradient")->value==Json("least-squares-cell-based"));CHECK(by.at("methods/gradient")->value_label=="Least Squares Cell Based");
  CHECK(by.at("methods/pressure-scheme")->value_label=="Second Order");CHECK_FALSE(by.at("methods/pressure-scheme")->editable);
  // SST solves k and omega; the epsilon scheme is listed but unused.
  CHECK(by.at("methods/omega-scheme")->active);CHECK(by.at("methods/omega-scheme")->editable);
  CHECK_FALSE(by.at("methods/epsilon-scheme")->active);CHECK_FALSE(by.at("methods/epsilon-scheme")->editable);
  MeshTemp t;
  mesh_text(t.path("simplec.json"),"{\"v\":1,\"overrides\":{\"methods/pressure-velocity-coupling\":\"simplec\"}}");
  CHECK(kutta::setup::apply_setup(t.path("simplec.json"),settings).size()==1);
  const auto switched=kutta::setup::fluent_setup(c,mesh,wfgc);
  CHECK(switched.settings.coupling==kutta::fv::PressureVelocityCoupling::Simplec);
  // The saved Coupled pseudo-time controls still apply; only the path changes.
  CHECK(switched.settings.momentum_relaxation==.5);CHECK(switched.settings.pseudo_time_scale==saved.settings.pseudo_time_scale);
  CHECK(notes(switched).find("SIMPLEC selected by setup override for a Coupled case")!=std::string::npos);
  mesh_text(t.path("simple.json"),"{\"v\":1,\"overrides\":{\"methods/pressure-velocity-coupling\":\"simple\"}}");
  CHECK(kutta::setup::apply_setup(t.path("simple.json"),settings).size()==1);
  CHECK(kutta::setup::fluent_setup(c,mesh,wfgc).settings.coupling==kutta::fv::PressureVelocityCoupling::Simple);
#else
  MESSAGE("SKIP: HDF5 disabled");
#endif
}
TEST_CASE("Methods settings follow the case's model, gradient and pseudo-time formulation") {
#if KUTTA_HAS_HDF5
  using kutta::Json;
  const auto param=fluent_examples()/"mixing_elbow/elbow_param.cas.h5";
  const auto cylinder=fluent_examples()/"adjoint_cylinder_tutorial/cylinder_tutorial.cas";
  if(!std::filesystem::exists(param) || !std::filesystem::exists(cylinder)){MESSAGE("SKIP: optional Fluent elbow_param or cylinder absent");return;}
  const auto index=[](const std::vector<kutta::setup::Setting>& settings) {
    std::map<std::string,const kutta::setup::Setting*> by;for(const auto& s:settings)by[s.path]=&s;return by;
  };
  {
    // Standard k-epsilon, SIMPLE without pseudo time: Coupled is not offered,
    // and the epsilon scheme is the solved dissipation equation's.
    auto c=kutta::fluent::read_case(param);kutta::setup::RunControls run;
    const auto settings=kutta::setup::settings_inventory(c,run);const auto by=index(settings);
    const auto& coupling=*by.at("methods/pressure-velocity-coupling");
    CHECK(coupling.value==Json("simple"));CHECK(coupling.options==Json::Array{"simple","simplec"});CHECK(coupling.option_labels.size()==2);
    CHECK(by.at("methods/pseudo-transient")->value==Json(false));
    CHECK(by.at("methods/pressure-scheme")->value==Json("standard"));
    CHECK(by.at("methods/epsilon-scheme")->active);CHECK(by.at("methods/epsilon-scheme")->editable);
    CHECK(by.at("methods/epsilon-scheme")->value==Json("first-order-upwind"));
    CHECK_FALSE(by.at("methods/omega-scheme")->active);CHECK(by.at("methods/k-scheme")->active);
    MeshTemp t;
    mesh_text(t.path("coupled.json"),"{\"v\":1,\"overrides\":{\"methods/pressure-velocity-coupling\":\"coupled\"}}");
    CHECK_THROWS_AS(kutta::setup::apply_setup(t.path("coupled.json"),settings),kutta::setup::RejectedSettings);
    mesh_text(t.path("epsilon.json"),"{\"v\":1,\"overrides\":{\"methods/epsilon-scheme\":\"second-order-upwind\"}}");
    CHECK(kutta::setup::apply_setup(t.path("epsilon.json"),settings).size()==1);
    CHECK(c.numerics.epsilon==kutta::fluent::Scheme::second_upwind);
  }
  {
    // 2D laminar, Green-Gauss cell based, saved Coupled without pseudo time
    // (which runs SIMPLEC): no turbulence scheme is active.
    auto c=kutta::fluent::read_case(cylinder);kutta::setup::RunControls run;
    const auto settings=kutta::setup::settings_inventory(c,run);const auto by=index(settings);
    CHECK(by.at("general/dimension")->value==Json(2));
    CHECK(by.at("methods/gradient")->value==Json("green-gauss-cell-based"));CHECK(by.at("methods/gradient")->value_label=="Green-Gauss Cell Based");
    CHECK(by.at("methods/pressure-velocity-coupling")->value==Json("coupled"));
    CHECK(by.at("methods/pressure-velocity-coupling")->options==Json::Array{"simple","simplec","coupled"});
    CHECK(by.at("methods/pseudo-transient")->value==Json(false));
    for(auto p:{"methods/k-scheme","methods/omega-scheme","methods/epsilon-scheme"}){CAPTURE(p);CHECK_FALSE(by.at(p)->active);CHECK_FALSE(by.at(p)->editable);}
    CHECK(by.at("methods/momentum-scheme")->active);
    // Its inlet's saved turbulence inputs are unused and cannot be changed.
    CHECK_FALSE(by.at("boundary/inlet/turbulence-specification")->active);
    for(auto p:{"boundary/inlet/k","boundary/inlet/turb-intensity"}){CAPTURE(p);CHECK_FALSE(by.at(p)->active);CHECK_FALSE(by.at(p)->editable);}
    CHECK(by.at("boundary/inlet/velocity-magnitude")->editable);
    CHECK_FALSE(by.at("controls/k-relaxation")->active);CHECK_FALSE(by.at("controls/omega-relaxation")->active);
    CHECK(by.at("controls/mom-relaxation")->active);
  }
#else
  MESSAGE("SKIP: HDF5 disabled");
#endif
}
