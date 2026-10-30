#include "kutta/post/aero.hpp"
#include "kutta/core/json.hpp"
#include <doctest/doctest.h>
#include <algorithm>

namespace {
using kutta::Json;
using kutta::Vec3d;
using kutta::parse_json;
using namespace kutta::post;
const Json& at(const Json& j,const std::string& key){return std::get<Json::Object>(j.value).at(key);}
double num(const Json& j){return std::get<double>(j.value);}
bool null(const Json& j){return std::holds_alternative<std::nullptr_t>(j.value);}
const Json& row(const Json& result,const std::string& name) {
  const auto& rows=std::get<Json::Array>(at(result,"rows").value);
  const auto it=std::find_if(rows.begin(),rows.end(),[&](const Json& r){return at(r,"name")==Json(name);});
  REQUIRE(it!=rows.end());return *it;
}
Loads reference_loads() {
  Loads l;l.reference={1.0,1.2,20.0,2.5,0.0};l.iterations=600;l.provenance="APPROXIMATE RESULT: test";l.approximate=true;
  return l;
}
ZoneLoads zone(std::string name,Vec3d pressure,Vec3d viscous,Vec3d at) {
  // Point loads: moments about the origin are r x F.
  const auto cross=[](Vec3d a,Vec3d b){return Vec3d{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};};
  return {std::move(name),7,1.0,pressure,viscous,cross(at,pressure),cross(at,viscous)};
}
Json spec(std::string_view text){return parse_json(text);}
}

TEST_CASE("Aero coefficients, moments and balance match a hand-computed single load") {
  // F=(100,0,-1000) N at (1,0,0.3): M=(0,1030,0); q=240 Pa, A=1, L=2.5.
  auto loads=reference_loads();loads.zones.push_back(zone("car",{100,0,-1000},{},{1,0,.3}));
  const auto result=aero(loads,parse_aero_spec(spec(R"({"v":1,"vehicle":{"front_axle":[0,0,0],"wheelbase":2.5}})")));
  CHECK(num(at(at(result,"reference"),"dynamic_pressure"))==doctest::Approx(240));
  CHECK(at(at(result,"axes"),"side")==Json(Json::Array{0,1,0}));
  const auto& total=row(result,"Total");CHECK(at(total,"kind")==Json("total"));
  CHECK(num(at(at(at(total,"force"),"total"),"drag"))==doctest::Approx(100));
  CHECK(num(at(total,"downforce"))==doctest::Approx(1000));
  const auto& c=at(total,"coefficients");
  CHECK(num(at(c,"cd"))==doctest::Approx(0.416667).epsilon(1e-6));
  CHECK(num(at(c,"cl"))==doctest::Approx(-4.166667).epsilon(1e-6));
  CHECK(num(at(c,"cla"))==doctest::Approx(-4.166667).epsilon(1e-6));
  CHECK(num(at(c,"cda"))==doctest::Approx(0.416667).epsilon(1e-6));
  CHECK(num(at(at(total,"moment"),"pitch"))==doctest::Approx(1030));
  CHECK(num(at(c,"cm_pitch"))==doctest::Approx(1.716667).epsilon(1e-6));
  CHECK(num(at(total,"downforce_share"))==doctest::Approx(1));
  CHECK(num(at(result,"lift_to_drag"))==doctest::Approx(10));
  const auto& b=at(result,"balance");
  CHECK(at(b,"status")==Json("defined"));
  CHECK(num(at(b,"rear_load"))==doctest::Approx(412));CHECK(num(at(b,"front_load"))==doctest::Approx(588));
  CHECK(num(at(b,"front_percent"))==doctest::Approx(58.8));CHECK(num(at(b,"cop_from_front_axle"))==doctest::Approx(1.03));
  CHECK(at(result,"provenance")==Json("APPROXIMATE RESULT: test"));CHECK(at(result,"approximate")==Json(true));CHECK(at(result,"iterations")==Json(600));
  // The pressure/viscous split of an all-pressure load.
  CHECK(num(at(c,"cd_viscous"))==0);CHECK(num(at(c,"cl_pressure"))==doctest::Approx(-4.166667).epsilon(1e-6));
}

TEST_CASE("Aero groups split pressure and viscous loads and transfer moments to the moment centre") {
  auto loads=reference_loads();
  loads.zones={zone("front-wing",{2,0,-300},{1,0,0},{-.5,0,.2}),zone("rear-wing",{5,0,-400},{.5,0,0},{2,0,.8}),zone("floor",{0,0,-200},{3,0,0},{1,0,0})};
  // A load file round trip keeps every value the math reads.
  loads=parse_loads(loads_json(loads));
  const auto result=aero(loads,parse_aero_spec(spec(R"({"moment_centre":[1,0,0],"groups":[{"name":"Wings","zones":["front-wing","rear-wing"]}],
    "vehicle":{"front_axle":[0,0,0],"wheelbase":2.5}})")));
  std::vector<std::string> names;for(const auto& r:std::get<Json::Array>(at(result,"rows").value))names.push_back(std::get<std::string>(at(r,"name").value));
  CHECK(names==std::vector<std::string>{"front-wing","rear-wing","floor","Wings","Total"});
  CHECK(at(result,"ungrouped")==Json(Json::Array{"floor"}));
  const auto& wings=row(result,"Wings");CHECK(at(wings,"kind")==Json("group"));
  CHECK(at(wings,"zones")==Json(Json::Array{"front-wing","rear-wing"}));
  const auto& c=at(wings,"coefficients");
  CHECK(num(at(c,"cd_pressure"))==doctest::Approx(7./240));CHECK(num(at(c,"cd_viscous"))==doctest::Approx(1.5/240));
  CHECK(num(at(c,"cl_pressure"))==doctest::Approx(-700./240));CHECK(num(at(c,"cl_viscous"))==0);
  CHECK(num(at(at(at(wings,"force"),"viscous"),"drag"))==doctest::Approx(1.5));
  CHECK(num(at(wings,"downforce_share"))==doctest::Approx(7./9));
  CHECK(num(at(row(result,"floor"),"downforce_share"))==doctest::Approx(2./9));
  // About C=(1,0,0): front -449.4, rear 404.4, floor 0 N m.
  CHECK(num(at(at(row(result,"front-wing"),"moment"),"pitch"))==doctest::Approx(-449.4));
  CHECK(num(at(at(wings,"moment"),"pitch"))==doctest::Approx(-45));
  CHECK(num(at(at(row(result,"Total"),"moment"),"pitch"))==doctest::Approx(-45));
  CHECK(num(at(at(row(result,"Total"),"coefficients"),"cm_pitch"))==doctest::Approx(-45./600));
  CHECK(num(at(result,"lift_to_drag"))==doctest::Approx(900/11.5));
  // The balance is independent of the moment centre.
  const auto& b=at(result,"balance");
  CHECK(num(at(b,"rear_load"))==doctest::Approx(342));CHECK(num(at(b,"front_load"))==doctest::Approx(558));
  CHECK(num(at(b,"front_percent"))==doctest::Approx(62));CHECK(num(at(b,"cop_from_front_axle"))==doctest::Approx(.95));
  // Car zones default to the case's force-report zones.
  loads.report_zones={"rear-wing","floor"};
  const auto reported=aero(loads,{});
  CHECK(num(at(row(reported,"Total"),"downforce"))==doctest::Approx(600));
  CHECK(at(reported,"ungrouped")==Json(Json::Array{"rear-wing","floor"}));
}

TEST_CASE("Aero axes, missing references and balance statuses") {
  auto loads=reference_loads();loads.reference.length.reset();
  ZoneLoads z{"car",3,2.0,{10,-50,3},{},{1,2,3},{}};loads.zones.push_back(z);
  // Drag +x, lift +y: side is lift x drag = -z. Inputs are normalized.
  auto result=aero(loads,parse_aero_spec(spec(R"({"axes":{"drag":[2,0,0],"lift":[0,3,0]}})")));
  CHECK(at(at(result,"axes"),"side")==Json(Json::Array{0,0,-1}));CHECK(at(at(result,"axes"),"lift")==Json(Json::Array{0,1,0}));
  const auto& total=row(result,"Total");const auto& f=at(at(total,"force"),"total");
  CHECK(num(at(f,"drag"))==10);CHECK(num(at(f,"lift"))==-50);CHECK(num(at(f,"side"))==-3);CHECK(num(at(total,"downforce"))==50);
  const auto& m=at(total,"moment");CHECK(num(at(m,"pitch"))==-3);CHECK(num(at(m,"roll"))==1);CHECK(num(at(m,"yaw"))==2);
  // No length: moment coefficients are null, force coefficients are not.
  CHECK(null(at(at(total,"coefficients"),"cm_pitch")));CHECK(num(at(at(total,"coefficients"),"cd"))==doctest::Approx(10./240));
  const auto& b=at(result,"balance");CHECK(at(b,"status")==Json("needs_vehicle"));
  for(auto key:{"front_axle","wheelbase","front_load","rear_load","front_percent","cop_from_front_axle"})CHECK(null(at(b,key)));
  // No density: every coefficient and the dynamic pressure are null.
  loads.reference.density.reset();result=aero(loads,{});
  CHECK(null(at(at(result,"reference"),"dynamic_pressure")));
  for(const auto& [key,value]:std::get<Json::Object>(at(row(result,"Total"),"coefficients").value)){CAPTURE(key);CHECK(null(value));}
  // Lift instead of downforce: loads are given, the split is not.
  loads=reference_loads();loads.zones={zone("car",{10,0,100},{},{1,0,0})};
  result=aero(loads,parse_aero_spec(spec(R"({"vehicle":{"front_axle":[0,0,0],"wheelbase":2}})")));
  const auto& lifted=at(result,"balance");CHECK(at(lifted,"status")==Json("no_downforce"));
  CHECK(num(at(lifted,"rear_load"))==doctest::Approx(-50));CHECK(num(at(lifted,"front_load"))==doctest::Approx(-50));
  CHECK(null(at(lifted,"front_percent")));CHECK(null(at(lifted,"cop_from_front_axle")));
  CHECK(num(at(result,"lift_to_drag"))==doctest::Approx(-10));
  loads.zones={zone("car",{0,0,-10},{},{})};CHECK(null(at(aero(loads,{}),"lift_to_drag")));
}

TEST_CASE("Aero specs are strict and zone names are checked against the run") {
  auto loads=reference_loads();loads.zones={zone("a",{0,0,-1},{},{}),zone("b",{0,0,-1},{},{})};
  const auto rejects=[](std::string_view text,const char* why) {
    CAPTURE(text);CHECK_THROWS_WITH_AS(parse_aero_spec(spec(text)),doctest::Contains(why),InvalidAeroSpec);
  };
  rejects(R"({"axes":{"drag":[1,0,0],"lift":[1,1,0]}})","orthogonal");
  rejects(R"({"axes":{"drag":[0,0,1]}})","orthogonal");
  rejects(R"({"axes":{"drag":[0,0,0]}})","nonzero");
  rejects(R"({"groups":[{"name":"x","zones":["a"]},{"name":"y","zones":["a","b"]}]})","more than one group");
  rejects(R"({"groups":[{"name":"x","zones":["a"]},{"name":"x","zones":["b"]}]})","duplicate group");
  rejects(R"({"groups":[{"name":"","zones":["a"]}]})","nonempty name");
  rejects(R"({"vehicle":{"front_axle":[0,0,0],"wheelbase":0}})","positive");
  rejects(R"({"vehicle":{"wheelbase":1}})","front_axle");
  rejects(R"({"drag":[1,0,0]})","unknown property");
  rejects(R"({"axes":{"side":[0,1,0]}})","axes.side: unknown property");
  rejects(R"({"groups":[{"name":"x","zones":["a"],"colour":"red"}]})","unknown property");
  rejects(R"({"v":2})","version");
  rejects(R"({"car_zones":["a","a"]})","twice");
  rejects(R"({"moment_centre":[0,0]})","[x,y,z]");
  rejects(R"([])","object");
  CHECK_THROWS_WITH_AS(aero(loads,parse_aero_spec(spec(R"({"car_zones":["a","missing"]})"))),doctest::Contains("missing"),InvalidAeroSpec);
  CHECK_THROWS_WITH_AS(aero(loads,parse_aero_spec(spec(R"({"car_zones":["a"],"groups":[{"name":"x","zones":["b"]}]})"))),doctest::Contains("not a car zone"),InvalidAeroSpec);
  CHECK(std::get<Json::Array>(at(aero(loads,parse_aero_spec(spec("{}"))),"rows").value).size()==3);
  // Strict JSON duplicates are rejected by the parser itself.
  CHECK_THROWS(parse_json(R"({"v":1,"v":1})"));
}
