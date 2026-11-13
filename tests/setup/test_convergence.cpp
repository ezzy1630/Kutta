#include "kutta/setup/convergence.hpp"
#include "kutta/core/json.hpp"
#include <doctest/doctest.h>

namespace {
using namespace kutta::setup;
using kutta::Json;
const Json& at(const Json& j,const std::string& key){return std::get<Json::Object>(j.value).at(key);}
std::map<std::string,double> residuals(double continuity,double k) {
  std::map<std::string,double> out;for(auto eq:residual_equations)out.emplace(eq,1.0);
  out["continuity"]=continuity;out["k"]=k;return out;
}
ConvergenceControls controls() {
  ConvergenceControls c;c.check=true;
  for(auto eq:residual_equations)c.residuals.push_back({std::string(eq),false,1e-3});
  return c;
}
// Report values at consecutive monitor samples; null entries mean no value.
std::vector<std::optional<double>> report_values(const ReportCriterion& r,const std::vector<double>& m,int first_iteration=1) {
  auto c=controls();c.reports={r};ConvergenceMonitor monitor(c);std::vector<std::optional<double>> out;
  for(std::size_t i=0;i<m.size();++i) {
    const auto s=monitor.sample(first_iteration+int(i),residuals(1,1),{{r.report,m[i]}});
    const Json state=s.json;const auto& condition=std::get<Json::Array>(at(state,"conditions").value).at(0);
    const auto& value=at(condition,"value");
    out.push_back(std::holds_alternative<double>(value.value)?std::optional(std::get<double>(value.value)):std::nullopt);
    CHECK(s.met==(out.back() && *out.back()<r.stop));CHECK(at(condition,"met")==Json(s.met));
  }
  return out;
}
}

TEST_CASE("Residual convergence is one condition over every checked equation") {
  auto c=controls();
  CHECK_FALSE(c.has_conditions());
  CHECK(ConvergenceMonitor(c).sample(20,residuals(0,0),{}).met==false); // nothing to meet
  c.residuals[0].check=true;c.residuals[4].check=true;c.residuals[4].absolute=2e-3; // continuity 1e-3, k 2e-3
  CHECK(c.has_conditions());
  ConvergenceMonitor monitor(c);
  auto s=monitor.sample(20,residuals(5e-4,3e-3),{});
  CHECK_FALSE(s.met);CHECK(s.met_by.empty());
  const Json state=s.json;const auto& condition=std::get<Json::Array>(at(state,"conditions").value).at(0);
  CHECK(at(condition,"kind")==Json("residuals"));
  const auto& equations=std::get<Json::Object>(at(condition,"equations").value);
  CHECK(equations.size()==2); // unchecked equations, all at 1.0 here, are not listed
  CHECK(at(equations.at("continuity"),"met")==Json(true));CHECK(at(equations.at("k"),"met")==Json(false));
  CHECK(at(equations.at("k"),"criterion")==Json(2e-3));CHECK(at(equations.at("k"),"value")==Json(3e-3));
  // At the criterion counts as met.
  s=monitor.sample(40,residuals(1e-3,2e-3),{});
  CHECK(s.met);CHECK(s.met_by==std::vector<std::string>{"residuals"});
  CHECK(at(Json(s.json),"met")==Json(true));
}

TEST_CASE("Report convergence uses Fluent lag residuals after ignore-before") {
  // Np=2: values at the 3rd..5th samples are max of the last two differences
  // over the newest value: 10/105, 5/104, 1/104.5.
  ReportCriterion r{"drag",true,1e-2,2,0};
  const auto v=report_values(r,{100,110,105,104,104.5});
  CHECK_FALSE(v[0]);CHECK_FALSE(v[1]);
  CHECK(*v[2]==doctest::Approx(10/105.));CHECK(*v[3]==doctest::Approx(5/104.));CHECK(*v[4]==doctest::Approx(1/104.5));
  // Np=3 meets at the 6th sample (1/104.6); ignoring iterations up to 3 drops
  // the large early changes but needs four later samples: null at 6, met at 7.
  r={"drag",true,1e-2,3,0};
  auto w=report_values(r,{100,110,105,104,104.5,104.6,104.55});
  CHECK(*w[4]==doctest::Approx(5/104.5));CHECK(*w[5]==doctest::Approx(1/104.6));
  r.ignore_before=3;w=report_values(r,{100,110,105,104,104.5,104.6,104.55});
  for(int i=0;i<6;++i){CAPTURE(i);CHECK_FALSE(w[i]);}
  CHECK(*w[6]==doctest::Approx(.5/104.55));
  // A zero newest value is never met.
  r={"drag",true,1e30,1,0};w=report_values(r,{0,0,1,0});
  CHECK_FALSE(w[0]);CHECK_FALSE(w[1]);CHECK(*w[2]==doctest::Approx(1));CHECK_FALSE(w[3]);
  CHECK_FALSE(lag_residual(std::vector<double>{1,2},0));
}

TEST_CASE("All needs every condition at one sample and any needs one") {
  auto c=controls();c.residuals[0].check=true;c.reports={{"drag",true,1e-3,1,0},{"lift",false}};
  const auto run=[&](ConvergenceControls::Combine combine) {
    c.combine=combine;ConvergenceMonitor monitor(c);
    monitor.sample(1,residuals(1e-4,1),{{"drag",100},{"lift",1}});
    return monitor.sample(2,residuals(1e-4,1),{{"drag",110},{"lift",1}}); // residuals met, drag not
  };
  auto s=run(ConvergenceControls::Combine::all);
  CHECK_FALSE(s.met);CHECK(s.met_by.empty());
  const Json state=s.json;CHECK(std::get<Json::Array>(at(state,"conditions").value).size()==2); // inactive lift is absent
  s=run(ConvergenceControls::Combine::any);
  CHECK(s.met);CHECK(s.met_by==std::vector<std::string>{"residuals"});
  // The setup object lists checked residuals and active reports only.
  const Json setup=c.json();
  CHECK(dump_json(setup)==R"({"condition":"any","enabled":true,"reports":[{"ignore_before":0,"previous_values":1,"report":"drag","stop_criterion":0.001}],"residuals":{"continuity":0.001}})");
  // A checked residual needs a criterion before a solve starts.
  c.residuals[2].absolute.reset();CHECK(criteria_issues(c).empty());
  c.residuals[2].check=true;
  const auto issues=criteria_issues(c);REQUIRE(issues.size()==1);CHECK(issues[0].path=="monitors/residual/y-velocity/absolute-criterion");
}

namespace {
ConvergenceControls settling(double tolerance=1e-3,int window=100) {
  auto c=controls();c.residuals[0].check=true; // continuity at 1e-3
  c.reports={{"drag",false},{"lift",false},{"pressure",false,1e-3,5,100,false}}; // the last is not a force report
  c.settle={true,tolerance,window};
  return c;
}
std::map<std::string,double> forces(double drag,double lift){return {{"drag",drag},{"lift",lift},{"pressure",0}};}
Json condition(const ConvergenceSample& s,const std::string& kind) {
  const Json state=s.json;
  for(const auto& c:std::get<Json::Array>(at(state,"conditions").value))if(at(c,"kind")==Json(kind))return c;
  throw std::runtime_error("no condition "+kind);
}
}

TEST_CASE("Settled forces need a full window and a confirming sample one iteration later") {
  const auto c=settling();
  CHECK(c.settle_condition());CHECK(c.has_conditions());
  CHECK(dump_json(Json(c.json()))==R"({"condition":"all","enabled":true,"forces":{"reports":["drag","lift"],"tolerance":0.001,"window":100},"reports":[],"residuals":{"continuity":0.001}})");
  ConvergenceMonitor monitor(c);
  // Samples every 20 iterations; the window needs one at or before n-100.
  for(int it=20;it<=100;it+=20) {
    const auto s=monitor.sample(it,residuals(1e-4,1),forces(2,1));
    CAPTURE(it);CHECK_FALSE(s.met);CHECK_FALSE(s.confirm);CHECK(at(at(condition(s,"forces"),"reports"),"drag")==Json(nullptr));
  }
  auto s=monitor.sample(120,residuals(1e-4,1),forces(2,1));
  CHECK_FALSE(s.met);CHECK(s.confirm);CHECK(at(condition(s,"forces"),"settled")==Json(true));CHECK(at(condition(s,"forces"),"met")==Json(false));
  CHECK(at(at(condition(s,"forces"),"reports"),"drag")==Json(0.0));
  s=monitor.sample(121,residuals(1e-4,1),forces(2.0015,1.0005)); // within 1e-3 of the means
  CHECK(s.met);CHECK_FALSE(s.confirm);CHECK(s.met_by==std::vector<std::string>{"residuals","forces"});CHECK_FALSE(s.residuals_waived);
}

TEST_CASE("A cycle aliased by the monitor stride fails the off-stride confirmation") {
  // Period two: every even iteration reads 2, every odd one 2.1. The stride
  // samples all read 2 and look settled; the confirming sample does not.
  ConvergenceMonitor monitor(settling());
  for(int it=20;it<=200;it+=20) {
    const auto s=monitor.sample(it,residuals(1e-4,1),forces(2,1));
    CHECK_FALSE(s.met);
    if(s.confirm) {
      const auto odd=monitor.sample(it+1,residuals(1e-4,1),forces(2.1,1));
      CAPTURE(it);CHECK_FALSE(odd.met);CHECK(at(condition(odd,"forces"),"settled")==Json(false));
    }
  }
}

TEST_CASE("Drifting forces never settle and stalled residuals are waived only after twice the window") {
  {
    ConvergenceMonitor monitor(settling());
    for(int it=20;it<=1000;it+=20) {
      const auto s=monitor.sample(it,residuals(1e-4,1),forces(2+it*4e-5,1)); // 0.2% per 100 iterations
      CAPTURE(it);CHECK_FALSE(s.met);CHECK_FALSE(s.confirm);
    }
  }
  ConvergenceMonitor monitor(settling());
  bool waived=false;int at_iteration=0;
  for(int it=20;it<=400 && !waived;it+=20) {
    auto s=monitor.sample(it,residuals(5e-3,1),forces(2,1)); // continuity stalled at 5e-3
    if(!s.confirm){CHECK_FALSE(s.met);continue;}
    s=monitor.sample(it+1,residuals(5e-3,1),forces(2,1));
    // Forces alone settle at 120 but the residual condition holds the run
    // until they have held for 200 iterations.
    CAPTURE(it);CHECK(it>=220);
    waived=s.met;at_iteration=it+1;
    if(waived){CHECK(s.residuals_waived);CHECK(s.met_by==std::vector<std::string>{"forces"});CHECK(at(condition(s,"residuals"),"waived")==Json(true));}
  }
  CHECK(waived);CHECK(at_iteration==221);
}

TEST_CASE("A report near zero is judged against one percent of the largest report") {
  // Lift is zero on average and wanders by 1e-4 of the drag: 1e-2 of 1% of
  // drag. Relative to its own mean it would never settle.
  auto c=settling(3e-2,100);c.residuals[0].check=false;
  ConvergenceMonitor monitor(c);ConvergenceSample s;
  for(int it=20;it<=140;it+=20)s=monitor.sample(it,residuals(1,1),forces(1,it%40?1e-4:-1e-4));
  CHECK(s.confirm);
  CHECK(std::get<double>(at(at(condition(s,"forces"),"reports"),"lift").value)==doctest::Approx(2e-4/1e-2));
}
