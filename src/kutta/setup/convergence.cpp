#include "convergence.hpp"
#include <algorithm>
#include <cmath>

namespace kutta::setup {
bool ConvergenceControls::has_conditions() const {
  return std::any_of(residuals.begin(),residuals.end(),[](const auto& r){return r.check;}) ||
         std::any_of(reports.begin(),reports.end(),[](const auto& r){return r.active;});
}
Json::Object ConvergenceControls::json() const {
  Json::Object checked;Json::Array active;
  for(const auto& r:residuals)if(r.check)checked.emplace(r.equation,r.absolute?Json(*r.absolute):Json(nullptr));
  for(const auto& r:reports)if(r.active)active.push_back(Json::Object{{"report",r.report},{"stop_criterion",r.stop},{"previous_values",r.previous_values},{"ignore_before",r.ignore_before}});
  return {{"enabled",check},{"condition",combine==Combine::all?"all":"any"},{"residuals",std::move(checked)},{"reports",std::move(active)}};
}
std::vector<fluent::SettingIssue> criteria_issues(const ConvergenceControls& c) {
  std::vector<fluent::SettingIssue> out;
  for(const auto& r:c.residuals)if(r.check && !r.absolute)out.push_back({"monitors/residual/"+r.equation+"/absolute-criterion","a checked residual needs an absolute criterion"});
  return out;
}
std::optional<double> lag_residual(std::span<const double> m,int previous_values) {
  if(previous_values<1 || m.size()<std::size_t(previous_values)+1 || m.back()==0)return std::nullopt;
  const auto n=m.size()-1;double worst=0;
  for(std::size_t k=1;k<=std::size_t(previous_values);++k)worst=std::max(worst,std::abs(m[n-k+1]-m[n-k]));
  const auto value=worst/std::abs(m.back());
  return std::isfinite(value)?std::optional(value):std::nullopt;
}
ConvergenceMonitor::ConvergenceMonitor(ConvergenceControls controls):controls_(std::move(controls)),history_(controls_.reports.size()) {}
ConvergenceSample ConvergenceMonitor::sample(int iteration,const std::map<std::string,double>& residuals,const std::map<std::string,double>& reports) {
  ConvergenceSample out;Json::Array conditions;std::vector<bool> met;
  if(std::any_of(controls_.residuals.begin(),controls_.residuals.end(),[](const auto& r){return r.check;})) {
    // All checked residuals together form one condition.
    Json::Object equations;bool all=true;
    for(const auto& r:controls_.residuals)if(r.check) {
      const auto value=residuals.at(r.equation);const bool ok=r.absolute && value<=*r.absolute;all&=ok;
      equations.emplace(r.equation,Json::Object{{"value",value},{"criterion",r.absolute?Json(*r.absolute):Json(nullptr)},{"met",ok}});
    }
    conditions.push_back(Json::Object{{"kind","residuals"},{"met",all},{"equations",std::move(equations)}});
    met.push_back(all);if(all)out.met_by.push_back("residuals");
  }
  for(std::size_t i=0;i<controls_.reports.size();++i) {
    const auto& r=controls_.reports[i];if(!r.active)continue;
    // Only samples after ignore-before count; a window of Np differences
    // needs Np+1 of them.
    auto& h=history_[i];
    if(iteration>r.ignore_before) {
      h.push_back(reports.at(r.report));
      if(h.size()>std::size_t(r.previous_values)+1)h.erase(h.begin());
    }
    const auto value=iteration>r.ignore_before?lag_residual(h,r.previous_values):std::nullopt;
    const bool ok=value && *value<r.stop;
    conditions.push_back(Json::Object{{"kind","report"},{"report",r.report},{"value",value?Json(*value):Json(nullptr)},{"criterion",r.stop},{"met",ok}});
    met.push_back(ok);if(ok)out.met_by.push_back("report:"+r.report);
  }
  out.met=!met.empty() && (controls_.combine==ConvergenceControls::Combine::all?std::all_of(met.begin(),met.end(),[](bool b){return b;}):std::any_of(met.begin(),met.end(),[](bool b){return b;}));
  if(!out.met)out.met_by.clear();
  out.json={{"met",out.met},{"conditions",std::move(conditions)}};
  return out;
}
std::vector<std::string_view> solved_residuals(std::optional<fluent::ViscousModel> model,bool transition,bool planar) {
  std::vector<std::string_view> out(residual_equations.begin(),residual_equations.begin()+(planar?3:4));
  if(model==fluent::ViscousModel::ke_standard) out.insert(out.end(),{"k","epsilon"});
  else if(model!=fluent::ViscousModel::laminar) out.insert(out.end(),{"k","omega"});
  // Fluent's residual monitor names the intermittency equation "intermit".
  if(transition) out.push_back("intermit");
  return out;
}
} // namespace kutta::setup
