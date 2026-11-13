#include "convergence.hpp"
#include <algorithm>
#include <cmath>

namespace kutta::setup {
bool ConvergenceControls::settle_condition() const {
  return settle.active && std::any_of(reports.begin(),reports.end(),[](const auto& r){return r.force;});
}
bool ConvergenceControls::has_conditions() const {
  return std::any_of(residuals.begin(),residuals.end(),[](const auto& r){return r.check;}) ||
         std::any_of(reports.begin(),reports.end(),[](const auto& r){return r.active;}) || settle_condition();
}
Json::Object ConvergenceControls::json() const {
  Json::Object checked;Json::Array active;
  for(const auto& r:residuals)if(r.check)checked.emplace(r.equation,r.absolute?Json(*r.absolute):Json(nullptr));
  for(const auto& r:reports)if(r.active)active.push_back(Json::Object{{"report",r.report},{"stop_criterion",r.stop},{"previous_values",r.previous_values},{"ignore_before",r.ignore_before}});
  Json::Object out{{"enabled",check},{"condition",combine==Combine::all?"all":"any"},{"residuals",std::move(checked)},{"reports",std::move(active)}};
  if(settle_condition()) {
    Json::Array watched;for(const auto& r:reports)if(r.force)watched.push_back(r.report);
    out.emplace("forces",Json::Object{{"tolerance",settle.tolerance},{"window",settle.window},{"reports",std::move(watched)}});
  }
  return out;
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
ConvergenceMonitor::ConvergenceMonitor(ConvergenceControls controls):controls_(std::move(controls)),history_(controls_.reports.size()) {
  if(controls_.settle_condition())for(const auto& r:controls_.reports)if(r.force)force_reports_.push_back(r.report);
}
std::vector<std::optional<double>> ConvergenceMonitor::ranges(int iteration,int window) const {
  std::vector<std::optional<double>> out(force_reports_.size());
  const auto from=iteration-window;
  if(forces_.empty() || forces_.front().iteration>from)return out;
  std::vector<double> lo(out.size(),INFINITY),hi(out.size(),-INFINITY),mean(out.size(),0);
  int n=0;
  for(auto it=forces_.rbegin();it!=forces_.rend() && it->iteration>=from;++it,++n)
    for(std::size_t i=0;i<out.size();++i){lo[i]=std::min(lo[i],it->values[i]);hi[i]=std::max(hi[i],it->values[i]);mean[i]+=it->values[i];}
  double largest=0;
  for(auto& m:mean){m=std::abs(m/n);largest=std::max(largest,m);}
  for(std::size_t i=0;i<out.size();++i) {
    const auto scale=std::max(mean[i],.01*largest);
    if(scale>0 && std::isfinite(hi[i]-lo[i]))out[i]=(hi[i]-lo[i])/scale;
  }
  return out;
}
ConvergenceSample ConvergenceMonitor::sample(int iteration,const std::map<std::string,double>& residuals,const std::map<std::string,double>& reports) {
  ConvergenceSample out;Json::Array conditions;
  // Each condition's state; the residual condition may be waived below.
  struct State {bool met;std::string name;};std::vector<State> states;
  std::optional<std::size_t> residual_index;
  if(std::any_of(controls_.residuals.begin(),controls_.residuals.end(),[](const auto& r){return r.check;})) {
    // All checked residuals together form one condition.
    Json::Object equations;bool all=true;
    for(const auto& r:controls_.residuals)if(r.check) {
      const auto value=residuals.at(r.equation);const bool ok=r.absolute && value<=*r.absolute;all&=ok;
      equations.emplace(r.equation,Json::Object{{"value",value},{"criterion",r.absolute?Json(*r.absolute):Json(nullptr)},{"met",ok}});
    }
    residual_index=conditions.size();
    conditions.push_back(Json::Object{{"kind","residuals"},{"met",all},{"equations",std::move(equations)}});
    states.push_back({all,"residuals"});
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
    states.push_back({ok,"report:"+r.report});
  }
  bool settled=false,settled_twice=false;
  if(!force_reports_.empty()) {
    const auto& s=controls_.settle;
    Forces f{iteration,{}};
    for(const auto& name:force_reports_)f.values.push_back(reports.at(name));
    forces_.push_back(std::move(f));
    // Keep one sample at or before the start of the longer (2x) window.
    while(forces_.size()>1 && forces_[1].iteration<=iteration-2*s.window)forces_.erase(forces_.begin());
    const auto within=[&](const std::vector<std::optional<double>>& r) {
      return std::all_of(r.begin(),r.end(),[&](const auto& v){return v && *v<=s.tolerance;});
    };
    const auto window=ranges(iteration,s.window);
    settled=within(window);settled_twice=settled && within(ranges(iteration,2*s.window));
    // Two samples one iteration apart: the first is on the monitor stride,
    // the second off it.
    const bool confirmed=settled && last_iteration_==iteration-1 && last_settled_;
    Json::Object values;
    for(std::size_t i=0;i<force_reports_.size();++i)values.emplace(force_reports_[i],window[i]?Json(*window[i]):Json(nullptr));
    conditions.push_back(Json::Object{{"kind","forces"},{"met",confirmed},{"settled",settled},{"criterion",s.tolerance},
                                      {"window",s.window},{"reports",std::move(values)}});
    states.push_back({confirmed,"forces"});
    // A stalled residual must not hold a run whose forces stopped changing
    // long ago; the waiver needs the same confirmation.
    if(residual_index && !states[*residual_index].met && settled_twice && last_iteration_==iteration-1 && last_settled_twice_) {
      states[*residual_index].met=true;out.residuals_waived=true;
      std::get<Json::Object>(conditions[*residual_index].value).emplace("waived",true);
    }
    const auto combined=[&](bool forces_met) {
      std::vector<bool> m;for(const auto& st:states)m.push_back(st.name=="forces"?forces_met:st.met);
      return controls_.combine==ConvergenceControls::Combine::all?std::all_of(m.begin(),m.end(),[](bool b){return b;})
                                                                  :std::any_of(m.begin(),m.end(),[](bool b){return b;});
    };
    // Settled but not yet confirmed: ask for the next iteration's sample when
    // it is what stands between this sample and convergence. A residual
    // waiver can only come with that confirmation too.
    if(settled && !confirmed && !combined(false)) {
      auto waived=states;
      if(residual_index && settled_twice)waived[*residual_index].met=true;
      std::vector<bool> m;for(const auto& st:waived)m.push_back(st.name=="forces" || st.met);
      out.confirm=controls_.combine==ConvergenceControls::Combine::all?std::all_of(m.begin(),m.end(),[](bool b){return b;}):true;
    }
    last_iteration_=iteration;last_settled_=settled;last_settled_twice_=settled_twice;
  }
  // A waived residual condition counts toward the combination but was not met.
  for(const auto& st:states)if(st.met && !(out.residuals_waived && st.name=="residuals"))out.met_by.push_back(st.name);
  out.met=!states.empty() && (controls_.combine==ConvergenceControls::Combine::all?std::all_of(states.begin(),states.end(),[](const State& st){return st.met;})
                                                                                   :std::any_of(states.begin(),states.end(),[](const State& st){return st.met;}));
  if(!out.met){out.met_by.clear();out.residuals_waived=false;out.confirm=out.confirm && !out.met;}
  else out.confirm=false;
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
