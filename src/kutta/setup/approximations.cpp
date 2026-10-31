#include "approximations.hpp"
#include "kutta/core/types.hpp"
#include <algorithm>
#include <array>

namespace kutta::setup {
namespace {
// This is the complete opt-in allowlist. Evidence is a case comparison or a
// case-specific bound, not an isolated error estimate for other meshes.
constexpr std::array allowlist{ApproximationInfo{
  Approximation::wfgc,"wfgc","recon/cell-lsf?","#t",
  "run without warped-face gradient correction (WFGC off, #f)",
  "src/kutta/fv/VERIFICATION.md#v6-initialization-comparison",
  "V6 CPU hybrid diagnostic versus Fluent: velocity L2 3.601459%, pressure L2 2.961922%, "
  "k L2 15.478919%, omega L2 8.474918%, pressure-drop error 1.185077%; "
  "combined solver discrepancy, not isolated WFGC error or a general bound",
  "Warped-face gradient correction",
  "Fluent corrects gradients for nonplanar faces. Kutta does not yet implement Fluent's warped-face correction and will run with WFGC off only when you opt in with --approximate wfgc."},
  ApproximationInfo{
  Approximation::compressibility,"compressibility","turb-compress-mod?","#t",
  "run without the turbulence compressibility (dilatation dissipation) correction",
  "validation/results/case_coverage/REPORT.md#compressibility-bound",
  "the omitted dilatation dissipation relative to dissipation is 2k/a^2 (Sarkar, TG 4.3.8); "
  "for the elbow_param and elbow-2 air elbows, Fluent's saved k_max 9.99 m^2/s^2 and the ideal-gas sound "
  "speed 346.9 m/s at the coldest saved 299.4 K bound it by 1.7e-4; other cases need their own bound",
  "Turbulence compressibility correction",
  "Fluent adds dilatation dissipation to the turbulence equations. Kutta solves constant-density flow without it and will run only when you opt in with --approximate compressibility."},
  ApproximationInfo{
  Approximation::ke_options,"ke-options","ke-source-mod/ke-time-marching?","#t",
  "run the documented standard k-epsilon (TG 4.3.1) without these undocumented Fluent k-epsilon solver switches",
  "validation/results/case_coverage/REPORT.md#standard-k-epsilon-elbows",
  "not isolated, because the switches are undocumented; elbow_param CPU cold start versus Fluent's converged "
  "saved fields, with compressibility also omitted: velocity L2 0.93%, pressure L2 6.37%, k L2 9.08%, "
  "epsilon L2 17.59%, pressure-drop error 0.91%; combined solver discrepancy, not an isolated error or a general bound",
  "k-epsilon solver options",
  "Fluent enables k-epsilon solver switches whose equations are not documented. Kutta solves the documented standard k-epsilon model and will run only when you opt in with --approximate ke-options."},
  ApproximationInfo{
  Approximation::reverse_flow,"reverse-flow","prevent-reverse-flow?","absent (Fluent 18.0 or earlier writer)",
  "run pressure outlets without reverse-flow suppression",
  "validation/results/cylinder_2d/REPORT.md#reverse-flow",
  "suppression acts only on outlet faces with inflow; for the Fluent 18.0 Re 40 cylinder, none of the outlet's 30 faces "
  "has inflow in Fluent's saved converged fields or in Kutta's converged CPU and Metal solutions, so the effect there is zero; "
  "other cases need their own check",
  "Pressure-outlet reverse flow",
  "Fluent can suppress flow entering through a pressure outlet. Cases written by Fluent 18.0 and earlier do not record whether that was on, and the earlier behaviour is undocumented. Kutta runs without suppression only when you opt in with --approximate reverse-flow."}};
}
std::span<const ApproximationInfo> approximations() {return allowlist;}
const ApproximationInfo& approximation_info(Approximation id) {
  const auto it=std::find_if(allowlist.begin(),allowlist.end(),[&](const auto& info){return info.id==id;});
  KUTTA_REQUIRE(it!=allowlist.end(),"unknown approximation id");
  return *it;
}
std::vector<Approximation> parse_approximations(std::string_view names) {
  std::vector<Approximation> out;
  for(;;) {
    const auto end=names.find(',');const auto name=names.substr(0,end);
    const auto it=std::find_if(allowlist.begin(),allowlist.end(),[&](const auto& info){return info.name==name;});
    std::string choices;
    for(const auto& info:allowlist) {if(!choices.empty()) choices+=',';choices+=info.name;}
    KUTTA_REQUIRE(it!=allowlist.end(),"unknown --approximate setting '"+std::string(name)+"'; allowed: "+choices);
    KUTTA_REQUIRE(std::find(out.begin(),out.end(),it->id)==out.end(),"duplicate --approximate setting: "+std::string(name));
    out.push_back(it->id);
    if(end==std::string_view::npos) return out;
    names.remove_prefix(end+1);
  }
}
std::string approximation_notice(std::span<const Approximation> used) {
  if(used.empty()) return {};
  std::string out="APPROXIMATE RESULT: not a faithful Fluent reproduction";
  for(auto id:used) {
    const auto& info=approximation_info(id);
    out+="; "+std::string(info.name)+": "+std::string(info.setting)+" Fluent "+std::string(info.fluent_value)+
      " -> Kutta "+std::string(info.kutta_action)+"; "+std::string(info.measured_effect)+"; evidence: "+std::string(info.effect_documentation);
  }
  return out;
}
} // namespace kutta::setup
