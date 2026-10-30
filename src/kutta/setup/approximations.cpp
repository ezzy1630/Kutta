#include "approximations.hpp"
#include "kutta/core/types.hpp"
#include <algorithm>
#include <array>

namespace kutta::setup {
namespace {
// This is the complete opt-in allowlist. Evidence is a case comparison, not
// an isolated WFGC error estimate or a bound for other meshes.
constexpr std::array allowlist{ApproximationInfo{
  Approximation::wfgc,"wfgc","recon/cell-lsf?","#t",
  "run without warped-face gradient correction (WFGC off, #f)",
  "src/kutta/fv/VERIFICATION.md#v6-initialization-comparison",
  "V6 CPU hybrid diagnostic versus Fluent: velocity L2 3.601459%, pressure L2 2.961922%, "
  "k L2 15.478919%, omega L2 8.474918%, pressure-drop error 1.185077%; "
  "combined solver discrepancy, not isolated WFGC error or a general bound"}};
}
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
