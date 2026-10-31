#include "settings.hpp"
#include "kutta/mesh/detail.hpp"
#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <limits>
#include <set>

namespace kutta::setup {
namespace {
using namespace fluent;
std::string segment(std::string_view text) {
  // Percent-encode rather than slugify: Fluent names are exact identities.
  constexpr char hex[]="0123456789abcdef";std::string out;
  for(unsigned char c:text)if((c>='a' && c<='z') || (c>='0' && c<='9') || c=='-' || c=='_' || c=='.')out+=char(c);else {out+='%';out+=hex[c>>4];out+=hex[c&15];}
  return out.empty()?"%00":out;
}
std::string key_path(std::string_view key) {
  std::string out;for(;;){const auto end=key.find('/');if(!out.empty())out+='/';out+=segment(key.substr(0,end));if(end==std::string_view::npos)break;key.remove_prefix(end+1);}return out;
}
Json value(const SExpr* v) {
  if(!v)return nullptr;
  if(v->kind==SExpr::Kind::boolean)return v->boolean;
  if(const auto n=v->number())return std::isfinite(*n)?Json(*n):Json(nullptr);
  if(auto text=v->atom_text())return std::string(*text);
  if(v->kind==SExpr::Kind::list && !v->tail && v->elements.size()==3 && std::all_of(v->elements.begin(),v->elements.end(),[](const auto& x){return x.number().has_value();}))return Json::Array{*v->elements[0].number(),*v->elements[1].number(),*v->elements[2].number()};
  return to_string(*v);
}
Json value(const std::optional<ScalarValue>& v) {if(!v)return nullptr;const auto n=v->constant?v->constant:v->resolved;return n?Json(*n):Json(v->binding);}
Json vector_value(const ScalarVector& v) {Json::Array a;for(const auto& x:v)a.push_back(value(x));return a;}
std::string kind(const Json& v) {if(std::holds_alternative<bool>(v.value))return "bool";if(std::holds_alternative<double>(v.value))return "number";if(std::holds_alternative<Json::Array>(v.value))return "vec3";return "string";}
std::string group(std::string_view key) {
  if(key.starts_with("reference"))return "Reference Values";
  if(key.starts_with("monitor") || key.find("residual")!=std::string_view::npos)return "Monitors";
  if(key.find("niter")!=std::string_view::npos || key=="number-of-iterations")return key.starts_with("hyb")?"Initialization":"Run Calculation";
  if(key.starts_with("hyb") || key.starts_with("fmg") || key.find("/default")!=std::string_view::npos)return "Initialization";
  if(key.find("relax")!=std::string_view::npos || key.find("pseudo")!=std::string_view::npos || key.find("courant")!=std::string_view::npos)return "Controls";
  if(key.starts_with("recon") || key.find("/scheme")!=std::string_view::npos || key.starts_with("pbs/") || key.starts_with("wall-distance"))return "Methods";
  if(key.starts_with("kw") || key.starts_with("ke") || key.starts_with("sst") || key.starts_with("turb") || key.find("energy")!=std::string_view::npos)return "Models";
  return "General";
}
std::string scheme_label(Scheme s) {
  switch(s){case Scheme::first_upwind:return "First Order Upwind";case Scheme::second_upwind:return "Second Order Upwind";case Scheme::power_law:return "Power Law";case Scheme::central:return "Central Differencing";
    case Scheme::quick:return "QUICK";case Scheme::muscl:return "Third-Order MUSCL";case Scheme::bounded_central:return "Bounded Central Differencing";case Scheme::standard_pressure:return "Standard";
    case Scheme::linear_pressure:return "Linear";case Scheme::second_pressure:return "Second Order";case Scheme::body_force_weighted:return "Body Force Weighted";case Scheme::presto:return "PRESTO!";}
  return "";
}
std::string viscous_label(ViscousModel m) {
  switch(m){case ViscousModel::kw_sst:return "k-omega SST";case ViscousModel::kw_bsl:return "k-omega BSL";case ViscousModel::kw_standard:return "k-omega Standard";case ViscousModel::ke_standard:return "k-epsilon Standard";
    case ViscousModel::ke_rng:return "k-epsilon RNG";case ViscousModel::ke_realizable:return "k-epsilon Realizable";case ViscousModel::spalart_allmaras:return "Spalart-Allmaras";case ViscousModel::laminar:return "Laminar";}
  return "";
}
// residuals/settings rows as Fluent 22.2 saves them (the mixing elbow):
// (equation monitor? <number> check-convergence? absolute-criterion <number> <number>).
// Only that layout is read; any other shape leaves every check off and every
// criterion unknown rather than guessing. Returns whether every row maps to a
// residual Kutta reports, i.e. whether the raw value is fully represented.
bool residual_criteria(const SExpr* v,std::vector<ResidualCriterion>& out) {
  if(!v || v->kind!=SExpr::Kind::list || v->tail)return false;
  std::map<std::string,std::pair<bool,double>> rows;
  for(const auto& e:v->elements) {
    if(e.kind!=SExpr::Kind::list || e.tail || e.elements.size()!=7 || !e.elements[0].atom_text() || e.elements[3].kind!=SExpr::Kind::boolean)return false;
    const auto criterion=e.elements[4].number();
    if(!criterion || !std::isfinite(*criterion) || *criterion<0 || !rows.emplace(std::string(*e.elements[0].atom_text()),std::pair{e.elements[3].boolean,*criterion}).second)return false;
  }
  bool complete=true;
  for(const auto& [name,row]:rows) {
    const auto it=std::find_if(out.begin(),out.end(),[&](const auto& r){return r.equation==name;});
    if(it==out.end())complete=false;else{it->check=row.first;it->absolute=row.second;}
  }
  return complete;
}
// Label the current value from the parallel option lists.
void label_value(Setting& s) {for(std::size_t i=0;i<s.options.size() && i<s.option_labels.size();++i)if(s.options[i]==s.value)s.value_label=std::get<std::string>(s.option_labels[i].value);}
class Inventory {
  // A deque keeps references from add() valid while later entries are added.
  std::deque<Setting> out;
  std::map<const RawSettings*,std::set<std::string>> consumed;
 public:
  void consume(const RawSettings& raw,std::string key){consumed[&raw].insert(std::move(key));}
  Setting& add(std::string path,std::string label,std::string group,Json val,std::string unit="",std::string type="") {
    Setting s;s.path=std::move(path);s.label=std::move(label);s.group=std::move(group);s.unit=std::move(unit);s.value=std::move(val);s.kind=type.empty()?kind(s.value):std::move(type);out.push_back(std::move(s));return out.back();
  }
  void raw(const RawSettings& raw,std::string prefix,std::string fixed_group="") {
    std::set<std::string> names;
    for(const auto& entry:raw.alist.elements) {
      if(entry.elements.empty() || !entry.elements[0].atom_text())continue;
      const std::string key(*entry.elements[0].atom_text());if(consumed[&raw].contains(key) || !names.insert(key).second)continue;
      auto& s=add(prefix+key_path(key),key,fixed_group.empty()?group(key):fixed_group,value(raw.value(key)));
      if(const auto* v=raw.value(key);v && v->kind==SExpr::Kind::integer)s.kind="integer";
    }
  }
  void editable(Setting& s,std::function<void(const Json&)> setter,std::optional<double> min={},std::optional<double> max={}) {s.editable=true;s.reason.clear();s.apply=std::move(setter);s.min=min;s.max=max;}
  void depends(Setting& s,const std::string& control,Json::Array when){s.depends_on=control;s.active_when=std::move(when);}
  Setting& scalar(std::string path,std::string label,std::string group,std::optional<ScalarValue>& v,std::string unit,bool active,std::optional<double> min={},std::optional<double> max={}) {
    auto& s=add(std::move(path),std::move(label),std::move(group),value(v),std::move(unit),"number");
    if(active && v && (v->constant || v->resolved))editable(s,[target=&v](const Json& j){ScalarValue v;v.kind=ScalarValue::Kind::constant;v.constant=std::get<double>(j.value);*target=std::move(v);},min,max);
    return s;
  }
  void number(std::string path,std::string label,std::string group,std::optional<double>& v,std::string unit,bool active,std::optional<double> min={},std::optional<double> max={}) {
    auto& s=add(std::move(path),std::move(label),std::move(group),v?Json(*v):Json(nullptr),std::move(unit),"number");
    if(active && v)editable(s,[target=&v](const Json& j){*target=std::get<double>(j.value);},min,max);
  }
  Setting& vector(std::string path,std::string label,ScalarVector& v,std::string unit,bool active,bool nonzero=false) {
    auto& s=add(std::move(path),std::move(label),"Boundary Conditions",vector_value(v),std::move(unit),"vec3");
    if(active && std::all_of(v.begin(),v.end(),[](const auto& x){return x && (x->constant || x->resolved);})) {
      editable(s,[target=&v](const Json& j){const auto& a=std::get<Json::Array>(j.value);for(int i=0;i<3;++i){ScalarValue x;x.kind=ScalarValue::Kind::constant;x.constant=std::get<double>(a[i].value);(*target)[i]=x;}},-1e30,1e30);
      if(nonzero)s.validate=[](const Json& j){const auto& a=std::get<Json::Array>(j.value);const double n=std::hypot(std::get<double>(a[0].value),std::get<double>(a[1].value),std::get<double>(a[2].value));return n>0?std::string{}:"direction must be nonzero";};
    }
    return s;
  }
  std::vector<Setting> finish() {
    std::vector<Setting> sorted(std::make_move_iterator(out.begin()),std::make_move_iterator(out.end()));
    std::sort(sorted.begin(),sorted.end(),[](const auto& a,const auto& b){return a.path<b.path;});
    std::map<std::string,const Json*> values;
    for(std::size_t i=0;i<sorted.size();++i){KUTTA_REQUIRE(i==0 || sorted[i-1].path!=sorted[i].path,"duplicate settings identity: "+sorted[i].path);values.emplace(sorted[i].path,&sorted[i].value);}
    // Activity follows the controlling value as saved, so it cannot disagree
    // with the inventory the app displays.
    for(auto& s:sorted)if(s.depends_on){const auto it=values.find(*s.depends_on);KUTTA_REQUIRE(it!=values.end(),"missing controlling setting: "+*s.depends_on);s.active=std::find(s.active_when.begin(),s.active_when.end(),*it->second)!=s.active_when.end();}
    return sorted;
  }
};
}
std::vector<Setting> settings_inventory(FluentCase& c,RunControls& run,bool seeded) {
  Inventory b;
  constexpr double positive=double(std::numeric_limits<float>::min());
  auto& general=b.add("general/operating-pressure","Operating pressure","General",c.operating_pressure?Json(*c.operating_pressure):Json(nullptr),"Pa","number");
  general.reason="not yet configurable in Kutta";b.consume(c.variables,"operating-pressure");
  b.add("general/dimension","Dimension","General",c.dimension?Json(*c.dimension):Json(nullptr),"","integer");
  b.add("models/viscous","Viscous model","Models",c.viscous.model?Json(std::string(name(*c.viscous.model))):Json(nullptr),"","enum").value_label=c.viscous.model?viscous_label(*c.viscous.model):"";
  b.add("models/energy","Energy equation","Models",c.energy?Json(*c.energy):Json(nullptr),"","bool");
  {
    auto& s=b.add("models/transition","Transition model","Models",c.viscous.transition?Json(std::string(name(*c.viscous.transition))):Json(nullptr),"","enum");
    s.options={"none","gamma-transport-eqn","gamma-algebraic","transition-sst"};
    s.option_labels={"None","Intermittency Transition Model","Algebraic Transition Model","Transition SST"};label_value(s);
  }
  const bool transition=c.viscous.transition==TransitionModel::intermittency;
  b.consume(c.variables,"case-config");b.raw(c.config,"general/case-config/");
  b.consume(c.variables,"materials");
  for(auto& m:c.materials) {
    const auto prefix="materials/"+segment(m.name)+"/";
    const bool active=std::any_of(c.cell_zones.begin(),c.cell_zones.end(),[&](const auto& z){return z.material==m.name && z.type==ZoneType::fluid;});
    b.scalar(prefix+"density","Density","Materials",m.density,"kg/m^3",active && m.nonconstant_density==false,positive,1e30);
    b.scalar(prefix+"viscosity","Dynamic viscosity","Materials",m.viscosity,"Pa s",active && m.nonconstant_viscosity==false,positive,1e30);
    b.consume(m.raw,"density");b.consume(m.raw,"viscosity");b.raw(m.raw,prefix+"raw/","Materials");
    b.add(prefix+"type","Material type","Materials",m.type,"","enum").value_label=m.type;
  }
  for(auto& z:c.cell_zones){const auto prefix="cell-zone/"+segment(z.name)+"/";b.add(prefix+"type","Zone type","Cell Zone Conditions",z.fluent_type,"","enum").value_label=z.fluent_type;b.raw(z.raw,prefix+"raw/","Cell Zone Conditions");}
  const bool omega_model=c.viscous.model==ViscousModel::kw_sst || c.viscous.model==ViscousModel::kw_bsl || c.viscous.model==ViscousModel::kw_standard;
  const bool epsilon_model=c.viscous.model==ViscousModel::ke_standard || c.viscous.model==ViscousModel::ke_rng || c.viscous.model==ViscousModel::ke_realizable;
  // Laminar (and one-equation) cases read no two-equation turbulence inputs.
  const bool two_equation=omega_model || epsilon_model;
  for(auto& z:c.boundaries) {
    const auto prefix="boundary/"+segment(z.name)+"/";
    b.add(prefix+"type","Boundary type","Boundary Conditions",z.fluent_type,"","enum").value_label=z.fluent_type;
    const bool inlet=z.type==ZoneType::velocity_inlet,outlet=z.type==ZoneType::pressure_outlet;
    const auto velocity_spec=prefix+"velocity-specification";
    if(inlet) {
      auto& s=b.add(velocity_spec,"Velocity specification method","Boundary Conditions",z.velocity_spec?Json(int(*z.velocity_spec)):Json(nullptr),"","enum");
      s.options={0,1,2};s.option_labels={"Magnitude and Direction","Components","Magnitude, Normal to Boundary"};label_value(s);b.consume(z.raw,"velocity-spec");
    }
    // Inlet velocity inputs follow Fluent's specification method; other zone
    // types never use them as inlet velocity.
    const auto velocity=[&](Setting& s,Json::Array when){if(inlet)b.depends(s,velocity_spec,std::move(when));else s.active=false;};
    velocity(b.scalar(prefix+"velocity-magnitude","Velocity magnitude","Boundary Conditions",z.velocity_magnitude,"m/s",inlet && (z.velocity_spec==VelocitySpec::normal || z.velocity_spec==VelocitySpec::magnitude_direction),0,1e15),{0,2});
    velocity(b.vector(prefix+"velocity","Velocity components",z.velocity,"m/s",inlet && z.velocity_spec==VelocitySpec::components),{1});
    velocity(b.vector(prefix+"direction","Flow direction",z.direction,"1",inlet && z.velocity_spec==VelocitySpec::magnitude_direction,true),{0});
    b.scalar(prefix+"gauge-pressure","Gauge pressure","Boundary Conditions",z.gauge_pressure,"Pa",outlet,-1e30,1e30).active=outlet;
    for(auto key:{"vmag","p","u","v","w","ni","nj","nk","direction-0","direction-1","direction-2","velocity-vector","flow-direction-component","reference-direction-components","flow-direction"})b.consume(z.raw,key);
    const auto spec=z.turbulence_spec;const auto spec_path=prefix+"turbulence-specification";
    // The values the adapter reads for each spec (setup/fluent.cpp
    // turbulence()). A switch is offered only to specs whose values were saved,
    // so it can never run on a guessed value.
    const std::array<std::vector<std::string>,4> required{{{"k","omega"},{"turb-intensity","turb-length-scale"},{"turb-intensity","turb-viscosity-ratio"},{"turb-intensity","turb-hydraulic-diam"}}};
    const Json::Array labels{epsilon_model?"K and Epsilon":"K and Omega","Intensity and Length Scale","Intensity and Viscosity Ratio","Intensity and Hydraulic Diameter"};
    const auto saved=[&](const std::string& key){const auto it=z.turbulence.find(key);if(it==z.turbulence.end())return false;const auto n=it->second.constant?it->second.constant:it->second.resolved;return n && std::isfinite(*n);};
    auto& specification=b.add(spec_path,"Turbulence specification method","Boundary Conditions",spec?Json(int(*spec)):Json(nullptr),"","enum");
    specification.options={0,1,2,3};specification.option_labels=labels;label_value(specification);specification.active=spec.has_value() && two_equation;b.consume(z.raw,"ke-spec");
    if((inlet || outlet) && omega_model && spec) {
      specification.options.clear();specification.option_labels.clear();
      for(int i=0;i<4;++i)if(i==int(*spec) || std::all_of(required[i].begin(),required[i].end(),saved)){specification.options.push_back(i);specification.option_labels.push_back(labels[i]);}
      b.editable(specification,[target=&z.turbulence_spec](const Json& j){*target=TurbulenceSpec(int(std::get<double>(j.value)));});
    }
    for(auto& [key,v]:z.turbulence) {
      static const std::map<std::string,std::tuple<std::string,std::string,std::vector<int>>> known{
        {"turb-intensity",{"Turbulent intensity","1",{1,2,3}}},{"turb-length-scale",{"Turbulent length scale","m",{1}}},{"turb-viscosity-ratio",{"Turbulent viscosity ratio","1",{2}}},
        {"turb-hydraulic-diam",{"Hydraulic diameter","m",{3}}},{"k",{"Turbulent kinetic energy","m^2/s^2",{0}}},{"omega",{"Specific dissipation rate","1/s",{0}}},{"epsilon",{"Turbulent dissipation rate","m^2/s^3",{0}}}};
      const auto it=known.find(key);
      std::optional<ScalarValue> wrapped=v;
      auto& s=b.add(prefix+key,it!=known.end()?std::get<0>(it->second):key,"Boundary Conditions",value(wrapped),it!=known.end()?std::get<1>(it->second):"1","number");
      Json::Array when;
      if(it!=known.end() && two_equation && (key!="omega" || omega_model) && (key!="epsilon" || epsilon_model))for(auto i:std::get<2>(it->second))when.push_back(i);
      if(it!=known.end())b.depends(s,spec_path,when);
      else s.active=false; // e.g. SA nut: its selector mapping is unsupported
      // Inactive values stay editable, so a user can switch method and then
      // set its values; Fluent likewise keeps values the method does not use.
      // A value no method of this model reads (epsilon for k-omega) is not.
      if((inlet || outlet) && !when.empty() && (v.constant || v.resolved))b.editable(s,[target=&v](const Json& j){ScalarValue x;x.kind=ScalarValue::Kind::constant;x.constant=std::get<double>(j.value);*target=x;},key=="k" || key=="turb-intensity"?0:positive,key=="turb-intensity"?1:1e30);
      b.consume(z.raw,key=="omega"?"o":key=="epsilon"?"e":key);
    }
    if(transition && (inlet || outlet)) {
      // Fluent's inlet (outlet: backflow) intermittency, used on inflow faces.
      b.scalar(prefix+"intermittency","Intermittency","Boundary Conditions",z.intermittency,"1",true,0,1);
      b.consume(z.raw,"intermit");
    }
    b.raw(z.raw,prefix+"raw/","Boundary Conditions");
  }
  for(auto [key,label,unit,v,editable]:std::vector<std::tuple<std::string,std::string,std::string,std::optional<double>*,bool>>{
      {"area","Reference area","m^2",&c.reference.area,true},{"length","Reference length","m",&c.reference.length,true},
      {"density","Reference density","kg/m^3",&c.reference.density,true},{"velocity","Reference velocity","m/s",&c.reference.velocity,true},
      {"pressure","Reference pressure","Pa",&c.reference.pressure,true},{"viscosity","Reference viscosity","Pa s",&c.reference.viscosity,false}}) {
    b.number("reference/"+key,label,"Reference Values",*v,unit,editable,key=="pressure"?-1e30:positive,1e30);b.consume(c.variables,"reference-"+key);
  }
  // Retain uncommon nested reference values, without duplicating the six
  // canonical values above.
  if(const auto* r=c.variables.value("reference-values")){RawSettings nested;nested.alist=*r;for(auto key:{"area","length","density","velocity","pressure","viscosity"})b.consume(nested,key);b.raw(nested,"reference/raw/","Reference Values");}
  b.consume(c.variables,"reference-values");
  // Fluent's Methods panel. Each transport equation's scheme is listed for
  // every case; only the equations the case's model solves are active, and
  // only active saved schemes are editable.
  for(auto [path,label,key,target,active]:std::vector<std::tuple<std::string,std::string,std::string,std::optional<Scheme>*,bool>>{
      {"momentum-scheme","Momentum","mom/scheme",&c.numerics.momentum,true},{"k-scheme","Turbulent kinetic energy","k/scheme",&c.numerics.k,two_equation},
      {"omega-scheme","Specific dissipation rate","omega/scheme",&c.numerics.omega,omega_model},{"epsilon-scheme","Turbulent dissipation rate","epsilon/scheme",&c.numerics.epsilon,epsilon_model}}) {
    auto& s=b.add("methods/"+path,label,"Methods",*target?Json(**target==Scheme::first_upwind?"first-order-upwind":**target==Scheme::second_upwind?"second-order-upwind":std::string(name(**target))):Json(nullptr),"","enum");
    s.options={"first-order-upwind","second-order-upwind"};s.option_labels={"First Order Upwind","Second Order Upwind"};if(*target)s.value_label=scheme_label(**target);
    s.active=active;
    if(*target && active)b.editable(s,[target](const Json& j){*target=std::get<std::string>(j.value)=="first-order-upwind"?Scheme::first_upwind:Scheme::second_upwind;});
    else if(!active)s.reason="not solved by the case's viscous model";
    b.consume(c.variables,key);
  }
  {
    const auto& g=c.numerics.gradient;
    auto& s=b.add("methods/gradient","Gradient","Methods",g?Json(*g==Gradient::least_squares_cell?"least-squares-cell-based":*g==Gradient::green_gauss_cell?"green-gauss-cell-based":"green-gauss-node-based"):Json(nullptr),"","enum");
    s.options={"least-squares-cell-based","green-gauss-cell-based","green-gauss-node-based"};s.option_labels={"Least Squares Cell Based","Green-Gauss Cell Based","Green-Gauss Node Based"};label_value(s);
    b.consume(c.variables,"recon/cell-lsq?");b.consume(c.variables,"recon/node-lsq?");
  }
  {
    const auto& p=c.numerics.pressure;
    auto& s=b.add("methods/pressure-scheme","Pressure","Methods",p?Json(*p==Scheme::standard_pressure?"standard":*p==Scheme::second_pressure?"second-order":std::string(name(*p))):Json(nullptr),"","enum");
    s.options={"standard","second-order"};s.option_labels={"Standard","Second Order"};if(p)s.value_label=scheme_label(*p);
    // Fluent's pressure/scheme: the two interpolations the solve adapter runs.
    if(p)b.editable(s,[target=&c.numerics.pressure](const Json& j){
      *target=std::get<std::string>(j.value)=="standard"?Scheme::standard_pressure:Scheme::second_pressure;});
    b.consume(c.variables,"pressure/scheme");
  }
  if(transition) {
    auto& target=c.numerics.intermittency;
    auto& s=b.add("methods/intermittency-scheme","Intermittency","Methods",target?Json(*target==Scheme::first_upwind?"first-order-upwind":*target==Scheme::second_upwind?"second-order-upwind":std::string(name(*target))):Json(nullptr),"","enum");
    s.options={"first-order-upwind","second-order-upwind"};s.option_labels={"First Order Upwind","Second Order Upwind"};if(target)s.value_label=scheme_label(*target);
    if(target)b.editable(s,[target=&target](const Json& j){*target=std::get<std::string>(j.value)=="first-order-upwind"?Scheme::first_upwind:Scheme::second_upwind;});
    b.consume(c.variables,"intermit/scheme");
  }
  b.add("methods/pseudo-transient","Pseudo-transient formulation","Methods",Json(uses_pseudo_time(c)),"","bool");
  b.consume(c.variables,"pseudo-transient-formulation?");
  auto& limit=b.add("methods/gradient-limiting","Gradient limiting","Methods",c.numerics.limiting?Json(*c.numerics.limiting):Json(nullptr),"","bool");
  if(c.numerics.limiting)b.editable(limit,[target=&c.numerics.limiting](const Json& j){*target=std::get<bool>(j.value);});
  b.consume(c.variables,"recon/limit?");
  // Fluent's flow/scheme: SIMPLE, SIMPLEC and Coupled each run natively and may
  // be switched to one another; PISO is rejected by the solve adapter. Kutta's
  // Coupled solver is pseudo-transient, so a case without pseudo time offers
  // Coupled only when it was saved with it (and then runs SIMPLEC).
  {
    const auto& scheme=c.numerics.run_coupling?c.numerics.run_coupling:c.numerics.coupling;
    const auto id=[](Coupling v) {return v==Coupling::coupled?std::string("coupled"):v==Coupling::simplec?std::string("simplec"):v==Coupling::simple?std::string("simple"):std::string("piso");};
    auto& s=b.add("methods/pressure-velocity-coupling","Pressure-velocity coupling","Methods",scheme?Json(id(*scheme)):Json(nullptr),"","enum");
    s.options={"simple","simplec"};s.option_labels={"SIMPLE","SIMPLEC"};
    if(uses_pseudo_time(c) || c.numerics.coupling==Coupling::coupled){s.options.push_back("coupled");s.option_labels.push_back("Coupled");}
    label_value(s);
    if(scheme && s.value_label.empty())s.value_label=std::string(name(*scheme));
    if(scheme)b.editable(s,[target=&c.numerics.run_coupling](const Json& j){
      const auto& v=std::get<std::string>(j.value);*target=v=="coupled"?Coupling::coupled:v=="simple"?Coupling::simple:Coupling::simplec;});
    b.consume(c.variables,"flow/scheme");
  }
  b.number("controls/pseudo-time-scale-factor","Pseudo-time scale factor","Controls",c.numerics.pseudo.scale_factor,"1",true,positive,1e15);b.consume(c.variables,"pseudo-auto-time-step-scale-factor");
  const bool coupled=c.numerics.coupling==Coupling::coupled,pseudo=uses_pseudo_time(c);
  const auto dissipation=c.viscous.model==ViscousModel::ke_standard?"epsilon":"omega";
  std::vector<const char*> relaxed{"mom","pressure","k",dissipation};
  if(transition)relaxed.push_back("intermit");
  for(auto eq:relaxed) {
    const bool flow=std::string_view(eq)=="mom" || std::string_view(eq)=="pressure";
    // The factors the solve adapter uses: pseudo-time factors, classic
    // under-relaxation, or (Coupled without pseudo time) the stored explicit
    // relaxation, which the adapter reads directly and so stays read-only.
    if(!pseudo && coupled && flow) {
      const auto key="pressure-coupled/"+std::string(eq)+"/explicit-relax";const auto* v=c.variables.value(key);
      auto& s=b.add("controls/"+std::string(eq)+"-relaxation",std::string(eq)+" relaxation","Controls",v && v->number()?Json(*v->number()):Json(nullptr),"1","number");
      s.reason="Coupled explicit relaxation is not yet configurable in Kutta";b.consume(c.variables,key);
      continue;
    }
    if(!flow && !two_equation) {
      // Laminar: no turbulence equation is relaxed.
      auto& s=b.add("controls/"+std::string(eq)+"-relaxation",std::string(eq)+" relaxation","Controls",Json(nullptr),"1","number");
      s.active=false;s.reason="not solved by the case's viscous model";
      continue;
    }
    if(!pseudo) {
      const auto it=c.numerics.relaxation.find(eq);auto& s=b.add("controls/"+std::string(eq)+"-relaxation",std::string(eq)+" relaxation","Controls",it!=c.numerics.relaxation.end()?Json(it->second):Json(nullptr),"1","number");
      if(it!=c.numerics.relaxation.end())b.editable(s,[target=&it->second](const Json& j){*target=std::get<double>(j.value);},positive,1);
      b.consume(c.variables,std::string(eq)+"/relax");
      continue;
    }
    const std::string key=coupled && flow?"pressure-coupled/"+std::string(eq)+"/pseudo-explicit-relax":std::string(eq)+"/pseudo-relax";
    const auto it=c.numerics.pseudo_relaxation.find(key);auto& s=b.add("controls/"+std::string(eq)+"-relaxation",std::string(eq)+" relaxation","Controls",it!=c.numerics.pseudo_relaxation.end()?Json(it->second):Json(nullptr),"1","number");
    if(it!=c.numerics.pseudo_relaxation.end())b.editable(s,[target=&it->second](const Json& j){*target=std::get<double>(j.value);},positive,1);
    b.consume(c.variables,key);
  }
  for(auto [key,label,target,max]:std::vector<std::tuple<std::string,std::string,int*,double>>{{"iterations","Number of iterations",&run.iterations,1e9},{"monitor-interval","Monitor interval",&run.monitor_interval,100000}}) {
    auto& s=b.add("run/"+key,label,"Run Calculation",*target,"","integer");b.editable(s,[target](const Json& j){*target=int(std::get<double>(j.value));},1,max);
  }
  // Convergence conditions are run controls: residual criteria come from the
  // case, report conditions use Kutta's defaults (docs/app-protocol.md).
  auto& conv=run.convergence;conv=ConvergenceControls{};
  for(auto eq:solved_residuals(c.viscous.model,transition))conv.residuals.push_back({std::string(eq)});
  if(residual_criteria(c.variables.value("residuals/settings"),conv.residuals))b.consume(c.variables,"residuals/settings");
  for(const auto& r:c.reports)conv.reports.push_back({r.name});
  b.editable(b.add("monitors/check-convergence","Check convergence","Monitors",conv.check),[t=&conv.check](const Json& j){*t=std::get<bool>(j.value);});
  auto& combine=b.add("monitors/condition","Condition","Monitors","all","","enum");
  combine.options={"all","any"};combine.option_labels={"All Conditions are Met","Any Condition is Met"};label_value(combine);
  b.editable(combine,[t=&conv.combine](const Json& j){*t=std::get<std::string>(j.value)=="any"?ConvergenceControls::Combine::any:ConvergenceControls::Combine::all;});
  for(auto& r:conv.residuals) {
    const auto prefix="monitors/residual/"+segment(r.equation)+"/";
    b.editable(b.add(prefix+"check","Check Convergence ("+r.equation+")","Monitors",r.check),[t=&r.check](const Json& j){*t=std::get<bool>(j.value);});
    b.editable(b.add(prefix+"absolute-criterion","Absolute Criteria ("+r.equation+")","Monitors",r.absolute?Json(*r.absolute):Json(nullptr),"1","number"),[t=&r.absolute](const Json& j){*t=std::get<double>(j.value);},0,1e30);
  }
  for(std::size_t i=0;i<c.reports.size();++i) {
    using Kind=ReportDefinition::Kind;auto& r=conv.reports[i];const auto kind=c.reports[i].kind;
    const bool force=kind==Kind::force || kind==Kind::drag || kind==Kind::lift || kind==Kind::downforce;
    const auto prefix="monitors/report/"+segment(r.report)+"/";
    const auto bind=[&](Setting& s,std::function<void(const Json&)> setter,std::optional<double> min={},std::optional<double> max={}) {
      if(force)b.editable(s,std::move(setter),min,max);else s.reason="only force, drag, lift and downforce reports can be convergence conditions in Kutta";
    };
    bind(b.add(prefix+"active","Active ("+r.report+")","Monitors",r.active),[t=&r.active](const Json& j){*t=std::get<bool>(j.value);});
    bind(b.add(prefix+"stop-criterion","Stop Criterion ("+r.report+")","Monitors",r.stop,"1","number"),[t=&r.stop](const Json& j){*t=std::get<double>(j.value);},positive,1e30);
    bind(b.add(prefix+"previous-values","Use Iterations ("+r.report+")","Monitors",r.previous_values,"","integer"),[t=&r.previous_values](const Json& j){*t=int(std::get<double>(j.value));},1,100000);
    bind(b.add(prefix+"ignore-before","Ignore Iterations Before ("+r.report+")","Monitors",r.ignore_before,"","integer"),[t=&r.ignore_before](const Json& j){*t=int(std::get<double>(j.value));},0,1e9);
  }
  const auto* hybrid=c.variables.value("hyb-init?");
  auto& init=b.add("initialization/method","Initialization method","Initialization",seeded?"fluent-data":hybrid && hybrid->kind==SExpr::Kind::boolean?(hybrid->boolean?"hybrid":"standard"):"unparsed","","enum");
  init.options={"standard","hybrid","fluent-data"};init.option_labels={"Standard Initialization","Hybrid Initialization","From Fluent data file"};label_value(init);b.consume(c.variables,"hyb-init?");
  // The method itself stays read-only: changing it would require values which
  // a saved case may not contain. Standard numeric defaults can be overridden.
  if(!seeded && hybrid && hybrid->kind==SExpr::Kind::boolean && !hybrid->boolean)for(auto [key,label,unit,min]:std::vector<std::tuple<std::string,std::string,std::string,double>>{
      {"x-velocity/default","Initial x velocity","m/s",-1e30},{"y-velocity/default","Initial y velocity","m/s",-1e30},{"z-velocity/default","Initial z velocity","m/s",-1e30},
      {"pressure/default","Initial gauge pressure","Pa",-1e30},{"k/default","Initial turbulent kinetic energy","m^2/s^2",0},{"omega/default","Initial specific dissipation","1/s",positive}}) {
    const auto* entry=c.variables.value(key);auto& setting=b.add("initialization/"+key.substr(0,key.find('/')),label,"Initialization",value(entry),unit,"number");
    if(entry && entry->number()) {
      // Locate the actual data-bearing entry rather than registration metadata.
      b.editable(setting,[raw=&c.variables,key](const Json& j){
        const auto* selected=raw->entry(key);
        for(auto& e:raw->alist.elements)if(&e==selected){
          SExpr v;v.kind=SExpr::Kind::real;v.real=std::get<double>(j.value);
          e.elements={SExpr::symbol(key),std::move(v)};e.tail.reset();return;
        }
        throw Error("initialization setting disappeared: "+key);
      },min,1e30);
    }
    b.consume(c.variables,key);
  }
  b.consume(c.variables,"monitor/report-definitions");
  for(auto& r:c.reports) {
    const auto prefix="reports/"+segment(r.name)+"/";b.add(prefix+"type","Report type","Report Definitions",r.fluent_type,"","enum").value_label=r.fluent_type;
    auto& direction=b.add(prefix+"direction","Force direction","Report Definitions",r.direction?Json(Json::Array{(*r.direction)[0],(*r.direction)[1],(*r.direction)[2]}):Json(nullptr),"1","vec3");
    if(r.direction && (r.kind==ReportDefinition::Kind::force || r.kind==ReportDefinition::Kind::drag || r.kind==ReportDefinition::Kind::lift || r.kind==ReportDefinition::Kind::downforce)) {
      b.editable(direction,[target=&r.direction](const Json& j){const auto& a=std::get<Json::Array>(j.value);*target=std::array{std::get<double>(a[0].value),std::get<double>(a[1].value),std::get<double>(a[2].value)};},-1e15,1e15);
      direction.validate=[](const Json& j){const auto& a=std::get<Json::Array>(j.value);return std::hypot(std::get<double>(a[0].value),std::get<double>(a[1].value),std::get<double>(a[2].value))>0?std::string{}:"direction must be nonzero";};
    }
    b.raw(r.raw,prefix+"raw/","Report Definitions");
  }
  for(const auto& [id,raw]:c.domains)b.raw(raw,"cell-zone/domain-"+std::to_string(id)+"/","Cell Zone Conditions");
  b.raw(c.variables,"case/");
  for(const auto& [name,forms]:c.raw_sections)if(name!="Rampant Variables" && name!="Thread Variables" && name!="Domain Variables")for(std::size_t i=0;i<forms.size();++i)b.add("general/source/"+segment(name)+"/"+std::to_string(i),name,"General",to_string(forms[i]));
  return b.finish();
}
Json settings_json(std::span<const Setting> settings) {
  Json::Array out;out.reserve(settings.size());
  for(const auto& s:settings){Json::Object row{{"path",s.path},{"label",s.label},{"group",s.group},{"value",s.value},{"unit",s.unit},{"kind",s.kind},{"options",s.options},{"option_labels",s.option_labels},{"value_label",s.value_label},{"active",s.active},{"depends_on",s.depends_on?Json(*s.depends_on):Json(nullptr)},{"active_when",s.active_when},{"min",s.min?Json(*s.min):Json(nullptr)},{"max",s.max?Json(*s.max):Json(nullptr)},{"editable",s.editable},{"reason",s.reason}};out.emplace_back(std::move(row));}
  return out;
}
Json::Object apply_setup(const std::filesystem::path& path,std::span<const Setting> settings) {
  KUTTA_REQUIRE(std::filesystem::file_size(path)<=1024*1024,"setup JSON exceeds 1 MiB");std::ifstream in(path,std::ios::binary);KUTTA_REQUIRE(in,"cannot open setup JSON");
  const std::string text{std::istreambuf_iterator<char>(in),{}};Json parsed;
  try {parsed=parse_json(text);}catch(const std::exception& e){throw RejectedSettings({{"setup",e.what()}});}
  const auto* root=std::get_if<Json::Object>(&parsed.value);
  std::vector<SettingIssue> errors;
  const auto fail=[&](std::string path,std::string reason){errors.push_back({std::move(path),std::move(reason)});};
  if(!root)throw RejectedSettings({{"setup","expected object with v and overrides"}});
  const auto version=root->find("v"),overrides=root->find("overrides");
  if(version==root->end() || version->second!=Json(1))fail("setup/v","expected protocol version 1");
  for(const auto& [key,v]:*root){(void)v;if(key!="v" && key!="overrides")fail("setup/"+key,"unknown setup property");}
  const auto* object=overrides==root->end()?nullptr:std::get_if<Json::Object>(&overrides->second.value);
  if(!object)fail("setup/overrides","expected an object keyed by setting path");
  if(object)for(const auto& [path,value]:*object) {
    auto it=std::find_if(settings.begin(),settings.end(),[&](const auto& s){return s.path==path;});
    if(it==settings.end()){fail(path,"unknown setting path");continue;}
    const auto& s=*it;if(!s.editable){fail(path,s.reason);continue;}
    std::vector<double> numbers;bool valid=false;
    if(s.kind=="number" || s.kind=="integer"){if(const auto* n=std::get_if<double>(&value.value)){valid=std::isfinite(*n) && (s.kind!="integer" || std::floor(*n)==*n);numbers.push_back(*n);}}
    else if(s.kind=="bool")valid=std::holds_alternative<bool>(value.value);
    else if(s.kind=="vec3"){if(const auto* a=std::get_if<Json::Array>(&value.value);a && a->size()==3){valid=true;for(const auto& x:*a){if(const auto* n=std::get_if<double>(&x.value);n && std::isfinite(*n))numbers.push_back(*n);else valid=false;}}}
    else if(s.kind=="enum")valid=std::holds_alternative<std::string>(value.value) || std::holds_alternative<double>(value.value);
    else valid=std::holds_alternative<std::string>(value.value);
    if(!valid){fail(path,"expected "+s.kind);continue;}
    if(!s.options.empty() && std::find(s.options.begin(),s.options.end(),value)==s.options.end()){fail(path,"value is not one of the supported options");continue;}
    if(std::any_of(numbers.begin(),numbers.end(),[&](double n){return (s.min && n<*s.min) || (s.max && n>*s.max);})){fail(path,"value is outside the advertised range");continue;}
    if(s.validate){const auto reason=s.validate(value);if(!reason.empty())fail(path,reason);}
  }
  if(!errors.empty())throw RejectedSettings(std::move(errors));
  for(const auto& [path,v]:*object){const auto it=std::find_if(settings.begin(),settings.end(),[&](const auto& s){return s.path==path;});it->apply(v);}
  return *object;
}
} // namespace kutta::setup
