#include "settings.hpp"
#include "kutta/mesh/detail.hpp"
#include <algorithm>
#include <cmath>
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
class Inventory {
  std::vector<Setting> out;
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
  void scalar(std::string path,std::string label,std::string group,std::optional<ScalarValue>& v,std::string unit,bool active,std::optional<double> min={},std::optional<double> max={}) {
    auto& s=add(std::move(path),std::move(label),std::move(group),value(v),std::move(unit),"number");
    if(active && v && (v->constant || v->resolved))editable(s,[target=&v](const Json& j){ScalarValue v;v.kind=ScalarValue::Kind::constant;v.constant=std::get<double>(j.value);*target=std::move(v);},min,max);
  }
  void number(std::string path,std::string label,std::string group,std::optional<double>& v,std::string unit,bool active,std::optional<double> min={},std::optional<double> max={}) {
    auto& s=add(std::move(path),std::move(label),std::move(group),v?Json(*v):Json(nullptr),std::move(unit),"number");
    if(active && v)editable(s,[target=&v](const Json& j){*target=std::get<double>(j.value);},min,max);
  }
  void vector(std::string path,std::string label,ScalarVector& v,std::string unit,bool active,bool nonzero=false) {
    auto& s=add(std::move(path),std::move(label),"Boundary Conditions",vector_value(v),std::move(unit),"vec3");
    if(active && std::all_of(v.begin(),v.end(),[](const auto& x){return x && (x->constant || x->resolved);})) {
      editable(s,[target=&v](const Json& j){const auto& a=std::get<Json::Array>(j.value);for(int i=0;i<3;++i){ScalarValue x;x.kind=ScalarValue::Kind::constant;x.constant=std::get<double>(a[i].value);(*target)[i]=x;}},-1e30,1e30);
      if(nonzero)s.validate=[](const Json& j){const auto& a=std::get<Json::Array>(j.value);const double n=std::hypot(std::get<double>(a[0].value),std::get<double>(a[1].value),std::get<double>(a[2].value));return n>0?std::string{}:"direction must be nonzero";};
    }
  }
  std::vector<Setting> finish(){std::sort(out.begin(),out.end(),[](const auto& a,const auto& b){return a.path<b.path;});for(std::size_t i=1;i<out.size();++i)KUTTA_REQUIRE(out[i-1].path!=out[i].path,"duplicate settings identity: "+out[i].path);return std::move(out);}
};
}
std::vector<Setting> settings_inventory(FluentCase& c,RunControls& run,bool seeded) {
  Inventory b;
  constexpr double positive=double(std::numeric_limits<float>::min());
  auto& general=b.add("general/operating-pressure","Operating pressure","General",c.operating_pressure?Json(*c.operating_pressure):Json(nullptr),"Pa","number");
  general.reason="not yet configurable in Kutta";b.consume(c.variables,"operating-pressure");
  b.add("general/dimension","Dimension","General",c.dimension?Json(*c.dimension):Json(nullptr),"","integer");
  b.add("models/viscous","Viscous model","Models",c.viscous.model?Json(std::string(name(*c.viscous.model))):Json(nullptr),"","enum");
  b.add("models/energy","Energy equation","Models",c.energy?Json(*c.energy):Json(nullptr),"","bool");
  b.consume(c.variables,"case-config");b.raw(c.config,"general/case-config/");
  b.consume(c.variables,"materials");
  for(auto& m:c.materials) {
    const auto prefix="materials/"+segment(m.name)+"/";
    const bool active=std::any_of(c.cell_zones.begin(),c.cell_zones.end(),[&](const auto& z){return z.material==m.name && z.type==ZoneType::fluid;});
    b.scalar(prefix+"density","Density","Materials",m.density,"kg/m^3",active && m.nonconstant_density==false,positive,1e30);
    b.scalar(prefix+"viscosity","Dynamic viscosity","Materials",m.viscosity,"Pa s",active && m.nonconstant_viscosity==false,positive,1e30);
    b.consume(m.raw,"density");b.consume(m.raw,"viscosity");b.raw(m.raw,prefix+"raw/","Materials");
    b.add(prefix+"type","Material type","Materials",m.type,"","enum");
  }
  for(auto& z:c.cell_zones){const auto prefix="cell-zone/"+segment(z.name)+"/";b.add(prefix+"type","Zone type","Cell Zone Conditions",z.fluent_type,"","enum");b.raw(z.raw,prefix+"raw/","Cell Zone Conditions");}
  for(auto& z:c.boundaries) {
    const auto prefix="boundary/"+segment(z.name)+"/";
    b.add(prefix+"type","Boundary type","Boundary Conditions",z.fluent_type,"","enum");
    const bool inlet=z.type==ZoneType::velocity_inlet,outlet=z.type==ZoneType::pressure_outlet;
    b.scalar(prefix+"velocity-magnitude","Velocity magnitude","Boundary Conditions",z.velocity_magnitude,"m/s",inlet && (z.velocity_spec==VelocitySpec::normal || z.velocity_spec==VelocitySpec::magnitude_direction),0,1e15);
    b.vector(prefix+"velocity","Velocity components",z.velocity,"m/s",inlet && z.velocity_spec==VelocitySpec::components);
    b.vector(prefix+"direction","Flow direction",z.direction,"1",inlet && z.velocity_spec==VelocitySpec::magnitude_direction,true);
    b.scalar(prefix+"gauge-pressure","Gauge pressure","Boundary Conditions",z.gauge_pressure,"Pa",outlet,-1e30,1e30);
    for(auto key:{"vmag","p","u","v","w","ni","nj","nk","direction-0","direction-1","direction-2","velocity-vector","flow-direction-component","reference-direction-components","flow-direction"})b.consume(z.raw,key);
    const auto spec=z.turbulence_spec;
    b.add(prefix+"turbulence-specification","Turbulence specification","Boundary Conditions",spec?Json(int(*spec)):Json(nullptr),"","enum");b.consume(z.raw,"ke-spec");
    for(auto& [key,v]:z.turbulence) {
      const bool active=(inlet || outlet) && ((spec==TurbulenceSpec::k_dissipation && (key=="k" || key=="omega")) ||
        (spec && spec!=TurbulenceSpec::k_dissipation && (key=="turb-intensity" || (key=="turb-viscosity-ratio" && spec==TurbulenceSpec::intensity_viscosity_ratio) ||
        (key=="turb-length-scale" && spec==TurbulenceSpec::intensity_length) || (key=="turb-hydraulic-diam" && spec==TurbulenceSpec::intensity_diameter))));
      std::optional<ScalarValue> wrapped=v;
      const auto unit=key=="k"?"m^2/s^2":key=="omega"?"1/s":key=="turb-length-scale" || key=="turb-hydraulic-diam"?"m":"1";
      auto& s=b.add(prefix+key,key,"Boundary Conditions",value(wrapped),unit,"number");
      if(active && (v.constant || v.resolved))b.editable(s,[target=&v](const Json& j){ScalarValue x;x.kind=ScalarValue::Kind::constant;x.constant=std::get<double>(j.value);*target=x;},key=="k" || key=="turb-intensity"?0:positive,key=="turb-intensity"?1:1e30);
      b.consume(z.raw,key=="omega"?"o":key=="epsilon"?"e":key);
    }
    b.raw(z.raw,prefix+"raw/","Boundary Conditions");
  }
  for(auto [key,label,unit,v,editable]:std::vector<std::tuple<std::string,std::string,std::string,std::optional<double>*,bool>>{
      {"area","Reference area","m^2",&c.reference.area,true},{"length","Reference length","m",&c.reference.length,false},
      {"density","Reference density","kg/m^3",&c.reference.density,true},{"velocity","Reference velocity","m/s",&c.reference.velocity,true},
      {"pressure","Reference pressure","Pa",&c.reference.pressure,true},{"viscosity","Reference viscosity","Pa s",&c.reference.viscosity,false}}) {
    b.number("reference/"+key,label,"Reference Values",*v,unit,editable,key=="pressure"?-1e30:positive,1e30);b.consume(c.variables,"reference-"+key);
  }
  // Retain uncommon nested reference values, without duplicating the six
  // canonical values above.
  if(const auto* r=c.variables.value("reference-values")){RawSettings nested;nested.alist=*r;for(auto key:{"area","length","density","velocity","pressure","viscosity"})b.consume(nested,key);b.raw(nested,"reference/raw/","Reference Values");}
  b.consume(c.variables,"reference-values");
  for(auto [path,label,key,target]:std::vector<std::tuple<std::string,std::string,std::string,std::optional<Scheme>*>>{
      {"momentum-scheme","Momentum","mom/scheme",&c.numerics.momentum},{"k-scheme","Turbulent kinetic energy","k/scheme",&c.numerics.k},{"omega-scheme","Specific dissipation rate","omega/scheme",&c.numerics.omega}}) {
    auto& s=b.add("methods/"+path,label,"Methods",*target?Json(**target==Scheme::first_upwind?"first-order-upwind":**target==Scheme::second_upwind?"second-order-upwind":std::string(name(**target))):Json(nullptr),"","enum");
    s.options={"first-order-upwind","second-order-upwind"};
    if(*target)b.editable(s,[target](const Json& j){*target=std::get<std::string>(j.value)=="first-order-upwind"?Scheme::first_upwind:Scheme::second_upwind;});
    b.consume(c.variables,key);
  }
  auto& limit=b.add("methods/gradient-limiting","Gradient limiting","Methods",c.numerics.limiting?Json(*c.numerics.limiting):Json(nullptr),"","bool");
  if(c.numerics.limiting)b.editable(limit,[target=&c.numerics.limiting](const Json& j){*target=std::get<bool>(j.value);});
  b.consume(c.variables,"recon/limit?");
  b.number("controls/pseudo-time-scale-factor","Pseudo-time scale factor","Controls",c.numerics.pseudo.scale_factor,"1",true,positive,1e15);b.consume(c.variables,"pseudo-auto-time-step-scale-factor");
  const bool coupled=c.numerics.coupling==Coupling::coupled;
  for(auto eq:{"mom","pressure","k","omega"}) {
    const std::string key=coupled && (std::string_view(eq)=="mom" || std::string_view(eq)=="pressure")?"pressure-coupled/"+std::string(eq)+"/pseudo-explicit-relax":std::string(eq)+"/pseudo-relax";
    const auto it=c.numerics.pseudo_relaxation.find(key);auto& s=b.add("controls/"+std::string(eq)+"-relaxation",std::string(eq)+" relaxation","Controls",it!=c.numerics.pseudo_relaxation.end()?Json(it->second):Json(nullptr),"1","number");
    if(it!=c.numerics.pseudo_relaxation.end())b.editable(s,[target=&it->second](const Json& j){*target=std::get<double>(j.value);},positive,1);
    b.consume(c.variables,key);
  }
  for(auto [key,label,target,max]:std::vector<std::tuple<std::string,std::string,int*,double>>{{"iterations","Number of iterations",&run.iterations,1e9},{"monitor-interval","Monitor interval",&run.monitor_interval,100000}}) {
    auto& s=b.add("run/"+key,label,"Run Calculation",*target,"","integer");b.editable(s,[target](const Json& j){*target=int(std::get<double>(j.value));},1,max);
  }
  b.add("monitors/convergence-criteria","Convergence criteria","Monitors",nullptr,"","string");
  const auto* hybrid=c.variables.value("hyb-init?");
  auto& init=b.add("initialization/method","Initialization method","Initialization",seeded?"fluent-data":hybrid && hybrid->kind==SExpr::Kind::boolean?(hybrid->boolean?"hybrid":"standard"):"unparsed","","enum");
  init.options={"standard","hybrid","fluent-data"};b.consume(c.variables,"hyb-init?");
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
    const auto prefix="reports/"+segment(r.name)+"/";b.add(prefix+"type","Report type","Report Definitions",r.fluent_type,"","enum");
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
  for(const auto& s:settings){Json::Object row{{"path",s.path},{"label",s.label},{"group",s.group},{"value",s.value},{"unit",s.unit},{"kind",s.kind},{"options",s.options},{"min",s.min?Json(*s.min):Json(nullptr)},{"max",s.max?Json(*s.max):Json(nullptr)},{"editable",s.editable},{"reason",s.reason}};out.emplace_back(std::move(row));}
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
