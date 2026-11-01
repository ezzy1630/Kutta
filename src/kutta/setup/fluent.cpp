#include "kutta/setup/fluent.hpp"
#include "kutta/setup/fluent_physics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <type_traits>

#include "kutta/mesh/detail.hpp"

namespace kutta::setup {
using namespace fluent;
using namespace fv;
using namespace kutta::detail;
namespace {
void require(bool ok,const std::string& key) {KUTTA_REQUIRE(ok,"Fluent setup: unsupported or absent "+key);}
// Accumulate independent setting failures before returning any setup. Placeholder
// values below exist only during validation and never reach a solver.
class SetupValidation {
  std::vector<SettingIssue> errors;
  bool epsilon_model=false,transition=false;
 public:
  void require(bool ok,const std::string& key,const std::string& reason="value is unsupported or required setting is absent") {
    if(!ok) errors.push_back({key,reason});
  }
  template<class T> T required(const std::optional<T>& value,const std::string& key) {
    require(value.has_value(),key,"required setting absent or unparsed");return value.value_or(T{});
  }
  void finish() const {
    if(errors.empty()) return;
    throw RejectedSettings(errors);
  }
double scalar(const std::optional<ScalarValue>& value,const std::string& key) {
  if(!value) {require(false,key,"required scalar absent or unparsed");return 0;}
  const auto& v=*value;const auto n=v.constant?v.constant:v.resolved;
  require(n && std::isfinite(*n),key,"requires a finite constant or resolved expression");return n.value_or(0);
}
double number(const RawSettings& raw,const std::string& key) {
  const auto* value=raw.value(key);require(value && value->number() && std::isfinite(*value->number()),key,"requires a finite number");
  return value && value->number()?*value->number():0;
}
bool flag(const RawSettings& raw,const std::string& key) {
  const auto* value=raw.value(key);require(value && value->kind==SExpr::Kind::boolean,key,"requires an explicit boolean");return value && value->kind==SExpr::Kind::boolean?value->boolean:false;
}
Vec3d vector(const ScalarVector& value,const std::string& key) {
  return {scalar(value[0],key+"/x"),scalar(value[1],key+"/y"),scalar(value[2],key+"/z")};
}
ConvectionScheme convection(std::optional<Scheme> scheme,const std::string& key) {
  require(scheme==Scheme::first_upwind || scheme==Scheme::second_upwind,key,"only first-order and second-order upwind are supported");
  return scheme==Scheme::first_upwind?ConvectionScheme::FirstOrderUpwind:ConvectionScheme::SecondOrderUpwind;
}
TurbulenceCondition turbulence(const BoundaryZone& z) {
  const auto key="zone/"+std::to_string(z.id)+"/";
  const auto get=[&](const std::string& name) {
    auto it=z.turbulence.find(name);require(it!=z.turbulence.end(),key+name);
    return it==z.turbulence.end()?0:scalar(it->second,key+name);
  };
  switch(required(z.turbulence_spec,key+"ke-spec")) {
    case TurbulenceSpec::k_dissipation:
      if(epsilon_model) return KeValues{get("k"),get("epsilon")};
      return TurbulenceValues{get("k"),get("omega")};
    case TurbulenceSpec::intensity_length:return IntensityLength{get("turb-intensity"),get("turb-length-scale")};
    case TurbulenceSpec::intensity_diameter:return IntensityLength{get("turb-intensity"),.07*get("turb-hydraulic-diam")};
    case TurbulenceSpec::intensity_viscosity_ratio:return IntensityViscosityRatio{get("turb-intensity"),get("turb-viscosity-ratio")};
  }
  require(false,key+"ke-spec");return TurbulenceValues{};
}
// Inflow gamma of the intermittency model; unused (1) without it.
double intermittency(const BoundaryZone& z) {
  if(!transition) return 1;
  const auto key="zone/"+std::to_string(z.id)+"/intermit";
  const auto value=scalar(z.intermittency,key);
  require(value>=0 && value<=1,key,"intermittency must lie in [0, 1]");return value;
}
FluentSetup build(const FluentCase& c,const PolyMesh& mesh,std::span<const Approximation> allowed);
};
bool output_only(const ReportDefinition& r) {
  return r.kind==ReportDefinition::Kind::moment || r.kind==ReportDefinition::Kind::other;
}
// Issues of a report Kutta does not evaluate. Such reports are outputs only:
// Kutta reads no report convergence condition from the case, and only force
// reports can be activated as one (setup/settings.cpp).
bool unevaluated_report(const FluentCase& c,const std::string& path) {
  return std::any_of(c.reports.begin(),c.reports.end(),[&](const auto& r) {
    return output_only(r) && (path=="report/"+r.name || path.starts_with("report/"+r.name+"/"));
  });
}
int version_major(const std::string& version) {
  int major=0;for(char ch:version) {if(ch<'0' || ch>'9') break;major=10*major+(ch-'0');}
  return major;
}
// A 2D Fluent case is admitted only as read_fluent_msh's one-cell extrusion:
// exactly the two zone-less symmetry caps front/back, one face per cell each,
// and a thickness equal to Fluent's reference depth, so wall forces are the
// forces Fluent reports for that depth.
bool planar_extrusion(const FluentCase& c,const PolyMesh& mesh) {
  if(c.mesh.dimension!=2) return false;
  std::vector<const Patch*> caps;
  for(const auto& p:mesh.patches) if(p.source_zone_id<0) caps.push_back(&p);
  if(caps.size()!=2) return false;
  for(const auto* p:caps)
    if(p->kind!=PatchKind::Symmetry || p->size!=mesh.n_cells || (p->name!="front" && p->name!="back") ||
       std::any_of(c.boundaries.begin(),c.boundaries.end(),[&](const auto& z){return z.name==p->name;})) return false;
  const auto* depth=c.variables.value("reference-depth");
  if(!depth || !depth->number() || !(*depth->number()>0)) return false;
  double lo=std::numeric_limits<double>::infinity(),hi=-lo;
  for(const auto& point:mesh.points) {lo=std::min(lo,point.z);hi=std::max(hi,point.z);}
  return std::abs((hi-lo)-*depth->number())<=1e-12*std::max(1.,*depth->number());
}
std::vector<std::uint64_t> cell_ids(const FluentCase& c,const PolyMesh& mesh) {
  std::vector<std::uint64_t> ids(mesh.n_cells,0);
  for(const auto& zone:mesh.cell_zones) {
    const auto it=std::find_if(c.mesh.cell_zones.begin(),c.mesh.cell_zones.end(),[&](const auto& z){return z.zone_id==zone.source_zone_id;});
    require(it!=c.mesh.cell_zones.end() && it->size()==zone.cells.size(),"mesh/cell-zone identity");
    for(std::size_t i=0;i<zone.cells.size();++i) {
      const auto id=zone.cells[i];require(id>=0 && id<mesh.n_cells && ids[id]==0,"mesh/cell-zone mapping");
      ids[id]=it->first+i;
    }
  }
  require(std::none_of(ids.begin(),ids.end(),[](auto id){return id==0;}),"mesh/cell-zone coverage");
  return ids;
}
}
RejectedSettings::RejectedSettings(std::vector<SettingIssue> issues)
  : Error([&] {
      std::string message="Fluent setup rejected:";
      for(const auto& i:issues) message+="\n  "+i.path+": "+i.reason;
      return message;
    }()), issues_(std::move(issues)) {}
TurbulenceValues turbulence_values(const TurbulenceCondition& spec,double speed,double rho,double mu) {
  require(std::isfinite(speed) && speed>=0 && rho>0 && mu>0,"turbulence inputs");
  TurbulenceValues out;
  if(const auto* v=std::get_if<TurbulenceValues>(&spec)) out=*v;
  else if(const auto* v=std::get_if<IntensityLength>(&spec)) {
    require(v->intensity>=0 && v->length>0,"intensity/length");
    out.k=1.5*std::pow(speed*v->intensity,2);out.omega=std::sqrt(out.k)/(std::pow(.09,.25)*v->length);
  } else {
    const auto& vr=std::get<IntensityViscosityRatio>(spec);
    require(vr.intensity>=0 && vr.ratio>0,"intensity/viscosity-ratio");
    out.k=1.5*std::pow(speed*vr.intensity,2);out.omega=rho*out.k/(mu*vr.ratio);
  }
  require(std::isfinite(out.k) && out.k>=0 && std::isfinite(out.omega) && out.omega>=0,"k/omega");return out;
}
FluentSetup SetupValidation::build(const FluentCase& c,const PolyMesh& mesh,std::span<const Approximation> allowed) {
  const bool planar=c.dimension==2 && planar_extrusion(c,mesh);
  require(c.dimension==3 || planar,"dimension","requires a 3D case, or a 2D case read as a one-cell extrusion of the reference depth with front/back symmetry");
  require(c.pressure_based==true,"rp-seg?","requires the pressure-based solver");require(c.unsteady==false,"rp-unsteady?","only steady cases are supported");
  const bool laminar=c.viscous.model==ViscousModel::laminar,sst=c.viscous.model==ViscousModel::kw_sst;
  epsilon_model=c.viscous.model==ViscousModel::ke_standard;
  require(laminar || sst || epsilon_model,"viscous-model","requires laminar, k-omega SST or standard k-epsilon");
  require(c.mesh.cells && *c.mesh.cells==std::uint64_t(mesh.n_cells),"mesh/cell-count");
  FluentSetup out;auto& s=out.settings;const auto& n=c.numerics;const auto& t=c.viscous;
  try {s.cell_zone_physics=fluent_zone_physics(c);} catch(const Error& e) {errors.push_back({"cell-zone-physics",e.what()});}
  require(!c.cell_zones.empty(),"cell zones");
  std::optional<double> rho,mu;
  for(const auto& z:c.cell_zones) {
    const auto prefix="zone/"+std::to_string(z.id)+"/";
    require(z.type==ZoneType::fluid,prefix+"type");require(laminar || z.laminar==false,prefix+"laminar?");
    const auto it=std::find_if(c.materials.begin(),c.materials.end(),[&](const auto& m){return z.material==m.name;});
    require(it!=c.materials.end(),prefix+"material");
    if(it==c.materials.end()) continue;
    const auto r=scalar(it->density,"material/"+it->name+"/density");
    const auto v=scalar(it->viscosity,"material/"+it->name+"/viscosity");
    require(r>0 && v>0 && (!rho || *rho==r) && (!mu || *mu==v),"constant uniform density/viscosity");rho=r;mu=v;
  }
  s.density=rho.value_or(0);s.viscosity=mu.value_or(0);
  require(c.gravity==false,"gravity? (body force not adapted)");
  if(c.energy==true) out.convergence_notes.push_back("energy is decoupled: active material density and viscosity are constant");
  if(n.gradient==Gradient::green_gauss_cell) {
    s.gradient=GradientScheme::GreenGaussCell;
    require(!sst,"recon/cell-lsq?","Green-Gauss cell gradients are adapted for laminar and standard k-epsilon flow only");
  } else require(n.gradient==Gradient::least_squares_cell,"recon/cell-lsq?","requires least-squares or Green-Gauss cell-based gradients");
  if(n.pressure==Scheme::standard_pressure) s.pressure_interpolation=PressureInterpolation::Standard;
  else require(n.pressure==Scheme::second_pressure,"pressure/scheme","requires second-order or standard pressure interpolation");
  s.convection=convection(n.momentum,"mom/scheme");
  require(n.limiter==Limiter::standard,"recon/lim/type","only the standard limiter (0) is supported");
  s.limit=required(n.limiting,"recon/limit?");
  require(n.limit_to_face==true,"recon/lim/dx-to-face?","requires face-based limiting");require(n.limiter_filter==false,"recon/lim/filter/on?","limiter filtering is not implemented");
  const auto& wfgc=approximation_info(Approximation::wfgc);
  if(n.warped_face_correction==true) {
    if(std::find(allowed.begin(),allowed.end(),wfgc.id)!=allowed.end()) out.approximations.push_back(wfgc.id);
    else require(false,std::string(wfgc.setting),"warped-face gradient correction is not implemented; opt in with --approximate "+
      std::string(wfgc.name)+" to "+std::string(wfgc.kutta_action)+"; effect: "+std::string(wfgc.effect_documentation));
  } else require(n.warped_face_correction==false,std::string(wfgc.setting),"requires an explicit boolean; approximation cannot replace absent or unparsed settings");
  // Rhie-Chow face velocity. Cases written before the flux-type selector
  // existed store neither key; the Fluent 12.0 Theory Guide, section 18.4.1,
  // documents only momentum-weighted (a_P) averaging for them.
  std::set<std::string> handled{"rp/turbulence/wall-omega-treatment","rp/kw-dissipation-newton?",
    "rp/kw-wall-omega-treatment-r13?","rp/sst-des-iddes-simplify_omit_fe?","rp/sst-linearize-neg-cross-diff?"};
  if(!c.variables.entry("pbs/flux-index") && !c.variables.entry("pbs/flux-auto-select?")) {
    handled.insert({"rp/pbs/flux-index","rp/pbs/flux-auto-select?"});
    out.convergence_notes.push_back("pbs/flux-index: absent (case predates the flux-type selector); momentum-weighted Rhie-Chow, Fluent 12.0 TG 18.4.1");
  } else {
    const auto automatic=required(n.flux_auto_select,"pbs/flux-auto-select?");
    // Auto-select picks momentum-based for steady incompressible flow; a stored
    // distance index under auto-select would be contradictory evidence.
    require(n.flux_index==1 || (n.flux_index==0 && !automatic),"pbs/flux-index",
      "requires momentum-based (1), or distance-based (0) with auto-select off");
    if(n.flux_index==0) s.face_velocity=FaceVelocityWeights::Distance;
  }
  const auto coupling=required(n.coupling,"flow/scheme");
  require(coupling!=Coupling::piso,"flow/scheme","PISO is not implemented");
  // An override changes the algorithm only: the saved scheme's pseudo-time
  // controls and relaxation factors stay in force, so the fixed point is the
  // same and only the convergence path changes.
  const auto run=n.run_coupling.value_or(coupling);
  if(run!=coupling) out.convergence_notes.push_back("flow/scheme: "+std::string(name(run))+
    " selected by setup override for a "+std::string(name(coupling))+" case; its pseudo-time controls and relaxation factors apply; convergence path only");
  if(run==Coupling::simple) s.coupling=fv::PressureVelocityCoupling::Simple;
  if(required(n.high_order_relaxation,"recon/relax/relax?")) out.convergence_notes.push_back("high-order term relaxation omitted; same fixed point");
  const auto method=coupling==Coupling::coupled?n.pseudo.coupled_method:n.pseudo.segregated_method;
  if(method) require(*method==0 || *method==1,"pseudo-time-method/dt-method","only off (0) and global time step (1) are implemented");
  else handled.insert({"rp/pseudo-time-method/coupled-pbns/dt-method","rp/pseudo-time-method/segregated-pbns/dt-method"});
  const bool pseudo=uses_pseudo_time(c);
  if(!n.pseudo.formulation) required(method,"pseudo-time-method/dt-method");
  else if(method && (*method==1)!=pseudo)
    out.convergence_notes.push_back("pseudo-time-method/dt-method disagrees with pseudo-transient-formulation?; the formulation flag governs");
  if(pseudo && !method) out.convergence_notes.push_back("pseudo-time-method: absent (case predates the selector); pseudo-transient-formulation? #t runs global automatic pseudo time");
  s.pseudo_transient=pseudo;
  // Kutta's Coupled solver is pseudo-transient. A saved Coupled case without
  // pseudo time falls back to SIMPLEC below; an override cannot ask for it.
  if(run==Coupling::coupled && pseudo) {
    s.coupling=fv::PressureVelocityCoupling::Coupled;
    // The saved Flow Courant Number sets the Coupled implicit relaxation;
    // absent, Fluent's default 200 stays (fv::SolverSettings::coupled_courant).
    if(n.courant) {
      require(std::isfinite(*n.courant) && *n.courant>0,"pressure-coupled/courant-number","requires a finite positive Courant number");
      if(std::isfinite(*n.courant) && *n.courant>0) s.coupled_courant=*n.courant;
    }
  }
  else if(run==Coupling::coupled && run!=coupling)
    require(false,"flow/scheme","the Coupled override requires the pseudo-transient formulation");
  if(pseudo) {
    require(n.pseudo.mode==0,"pseudo-transient-mode");
    require(n.pseudo.length_scale_method==1,"pseudo/autotime-lscale-option");
    s.pseudo_time_scale=required(n.pseudo.scale_factor,"pseudo-auto-time-step-scale-factor");
    if(!required(n.pseudo.automatic,"pseudo-auto-time-step?")) {
      s.pseudo_time_step=required(n.pseudo.time_step,"pseudo-time-step");
      require(s.pseudo_time_step>0,"pseudo-time-step","an explicit step must be positive; zero would select automatic stepping");
    }
  } else {
    handled.insert("rp/pseudo-transient-mode");
    if(coupling==Coupling::coupled) out.convergence_notes.push_back("Coupled without pseudo time: its explicit momentum and pressure relaxation become SIMPLEC under-relaxation; the Courant-number implicit relaxation is not reproduced; convergence path only");
  }
  // Pseudo-transient factors (TG 23-65), or classic implicit under-relaxation
  // (TG 23-66). Coupled without pseudo time stores explicit relaxation for the
  // coupled pressure and momentum system.
  const auto relax=[&](const std::string& eq) {
    double value=0;std::string key;
    if(pseudo) {
      key=coupling==Coupling::coupled && (eq=="mom" || eq=="pressure")?
        "pressure-coupled/"+eq+"/pseudo-explicit-relax":eq+"/pseudo-relax";
      const auto it=n.pseudo_relaxation.find(key);require(it!=n.pseudo_relaxation.end(),key);
      if(it==n.pseudo_relaxation.end()) return 0.0;
      value=it->second;
    } else if(coupling==Coupling::coupled && (eq=="mom" || eq=="pressure")) {
      key="pressure-coupled/"+eq+"/explicit-relax";value=number(c.variables,key);
    } else {
      key=eq+"/relax";
      const auto it=n.relaxation.find(eq);require(it!=n.relaxation.end(),key);
      if(it==n.relaxation.end()) return 0.0;
      value=it->second;
    }
    require(value>0 && value<=1,key);return value;
  };
  require(std::isfinite(s.pseudo_time_scale) && s.pseudo_time_scale>0,"pseudo-auto-time-step-scale-factor","requires a finite positive scale");
  require(std::isfinite(s.pseudo_time_step) && s.pseudo_time_step>=0,"pseudo-time-step","requires a finite nonnegative time step");
  s.momentum_relaxation=relax("mom");s.pressure_relaxation=relax("pressure");
  require(relax("density")==1,"density/relax");require(relax("body-force")==1,"body-force/relax");
  // Compressibility corrections need a thermodynamic speed of sound; with a
  // constant density Fluent's value is undocumented, so dropping it is named.
  if(!laminar && c.viscous.compressibility==true) {
    const auto& info=approximation_info(Approximation::compressibility);
    if(std::find(allowed.begin(),allowed.end(),info.id)!=allowed.end()) out.approximations.push_back(info.id);
    else require(false,std::string(info.setting),"compressibility correction is not implemented; opt in with --approximate "+
      std::string(info.name)+" to "+std::string(info.kutta_action)+"; effect: "+std::string(info.effect_documentation));
  }
  if(sst) {
  require(n.wall_distance_method==2,"wall-distance/method","only wall-distance method 2 is supported");
  require(n.wall_distance_quad_correction==true,"wall-distance/quad-correction?","requires quadrilateral wall-distance correction");
  require(relax("turb-viscosity")==1,"turb-viscosity/relax");
  s.sst=SstSettings{};auto& opt=*s.sst;auto& model=opt.model;
  opt.k_convection=convection(n.k,"k/scheme");opt.omega_convection=convection(n.omega,"omega/scheme");
  opt.k_relaxation=relax("k");opt.omega_relaxation=relax("omega");
  model.production_limiter=required(t.production_limiter,"kw-not-trans-sst-production-limiter?");
  model.production_clip=float(required(t.production_factor,"sst-pk-factor"));
  model.low_re=required(t.low_re,"kw-low-re-mod?");
  require(t.curvature_correction==false,"curvature-correction?","curvature correction is not implemented");
  require(t.sst_strain==false && t.shear_modification==true,"kw-sst-strain?/kw-shear-mod?");
  model.kato_launder=required(t.kato_launder,"kw-not-trans-sst-kato-launder?");
  model.vorticity_production=required(t.vorticity_production,"kw-vorticity-pk?");
  if(t.wall_omega) require(t.wall_omega==WallOmega::correlation,"turbulence/wall-omega-treatment");
  else require(flag(c.variables,"kw-wall-omega-treatment-r13?"),"kw-wall-omega-treatment-r13?");
  const std::map<std::string,double> coefficients{{"sst-a1",.31},{"sst-beta-i1",.075},{"sst-beta-i2",.0828},
    {"sst-sig-k1",1.176},{"sst-sig-k2",1},{"sst-sig-w1",2},{"sst-sig-w2",1.168},{"sst-kappa",.41},
    {"kw-beta-star-inf",.09},{"kw-alpha-star-inf",1},{"kw-beta-i",.072},{"kw-alpha-0",.11111},
    {"kw-r-w",2.95},{"kw-r-k",6},{"kw-r-beta",8},{"kw-mt0",.25},{"kw-zeta-star",1.5},
    {"kw-alpha-inf",.52},{"kw-sig-k",2},{"kw-sig-w",2}};
  for(const auto& [key,value]:t.coefficients) {const auto it=coefficients.find(key);require(it!=coefficients.end() && std::abs(value-it->second)<1e-10,key);}
  model.omega_floor=float(required(t.omega_floor,"omega/min-limit"));model.k_floor=float(number(c.variables,"k/min-limit"));
  model.max_viscosity_ratio=float(number(c.variables,"turb-viscosity/max-limit"));
  model.cross_diffusion_clip=float(number(c.variables,"sst-limit-cross-diff-to-n-times-w-destruct"));
  for(const auto& [key,value]:std::map<std::string,float>{
      {"sst-pk-factor",model.production_clip},{"omega/min-limit",model.omega_floor},
      {"k/min-limit",model.k_floor},{"turb-viscosity/max-limit",model.max_viscosity_ratio},
      {"sst-limit-cross-diff-to-n-times-w-destruct",model.cross_diffusion_clip}})
    require(std::isfinite(value) && value>0,key,"requires a positive value representable in float storage");
  require(!(model.kato_launder && model.vorticity_production),"kw-not-trans-sst-kato-launder?/kw-vorticity-pk?","production options are mutually exclusive");
  require(flag(c.variables,"sst-linearize-neg-cross-diff?"),"sst-linearize-neg-cross-diff?");
  out.convergence_notes.push_back("kw-dissipation-newton?: Newton implicit destruction with matching RHS correction");
  }
  // Fluent's Intermittency Transition Model (TG 4.8) on SST. Every constant
  // Fluent lets a user change is read and used; only the switches whose
  // physics Kutta lacks are rejected. Transition SST and the algebraic model
  // are unsupported issues of the case itself (fluent_case.cpp).
  transition=t.transition==TransitionModel::intermittency;
  if(transition) require(sst,"case-config/rp-trans-one-eqn?","the intermittency transition model is implemented with k-omega SST only");
  if(transition && sst) {
    auto& tr=s.sst->transition.emplace();auto& m=tr.model;
    const auto coefficient=[&](const std::string& key) {
      const auto it=t.transition_coefficients.find(key);
      require(it!=t.transition_coefficients.end(),"rp/"+key,"required setting absent or unparsed");
      const double v=it==t.transition_coefficients.end()?0:it->second;
      require(std::isfinite(float(v)),"rp/"+key,"requires a value representable in float storage");return float(v);
    };
    for(auto [key,target]:std::vector<std::pair<std::string,float*>>{{"trans-one-eqn-ctu1",&m.ctu1},{"trans-one-eqn-ctu2",&m.ctu2},
        {"trans-one-eqn-ctu3",&m.ctu3},{"trans-one-eqn-cpg1",&m.cpg1},{"trans-one-eqn-cpg2",&m.cpg2},{"trans-one-eqn-cpg3",&m.cpg3},
        {"trans-one-eqn-lim-cpg1",&m.cpg1_lim},{"trans-one-eqn-lim-cpg2",&m.cpg2_lim},{"trans-one-eqn-csep",&m.csep},
        {"trans-one-eqn-tke-term-coef",&m.ck},{"intermit-sig",&m.sigma},{"trans-intermit/min-limit",&m.minimum},
        {"trans-intermit/max-limit",&m.maximum}}) *target=coefficient(key);
    require(m.csep>=0 && m.ck>=0,"rp/trans-one-eqn-csep/tke-term-coef","requires nonnegative values");
    require(m.sigma>0,"rp/intermit-sig","requires a positive Prandtl number");
    require(m.minimum>0 && m.minimum<m.maximum && m.maximum<=1,"rp/trans-intermit/min-limit/max-limit","requires 0 < min < max <= 1");
    require(t.transition_crossflow==false,"rp/trans-one-eqn-crossflow?","crossflow transition (TG 4.8.2, local C1 criterion) is not implemented");
    require(t.transition_smoothing==false,"rp/trans-intermit-smoothing?","intermittency smoothing is undocumented and not implemented");
    // Fluent documents dV/dy but not its selector; 1 is the stored default.
    require(t.transition_normal_derivative==1,"rp/trans-one-eqn-dudn-method","only the default method (1) is admitted");
    out.convergence_notes.push_back("trans-one-eqn-dudn-method 1: dV/dy = n.grad(n.U) with n from the nearest wall point "
      "(Menter et al. 2015); Fluent does not document the method it selects");
    tr.convection=convection(n.intermittency,"intermit/scheme");tr.relaxation=relax("intermit");
  }
  if(epsilon_model) {
    // Standard k-epsilon, TG 4.3.1, with standard wall functions.
    auto o=standard_ke_options();
    const std::map<std::string,float*> constants{{"kecmu",&o.cmu_standard},{"kec1",&o.c1_standard},{"kec2",&o.c2_standard},
      {"kesigk",&o.sigma_k},{"kesige",&o.sigma_epsilon}};
    for(const auto& [key,target]:constants) {
      const double value=number(c.variables,key);
      // Fluent stores these in single precision (kec1 1.4400001).
      require(std::abs(value-*target)<=2e-7*std::abs(*target)+1e-12,key,"requires the standard k-epsilon constant "+std::to_string(*target));
    }
    const auto option=[&](const std::string& key) {
      const auto it=t.unsupported_options.find(key);return it==t.unsupported_options.end()?std::optional<bool>{}:std::optional<bool>{it->second};
    };
    require(option("ke-std-wall?")==true,"ke-std-wall?","requires standard wall functions");
    for(auto key:{"ke-enh-wall?","ke-scalable-wall?","ke-non-eqm","ke-nwt-ml?","ke-udf-wall?"})
      require(option(key)!=true,key,"only standard wall functions are implemented");
    require(t.low_re==false,"ke-low-re?","low-Reynolds k-epsilon is not implemented");
    require(t.kato_launder==false,"ke-kato-launder?","Kato-Launder production is not implemented");
    require(t.production_limiter==false,"ke-production-limiter?","the production limiter is not implemented");
    require(t.curvature_correction!=true,"curvature-correction?","curvature correction is not implemented");
    // ke-rough-wall-treatment-r14? selects a roughness law; walls are smooth
    // (rough-bc 0, zero height, checked per wall). ke-uncoupled? solves k and
    // epsilon as separate scalar equations, as Kutta does.
    handled.insert({"rp/ke-uncoupled?","rp/ke-rough-wall-treatment-r14?"});
    const auto& ke=approximation_info(Approximation::ke_options);
    if(option("ke-source-mod")==true || option("ke-time-marching?")==true) {
      if(std::find(allowed.begin(),allowed.end(),ke.id)!=allowed.end()) {
        out.approximations.push_back(ke.id);handled.insert({"rp/ke-source-mod","rp/ke-time-marching?"});
      } else require(false,std::string(ke.setting),"undocumented Fluent k-epsilon solver switches; opt in with --approximate "+
        std::string(ke.name)+" to "+std::string(ke.kutta_action)+"; effect: "+std::string(ke.effect_documentation));
    }
    o.k_floor=float(number(c.variables,"k/min-limit"));o.epsilon_floor=float(number(c.variables,"epsilon/min-limit"));
    o.max_viscosity_ratio=float(number(c.variables,"turb-viscosity/max-limit"));
    for(const auto& [key,value]:std::map<std::string,float>{{"k/min-limit",o.k_floor},{"epsilon/min-limit",o.epsilon_floor},{"turb-viscosity/max-limit",o.max_viscosity_ratio}})
      require(std::isfinite(value) && value>0,key,"requires a positive value representable in float storage");
    s.realizable_ke=RealizableKeSettings{};auto& ke_settings=*s.realizable_ke;ke_settings.model=o;
    ke_settings.k_convection=convection(n.k,"k/scheme");ke_settings.epsilon_convection=convection(n.epsilon,"epsilon/scheme");
    ke_settings.k_relaxation=relax("k");ke_settings.epsilon_relaxation=relax("epsilon");
    require(relax("turb-viscosity")==1,"turb-viscosity/relax");
  }
  for(const auto& issue:c.issues) if(issue.classification==IssueClass::unsupported) {
    // Report definitions are optional. An absent list creates no reports;
    // a present but unparsed list must still be rejected.
    const bool absent_reports=issue.path=="rp/monitor/report-definitions" &&
      !c.variables.entry("monitor/report-definitions");
    if(!handled.contains(issue.path) && !absent_reports && !unevaluated_report(c,issue.path)) errors.push_back(issue);
  }
  for(const auto& r:c.reports) if(output_only(r))
    out.convergence_notes.push_back("report "+r.name+" ("+r.fluent_type+"): not evaluated; Kutta evaluates force, drag, lift and downforce reports only");
  require(!flag(c.config,"rp-des?"),"rp-des?");
  s.reference_point={number(c.variables,"pressure-reference/x"),number(c.variables,"pressure-reference/y"),number(c.variables,"pressure-reference/z")};
  for(const auto& z:c.boundaries) if(z.type==ZoneType::fan)
    require(false,"zone/"+std::to_string(z.id)+"/"+z.fluent_type,
      "the fan polynomial solver path is analytic-verified; native Fluent fan polynomial input conventions still require a reference fixture");
  for(const auto& z:c.boundaries) if(z.type==ZoneType::porous_jump) {
    const auto key="zone/"+std::to_string(z.id)+"/";
    out.settings.interior_conditions.push_back({int(z.id),PorousJump{
      scalar(z.jump_permeability,key+"alpha"),scalar(z.jump_thickness,key+"dm"),scalar(z.jump_inertial,key+"c2")}});
    const auto face_zone=std::find_if(mesh.internal_face_zones.begin(),mesh.internal_face_zones.end(),[&](const auto& f){return f.source_zone_id==z.id;});
    require(face_zone!=mesh.internal_face_zones.end(),key+"identity","porous jump requires retained internal face-zone metadata");
  }
  for(const auto& patch:mesh.patches) {
    // A planar case has no z flux, z derivative or z face: the extrusion caps
    // are symmetry planes and w stays zero (docs/research/fluent-reference-cases.md).
    if(planar && patch.source_zone_id<0) {out.boundaries.push_back({patch.name,Symmetry{}});continue;}
    const auto it=std::find_if(c.boundaries.begin(),c.boundaries.end(),[&](const auto& z){return z.id==patch.source_zone_id;});
    require(it!=c.boundaries.end() && it->name==patch.name,"mesh/patch identity: "+patch.name);
    if(it==c.boundaries.end()) continue;
    const auto& z=*it;
    const auto prefix="zone/"+std::to_string(z.id)+"/";
    require(patch.size>0,prefix+"empty patch");
    if(patch.size<=0) continue;
    if(z.type==ZoneType::wall) {
      require(z.specified_shear==false && z.shear_bc==0 && z.rough_bc==0,prefix+"shear-bc","requires specified-shear?=false, shear-bc=0 and rough-bc=0");
      require(scalar(z.roughness_height,prefix+"roughness-height")==0,prefix+"roughness-height");
      try {out.boundaries.push_back({patch.name,fluent_wall(z)});} catch(const Error& e) {errors.push_back({prefix+"motion",e.what()});}
    } else if(z.type==ZoneType::symmetry) out.boundaries.push_back({patch.name,Symmetry{}});
    else if(z.type==ZoneType::velocity_inlet) {
      require(z.frame_of_reference==0,prefix+"frame-of-reference");
      Vec3d velocity;
      if(z.velocity_spec==VelocitySpec::components) velocity=vector(z.velocity,prefix+"velocity");
      else if(z.velocity_spec==VelocitySpec::magnitude_direction) {
        const auto direction=vector(z.direction,prefix+"direction");require(norm(direction)>0,prefix+"direction");
        velocity=direction*(scalar(z.velocity_magnitude,prefix+"vmag")/norm(direction));
      } else {
        require(z.velocity_spec==VelocitySpec::normal,prefix+"velocity-spec");
        const auto area=[&](Index f) {return face_geometry(mesh,f).area;};
        const auto normal=area(patch.start)/norm(area(patch.start));
        for(Index f=patch.start;f<patch.start+patch.size;++f)
          require(norm(area(f)/norm(area(f))-normal)<1e-6,prefix+"normal velocity on nonplanar patch");
        velocity=normal*(-scalar(z.velocity_magnitude,prefix+"vmag"));
      }
      out.boundaries.push_back({patch.name,VelocityInlet{velocity,turbulence(z),intermittency(z)}});
    } else if(z.type==ZoneType::pressure_outlet) {
      // Writers up to Fluent 18.0 store neither key (the corpus's 17.2 and 18.0
      // cases; 21.2 and later store both). What such a writer did at backflow
      // faces is undocumented, so running without suppression is named.
      const bool early=c.solver_version && version_major(*c.solver_version)<=18;
      if(!z.raw.entry("prevent-reverse-flow?") && early) {
        const auto& info=approximation_info(Approximation::reverse_flow);
        if(std::find(allowed.begin(),allowed.end(),info.id)!=allowed.end()) {
          if(std::find(out.approximations.begin(),out.approximations.end(),info.id)==out.approximations.end()) out.approximations.push_back(info.id);
        } else require(false,prefix+std::string(info.setting),"absent (Fluent "+*c.solver_version+" writer); opt in with --approximate "+
          std::string(info.name)+" to "+std::string(info.kutta_action)+"; effect: "+std::string(info.effect_documentation));
      } else require(z.prevent_reverse_flow==false,prefix+"prevent-reverse-flow?");
      require(number(z.raw,"direction-spec")==1,prefix+"direction-spec (normal backflow)");
      // Any multiplier of a constant zero gauge pressure is zero.
      const auto gauge=scalar(z.gauge_pressure,prefix+"p");
      if(!z.raw.entry("p-profile-multiplier") && gauge==0 && z.gauge_pressure && z.gauge_pressure->kind==ScalarValue::Kind::constant)
        out.convergence_notes.push_back(prefix+"p-profile-multiplier: absent (Fluent "+c.solver_version.value_or("unknown")+" writer); immaterial for the constant zero gauge pressure");
      else require(number(z.raw,"p-profile-multiplier")==1,prefix+"p-profile-multiplier");
      const auto type=required(z.backflow_pressure_spec,prefix+"p-backflow-spec");require(type==0 || type==1,prefix+"p-backflow-spec");
      out.boundaries.push_back({patch.name,PressureOutlet{scalar(z.gauge_pressure,prefix+"p"),NormalBackflow{},
        type==0?BackflowPressure::Total:BackflowPressure::Static,turbulence(z),intermittency(z)}});
    } else if(z.type==ZoneType::pressure_inlet || z.type==ZoneType::mass_flow_inlet) {
      require(z.prevent_reverse_flow!=true,prefix+"prevent-reverse-flow?","reverse-flow suppression is unsupported");
      require(z.frame_of_reference==0,prefix+"frame-of-reference","only absolute inlet data is implemented");
      require(z.direction_spec==0 || z.direction_spec==1,prefix+"direction-spec","only Cartesian direction or inward normal is implemented");
      std::optional<Vec3d> direction;
      if(z.direction_spec==0) {direction=vector(z.direction,prefix+"direction");require(norm(*direction)>0,prefix+"direction");}
      if(z.type==ZoneType::pressure_inlet) out.boundaries.push_back({patch.name,PressureInlet{scalar(z.total_pressure,prefix+"p0"),direction,turbulence(z),intermittency(z)}});
      else {
        require(z.mass_flow_spec==0 || z.mass_flow_spec==1,prefix+"flow-spec","only uniform mass-flow rate or uniform mass flux are implemented");
        double rate=0;
        if(z.mass_flow_spec==0) rate=scalar(z.mass_flow,prefix+"mass-flow");
        else {double area=0;for(Index f=patch.start;f<patch.start+patch.size;++f) area+=norm(face_geometry(mesh,f).area);rate=area*scalar(z.mass_flux,prefix+"mass-flux");}
        require(rate>=0,prefix+"mass-flow","negative mass-flow inlet is unsupported");
        out.boundaries.push_back({patch.name,MassFlowInlet{rate,direction,turbulence(z),intermittency(z)}});
      }
    } else require(false,prefix+"type");
  }
  finish();return out;
}
std::vector<ForceDefinition> force_definitions(const fluent::FluentCase& c,const PolyMesh& mesh) {
  std::vector<ForceDefinition> out;std::vector<SettingIssue> errors;
  for(const auto& r:c.reports) {
    if(output_only(r)) continue;
    const auto fail=[&](const std::string& why){errors.push_back({"report/"+r.name,why});};
    ForceDefinition d;d.name=r.name;
    if(!r.zone_ids || r.zone_ids->empty()) fail("wall zone identities are unresolved");
    else for(auto id:*r.zone_ids) {
      auto p=std::find_if(mesh.patches.begin(),mesh.patches.end(),[&](const auto& p){return p.source_zone_id==id && p.kind==PatchKind::Wall;});
      if(p==mesh.patches.end()) fail("zone "+std::to_string(id)+" is not a wall patch");
      else if(std::find(d.patches.begin(),d.patches.end(),p->name)==d.patches.end()) d.patches.push_back(p->name);
    }
    if(!r.direction) fail("force direction is absent");
    else {
      const auto& v=*r.direction;d.direction=Vec3d{v[0],v[1],v[2]};
      if(!std::isfinite(norm(*d.direction)) || norm(*d.direction)==0) fail("force direction must be finite and nonzero");
    }
    if(!r.scaled) fail("scaled? is absent");else d.scaled=*r.scaled;
    if(d.scaled) for(auto [name,v]:{std::pair{"density",c.reference.density},{"velocity",c.reference.velocity},{"area",c.reference.area}})
      if(!v || !std::isfinite(*v) || *v<=0) fail(std::string("reference/")+name+" must be finite and positive for coefficients");
    if(std::any_of(out.begin(),out.end(),[&](const auto& prior){return prior.name==d.name;})) fail("duplicate report name");
    out.push_back(std::move(d));
  }
  if(!errors.empty()) throw RejectedSettings(std::move(errors));
  return out;
}
bool uses_pseudo_time(const FluentCase& c) {
  const auto& p=c.numerics.pseudo;
  const auto method=c.numerics.coupling==Coupling::coupled?p.coupled_method:p.segregated_method;
  return p.formulation?*p.formulation:method==1;
}
FluentSetup fluent_setup(const FluentCase& c,const PolyMesh& mesh,std::span<const Approximation> allowed) {
  return SetupValidation{}.build(c,mesh,allowed);
}
CellFields fluent_fields(const FluentCase& c,const FluentData& d,const PolyMesh& mesh) {
  const auto ids=cell_ids(c,mesh);CellFields out;
  const auto get=[&](const std::string& key,std::uint64_t id) {
    const auto it=d.cell_fields.find(key);require(it!=d.cell_fields.end(),key);
    const auto value=it->second.at(id);require(value.size()==1 && std::isfinite(value[0]),key);return value[0];
  };
  const bool sst=c.viscous.model==ViscousModel::kw_sst,epsilon=c.viscous.model==ViscousModel::ke_standard;
  const bool transition=sst && c.viscous.transition==TransitionModel::intermittency;
  // A planar case stores no SV_W; its w is zero.
  const bool planar=c.dimension==2;
  require(!planar || !d.cell_fields.contains("SV_W"),"SV_W in a 2D case");
  for(auto id:ids) {
    out.velocity.push_back({get("SV_U",id),get("SV_V",id),planar?0.:get("SV_W",id)});out.pressure.push_back(get("SV_P",id));
    if(sst) {out.k.push_back(get("SV_K",id));out.omega.push_back(get("SV_O",id));}
    if(epsilon) {out.k.push_back(get("SV_K",id));out.epsilon.push_back(get("SV_D",id));}
    // Fluent's intermittency solver variable (UDF macro C_INTERMIT).
    if(transition) out.intermittency.push_back(get("SV_INTERMIT",id));
  }
  return out;
}
Initialization initial_fields(const FluentCase& c,const PolyMesh& mesh) {
  SetupValidation v;
  const bool hybrid=v.flag(c.variables,"hyb-init?");
  v.require(!v.flag(c.variables,"fmg-init/enable?"),"fmg-init/enable?");
  if(hybrid) {
    HybridInitialization h;
    const auto pair=[&](const std::string& key,auto& out) {
      const auto* value=c.variables.value(key);
      v.require(value && value->kind==SExpr::Kind::list && value->elements.size()==2,key,"requires two scalar equation controls");
      if(!value || value->elements.size()!=2) return;
      for(int j=0;j<2;++j) {
        using T=typename std::remove_reference_t<decltype(out)>::value_type;
        if constexpr(std::is_same_v<T,bool>) {
          v.require(value->elements[j].kind==SExpr::Kind::boolean,key,"requires booleans");
          out[j]=value->elements[j].boolean;
        } else {
          const auto n=value->elements[j].number();
          v.require(n && std::isfinite(*n),key,"requires finite numbers");out[j]=n.value_or(0);
        }
      }
    };
    const double count=v.number(c.variables,"hyb-init/niter");
    v.require(count>=1 && count<=100000 && count==std::floor(count),"hyb-init/niter","requires a positive integer <=100000");
    if(count>=1 && count<=100000) h.iterations=int(count);
    pair("hyb-init/explicit-relax",h.explicit_relaxation);pair("hyb-init/relax",h.implicit_relaxation);
    pair("hyb-init/residual-tolerance",h.tolerance);pair("hyb-init/solve?",h.solve);pair("hyb-init/patch",h.initial_scalar);
    for(int j=0;j<2;++j) {
      v.require(h.explicit_relaxation[j]>0 && h.explicit_relaxation[j]<=1,"hyb-init/explicit-relax");
      v.require(h.implicit_relaxation[j]>0 && h.implicit_relaxation[j]<=1,"hyb-init/relax");
      v.require(h.tolerance[j]>0 && h.tolerance[j]<1,"hyb-init/residual-tolerance");
    }
    h.averaged_turbulence=v.flag(c.variables,"hyb-init/turb-domain-averaged?");
    h.linear_potential=v.flag(c.variables,"hyb-init/patch-ads?");
    h.constant_speed=v.flag(c.variables,"hyb-init/constant-vel-vector?");
    h.turbulence_intensity=v.number(c.variables,"hyb-init/turb-intens");
    h.turbulence_viscosity_ratio=v.number(c.variables,"hyb-init/turb-vr");
    v.require(h.turbulence_intensity>=0 && h.turbulence_viscosity_ratio>0,"hyb-init/turb-intens/turb-vr");
    if(v.flag(c.variables,"hyb-init/pressure-from-bc?")) {
      for(const auto& z:c.boundaries) if(z.type==ZoneType::velocity_inlet)
        h.inlet_pressure.emplace_back(z.name,v.scalar(z.gauge_pressure,"zone/"+std::to_string(z.id)+"/p"));
    }
    v.finish();return h;
  }
  CellFields out;out.velocity.assign(mesh.n_cells,{v.number(c.variables,"x-velocity/default"),v.number(c.variables,"y-velocity/default"),v.number(c.variables,"z-velocity/default")});
  out.pressure.assign(mesh.n_cells,v.number(c.variables,"pressure/default"));
  const bool sst=c.viscous.model==ViscousModel::kw_sst;
  if(sst) {out.k.assign(mesh.n_cells,v.number(c.variables,"k/default"));out.omega.assign(mesh.n_cells,v.number(c.variables,"omega/default"));}
  if(sst && c.viscous.transition==TransitionModel::intermittency) {
    const double gamma=v.number(c.variables,"intermit/default");
    v.require(gamma>=0 && gamma<=1,"intermit/default","intermittency must lie in [0, 1]");
    out.intermittency.assign(mesh.n_cells,gamma);
  }
  if(c.viscous.model==ViscousModel::ke_standard) {
    out.k.assign(mesh.n_cells,v.number(c.variables,"k/default"));out.epsilon.assign(mesh.n_cells,v.number(c.variables,"epsilon/default"));
    v.require(out.epsilon.empty() || out.epsilon.front()>0,"epsilon/default","requires positive dissipation");
  }
  // BEGIN FSAE physics: Fluent standard initialization may be relative even
  // when the transported solution variables are absolute velocities.
  if(std::any_of(c.cell_zones.begin(),c.cell_zones.end(),[](const auto& z){return z.frame_motion==true;})) {
    const bool relative=v.flag(c.variables,"init-relative?");
    try {
      const auto physics=fluent_zone_physics(c);const auto frames=cell_frames(mesh,physics);const auto geometry=compute_geometry(mesh);
      if(relative) for(Index i=0;i<mesh.n_cells;++i) out.velocity[i]=out.velocity[i]+frames[i].velocity(geometry.cell_centroid[i]);
    } catch(const Error& e) {v.require(false,"initialization/frame",e.what());}
  }
  // END FSAE physics
  v.require(out.k.empty() || out.k.front()>=0,"k/default","requires nonnegative turbulent kinetic energy");
  v.require(out.omega.empty() || out.omega.front()>0,"omega/default","requires positive specific dissipation");
  v.finish();return out;
}
void initialize(SteadySolver& solver,const Initialization& init) {
  if(const auto* fields=std::get_if<CellFields>(&init)) {
    if(!fields->epsilon.empty()) solver.initialize_ke_fields(fields->velocity,fields->pressure,fields->k,fields->epsilon);
    else solver.initialize_fields(fields->velocity,fields->pressure,fields->k,fields->omega,fields->intermittency);
  }
  else solver.initialize_hybrid(std::get<HybridInitialization>(init));
}
} // namespace kutta::setup
