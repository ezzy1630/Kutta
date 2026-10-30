#include "kutta/setup/fluent.hpp"
#include "kutta/setup/fluent_physics.hpp"

#include <algorithm>
#include <cmath>
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
  std::vector<std::string> errors;
 public:
  void require(bool ok,const std::string& key,const std::string& reason="value is unsupported or required setting is absent") {
    if(!ok) errors.push_back(key+": "+reason);
  }
  template<class T> T required(const std::optional<T>& value,const std::string& key) {
    require(value.has_value(),key,"required setting absent or unparsed");return value.value_or(T{});
  }
  void finish() const {
    if(errors.empty()) return;
    std::string message="Fluent setup rejected:";
    for(const auto& e:errors) message+="\n  "+e;
    throw Error(message);
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
    case TurbulenceSpec::k_dissipation:return TurbulenceValues{get("k"),get("omega")};
    case TurbulenceSpec::intensity_length:return IntensityLength{get("turb-intensity"),get("turb-length-scale")};
    case TurbulenceSpec::intensity_diameter:return IntensityLength{get("turb-intensity"),.07*get("turb-hydraulic-diam")};
    case TurbulenceSpec::intensity_viscosity_ratio:return IntensityViscosityRatio{get("turb-intensity"),get("turb-viscosity-ratio")};
  }
  require(false,key+"ke-spec");return TurbulenceValues{};
}
FluentSetup build(const FluentCase& c,const PolyMesh& mesh,std::span<const Approximation> allowed);
};
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
  require(c.dimension==3,"dimension","requires a 3D case");require(c.pressure_based==true,"rp-seg?","requires the pressure-based solver");require(c.unsteady==false,"rp-unsteady?","only steady cases are supported");
  require(c.viscous.model==ViscousModel::kw_sst,"viscous-model (requires SST)");
  require(c.mesh.cells && *c.mesh.cells==std::uint64_t(mesh.n_cells),"mesh/cell-count");
  FluentSetup out;auto& s=out.settings;const auto& n=c.numerics;const auto& t=c.viscous;
  try {s.cell_zone_physics=fluent_zone_physics(c);} catch(const Error& e) {errors.push_back(e.what());}
  require(!c.cell_zones.empty(),"cell zones");
  std::optional<double> rho,mu;
  for(const auto& z:c.cell_zones) {
    const auto prefix="zone/"+std::to_string(z.id)+"/";
    require(z.type==ZoneType::fluid,prefix+"type");require(z.laminar==false,prefix+"laminar?");
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
  require(n.gradient==Gradient::least_squares_cell,"recon/cell-lsq?","requires least-squares cell gradients");
  require(n.pressure==Scheme::second_pressure,"pressure/scheme","requires second-order pressure interpolation");
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
  require(n.flux_index==1,"pbs/flux-index","only flux index 1 is supported");required(n.flux_auto_select,"pbs/flux-auto-select?");
  require(n.wall_distance_method==2,"wall-distance/method","only wall-distance method 2 is supported");
  require(n.wall_distance_quad_correction==true,"wall-distance/quad-correction?","requires quadrilateral wall-distance correction");
  const auto coupling=required(n.coupling,"flow/scheme");
  if(coupling!=Coupling::simplec) out.convergence_notes.push_back("flow/scheme: SIMPLEC replaces "+std::string(name(coupling))+"; convergence path only, Majumdar physical diagonal");
  if(required(n.high_order_relaxation,"recon/relax/relax?")) out.convergence_notes.push_back("high-order term relaxation omitted; same fixed point");
  const bool pseudo=required(n.pseudo.formulation,"pseudo-transient-formulation?");
  require(pseudo,"pseudo-transient-formulation? (zero-pseudo SIMPLEC not adapted)");
  require((coupling==Coupling::coupled?n.pseudo.coupled_method:n.pseudo.segregated_method)==1,"pseudo-time-method/dt-method");
  require(n.pseudo.mode==0,"pseudo-transient-mode");
  require(n.pseudo.length_scale_method==1,"pseudo/autotime-lscale-option");
  s.pseudo_time_scale=required(n.pseudo.scale_factor,"pseudo-auto-time-step-scale-factor");
  if(!required(n.pseudo.automatic,"pseudo-auto-time-step?")) {
    s.pseudo_time_step=required(n.pseudo.time_step,"pseudo-time-step");
    require(s.pseudo_time_step>0,"pseudo-time-step","an explicit step must be positive; zero would select automatic stepping");
  }
  const auto relax=[&](const std::string& eq) {
    const auto key=coupling==Coupling::coupled && (eq=="mom" || eq=="pressure")?
      "pressure-coupled/"+eq+"/pseudo-explicit-relax":eq+"/pseudo-relax";
    const auto it=n.pseudo_relaxation.find(key);require(it!=n.pseudo_relaxation.end(),key);
    if(it==n.pseudo_relaxation.end()) return 0.0;
    require(it->second>0 && it->second<=1,key);return it->second;
  };
  require(std::isfinite(s.pseudo_time_scale) && s.pseudo_time_scale>0,"pseudo-auto-time-step-scale-factor","requires a finite positive scale");
  require(std::isfinite(s.pseudo_time_step) && s.pseudo_time_step>=0,"pseudo-time-step","requires a finite nonnegative time step");
  s.momentum_relaxation=relax("mom");s.pressure_relaxation=relax("pressure");
  require(relax("density")==1,"density/pseudo-relax");require(relax("body-force")==1,"body-force/pseudo-relax");
  require(relax("turb-viscosity")==1,"turb-viscosity/pseudo-relax");
  s.sst=SstSettings{};auto& opt=*s.sst;auto& model=opt.model;
  opt.k_convection=convection(n.k,"k/scheme");opt.omega_convection=convection(n.omega,"omega/scheme");
  opt.k_relaxation=relax("k");opt.omega_relaxation=relax("omega");
  model.production_limiter=required(t.production_limiter,"kw-not-trans-sst-production-limiter?");
  model.production_clip=float(required(t.production_factor,"sst-pk-factor"));
  model.low_re=required(t.low_re,"kw-low-re-mod?");
  require(t.curvature_correction==false,"curvature-correction?","curvature correction is not implemented");require(t.compressibility==false,"turb-compress-mod?","compressibility correction is not implemented");
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
  const std::set<std::string> handled{"rp/turbulence/wall-omega-treatment","rp/kw-dissipation-newton?",
    "rp/kw-wall-omega-treatment-r13?","rp/sst-des-iddes-simplify_omit_fe?","rp/sst-linearize-neg-cross-diff?"};
  for(const auto& issue:c.issues) if(issue.classification==IssueClass::unsupported) {
    // Report definitions are optional. An absent list creates no reports;
    // a present but unparsed list must still be rejected.
    const bool absent_reports=issue.path=="rp/monitor/report-definitions" &&
      !c.variables.entry("monitor/report-definitions");
    if(!handled.contains(issue.path) && !absent_reports) errors.push_back(issue.path+": "+issue.reason);
  }
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
    const auto it=std::find_if(c.boundaries.begin(),c.boundaries.end(),[&](const auto& z){return z.id==patch.source_zone_id;});
    require(it!=c.boundaries.end() && it->name==patch.name,"mesh/patch identity: "+patch.name);
    if(it==c.boundaries.end()) continue;
    const auto& z=*it;
    const auto prefix="zone/"+std::to_string(z.id)+"/";
    require(patch.size>0,prefix+"empty patch");
    if(patch.size<=0) continue;
    if(z.type==ZoneType::wall) {
      require(z.wall_motion!=WallMotion::rotating,prefix+"motion-bc",
        "native rotating-wall support held: V7.2 polyhedral pressure accuracy still exceeds the 1% gate; the C++ verification path remains available");
      require(z.specified_shear==false && z.shear_bc==0 && z.rough_bc==0,prefix+"shear-bc","requires specified-shear?=false, shear-bc=0 and rough-bc=0");
      require(scalar(z.roughness_height,prefix+"roughness-height")==0,prefix+"roughness-height");
      try {out.boundaries.push_back({patch.name,fluent_wall(z)});} catch(const Error& e) {errors.push_back(e.what());}
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
        const auto geom=compute_geometry(mesh);const auto normal=geom.face_area[patch.start]/norm(geom.face_area[patch.start]);
        for(Index f=patch.start;f<patch.start+patch.size;++f)
          require(norm(geom.face_area[f]/norm(geom.face_area[f])-normal)<1e-6,prefix+"normal velocity on nonplanar patch");
        velocity=normal*(-scalar(z.velocity_magnitude,prefix+"vmag"));
      }
      out.boundaries.push_back({patch.name,VelocityInlet{velocity,turbulence(z)}});
    } else if(z.type==ZoneType::pressure_outlet) {
      require(z.prevent_reverse_flow==false,prefix+"prevent-reverse-flow?");
      require(number(z.raw,"direction-spec")==1,prefix+"direction-spec (normal backflow)");
      require(number(z.raw,"p-profile-multiplier")==1,prefix+"p-profile-multiplier");
      const auto type=required(z.backflow_pressure_spec,prefix+"p-backflow-spec");require(type==0 || type==1,prefix+"p-backflow-spec");
      out.boundaries.push_back({patch.name,PressureOutlet{scalar(z.gauge_pressure,prefix+"p"),NormalBackflow{},
        type==0?BackflowPressure::Total:BackflowPressure::Static,turbulence(z)}});
    } else if(z.type==ZoneType::pressure_inlet || z.type==ZoneType::mass_flow_inlet) {
      require(z.prevent_reverse_flow!=true,prefix+"prevent-reverse-flow?","reverse-flow suppression is unsupported");
      require(z.frame_of_reference==0,prefix+"frame-of-reference","only absolute inlet data is implemented");
      require(z.direction_spec==0 || z.direction_spec==1,prefix+"direction-spec","only Cartesian direction or inward normal is implemented");
      std::optional<Vec3d> direction;
      if(z.direction_spec==0) {direction=vector(z.direction,prefix+"direction");require(norm(*direction)>0,prefix+"direction");}
      if(z.type==ZoneType::pressure_inlet) out.boundaries.push_back({patch.name,PressureInlet{scalar(z.total_pressure,prefix+"p0"),direction,turbulence(z)}});
      else {
        require(z.mass_flow_spec==0 || z.mass_flow_spec==1,prefix+"flow-spec","only uniform mass-flow rate or uniform mass flux are implemented");
        double rate=0;
        if(z.mass_flow_spec==0) rate=scalar(z.mass_flow,prefix+"mass-flow");
        else {const auto geometry=compute_geometry(mesh);double area=0;for(Index f=patch.start;f<patch.start+patch.size;++f) area+=norm(geometry.face_area[f]);rate=area*scalar(z.mass_flux,prefix+"mass-flux");}
        require(rate>=0,prefix+"mass-flow","negative mass-flow inlet is unsupported");
        out.boundaries.push_back({patch.name,MassFlowInlet{rate,direction,turbulence(z)}});
      }
    } else require(false,prefix+"type");
  }
  finish();return out;
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
  for(auto id:ids) {
    out.velocity.push_back({get("SV_U",id),get("SV_V",id),get("SV_W",id)});out.pressure.push_back(get("SV_P",id));
    out.k.push_back(get("SV_K",id));out.omega.push_back(get("SV_O",id));
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
  out.k.assign(mesh.n_cells,v.number(c.variables,"k/default"));out.omega.assign(mesh.n_cells,v.number(c.variables,"omega/default"));
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
  if(const auto* fields=std::get_if<CellFields>(&init))
    solver.initialize_fields(fields->velocity,fields->pressure,fields->k,fields->omega);
  else solver.initialize_hybrid(std::get<HybridInitialization>(init));
}
} // namespace kutta::setup
