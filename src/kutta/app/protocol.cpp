#include "protocol.hpp"
#include "kutta/setup/settings.hpp"
#include "kutta/core/json.hpp"
#include "kutta/core/path.hpp"
#include "kutta/core/sha256.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>

namespace kutta::app {
namespace {
using namespace fluent;
Issue describe(const SettingIssue& i) {
  Issue out{i.path,i.path,i.classification==IssueClass::unsupported?"unsupported":"informational",{}, {},{}};
  struct Explanation { std::string_view key,label,text; };
  static constexpr Explanation descriptions[]={
    {"curvature", "Streamline curvature correction", "Fluent modifies SST production for rotating or curved strain. Kutta does not yet compute the strain material derivative required by that correction."},
    {"compress", "Compressibility correction", "Fluent modifies turbulence for compressible flow. Kutta currently advances constant-density incompressible momentum; the k-epsilon dilatation term can only be omitted by explicit approximation."},
    {"viscous-model", "Viscous model", "This selects Fluent's laminar or turbulence closure. The Fluent solve adapter implements laminar flow, k-omega SST and standard k-epsilon with standard wall functions; substituting another closure changes turbulent stresses and separation."},
    {"rp-unsteady", "Transient flow", "Fluent advances physical time for this setting. Kutta's case adapter currently solves a steady fixed point; pseudo time is a convergence aid and cannot reproduce physical transients."},
    {"rp-seg", "Pressure-based solver", "This selects Fluent's pressure-based or density-based formulation. Kutta implements the incompressible pressure-based equations and cannot reproduce a density-based solution."},
    {"dimension", "Mesh dimension", "This describes the case dimension. The Fluent solve adapter solves 3D cases, and planar 2D cases as a one-cell extrusion of Fluent's reference depth between two symmetry planes; axisymmetric and other 2D readings are not adapted."},
    {"rf-energy", "Energy equation", "Fluent solves temperature or enthalpy. Kutta omits energy; it is decoupled from momentum only when all active density and viscosity laws are constant. Temperature results are unavailable."},
    {"energy", "Energy equation", "Fluent solves temperature or enthalpy. Kutta omits energy; it is decoupled from momentum only when all active density and viscosity laws are constant. Temperature results are unavailable."},
    {"gravity", "Gravity and buoyancy", "Fluent applies a gravitational body force. The Fluent adapter does not yet transfer this source or its pressure convention, even when density is constant."},
    {"mrf", "Moving reference frame", "Fluent includes rotating-frame momentum sources and frame-relative velocities. The adapter currently implements stationary frames and cannot reproduce Coriolis and centrifugal effects."},
    {"porous", "Porous medium", "Fluent adds directional viscous and inertial resistance. Kutta has no porous-zone source adapter, so the pressure loss cannot be reproduced."},
    {"motion", "Wall or mesh motion", "Fluent prescribes moving walls or a moving frame. The Fluent adapter currently accepts stationary no-slip walls and stationary fluid zones only."},
    {"rough", "Wall roughness", "Fluent modifies the wall law using roughness height and constant. Kutta currently implements smooth-wall laws only, so roughness drag and y-plus cannot be reproduced."},
    {"shear-bc", "Wall shear condition", "Fluent can prescribe shear or slip at a wall. Kutta's adapter currently requires a smooth no-slip wall and computes shear with the active model's wall treatment."},
    {"wall-omega", "SST wall omega treatment", "Fluent selects a wall-omega correlation or a tabulated law. Kutta implements the correlation only; the alternative table is unavailable. Correlation face-omega offset was about 0.42% on the saved Fluent 22.2 elbow fields; this is not a flow-solution error bound."},
    {"recon/cell-lsq", "Gradient reconstruction", "Fluent selects its cell or node gradient method. Kutta implements least-squares and Green-Gauss cell-based gradients, Green-Gauss for laminar and standard k-epsilon flow; node-based reconstruction is not adapted."},
    {"recon/node", "Gradient reconstruction", "Fluent selects node-based gradient reconstruction. Kutta's Fluent adapter implements cell-based gradients only; a different reconstruction changes transport fluxes."},
    {"recon/lim", "Gradient limiter", "Fluent limits reconstructed values before evaluating face fluxes. Kutta implements the standard face-based limiter without filtering; other limiter settings alter the discrete transport equations."},
    {"pressure/scheme", "Pressure interpolation", "Fluent selects the pressure interpolation at faces. Kutta implements standard and second-order pressure; PRESTO, linear and body-force-weighted interpolation are not implemented."},
    {"/scheme", "Discretization scheme", "Fluent selects the face interpolation for this equation. Kutta's transport adapter implements first- and second-order upwind only; other schemes change numerical diffusion and the solution."},
    {"pseudo", "Pseudo-transient controls", "Fluent uses these controls to converge a steady solution. Kutta adapts global pseudo-time stepping and the classic relaxation-only iteration; local or other pseudo-time formulations are unavailable."},
    {"relax", "Equation relaxation", "Fluent damps iterative updates using this factor. Kutta uses the factors of the active formulation (pseudo-time or classic); missing or incompatible controls cannot be replaced by guessed defaults."},
    {"wall-distance", "Wall-distance calculation", "Fluent selects how distance to the wall is computed. Kutta uses exact nearest wall-surface distances with quadrilateral correction; other selections are not adapted."},
    {"hyb-init", "Hybrid initialization", "Fluent solves potential-flow initialization using saved scalar controls. Kutta uses its own AMG cycles; absent or invalid controls cannot be inferred. This affects startup, not the target steady equations."},
    {"fmg-init", "Full multigrid initialization", "Fluent initializes on a hierarchy of coarsened flow meshes. Kutta has no FMG initialization adapter and cannot reproduce that startup path."},
    {"/default", "Standard initialization", "Fluent uses this initial flow or turbulence value. Kutta requires finite values and positive omega or epsilon; it does not invent absent initialization settings."},
    {"material", "Material properties", "Fluent binds a material and its density and viscosity laws to a cell zone. Kutta currently requires one positive, constant, uniform density and viscosity; variable properties or unresolved profiles change the equations."},
    {"density/viscosity", "Density and viscosity", "Fluent assigns density and viscosity by material. Kutta's current incompressible adapter requires constant positive properties shared by every active fluid zone."},
    {"report/", "Report definition", "Fluent defines the selected surfaces, force direction and scaling. Kutta can evaluate wall forces, drag and lift only when names resolve exactly to wall zones and the reference values are valid. ID-only selections cannot safely be treated as thread IDs."},
    {"report-definitions", "Report definitions", "Fluent stores force and other monitor definitions here. Kutta cannot translate a present malformed list; an absent list simply requests no reports."},
    {"zone/", "Zone boundary or cell condition", "Fluent prescribes this zone's physics or boundary values. Kutta's adapter supports stationary fluid, smooth no-slip walls, symmetry, constant velocity inlets and pressure outlets only; unparsed or other active conditions cannot be reproduced."},
    {"sst-", "SST model controls", "Fluent uses this coefficient or switch in SST production, blending or destruction. Kutta requires the implemented coefficient set and finite positive limits; other model controls are not adapted."},
    {"kw-", "k-omega model controls", "Fluent uses this switch or coefficient in the turbulence closure. Kutta implements the explicitly validated SST variants only; the requested or unparsed variant cannot be reproduced."},
    {"/ke-", "k-epsilon model controls", "Fluent uses this switch or coefficient in the k-epsilon closure or its wall treatment. Kutta implements the standard model with standard wall functions and Fluent's default coefficients; other variants are not adapted."},
    {"pbs/flux", "Mass-flux interpolation", "Fluent selects how pressure and velocity form face mass flux. Kutta adapts the momentum-weighted and distance-weighted Rhie-Chow fluxes; automatic selection or other selectors cannot reproduce that coupling."},
    {"mesh/", "Mesh identity", "The case's Fluent zone IDs, counts and names must match the mesh used to place fields. Kutta refuses inconsistent mappings because assigning values to different cells or faces would corrupt the solution."},
    {"case-config", "Active physics model", "Fluent enables a physics model with this switch. Kutta implements steady, constant-density, single-phase laminar, SST and standard k-epsilon flow; this active or unrecognized configuration is outside the case adapter's allowlist."}
  };
  for(const auto& d:descriptions) if(i.path.find(d.key)!=std::string::npos) {
    out.label=d.label;out.explanation=d.text;break;
  }
  if(out.explanation.empty()) out.explanation="This saved Fluent setting controls "+i.path+". Kutta has no validated translation for the requested value, or the required value was not present; using a default could change the case being solved.";
  out.explanation+=" Adapter detail: "+i.reason+".";
  if(out.classification=="unsupported") out.explanation+=" A run is blocked; no isolated quantitative error estimate is available for this setting.";
  if(i.path.find("wall-omega")!=std::string::npos) out.evidence="src/kutta/fv/turbulence/README.md#verification";
  return out;
}
template<class T> void number(std::ostream& s,const std::optional<T>& n) {if(n) json_number(s,double(*n));else s<<"null";}
void boolean(std::ostream& s,const std::optional<bool>& n) {s<<(n?(*n?"true":"false"):"null");}
void scalar(std::ostream& s,const std::optional<ScalarValue>& v) {
  if(!v) s<<"null";
  else if(v->constant || v->resolved) number(s,v->constant?v->constant:v->resolved);
  else s<<"{\"binding\":"<<json_string(v->binding)<<",\"definition\":"<<json_string(v->definition)<<'}';
}
void vector(std::ostream& s,const ScalarVector& v) {s<<'[';for(int i=0;i<3;++i){if(i)s<<',';scalar(s,v[i]);}s<<']';}
std::string_view report_kind(ReportDefinition::Kind k) {
  switch(k) {case ReportDefinition::Kind::force:return "force";case ReportDefinition::Kind::drag:return "drag";case ReportDefinition::Kind::lift:return "lift";case ReportDefinition::Kind::downforce:return "downforce";case ReportDefinition::Kind::moment:return "moment";case ReportDefinition::Kind::other:return "other";}
  return "other";
}
}
std::vector<Issue> support_issues(const FluentCase& c,const PolyMesh& mesh,bool seeded) {
  std::vector<Issue> out;
  const auto append=[&](const SettingIssue& issue) {
    auto d=describe(issue);
    if(std::none_of(out.begin(),out.end(),[&](const auto& p){return p.setting==d.setting && p.classification==d.classification;})) out.push_back(std::move(d));
  };
  const auto check=[&](auto fn) {try {fn();}catch(const setup::RejectedSettings& e){for(const auto& i:e.issues()) append(i);}};
  std::vector<setup::Approximation> all;
  for(const auto& a:setup::approximations()) all.push_back(a.id);
  std::vector<setup::Approximation> used;
  check([&]{used=setup::fluent_setup(c,mesh,all).approximations;});
  if(!seeded) check([&]{(void)setup::initial_fields(c,mesh);});
  check([&]{(void)setup::force_definitions(c,mesh);});
  for(const auto& i:c.issues) if(i.classification==IssueClass::informational) append(i);
  // A rejected case still reports WFGC, which needs no other setting to decide.
  if(c.numerics.warped_face_correction==true && std::ranges::find(used,setup::Approximation::wfgc)==used.end())
    used.insert(used.begin(),setup::Approximation::wfgc);
  for(const auto id:used) {
    const auto& a=setup::approximation_info(id);
    out.push_back({std::string(a.setting),std::string(a.label),"approximable",std::string(a.name),
      std::string(a.explanation)+" "+std::string(a.measured_effect)+".",std::string(a.effect_documentation)});
  }
  return out;
}
void write_issues(std::ostream& s,std::span<const Issue> issues) {
  s<<'[';bool first=true;
  for(const auto& i:issues) {
    if(!first)s<<',';first=false;
    s<<"{\"path\":"<<json_string(i.setting)<<",\"setting\":"<<json_string(i.setting)<<",\"label\":"<<json_string(i.label)<<",\"class\":"<<json_string(i.classification)
      <<",\"explanation\":"<<json_string(i.explanation);
    if(!i.approximation.empty())s<<",\"approximation\":"<<json_string(i.approximation);
    if(!i.evidence.empty())s<<",\"evidence\":"<<json_string(i.evidence);
    s<<'}';
  }
  s<<']';
}
void error_event(std::string_view message,std::span<const Issue> issues) {
  std::cout<<"{\"v\":1,\"event\":\"error\",\"message\":"<<json_string(message)<<",\"issues\":";write_issues(std::cout,issues);std::cout<<'}'<<std::endl;
}
int devices(int argc,char** argv) {
  if(argc!=2 && !(argc==3 && std::string_view(argv[2])=="--json")) {std::cerr<<"Usage: kutta devices [--json]\n";return 2;}
  const bool json=argc==3;
  if(json) std::cout<<"{\"v\":1,\"devices\":[";
  bool first=true;
  for(const auto& d:Device::enumerate()) {
    if(json) {
      if(!first)std::cout<<',';first=false;
      std::cout<<"{\"backend\":"<<json_string(backend_name(d.backend))<<",\"name\":"<<json_string(d.name)<<",\"memory_bytes\":"<<d.memory_bytes<<",\"ordinal\":"<<d.ordinal<<",\"default\":"<<(d.is_default?"true":"false")<<'}';
    } else std::cout<<backend_name(d.backend)<<": "<<d.name<<" ("<<d.memory_bytes<<" bytes)"<<(d.is_default?" [default]":"")<<'\n';
  }
  if(json)std::cout<<"]}\n";
  return 0;
}
int inspect_json(int argc,char** argv) {
  std::filesystem::path input,data;bool json=false;
  if(argc<3) return 2;
  input=argv[2];
  for(int i=3;i<argc;++i) {
    const std::string_view arg=argv[i];
    if(arg=="--json" && !json)json=true;
    else if(arg=="--data" && data.empty() && i+1<argc)data=argv[++i];
    else {std::cerr<<"Usage: kutta inspect CASE [--data DAT] --json\n";return 2;}
  }
  try {
    auto c=read_case(input);const auto mesh=read_mesh(input);
    const auto d=data.empty()?std::optional<FluentData>{}:read_data(data,c);
    auto issues=support_issues(c,mesh,d.has_value());
    if(d) try {(void)setup::fluent_fields(c,*d,mesh);}catch(const std::exception& e){issues.push_back(describe({"data/seed",e.what()}));}
    const bool unsupported=std::any_of(issues.begin(),issues.end(),[](const auto& i){return i.classification=="unsupported";});
    const bool approximate=std::any_of(issues.begin(),issues.end(),[](const auto& i){return i.classification=="approximable";});
    std::size_t tet=0,hex=0;const auto cf=build_cell_faces(mesh);
    for(Index cell=0;cell<mesh.n_cells;++cell) {
      const auto begin=cf.offsets[cell],end=cf.offsets[cell+1];bool triangles=true,quads=true;
      for(auto i=begin;i<end;++i){const auto n=detail::nodes(mesh,cf.faces[i]).size();triangles &= n==3;quads &= n==4;}
      tet+=triangles && end-begin==4;hex+=quads && end-begin==6;
    }
    Vec3d lo=mesh.points.front(),hi=lo;
    for(auto p:mesh.points){lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
    auto& s=std::cout;s.imbue(std::locale::classic());s<<std::setprecision(17);
    s<<"{\"v\":1,\"file\":{\"path\":"<<json_string(utf8_path(std::filesystem::absolute(input)))<<",\"sha256\":"<<json_string(sha256_file(input))<<",\"fluent_version\":"<<(c.solver_version?json_string(*c.solver_version):"null")
      <<"},\"mesh\":{\"cells\":"<<mesh.n_cells<<",\"faces\":"<<mesh.n_faces()<<",\"nodes\":"<<mesh.n_points()<<",\"cell_types\":{\"tet\":"<<tet<<",\"hex\":"<<hex<<",\"poly\":"<<mesh.n_cells-tet-hex<<"},\"bounds\":[["<<lo.x<<','<<lo.y<<','<<lo.z<<"],["<<hi.x<<','<<hi.y<<','<<hi.z<<"]]},\"zones\":[";
    bool first=true;
    for(const auto& z:c.boundaries) {
      if(!first)s<<',';first=false;
      auto p=std::find_if(mesh.patches.begin(),mesh.patches.end(),[&](const auto& p){return p.source_zone_id==z.id;});
      const auto range=std::find_if(c.mesh.face_zones.begin(),c.mesh.face_zones.end(),[&](const auto& r){return r.zone_id==z.id;});
      const auto faces=p!=mesh.patches.end()?std::uint64_t(p->size):range!=c.mesh.face_zones.end()?range->size():0;
      s<<"{\"id\":"<<z.id<<",\"name\":"<<json_string(z.name)<<",\"type\":"<<json_string(z.fluent_type)<<",\"faces\":"<<faces<<",\"condition\":{\"velocity\":";vector(s,z.velocity);
      s<<",\"velocity_magnitude\":";scalar(s,z.velocity_magnitude);s<<",\"direction\":";vector(s,z.direction);
      s<<",\"pressure\":";scalar(s,z.gauge_pressure);s<<",\"roughness_height\":";scalar(s,z.roughness_height);
      s<<",\"wall_motion\":"<<(z.wall_motion?json_string(*z.wall_motion==WallMotion::stationary?"stationary":*z.wall_motion==WallMotion::rotating?"rotating":"translating"):"null")<<"}}";
    }
    s<<"],\"models\":{\"viscous\":"<<(c.viscous.model?json_string(*c.viscous.model==ViscousModel::kw_sst?std::string_view("k-omega-sst"):name(*c.viscous.model)):"null")<<",\"energy\":";boolean(s,c.energy);
    bool constant=!c.cell_zones.empty();for(const auto& z:c.cell_zones){auto m=std::find_if(c.materials.begin(),c.materials.end(),[&](const auto& m){return z.material==m.name;});constant &= z.type==ZoneType::fluid && m!=c.materials.end() && m->nonconstant_density==false && m->nonconstant_viscosity==false;}
    s<<",\"energy_coupled\":";if(c.energy==false || constant)s<<"false";else s<<"null";
    s<<"},\"materials\":[";first=true;
    for(const auto& m:c.materials){if(!first)s<<',';first=false;s<<"{\"name\":"<<json_string(m.name)<<",\"density\":";scalar(s,m.density);s<<",\"viscosity\":";scalar(s,m.viscosity);s<<'}';}
    s<<"],\"reports\":[";first=true;
    for(const auto& r:c.reports){if(!first)s<<',';first=false;s<<"{\"name\":"<<json_string(r.name)<<",\"kind\":"<<json_string(report_kind(r.kind))<<",\"zones\":[";if(r.zone_names)for(std::size_t i=0;i<r.zone_names->size();++i){if(i)s<<',';s<<json_string((*r.zone_names)[i]);}s<<"],\"direction\":";if(r.direction)s<<'['<<(*r.direction)[0]<<','<<(*r.direction)[1]<<','<<(*r.direction)[2]<<']';else s<<"null";s<<'}';}
    s<<"],\"support\":{\"status\":"<<json_string(unsupported?"unsupported":approximate?"approximate":"supported")<<",\"issues\":";write_issues(s,issues);
    setup::RunControls controls;
    const auto settings=setup::settings_inventory(c,controls,d.has_value());
    s<<"},\"settings\":"<<dump_json(setup::settings_json(settings));
    s<<",\"data\":{\"present\":"<<(d?"true":"false")<<",\"fields\":[";first=true;
    if(d){std::set<std::string> names;for(const auto* fields:{&d->cell_fields,&d->face_fields})for(const auto& [key,field]:*fields){(void)field;names.insert(key);}for(const auto& n:names){if(!first)s<<',';first=false;s<<json_string(n);}}
    s<<"],\"iterations_hint\":";number(s,d?d->iteration_count:std::optional<std::uint64_t>{});s<<"}}\n";
    return 0;
  }catch(const std::exception& e){error_event(e.what());std::cerr<<"kutta inspect: "<<e.what()<<'\n';return 1;}
}
} // namespace kutta::app
