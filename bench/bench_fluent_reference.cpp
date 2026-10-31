// Compare kutta solve's ASCII VTU with the matching Fluent saved fields without
// advancing either solution. Recreate the solver's deterministic RCM mapping;
// filenames and physical nearest-neighbour matching are not cell identities.
#include "kutta/core/sha256.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/setup/fluent.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

using namespace kutta;
using namespace kutta::detail;
namespace {
std::vector<double> field(const std::string& xml,const std::string& name,std::size_t size) {
  const auto cell=xml.find("<CellData>");
  KUTTA_REQUIRE(cell!=std::string::npos,"VTU lacks CellData");
  const auto start=xml.find("Name=\""+name+"\"",cell);
  KUTTA_REQUIRE(start!=std::string::npos,"VTU lacks field "+name);
  const auto begin=xml.find('>',start),end=xml.find("</DataArray>",begin);
  KUTTA_REQUIRE(begin!=std::string::npos && end!=std::string::npos,"malformed VTU field "+name);
  std::istringstream input(xml.substr(begin+1,end-begin-1));
  std::vector<double> out;double v;
  while(input>>v) {KUTTA_REQUIRE(std::isfinite(v),"nonfinite VTU field");out.push_back(v);}
  KUTTA_REQUIRE(input.eof() && out.size()==size,"VTU field size/encoding mismatch: "+name);
  return out;
}
double l2(std::span<const double> a,std::span<const double> b) {
  KUTTA_REQUIRE(a.size()==b.size(),"field sizes differ");double n=0,d=0;
  for(std::size_t i=0;i<a.size();++i) {n+=(a[i]-b[i])*(a[i]-b[i]);d+=b[i]*b[i];}
  KUTTA_REQUIRE(d>0,"zero reference L2 norm");return std::sqrt(n/d);
}
struct Zone {double pressure_sum=0,area=0;Vec3d pressure_force,shear_force;};
std::map<std::string,Zone> reference_zones(const fluent::FluentCase& c,const fluent::FluentData& d,const PolyMesh& mesh) {
  const auto g=compute_geometry(mesh);std::map<std::string,Zone> out;
  for(const auto& patch:mesh.patches) {
    if(patch.source_zone_id<0) continue;  // a 2D case's extrusion caps have no Fluent zone
    const auto z=std::find_if(c.mesh.face_zones.begin(),c.mesh.face_zones.end(),[&](const auto& z){return z.zone_id==patch.source_zone_id;});
    KUTTA_REQUIRE(z!=c.mesh.face_zones.end() && z->size()==std::uint64_t(patch.size),"source face identity mismatch");
    auto& r=out[patch.name];
    for(Index j=0;j<patch.size;++j) {
      const auto f=patch.start+j;const auto id=z->first+j;const double p=d.face_fields.at("SV_P").at(id)[0],a=norm(g.face_area[f]);
      r.area+=a;r.pressure_sum+=p*a;r.pressure_force=r.pressure_force+g.face_area[f]*(p-c.reference.pressure.value_or(0));
      if(patch.kind==PatchKind::Wall) {const auto s=d.face_fields.at("SV_WALL_SHEAR").at(id);KUTTA_REQUIRE(s.size()==3,"wall shear must have 3 components");r.shear_force=r.shear_force-Vec3d{s[0],s[1],s[2]};}
    }
  }
  return out;
}
void vector(const char* kind,const std::string& name,Vec3d a,Vec3d b) {
  std::cout<<"REFERENCE_FORCE "<<name<<' '<<kind<<" ours="<<a.x<<','<<a.y<<','<<a.z<<" Fluent="<<b.x<<','<<b.y<<','<<b.z;
  if(norm(b)>0) std::cout<<" relative="<<norm(a-b)/norm(b);else std::cout<<" absolute="<<norm(a-b);
  std::cout<<'\n';
}
}
int main(int argc,char** argv) {
  try {
    KUTTA_REQUIRE(argc==5 || argc==6,"Usage: bench_fluent_reference case.cas.h5 data.dat.h5 result.vtu cpu|metal [approximations]");
    const auto case_hash=sha256_file(argv[1]);
    const auto summary_path=std::filesystem::path(argv[3]).parent_path()/"summary.json";
    std::ifstream summary_input(summary_path);KUTTA_REQUIRE(summary_input,"result needs adjacent solve summary.json provenance");
    const std::string summary{std::istreambuf_iterator<char>(summary_input),{}};
    KUTTA_REQUIRE(summary.find("\"case_sha256\": \""+case_hash+"\"")!=std::string::npos,"result summary belongs to a different case");
    const auto c=fluent::read_case(argv[1]);const auto d=fluent::read_data(argv[2],c);const auto mesh=read_mesh(argv[1]);
    const auto ref=setup::fluent_fields(c,d,mesh);const auto zones=reference_zones(c,d,mesh);
    // The solve's --approximate list; the V6 elbow needs wfgc.
    const auto allowed=setup::parse_approximations(argc==6?argv[5]:"wfgc");const auto problem=setup::fluent_setup(c,mesh,allowed);
    KUTTA_REQUIRE(std::string_view(argv[4])=="cpu" || std::string_view(argv[4])=="metal","expected cpu|metal");
    auto dev=Device::create(std::string_view(argv[4])=="cpu"?Backend::Cpu:Backend::Metal);
    fv::SteadySolver solver(*dev,mesh,problem.boundaries,problem.settings);
    std::ifstream in(argv[3]);KUTTA_REQUIRE(in,"cannot read result VTU");const std::string xml{std::istreambuf_iterator<char>(in),{}};
    // Closure fields follow the case model: SST k/omega, standard k-epsilon
    // k/epsilon, laminar none (fluent_fields fills the same vectors).
    const std::size_t n=mesh.n_cells;auto u=field(xml,"u",3*n),p=field(xml,"p",n);
    const bool turbulent=!ref.k.empty(),ke=!ref.epsilon.empty();
    std::vector<double> k,second;
    if(turbulent) {k=field(xml,"k",n);second=field(xml,ke?"epsilon":"omega",n);}
    auto fields=ref;auto& reference_second=ke?ref.epsilon:ref.omega;auto& fields_second=ke?fields.epsilon:fields.omega;
    for(std::size_t i=0;i<n;++i) {
      const auto old=solver.permutation().cells_new_to_old[i];fields.velocity[old]={u[3*i],u[3*i+1],u[3*i+2]};fields.pressure[old]=p[i];
      if(turbulent) {fields.k[old]=k[i];fields_second[old]=second[i];}
    }
    std::vector<double> actual_velocity,reference_velocity;
    for(std::size_t i=0;i<n;++i) {
      const auto a=fields.velocity[i],b=ref.velocity[i];
      actual_velocity.insert(actual_velocity.end(),{a.x,a.y,a.z});
      reference_velocity.insert(reference_velocity.end(),{b.x,b.y,b.z});
    }
    std::cout<<std::setprecision(12)<<"REFERENCE_CASE sha256="<<sha256_file(argv[1])<<" data_sha256="<<sha256_file(argv[2])<<" result_sha256="<<sha256_file(argv[3])<<" backend="<<argv[4]<<" cells="<<n<<'\n';
    // How converged the saved reference is: Fluent's own scaled residual.
    if(const auto r=d.residuals.find("continuity");r!=d.residuals.end() && !r->second.rows.empty()) {
      const auto& last=r->second.rows.back();
      std::cout<<"REFERENCE_SAVED iteration="<<r->second.iterations.back()<<" continuity_scaled="<<last[0]/last[1]<<'\n';
    }
    // Saved data can outlive an edited case (same data file, different material).
    if(const auto rho=d.cell_fields.find("SV_DENSITY");rho!=d.cell_fields.end()) {
      double worst=0;
      for(const auto& section:rho->second.sections) for(const double v:section.values) worst=std::max(worst,std::abs(v/problem.settings.density-1));
      if(worst>1e-6) std::cout<<"REFERENCE_MISMATCH saved density differs from the case material by up to "<<worst<<" relative\n";
    }
    std::vector<std::pair<const char*,double>> errors{{"velocity",l2(actual_velocity,reference_velocity)},{"pressure",l2(fields.pressure,ref.pressure)}};
    if(turbulent) {errors.push_back({"k",l2(fields.k,ref.k)});errors.push_back({ke?"epsilon":"omega",l2(fields_second,reference_second)});}
    for(const auto& [name,error]:errors) std::cout<<"REFERENCE_FIELD "<<name<<" relative_L2="<<error<<'\n';
    if(ke) solver.initialize_ke_fields(fields.velocity,fields.pressure,fields.k,fields.epsilon);
    else solver.initialize_fields(fields.velocity,fields.pressure,fields.k,fields.omega);
    // Initialization rebuilds boundary reconstruction at these final fields.
    const auto bp=dev->download(solver.boundary_pressure());const auto& g=solver.geometry();
    double inlet=0,inlet_ref=0,ia=0,outlet=0,outlet_ref=0,oa=0;
    std::map<std::string,std::pair<Vec3d,Vec3d>> wall_forces;  // ours, Fluent
    for(const auto& patch:solver.mesh().patches) {
      if(!zones.contains(patch.name)) continue;
      const auto& r=zones.at(patch.name);double ps=0,a=0;
      for(Index f=patch.start;f<patch.start+patch.size;++f) {const double af=norm(g.face_area[f]);a+=af;ps+=af*bp[f-solver.mesh().n_internal_faces()];}
      std::cout<<"REFERENCE_ZONE "<<patch.name<<" ours="<<ps/a<<" Fluent="<<r.pressure_sum/r.area<<'\n';
      if(patch.kind==PatchKind::VelocityInlet) {inlet+=ps;inlet_ref+=r.pressure_sum;ia+=a;}
      if(patch.kind==PatchKind::PressureOutlet) {outlet+=ps;outlet_ref+=r.pressure_sum;oa+=a;}
      if(patch.kind==PatchKind::Wall) {const std::vector<std::string> names{patch.name};const auto f=solver.forces(names,c.reference.pressure.value_or(0));vector("pressure",patch.name,f.pressure,r.pressure_force);vector("shear",patch.name,f.viscous,r.shear_force);vector("total",patch.name,f.force(),r.pressure_force+r.shear_force);
        wall_forces[patch.name]={f.force(),r.pressure_force+r.shear_force};}
    }
    // Faces where flow enters through a pressure outlet: the only faces that
    // reverse-flow suppression acts on. Both fluxes are outward positive.
    const auto mass=dev->download(solver.mass_flux());
    for(const auto& patch:solver.mesh().patches) if(patch.kind==PatchKind::PressureOutlet) {
      const auto z=std::find_if(c.mesh.face_zones.begin(),c.mesh.face_zones.end(),[&](const auto& z){return z.zone_id==patch.source_zone_id;});
      Index ours=0,saved=0;
      for(Index j=0;j<patch.size;++j) {
        if(mass[patch.start+j]<0) ++ours;
        if(d.face_fields.at("SV_FLUX").at(z->first+j)[0]<0) ++saved;
      }
      std::cout<<"REFERENCE_BACKFLOW "<<patch.name<<" ours="<<ours<<" Fluent="<<saved<<" faces="<<patch.size<<'\n';
    }
    if(ia>0 && oa>0) {const double drop=inlet/ia-outlet/oa,reference=inlet_ref/ia-outlet_ref/oa;std::cout<<"REFERENCE_PRESSURE_DROP ours="<<drop<<" Fluent="<<reference<<" relative="<<std::abs(drop/reference-1)<<'\n';}
    for(const auto& r:c.reports) {
      std::cout<<"REFERENCE_REPORT "<<r.name<<" type="<<r.fluent_type<<" scaled="<<(r.scaled?(*r.scaled?"true":"false"):"absent")<<'\n';
      // Raw report definitions remain in inspect output. The solver's forces.csv
      // is authoritative for evaluated projections and coefficient reference values.
    }
    // Force reports evaluated from both sides' wall forces with the case's
    // reference values (Fluent's coefficient: F.d / (0.5 rho U^2 A)).
    for(const auto& d:setup::force_definitions(c,mesh)) {
      if(!d.direction) continue;
      Vec3d ours,saved;
      for(const auto& name:d.patches) {const auto& f=wall_forces.at(name);ours=ours+f.first;saved=saved+f.second;}
      const auto direction=*d.direction/norm(*d.direction);
      double a=dot(ours,direction),b=dot(saved,direction);
      if(d.scaled) {const double q=.5**c.reference.density**c.reference.velocity**c.reference.velocity**c.reference.area;a/=q;b/=q;}
      std::cout<<"REFERENCE_COEFFICIENT "<<d.name<<" scaled="<<(d.scaled?"true":"false")<<" ours="<<a<<" Fluent="<<b
        <<" relative="<<std::abs(a/b-1)<<'\n';
    }
    std::cout<<setup::approximation_notice(problem.approximations)<<'\n';return 0;
  } catch(const std::exception& e) {std::cerr<<"Fluent comparison: "<<e.what()<<'\n';return 1;}
}
