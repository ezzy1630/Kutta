#include "solve.hpp"
#include "kutta/core/sha256.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/setup/fluent.hpp"
#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <set>

namespace {
using namespace kutta;
using namespace kutta::detail;
using Clock=std::chrono::steady_clock;
struct Options {
  std::filesystem::path input,data,out="kutta-result";
  std::string backend="auto";
  int iterations=1000;
};
Options parse(int argc,char** argv) {
  KUTTA_REQUIRE(argc>=3,"Usage: kutta solve <case.cas.h5> [--data init.dat.h5] [--backend cuda|metal|cpu|auto] [--iterations N] [--out DIR]");
  Options o;o.input=argv[2];std::set<std::string> seen;
  for(int i=3;i<argc;i+=2) {
    const std::string key=argv[i];
    KUTTA_REQUIRE(i+1<argc && seen.insert(key).second,"missing value or duplicate option: "+key);
    const std::string v=argv[i+1];
    if(key=="--data") o.data=v;
    else if(key=="--out") o.out=v;
    else if(key=="--backend") o.backend=v;
    else if(key=="--iterations") {
      const auto [end,ec]=std::from_chars(v.data(),v.data()+v.size(),o.iterations);
      KUTTA_REQUIRE(ec==std::errc{} && end==v.data()+v.size() && o.iterations>0,"--iterations must be a positive integer");
    } else throw Error("unknown solve option: "+key);
  }
  KUTTA_REQUIRE(o.backend=="auto" || o.backend=="cpu" || o.backend=="cuda" || o.backend=="metal","unknown backend: "+o.backend);
  return o;
}
std::unique_ptr<Device> device(const std::string& backend) {
  if(backend=="cpu") return Device::create(Backend::Cpu);
  if(backend=="cuda") return Device::create(Backend::Cuda);
  if(backend=="metal") return Device::create(Backend::Metal);
  for(auto b:{Backend::Cuda,Backend::Metal}) {
    try {return Device::create(b);} catch(const Error& e) {std::cerr<<"auto: "<<e.what()<<'\n';}
  }
  return Device::create(Backend::Cpu);
}
std::ofstream output(const std::filesystem::path& path) {
  std::ofstream f;f.exceptions(std::ios::failbit|std::ios::badbit);f.open(path);
  f.imbue(std::locale::classic());f<<std::setprecision(17);return f;
}
std::string json(const std::string& text) {
  std::ostringstream s;s<<'"';
  for(unsigned char ch:text) {
    if(ch=='"' || ch=='\\') s<<'\\'<<ch;
    else if(ch<32) s<<"\\u"<<std::hex<<std::setfill('0')<<std::setw(4)<<int(ch)<<std::dec;
    else s<<ch;
  }
  return s.str()+'"';
}
std::string csv(const std::string& s) {
  std::string out="\"";for(auto c:s) out+=c=='"'?"\"\"":std::string(1,c);return out+'"';
}
double seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
struct ForceDefinition {
  std::string name;
  std::vector<std::string> patches;
  std::optional<Vec3d> direction;
  bool scaled=false;
};
std::vector<ForceDefinition> definitions(const fluent::FluentCase& c,const PolyMesh& mesh) {
  std::vector<ForceDefinition> out;std::string errors;
  for(const auto& r:c.reports) {
    using Kind=fluent::ReportDefinition::Kind;
    if(r.kind!=Kind::force && r.kind!=Kind::drag && r.kind!=Kind::lift) continue;
    const auto fail=[&](const std::string& why){errors+="\n  report/"+r.name+": "+why;};
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
    out.push_back(std::move(d));
  }
  if(!errors.empty()) throw Error("Force reports rejected:"+errors);
  return out;
}
void write_fields(Device& dev,const fv::SteadySolver& solver,const std::filesystem::path& path) {
  const auto download=[&](const Buffer<float>& b) {
    const auto v=dev.download(b);
    KUTTA_REQUIRE(std::all_of(v.begin(),v.end(),[](float q){return std::isfinite(q);}),"non-finite result field");
    return std::vector<double>(v.begin(),v.end());
  };
  const auto x=download(solver.velocity().x),y=download(solver.velocity().y),z=download(solver.velocity().z);
  std::vector<double> u;u.reserve(3*x.size());
  for(std::size_t i=0;i<x.size();++i) {u.push_back(x[i]);u.push_back(y[i]);u.push_back(z[i]);}
  const auto p=download(solver.pressure()),k=download(solver.kinetic_energy()),omega=download(solver.specific_dissipation()),mu=download(solver.turbulence().mu_t);
  const CellField fields[]={{"u",3,u},{"p",1,p},{"k",1,k},{"omega",1,omega},{"mu_t",1,mu}};
  // Fields and geometry both use the solver's owned RCM order.
  write_vtu(path,solver.mesh(),fields);
}
}
int solve(int argc,char** argv) {
  try {
    const auto o=parse(argc,argv);const auto start=Clock::now();
    const auto case_hash=sha256_file(o.input),data_hash=o.data.empty()?std::string{}:sha256_file(o.data);
    const auto c=fluent::read_case(o.input);auto mesh=read_mesh(o.input);
    std::optional<setup::FluentSetup> setup;std::optional<setup::CellFields> fields;
    std::vector<ForceDefinition> reports;std::string errors;
    const auto check=[&](auto fn) {try {fn();} catch(const std::exception& e) {errors+=std::string(e.what())+'\n';}};
    check([&]{setup=setup::fluent_setup(c,mesh);});
    check([&]{fields=o.data.empty()?setup::initial_fields(c,mesh):setup::fluent_fields(c,fluent::read_data(o.data,c),mesh);});
    check([&]{reports=definitions(c,mesh);});
    if(!errors.empty()) throw Error(errors);
    for(const auto& note:setup->convergence_notes) std::cout<<"Setup: "<<note<<'\n';
    std::cout<<"Case convergence criteria are not exposed by the adapter; using iteration limit "<<o.iterations<<".\n";
    auto dev=device(o.backend);
    fv::SteadySolver solver(*dev,std::move(mesh),setup->boundaries,setup->settings);
    solver.initialize_fields(fields->velocity,fields->pressure,fields->k,fields->omega);
    dev->synchronize();const auto setup_time=seconds(start);
    std::filesystem::create_directories(o.out);
    // A failed replacement run must not leave a previous success summary next
    // to partially replaced outputs. The summary is the completion marker.
    std::filesystem::remove(o.out/"summary.json");
    auto residuals=output(o.out/"residuals.csv");
    residuals<<"iteration,continuity,x-velocity,y-velocity,z-velocity,k,omega,seconds_per_iteration\n";
    std::cout<<"Backend: "<<backend_name(dev->backend())<<"; device: "<<dev->name()<<"\n"
             <<"  iter    continuity    x-velocity    y-velocity    z-velocity             k         omega       sec/iter\n";
    const auto run_start=Clock::now();
    while(solver.report().iteration<o.iterations) {
      const int count=std::min(setup->settings.monitor_interval,o.iterations-solver.report().iteration);
      const auto batch=Clock::now();solver.iterate(count);
      const auto r=solver.report().residual_iteration==solver.report().iteration?solver.report():solver.monitor();
      dev->synchronize();const auto per_iteration=seconds(batch)/count;
      const double values[]={r.continuity.fluent_scaled,r.u.fluent_scaled,r.v.fluent_scaled,r.w.fluent_scaled,r.k.fluent_scaled,r.omega.fluent_scaled};
      residuals<<r.iteration;std::cout<<std::setw(6)<<r.iteration<<std::scientific<<std::setprecision(6);
      for(double v:values) {
        KUTTA_REQUIRE(std::isfinite(v),"non-finite residual at iteration "+std::to_string(r.iteration));
        residuals<<','<<v;std::cout<<std::setw(14)<<v;
      }
      residuals<<','<<per_iteration<<'\n';std::cout<<std::setw(15)<<per_iteration<<std::endl;
    }
    const auto wall_time=seconds(run_start);residuals.close();
    write_fields(*dev,solver,o.out/"result.vtu");
    auto forces=output(o.out/"forces.csv");
    forces<<"kind,name,zone_id,pressure_x,pressure_y,pressure_z,viscous_x,viscous_y,viscous_z,total_x,total_y,total_z,direction_x,direction_y,direction_z,pressure_projection,viscous_projection,total_projection,coefficient\n";
    const auto row=[&](const std::string& kind,const ForceDefinition& d,std::optional<int> zone) {
      const auto f=solver.forces(d.patches,c.reference.pressure.value_or(0));const auto total=f.force();
      forces<<kind<<','<<csv(d.name)<<',';if(zone) forces<<*zone;
      for(auto v:{f.pressure,f.viscous,total}) {
        KUTTA_REQUIRE(std::isfinite(norm(v)),"non-finite force: "+d.name);
        forces<<','<<v.x<<','<<v.y<<','<<v.z;
      }
      if(d.direction) {
        const auto v=*d.direction;forces<<','<<v.x<<','<<v.y<<','<<v.z<<','<<dot(f.pressure,v)<<','<<dot(f.viscous,v)<<','<<dot(total,v)<<',';
        if(d.scaled) forces<<f.coefficient(v,{*c.reference.density,*c.reference.velocity,*c.reference.area,c.reference.pressure.value_or(0)});
      } else forces<<",,,,,,,";
      forces<<'\n';
    };
    for(const auto& p:solver.mesh().patches) if(p.kind==PatchKind::Wall) row("wall",{p.name,{p.name},{},false},p.source_zone_id);
    for(const auto& r:reports) row("report",r,{});
    forces.close();
    auto summary=output(o.out/"summary.json");
    summary<<"{\n  \"case_file\": "<<json(o.input.string())<<",\n  \"case_sha256\": "<<json(case_hash)
      <<",\n  \"data_file\": "<<(o.data.empty()?"null":json(o.data.string()))<<",\n  \"data_sha256\": "<<(o.data.empty()?"null":json(data_hash))
      <<",\n  \"kutta_version\": "<<json(KUTTA_VERSION)<<",\n  \"kutta_commit\": "<<json(KUTTA_COMMIT)
      <<",\n  \"backend\": "<<json(backend_name(dev->backend()))<<",\n  \"device\": "<<json(dev->name())
      <<",\n  \"cells\": "<<solver.mesh().n_cells<<",\n  \"iterations\": "<<solver.report().iteration
      <<",\n  \"wall_time_seconds\": "<<wall_time<<",\n  \"setup_time_seconds\": "<<setup_time
      <<",\n  \"peak_device_memory_bytes\": null,\n  \"converged\": false,\n  \"stop_reason\": \"iteration_limit\",\n  \"convergence_criteria_available\": false,\n  \"boundary_wall_shear_available\": false\n}\n";
    summary.close();std::cout<<"Wrote "<<o.out.string()<<" ("<<solver.report().iteration<<" iterations, "<<wall_time<<" s solve).\n";
    return 0;
  } catch(const std::exception& e) {std::cerr<<"kutta solve: "<<e.what()<<'\n';return 1;}
}
