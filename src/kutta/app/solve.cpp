#include "solve.hpp"
#include "protocol.hpp"
#include "export.hpp"
#include "stop.hpp"
#include "kutta/setup/settings.hpp"
#include "kutta/core/json.hpp"
#include "kutta/core/path.hpp"
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
  std::filesystem::path input,data,setup_file,out="kutta-result";
  std::string backend="auto";
  int iterations=1000;
  std::vector<setup::Approximation> approximations;
};
Options parse(int argc,char** argv) {
  KUTTA_REQUIRE(argc>=3,"Usage: kutta solve <case.cas.h5> [--data init.dat.h5] [--backend cuda|metal|cpu|auto] [--iterations N] [--out DIR] [--approximate SETTING[,SETTING...]]");
  Options o;o.input=argv[2];std::set<std::string> seen;
  for(int i=3;i<argc;i+=2) {
    const std::string key=argv[i];
    KUTTA_REQUIRE(i+1<argc && seen.insert(key).second,"missing value or duplicate option: "+key);
    const std::string v=argv[i+1];
    if(key=="--data") o.data=v;
    else if(key=="--setup") o.setup_file=v;
    else if(key=="--out") o.out=v;
    else if(key=="--backend") o.backend=v;
    else if(key=="--progress") {KUTTA_REQUIRE(v=="jsonl","--progress requires jsonl");}
    else if(key=="--approximate") o.approximations=setup::parse_approximations(v);
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
std::string csv(const std::string& s) {
  std::string out="\"";for(auto c:s) out+=c=='"'?"\"\"":std::string(1,c);return out+'"';
}
double seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
void write_fields(Device& dev,const fv::SteadySolver& solver,const std::filesystem::path& path,std::string_view notice) {
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
  write_vtu(path,solver.mesh(),fields,notice);
}
}
int solve(int argc,char** argv) {
  std::string notice;
  // Discover progress mode even if another option is malformed, so usage
  // errors still produce one machine-readable terminal error.
  bool jsonl=false;
  for(int i=3;i+1<argc;++i) if(std::string_view(argv[i])=="--progress" && std::string_view(argv[i+1])=="jsonl")jsonl=true;
  Options o;
  try {o=parse(argc,argv);}catch(const std::exception& e){if(jsonl)app::error_event(e.what());std::cerr<<e.what()<<'\n';return 2;}
  app::StopRequest stop(jsonl);
  std::vector<app::Issue> issues;
  try {
    const auto start=Clock::now();
    const auto case_hash=sha256_file(o.input),data_hash=o.data.empty()?std::string{}:sha256_file(o.data);
    auto c=fluent::read_case(o.input);
    setup::RunControls controls{o.iterations,20};
    auto inventory=setup::settings_inventory(c,controls,!o.data.empty());
    const auto applied=o.setup_file.empty()?Json::Object{}:setup::apply_setup(o.setup_file,inventory);
    o.iterations=controls.iterations;
    auto mesh=read_mesh(o.input);
    std::optional<setup::FluentSetup> setup;std::optional<setup::Initialization> fields;
    std::vector<setup::ForceDefinition> reports;std::string errors;std::optional<fluent::FluentData> reference;
    const auto check=[&](auto fn) {try {fn();} catch(const std::exception& e) {errors+=std::string(e.what())+'\n';}};
    check([&]{setup=setup::fluent_setup(c,mesh,o.approximations);});
    check([&]{if(!o.data.empty())reference=fluent::read_data(o.data,c);fields=reference?setup::Initialization{setup::fluent_fields(c,*reference,mesh)}:setup::initial_fields(c,mesh);});
    check([&]{reports=setup::force_definitions(c,mesh);});
    if(!errors.empty()) {
      issues=app::support_issues(c,mesh,!o.data.empty());
      throw Error(errors);
    }
    auto& human=jsonl?std::cerr:std::cout;
    const auto warning=[&](std::string_view message) {
      if(jsonl)std::cout<<"{\"v\":1,\"event\":\"warning\",\"message\":"<<json_string(message)<<'}'<<std::endl;
      human<<message<<'\n';
    };
    setup->settings.monitor_interval=controls.monitor_interval;
    notice=setup::approximation_notice(setup->approximations);
    if(!applied.empty()) {
      if(!notice.empty())notice+="; ";
      notice+="MODIFIED RESULT: setup overrides applied relative to the case file; "+dump_json(applied);
    }
    if(!notice.empty() && !jsonl) human<<notice<<std::endl;
    if(!jsonl)for(const auto& note:setup->convergence_notes) human<<"Setup: "<<note<<'\n';
    if(!jsonl)human<<"Case convergence criteria are not exposed by the adapter; using iteration limit "<<o.iterations<<".\n";
    auto dev=device(o.backend);
    fv::SteadySolver solver(*dev,std::move(mesh),setup->boundaries,setup->settings);
    setup::initialize(solver,*fields);
    if(std::holds_alternative<fv::HybridInitialization>(*fields))
      human<<"Initialization: hybrid (stored scalar controls; Kutta AMG replaces Fluent hybrid AMG cycles/methods).\n";
    dev->synchronize();const auto setup_time=seconds(start);
    std::cout.imbue(std::locale::classic());std::cout<<std::setprecision(17);
    if(jsonl) {
      std::cout<<"{\"v\":1,\"event\":\"setup\",\"cells\":"<<solver.mesh().n_cells<<",\"backend\":"<<json_string(backend_name(dev->backend()))
        <<",\"device\":"<<json_string(dev->name())<<",\"approximations\":[";
      for(std::size_t i=0;i<setup->approximations.size();++i){if(i)std::cout<<',';std::cout<<json_string(setup::approximation_info(setup->approximations[i]).name);}
      const auto init=!o.data.empty()?"fluent-data":std::holds_alternative<fv::HybridInitialization>(*fields)?"hybrid":"standard";
      std::cout<<"],\"init\":"<<json_string(init)<<",\"overrides\":"<<dump_json(applied)<<",\"modified\":"<<(applied.empty()?"false":"true")<<'}'<<std::endl;
      if(!notice.empty())warning(notice);
      for(const auto& note:setup->convergence_notes)warning(note);
      warning("Case convergence criteria are not exposed by the adapter; stopping at the requested iteration limit.");
    }
    std::filesystem::create_directories(o.out);
    // A failed replacement run must not leave a previous success summary next
    // to partially replaced outputs. The summary is the completion marker.
    std::filesystem::remove(o.out/"summary.json");
    auto residuals=output(o.out/"residuals.csv");
    if(!notice.empty()) residuals<<"# "<<notice<<'\n';
    residuals<<"iteration,continuity,x-velocity,y-velocity,z-velocity,k,omega,seconds_per_iteration\n";
    if(!jsonl)human<<"Backend: "<<backend_name(dev->backend())<<"; device: "<<dev->name()<<"\n"
             <<"  iter    continuity    x-velocity    y-velocity    z-velocity             k         omega       sec/iter\n";
    const auto run_start=Clock::now();
    bool stopped=false;
    auto batch=Clock::now();int batch_count=0;
    while(solver.report().iteration<o.iterations) {
      if(stop.requested()){stopped=true;if(batch_count==0)break;}
      else {
        // One outer iteration is the cancellation boundary on every backend.
        solver.iterate(1);++batch_count;stopped=stop.requested();
      }
      const auto iter=solver.report().iteration;
      if(iter%setup->settings.monitor_interval!=0 && iter!=o.iterations && !stopped)continue;
      const auto r=solver.report().residual_iteration==iter?solver.report():solver.monitor();
      dev->synchronize();const auto per_iteration=seconds(batch)/batch_count;
      const double values[]={r.continuity.fluent_scaled,r.u.fluent_scaled,r.v.fluent_scaled,r.w.fluent_scaled,r.k.fluent_scaled,r.omega.fluent_scaled};
      residuals<<r.iteration;
      if(!jsonl)human<<std::setw(6)<<r.iteration<<std::scientific<<std::setprecision(6);
      for(double v:values) {
        KUTTA_REQUIRE(std::isfinite(v),"non-finite residual at iteration "+std::to_string(r.iteration));
        residuals<<','<<v;if(!jsonl)human<<std::setw(14)<<v;
      }
      residuals<<','<<per_iteration<<'\n';
      if(!jsonl)human<<std::setw(15)<<per_iteration<<std::endl;
      else {
        std::ostringstream event;event.imbue(std::locale::classic());event<<std::setprecision(17);
        event<<"{\"v\":1,\"event\":\"iteration\",\"iter\":"<<r.iteration<<",\"residuals\":{";
        const char* names[]={"continuity","x-velocity","y-velocity","z-velocity","k","omega"};
        for(int i=0;i<6;++i){if(i)event<<',';event<<json_string(names[i])<<':'<<values[i];}
        event<<"},\"reports\":{";
        for(std::size_t i=0;i<reports.size();++i) {
          const auto& d=reports[i];if(i)event<<',';
          const auto f=solver.forces(d.patches,c.reference.pressure.value_or(0));
          const auto force=dot(f.force(),*d.direction);KUTTA_REQUIRE(std::isfinite(force),"non-finite report: "+d.name);
          event<<json_string(d.name)<<":{\"force\":"<<force;
          if(c.reference.density && *c.reference.density>0 && c.reference.velocity && *c.reference.velocity>0 && c.reference.area && *c.reference.area>0) {
            const auto coefficient=f.coefficient(*d.direction,{*c.reference.density,*c.reference.velocity,*c.reference.area,c.reference.pressure.value_or(0)});
            KUTTA_REQUIRE(std::isfinite(coefficient),"non-finite report coefficient: "+d.name);
            event<<",\"coefficient\":"<<coefficient;
          }
          event<<'}';
        }
        event<<"},\"sec_per_iter\":"<<per_iteration<<",\"elapsed\":"<<seconds(run_start)<<'}';
        std::cout<<event.str()<<std::endl;
      }
      if(stopped)break;
      batch=Clock::now();batch_count=0;
    }
    const auto wall_time=seconds(run_start);residuals.close();
    write_fields(*dev,solver,o.out/"result.vtu",notice);
    auto forces=output(o.out/"forces.csv");
    if(!notice.empty()) forces<<"# "<<notice<<'\n';
    forces<<"kind,name,zone_id,pressure_x,pressure_y,pressure_z,viscous_x,viscous_y,viscous_z,total_x,total_y,total_z,direction_x,direction_y,direction_z,pressure_projection,viscous_projection,total_projection,coefficient\n";
    const auto row=[&](const std::string& kind,const setup::ForceDefinition& d,std::optional<int> zone) {
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
    auto fluxes=output(o.out/"fluxes.csv");
    if(!notice.empty())fluxes<<"# "<<notice<<'\n';
    fluxes<<"name,zone_id,mass_flow_outward\n";
    const auto mass_flux=dev->download(solver.mass_flux());
    for(const auto& patch:solver.mesh().patches) {
      double total=0;for(Index f=patch.start;f<patch.start+patch.size;++f)total+=mass_flux[f];
      KUTTA_REQUIRE(std::isfinite(total),"non-finite boundary mass flux");
      fluxes<<csv(patch.name)<<','<<patch.source_zone_id<<','<<total<<'\n';
    }
    fluxes.close();
    app::save_view_input(o.out/"view-input.bin",*dev,solver,c,reference?&*reference:nullptr,*setup,notice);
    auto summary=output(o.out/"summary.json");
    summary<<"{\n  \"v\": 1,\n  \"modified\": "<<(applied.empty()?"false":"true")<<",\n  \"overrides\": "<<dump_json(applied)<<",\n  \"approximate\": "<<(setup->approximations.empty()?"false":"true")<<",\n  \"approximation_notice\": "<<json_string(setup::approximation_notice(setup->approximations))
      <<",\n  \"approximations\": [";
    for(std::size_t i=0;i<setup->approximations.size();++i) {
      const auto& info=setup::approximation_info(setup->approximations[i]);
      if(i) summary<<',';
      summary<<"\n    {\"name\": "<<json_string(info.name)<<", \"setting\": "<<json_string(info.setting)
        <<", \"fluent_value\": "<<json_string(info.fluent_value)<<", \"kutta_action\": "<<json_string(info.kutta_action)
        <<", \"effect_documentation\": "<<json_string(info.effect_documentation)<<", \"measured_effect\": "<<json_string(info.measured_effect)<<'}';
    }
    summary<<"],\n  \"case_file\": "<<json_string(o.input.string())<<",\n  \"case_sha256\": "<<json_string(case_hash)
      <<",\n  \"data_file\": "<<(o.data.empty()?"null":json_string(o.data.string()))<<",\n  \"data_sha256\": "<<(o.data.empty()?"null":json_string(data_hash))
      <<",\n  \"kutta_version\": "<<json_string(KUTTA_VERSION)<<",\n  \"kutta_commit\": "<<json_string(KUTTA_COMMIT)
      <<",\n  \"backend\": "<<json_string(backend_name(dev->backend()))<<",\n  \"device\": "<<json_string(dev->name())
      <<",\n  \"cells\": "<<solver.mesh().n_cells<<",\n  \"iterations\": "<<solver.report().iteration
      <<",\n  \"wall_time_seconds\": "<<wall_time<<",\n  \"setup_time_seconds\": "<<setup_time
      <<",\n  \"peak_device_memory_bytes\": null,\n  \"converged\": false,\n  \"stop_reason\": "<<json_string(stopped?"stopped":"iteration_limit")<<",\n  \"convergence_criteria_available\": false,\n  \"boundary_wall_shear_available\": true\n}\n";
    summary.close();
    if(jsonl)std::cout<<"{\"v\":1,\"event\":\"done\",\"converged\":false,\"iterations\":"<<solver.report().iteration<<",\"elapsed\":"<<wall_time<<",\"out\":"<<json_string(utf8_path(o.out))<<",\"stop_reason\":"<<json_string(stopped?"stopped":"iteration_limit")<<'}'<<std::endl;
    else human<<"Wrote "<<o.out.string()<<" ("<<solver.report().iteration<<" iterations, "<<wall_time<<" s solve).\n";
    if(!notice.empty() && !jsonl) human<<notice<<std::endl;
    return 0;
  } catch(const setup::RejectedSettings& e) {
    for(const auto& i:e.issues())issues.push_back({i.path,i.path,"unsupported",{},i.reason,{}});
    if(jsonl)app::error_event(e.what(),issues);
    std::cerr<<"kutta solve: "<<e.what()<<'\n';return 1;
  } catch(const std::exception& e) {
    if(jsonl)app::error_event(e.what(),issues);
    std::cerr<<"kutta solve: "<<e.what()<<'\n';
    if(!notice.empty()) std::cerr<<notice<<'\n';
    return 1;
  }
}
