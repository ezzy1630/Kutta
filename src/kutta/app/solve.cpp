#include "solve.hpp"
#include "protocol.hpp"
#include "export.hpp"
#include "stop.hpp"
#include "live.hpp"
#include "kutta/io/live_view.hpp"
#include "kutta/setup/settings.hpp"
#include "kutta/post/aero.hpp"
#include "kutta/core/json.hpp"
#include "kutta/core/memory.hpp"
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
#include <map>
#include <optional>
#include <sstream>
#include <set>

namespace {
using namespace kutta;
using namespace kutta::detail;
using Clock=std::chrono::steady_clock;
struct Options {
  std::filesystem::path input,data,setup_file,live_directory,out="kutta-result";
  std::string backend="auto";
  int iterations=1000,live_every=app::LiveCadence::default_every;
  std::vector<setup::Approximation> approximations;
};
Options parse(int argc,char** argv) {
  KUTTA_REQUIRE(argc>=3,"Usage: kutta solve <case.cas.h5> [--data init.dat.h5] [--backend cuda|metal|cpu|auto] [--iterations N] [--out DIR] [--live-view DIR] [--live-every N] [--approximate SETTING[,SETTING...]]");
  Options o;o.input=argv[2];std::set<std::string> seen;
  for(int i=3;i<argc;i+=2) {
    const std::string key=argv[i];
    KUTTA_REQUIRE(i+1<argc && seen.insert(key).second,"missing value or duplicate option: "+key);
    const std::string v=argv[i+1];
    if(key=="--data") o.data=v;
    else if(key=="--setup") o.setup_file=v;
    else if(key=="--out") o.out=v;
    else if(key=="--live-view") {KUTTA_REQUIRE(!v.empty(),"--live-view requires a directory");o.live_directory=v;}
    else if(key=="--backend") o.backend=v;
    else if(key=="--progress") {KUTTA_REQUIRE(v=="jsonl","--progress requires jsonl");}
    else if(key=="--approximate") o.approximations=setup::parse_approximations(v);
    else if(key=="--iterations" || key=="--live-every") {
      auto& value=key=="--iterations"?o.iterations:o.live_every;
      const auto [end,ec]=std::from_chars(v.data(),v.data()+v.size(),value);
      KUTTA_REQUIRE(ec==std::errc{} && end==v.data()+v.size() && value>0,key+" must be a positive integer");
    } else throw Error("unknown solve option: "+key);
  }
  KUTTA_REQUIRE(!seen.contains("--live-every") || !o.live_directory.empty(),"--live-every requires --live-view");
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
void write_fields(Device& dev,const fv::SteadySolver& solver,const fv::SolverSettings& settings,const std::filesystem::path& path,std::string_view notice) {
  // Each float download goes straight into the double array the writer reads,
  // one component at a time, so the export holds no second copy of a field.
  const auto download=[&](const Buffer<float>& b,std::vector<double>& to,std::size_t components=1,std::size_t component=0) {
    const auto v=dev.download(b);
    KUTTA_REQUIRE(std::all_of(v.begin(),v.end(),[](float q){return std::isfinite(q);}),"non-finite result field");
    to.resize(components*v.size());
    for(std::size_t i=0;i<v.size();++i) to[i*components+component]=v[i];
  };
  const auto field=[&](const Buffer<float>& b) {std::vector<double> to;download(b,to);return to;};
  std::vector<double> u;
  download(solver.velocity().x,u,3,0);download(solver.velocity().y,u,3,1);download(solver.velocity().z,u,3,2);
  const auto p=field(solver.pressure());
  std::vector<CellField> fields{{"u",3,u},{"p",1,p}};
  // Turbulence fields exist only for the active closure; laminar has none.
  std::vector<double> k,omega,mu;
  if(settings.sst) {
    k=field(solver.kinetic_energy());omega=field(solver.specific_dissipation());mu=field(solver.turbulence().mu_t);
    fields.insert(fields.end(),{{"k",1,k},{"omega",1,omega},{"mu_t",1,mu}});
  } else if(settings.realizable_ke) {
    k=field(solver.kinetic_energy());omega=field(solver.dissipation());mu=field(solver.realizable_turbulence().mu_t);
    fields.insert(fields.end(),{{"k",1,k},{"epsilon",1,omega},{"mu_t",1,mu}});
  }
  std::vector<double> gamma;
  if(settings.sst && settings.sst->transition) {gamma=field(solver.intermittency());fields.push_back({"intermittency",1,gamma});}
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
    memory::phase("start");
    const auto case_hash=sha256_file(o.input),data_hash=o.data.empty()?std::string{}:sha256_file(o.data);
    auto c=fluent::read_case(o.input);
    memory::phase("read-case");
    setup::RunControls controls{o.iterations,20};
    auto inventory=setup::settings_inventory(c,controls,!o.data.empty());
    const auto applied=o.setup_file.empty()?Json::Object{}:setup::apply_setup(o.setup_file,inventory);
    // Only case-changing overrides make a result MODIFIED; run controls are
    // still applied and echoed.
    Json::Object changes,run_controls;for(const auto& [path,v]:applied)(setup::run_control(path)?run_controls:changes).emplace(path,v);
    if(auto missing=setup::criteria_issues(controls.convergence);!missing.empty())throw setup::RejectedSettings(std::move(missing));
    const auto& convergence=controls.convergence;
    std::string convergence_warning;
    if(!convergence.check)convergence_warning="Convergence checking is off (monitors/check-convergence); stopping at the requested iteration limit.";
    else if(!convergence.has_conditions())convergence_warning="Convergence checking is on but no residual is checked and no report condition is active; stopping at the requested iteration limit.";
    const bool report_conditions=convergence.check && std::any_of(convergence.reports.begin(),convergence.reports.end(),[](const auto& r){return r.active;});
    o.iterations=controls.iterations;
    auto mesh=read_mesh(o.input);
    memory::phase("read-mesh");
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
    if(!changes.empty()) {
      if(!notice.empty())notice+="; ";
      notice+="MODIFIED RESULT: setup overrides applied relative to the case file; "+dump_json(changes);
    }
    if(!notice.empty() && !jsonl) human<<notice<<std::endl;
    if(!jsonl)for(const auto& note:setup->convergence_notes) human<<"Setup: "<<note<<'\n';
    if(!jsonl && !convergence_warning.empty())human<<convergence_warning<<'\n';
    using Coupling=fv::PressureVelocityCoupling;const auto scheme=setup->settings.coupling;
    const auto coupling=scheme==Coupling::Coupled?"coupled":scheme==Coupling::Simple?"simple":"simplec";
    const auto coupling_label=scheme==Coupling::Coupled?"Coupled":scheme==Coupling::Simple?"SIMPLE":"SIMPLEC";
    if(!jsonl)human<<"Pressure-velocity coupling: "<<coupling_label<<".\n";
    memory::phase("case-setup");
    auto dev=device(o.backend);
    memory::phase("device");
    fv::SteadySolver solver(*dev,std::move(mesh),setup->boundaries,setup->settings);
    setup::initialize(solver,*fields);
    if(std::holds_alternative<fv::HybridInitialization>(*fields))
      human<<"Initialization: hybrid (stored scalar controls; Kutta AMG replaces Fluent hybrid AMG cycles/methods).\n";
    dev->synchronize();
    memory::phase("initialize");
    std::optional<view::LiveView> live;
    if(!o.live_directory.empty())live.emplace(o.live_directory);
    app::LiveCadence cadence(o.live_every);
    const auto live_configuration=[&] {return Json::Object{{"directory",utf8_path(live->directory())},{"run_id",live->run_id()},{"every",cadence.every()},{"unit","iterations"},{"budget_fraction",app::LiveCadence::budget_fraction}};};
    const auto live_statistics=[&](double wall) {auto info=live_configuration();const auto compute=std::max(0.0,wall-cadence.total_seconds());info.emplace("snapshots",double(cadence.count()));info.emplace("snapshot_time_seconds",cadence.total_seconds());info.emplace("solve_time_seconds",compute);info.emplace("overhead_fraction",compute>0?Json(cadence.total_seconds()/compute):Json(nullptr));return info;};
    const auto setup_time=seconds(start);
    std::cout.imbue(std::locale::classic());std::cout<<std::setprecision(17);
    if(jsonl) {
      std::cout<<"{\"v\":1,\"event\":\"setup\",\"cells\":"<<solver.mesh().n_cells<<",\"backend\":"<<json_string(backend_name(dev->backend()))
        <<",\"device\":"<<json_string(dev->name())<<",\"approximations\":[";
      for(std::size_t i=0;i<setup->approximations.size();++i){if(i)std::cout<<',';std::cout<<json_string(setup::approximation_info(setup->approximations[i]).name);}
      const auto init=!o.data.empty()?"fluent-data":std::holds_alternative<fv::HybridInitialization>(*fields)?"hybrid":"standard";
      std::cout<<"],\"init\":"<<json_string(init)<<",\"overrides\":"<<dump_json(applied)<<",\"run_controls\":"<<dump_json(run_controls)<<",\"modified\":"<<(changes.empty()?"false":"true")
        <<",\"pressure_velocity_coupling\":"<<json_string(coupling)<<",\"convergence\":"<<dump_json(convergence.json());
      if(live)std::cout<<",\"live_view\":"<<dump_json(live_configuration());
      std::cout<<'}'<<std::endl;
      if(!notice.empty())warning(notice);
      for(const auto& note:setup->convergence_notes)warning(note);
      if(!convergence_warning.empty())warning(convergence_warning);
    }
    std::filesystem::create_directories(o.out);
    // A failed replacement run must not leave a previous success summary next
    // to partially replaced outputs. The summary is the completion marker.
    std::filesystem::remove(o.out/"summary.json");
    // Loads and aero results derived from a previous run would describe a
    // different solution.
    std::filesystem::remove(o.out/"loads.json");std::filesystem::remove(o.out/"aero.json");
    auto residuals=output(o.out/"residuals.csv");
    if(!notice.empty()) residuals<<"# "<<notice<<'\n';
    const bool transition=setup->settings.sst && setup->settings.sst->transition;
    const bool planar=c.dimension==2;
    const auto residual_names=setup::solved_residuals(c.viscous.model,transition,planar);
    residuals<<"iteration";for(auto name:residual_names)residuals<<','<<name;residuals<<",seconds_per_iteration\n";
    if(!jsonl) {
      human<<"Backend: "<<backend_name(dev->backend())<<"; device: "<<dev->name()<<"\n";
      if(solver.poor_cells())human<<"Poor mesh numerics: "<<solver.poor_cells()<<" cells at orthogonal quality <= "<<setup->settings.poor_cell_quality<<" use low-order face pressures\n";
      human<<"  iter";
      for(auto name:residual_names)human<<std::setw(14)<<name;
      human<<"       sec/iter\n";
    }
    const auto run_start=Clock::now();
    bool stopped=false,converged=false;
    setup::ConvergenceMonitor monitor(convergence);std::vector<std::string> met_by;
    auto batch=Clock::now();int batch_count=0;double batch_snapshot_time=0;
    while(solver.report().iteration<o.iterations) {
      if(stop.requested()){stopped=true;if(batch_count==0)break;}
      else {
        // One outer iteration is the cancellation boundary on every backend.
        solver.iterate(1);++batch_count;stopped=stop.requested();
        if(solver.report().iteration==1)memory::phase("first-iteration");
      }
      const auto iter=solver.report().iteration;
      const bool report_due=iter%setup->settings.monitor_interval==0 || iter==o.iterations || stopped;
      // Between monitors the boundary-pressure buffer is pressure-correction
      // scratch. Reuse the solver's scheduled refresh; never export that scratch
      // or add extra monitor work that could change the solve for live view.
      if(report_due && solver.report().residual_iteration!=iter)solver.monitor();
      if(live && solver.report().residual_iteration==iter && cadence.due(iter,seconds(run_start)-cadence.total_seconds())) {
        const auto snapshot_start=Clock::now();
        const auto snapshot_fields=app::collect_view_fields(*dev,solver,c,reference?&*reference:nullptr,*setup,notice,false);
        // Live view is a preview: a failed snapshot is reported and skipped,
        // never allowed to end the solve.
        std::optional<Json::Object> event;std::string failure;
        try {event=live->publish(solver.mesh(),snapshot_fields,iter,seconds(run_start));}
        catch(const std::exception& e) {failure=e.what();}
        const auto cost=seconds(snapshot_start);cadence.record(iter,cost);batch_snapshot_time+=cost;
        if(event && jsonl){event->emplace("event","snapshot");event->emplace("snapshot_seconds",cost);std::cout<<dump_json(*event)<<std::endl;}
        if(!failure.empty()) {
          const auto message="live snapshot at iteration "+std::to_string(iter)+" skipped: "+failure;
          if(jsonl)std::cout<<dump_json(Json::Object{{"v",1},{"event","warning"},{"message",message}})<<std::endl;
          std::cerr<<message<<std::endl;
        }
        stopped=stop.requested();
      }
      // Iterations 1-5 are already monitored (Fluent normalizes continuity by
      // their largest value), so their residuals are reported at no solver
      // cost, as Fluent's residual history starts at iteration 1. They are not
      // monitor samples: no force reports and no convergence evaluation.
      const bool sample=iter%setup->settings.monitor_interval==0 || iter==o.iterations || stopped;
      const bool early=!sample && iter<=5 && solver.report().residual_iteration==iter;
      if(!sample && !early)continue;
      const auto r=solver.report().residual_iteration==iter?solver.report():solver.monitor();
      dev->synchronize();const auto per_iteration=(seconds(batch)-batch_snapshot_time)/batch_count;
      std::vector<double> values{r.continuity.fluent_scaled,r.u.fluent_scaled,r.v.fluent_scaled};
      if(!planar) values.push_back(r.w.fluent_scaled);
      if(setup->settings.sst) values.insert(values.end(),{r.k.fluent_scaled,r.omega.fluent_scaled});
      if(setup->settings.realizable_ke) values.insert(values.end(),{r.k.fluent_scaled,r.epsilon.fluent_scaled});
      if(transition) values.push_back(r.intermittency.fluent_scaled);
      KUTTA_REQUIRE(values.size()==residual_names.size(),"residual names do not match the solved equations");
      residuals<<r.iteration;
      if(!jsonl)human<<std::setw(6)<<r.iteration<<std::scientific<<std::setprecision(6);
      for(double v:values) {
        KUTTA_REQUIRE(std::isfinite(v),"non-finite residual at iteration "+std::to_string(r.iteration));
        residuals<<','<<v;if(!jsonl)human<<std::setw(14)<<v;
      }
      residuals<<','<<per_iteration<<'\n';
      // One force evaluation per report and sample serves the event and the
      // convergence conditions.
      std::vector<fv::ForceReport> report_forces;
      if(sample && (jsonl || report_conditions))for(const auto& d:reports)report_forces.push_back(solver.forces(d.patches,c.reference.pressure.value_or(0)));
      std::optional<setup::ConvergenceSample> state;
      if(sample && convergence.check) {
        std::map<std::string,double> residual_values,report_values;
        for(std::size_t i=0;i<std::size(values);++i)residual_values.emplace(residual_names[i],values[i]);
        for(std::size_t i=0;i<report_forces.size();++i)report_values.emplace(reports[i].name,dot(report_forces[i].force(),*reports[i].direction));
        state=monitor.sample(r.iteration,residual_values,report_values);
      }
      if(!jsonl)human<<std::setw(15)<<per_iteration<<std::endl;
      else {
        std::ostringstream event;event.imbue(std::locale::classic());event<<std::setprecision(17);
        event<<"{\"v\":1,\"event\":\"iteration\",\"iter\":"<<r.iteration<<",\"residuals\":{";
        for(std::size_t i=0;i<std::size(values);++i){if(i)event<<',';event<<json_string(residual_names[i])<<':'<<values[i];}
        event<<"},\"reports\":{";
        for(std::size_t i=0;i<report_forces.size();++i) {
          const auto& d=reports[i];if(i)event<<',';
          const auto& f=report_forces[i];
          const auto force=dot(f.force(),*d.direction);KUTTA_REQUIRE(std::isfinite(force),"non-finite report: "+d.name);
          event<<json_string(d.name)<<":{\"force\":"<<force;
          if(c.reference.density && *c.reference.density>0 && c.reference.velocity && *c.reference.velocity>0 && c.reference.area && *c.reference.area>0) {
            const auto coefficient=f.coefficient(*d.direction,{*c.reference.density,*c.reference.velocity,*c.reference.area,c.reference.pressure.value_or(0)});
            KUTTA_REQUIRE(std::isfinite(coefficient),"non-finite report coefficient: "+d.name);
            event<<",\"coefficient\":"<<coefficient;
          }
          event<<'}';
        }
        event<<"},\"sec_per_iter\":"<<per_iteration<<",\"elapsed\":"<<seconds(run_start);
        if(state)event<<",\"convergence\":"<<dump_json(state->json);
        event<<'}';
        std::cout<<event.str()<<std::endl;
      }
      // Met conditions win over a stop requested during the same sample.
      if(state && state->met){converged=true;met_by=state->met_by;break;}
      if(stopped)break;
      batch=Clock::now();batch_count=0;batch_snapshot_time=0;
    }
    const auto wall_time=seconds(run_start);residuals.close();
    memory::phase("iterations");
    write_fields(*dev,solver,setup->settings,o.out/"result.vtu",notice);
    memory::phase("write-vtu");
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
    memory::phase("forces-fluxes");
    app::save_view_input(o.out/"view-input.bin",*dev,solver,c,reference?&*reference:nullptr,*setup,notice);
    memory::phase("view-input");
    {
      // The boundary pressure and wall shear the view snapshot holds, after
      // the forces.csv evaluation refreshed them; totals match its wall rows.
      post::Loads loads;loads.reference={c.reference.area,c.reference.density,c.reference.velocity,c.reference.length,c.reference.pressure};
      for(const auto& d:reports)for(const auto& p:d.patches)if(std::find(loads.report_zones.begin(),loads.report_zones.end(),p)==loads.report_zones.end())loads.report_zones.push_back(p);
      const auto pressure=dev->download(solver.boundary_pressure());const auto& shear=solver.wall_shear();
      const auto sx=dev->download(shear.x),sy=dev->download(shear.y),sz=dev->download(shear.z);
      const auto& g=solver.geometry();const auto& m=solver.mesh();const double p_ref=c.reference.pressure.value_or(0);
      for(const auto& p:m.patches)if(p.kind==PatchKind::Wall) {
        post::ZoneLoads z{p.name,p.source_zone_id};
        for(Index f=p.start;f<p.start+p.size;++f) {
          const auto b=std::size_t(f-m.n_internal_faces());const auto area=g.face_area[f];const double a=norm(area);
          const auto on_pressure=area*(double(pressure[b])-p_ref),on_shear=Vec3d{sx[b],sy[b],sz[b]}*(-a);
          z.area+=a;z.pressure=z.pressure+on_pressure;z.viscous=z.viscous+on_shear;
          z.pressure_moment=z.pressure_moment+cross(g.face_centroid[f],on_pressure);z.viscous_moment=z.viscous_moment+cross(g.face_centroid[f],on_shear);
        }
        KUTTA_REQUIRE(finite(z.pressure) && finite(z.viscous) && finite(z.pressure_moment) && finite(z.viscous_moment),"non-finite wall loads: "+p.name);
        loads.zones.push_back(std::move(z));
      }
      loads.provenance=notice;loads.approximate=!setup->approximations.empty();loads.modified=!changes.empty();loads.iterations=solver.report().iteration;
      auto file=output(o.out/"loads.json");file<<dump_json(post::loads_json(loads))<<'\n';file.close();
    }
    auto summary_convergence=convergence.json();
    summary_convergence.emplace("met_by",Json::Array(met_by.begin(),met_by.end()));
    summary_convergence.emplace("iteration",converged?Json(solver.report().iteration):Json(nullptr));
    const auto stop_reason=converged?"converged":stopped?"stopped":"iteration_limit";
    auto summary=output(o.out/"summary.json");
    summary<<"{\n  \"v\": 1,\n  \"modified\": "<<(changes.empty()?"false":"true")<<",\n  \"overrides\": "<<dump_json(applied)<<",\n  \"run_controls\": "<<dump_json(run_controls)<<",\n  \"approximate\": "<<(setup->approximations.empty()?"false":"true")<<",\n  \"approximation_notice\": "<<json_string(setup::approximation_notice(setup->approximations))
      <<",\n  \"approximations\": [";
    for(std::size_t i=0;i<setup->approximations.size();++i) {
      const auto& info=setup::approximation_info(setup->approximations[i]);
      if(i) summary<<',';
      summary<<"\n    {\"name\": "<<json_string(info.name)<<", \"setting\": "<<json_string(info.setting)
        <<", \"fluent_value\": "<<json_string(info.fluent_value)<<", \"kutta_action\": "<<json_string(info.kutta_action)
        <<", \"effect_documentation\": "<<json_string(info.effect_documentation)<<", \"measured_effect\": "<<json_string(info.measured_effect)<<", \"short\": "<<json_string(info.short_label)<<'}';
    }
    summary<<"],\n  \"case_file\": "<<json_string(o.input.string())<<",\n  \"case_sha256\": "<<json_string(case_hash)
      <<",\n  \"data_file\": "<<(o.data.empty()?"null":json_string(o.data.string()))<<",\n  \"data_sha256\": "<<(o.data.empty()?"null":json_string(data_hash))
      <<",\n  \"kutta_version\": "<<json_string(KUTTA_VERSION)<<",\n  \"kutta_commit\": "<<json_string(KUTTA_COMMIT)
      <<",\n  \"backend\": "<<json_string(backend_name(dev->backend()))<<",\n  \"device\": "<<json_string(dev->name())
      <<",\n  \"cells\": "<<solver.mesh().n_cells<<",\n  \"poor_quality_cells\": "<<solver.poor_cells()<<",\n  \"pressure_velocity_coupling\": "<<json_string(coupling)
      <<",\n  \"iterations\": "<<solver.report().iteration
      <<",\n  \"wall_time_seconds\": "<<wall_time<<",\n  \"setup_time_seconds\": "<<setup_time;
    if(live)summary<<",\n  \"live_view\": "<<dump_json(live_statistics(wall_time));
    summary<<",\n  \"peak_device_memory_bytes\": null,\n  \"converged\": "<<(converged?"true":"false")<<",\n  \"stop_reason\": "<<json_string(stop_reason)
      <<",\n  \"convergence\": "<<dump_json(summary_convergence)<<",\n  \"convergence_criteria_available\": true,\n  \"boundary_wall_shear_available\": true\n}\n";
    summary.close();
    memory::phase("loads-summary");
    if(jsonl)std::cout<<"{\"v\":1,\"event\":\"done\",\"converged\":"<<(converged?"true":"false")<<",\"iterations\":"<<solver.report().iteration<<",\"elapsed\":"<<wall_time<<",\"out\":"<<json_string(utf8_path(o.out))<<",\"stop_reason\":"<<json_string(stop_reason)<<",\"convergence\":"<<dump_json(summary_convergence);
    if(jsonl){if(live)std::cout<<",\"live_view\":"<<dump_json(live_statistics(wall_time));std::cout<<'}'<<std::endl;}
    else {
      if(converged) {
        human<<"Converged at iteration "<<solver.report().iteration<<" (met:";
        for(const auto& m:met_by)human<<' '<<m;
        human<<").\n";
      } else if(stopped)human<<"Stopped on request at iteration "<<solver.report().iteration<<".\n";
      else human<<"Reached the iteration limit ("<<solver.report().iteration<<" iterations)"<<(convergence.check && convergence.has_conditions()?" before the convergence conditions were met":"")<<".\n";
      human<<"Wrote "<<o.out.string()<<" ("<<solver.report().iteration<<" iterations, "<<wall_time<<" s solve).\n";
    }
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
