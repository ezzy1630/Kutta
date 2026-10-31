#include <map>
#include <deque>
#include <set>
#include <cstdlib>
#include "../validation/cases/tmr_sa.hpp"
#include "sst_channel.hpp"
#include "kutta/core/parse.hpp"
#include "kutta/setup/fluent.hpp"

// Independent convergence experiments; the frozen validation rigs are unchanged.
int main(int argc,char** argv) {
  using namespace kutta;using namespace fv;using namespace validation::tmr;
  try {
    std::map<std::string,std::string> args;
    const std::set<std::string> keys{"--case","--backend","--model","--level","--iterations","--data","--fixture",
      "--nx","--ny","--nz","--viscosity","--interval","--dt","--time-scale","--momentum-relaxation",
      "--pressure-relaxation","--momentum-tolerance","--pressure-tolerance","--absolute-tolerance",
      "--turbulence-tolerance","--turbulence-relaxation","--max-seconds","--scheme","--coupled-tolerance",
      "--coupled-iterations"};
    for(int i=1;i<argc;i+=2) {KUTTA_REQUIRE(i+1<argc,"expected --key value");KUTTA_REQUIRE(keys.contains(argv[i]),"unknown option: "+std::string(argv[i]));KUTTA_REQUIRE(args.emplace(argv[i],argv[i+1]).second,"duplicate option");}
    auto str=[&](const std::string& key,const std::string& fallback) {auto i=args.find(key);return i==args.end()?fallback:i->second;};
    auto number=[&](const std::string& key,double fallback) {auto i=args.find(key);if(i==args.end()) return fallback;auto q=parse_double(i->second);KUTTA_REQUIRE(q && std::isfinite(*q),"invalid number: "+key);return *q;};
    auto integer=[&](const std::string& key,int fallback,int minimum=1) {const auto q=number(key,fallback);KUTTA_REQUIRE(q>=minimum && q<=std::numeric_limits<int>::max() && std::floor(q)==q,"invalid integer: "+key);return int(q);};
    const auto name=str("--case","plate"),backend=str("--backend","cpu"),model=str("--model","sst");
    KUTTA_REQUIRE(backend=="cpu" || backend=="metal" || backend=="cuda","invalid backend");
    KUTTA_REQUIRE(model=="sst" || model=="sa","invalid model");
    const int level=integer("--level",0,0),iterations=integer("--iterations",3000);
    KUTTA_REQUIRE(level>=0 && level<=2 && iterations>0,"invalid level/iteration count");
    auto dev=Device::create(backend=="cpu"?Backend::Cpu:backend=="metal"?Backend::Metal:Backend::Cuda);
    PolyMesh mesh;std::vector<PatchCondition> bc;SolverSettings s;std::optional<setup::Initialization> initialization;
    // Fluent cases: the V6 elbow and the SST nozzle, each against its saved run.
    const bool plate=name=="plate",channel=name=="channel",nozzle=name=="nozzle",elbow=name=="elbow" || nozzle;
    KUTTA_REQUIRE(plate || name=="bump" || channel || elbow,"invalid case");
    KUTTA_REQUIRE(!(channel || elbow) || model=="sst","channel/elbow/nozzle require SST");
    double continuity_target=0;
    const auto c=plate?Case::FlatPlate:Case::Bump;
    if(channel) {
      mesh=bench::sst_channel(integer("--nx",24,2),integer("--ny",8,2),integer("--nz",4));
      bc=bench::channel_boundaries(mesh);s=bench::channel_settings();s.viscosity=number("--viscosity",.0001);
    } else if(elbow) {
      const auto* fixtures=std::getenv("KUTTA_FLUENT_EXAMPLES");
      const auto stem=str("--fixture",fixtures?fixtures:"validation/fluent_examples")+
        (nozzle?"/flprj_3891/nozzle_para_named_1.cffdb/nozzle_para_named_1-Solve/Base DP/nozzle_para_named_1":"/mixing_elbow/mixing_elbow");
      const auto path=stem+".cas.h5";
      auto fluent=fluent::read_case(path);mesh=read_mesh(path);
      // V6's criterion: the saved run's final raw continuity.
      continuity_target=fluent::read_data(stem+".dat.h5",fluent).residuals.at("continuity").rows.back()[0];
      auto problem=setup::fluent_setup(fluent,mesh,setup::parse_approximations("wfgc"));
      std::cout<<"DIAGNOSTIC "<<setup::approximation_notice(problem.approximations)<<std::endl;
      for(const auto& note:problem.convergence_notes) std::cout<<"CONVERGENCE_NOTE "<<note<<std::endl;
      s=problem.settings;bc=std::move(problem.boundaries);initialization=setup::initial_fields(fluent,mesh);
    } else {
      s.viscosity=plate?1./5e6:1./3e6;
      mesh=read_plot3d(grid(str("--data","validation/data/nasa_tmr"),c,level),span);
      for(const auto& p:mesh.patches) {
        const TurbulenceCondition t=model=="sa"?TurbulenceCondition{ModifiedViscosity{3*s.viscosity}}:TurbulenceCondition{TurbulenceValues{2.25e-7,25e-6/s.viscosity}};
        if(p.kind==PatchKind::Wall) bc.push_back({p.name,Wall{}});
        else if(p.source_type=="subsonic_inflow_pt") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0},t}});
        else if(p.source_type=="back_pressure" || p.source_type=="farfield_riem") bc.push_back({p.name,PressureOutlet{0,Vec3d{1,0,0},BackflowPressure::Static,t}});
        else if(p.kind==PatchKind::Symmetry) bc.push_back({p.name,Symmetry{}});
        else throw Error("unmapped patch "+p.name);
      }
      if(model=="sst") {s.sst=SstSettings{};s.sst->k_convection=s.sst->omega_convection=ConvectionScheme::FirstOrderUpwind;}
      else {s.sa=SaSettings{};s.sa->model.formulation=SaFormulation::TmrStandard;s.sa->convection=ConvectionScheme::FirstOrderUpwind;}
      if(!plate) {s.pseudo_time_step=.01;s.momentum_relaxation=.5;s.pressure_relaxation=.3;}
    }
    const int interval=integer("--interval",elbow?50:channel?20:100);s.monitor_interval=interval;s.pseudo_time_step=number("--dt",s.pseudo_time_step);s.pseudo_time_scale=number("--time-scale",s.pseudo_time_scale);
    const auto max_seconds=number("--max-seconds",600);KUTTA_REQUIRE(max_seconds>0,"max seconds must be positive");
    s.momentum_relaxation=number("--momentum-relaxation",s.momentum_relaxation);s.pressure_relaxation=number("--pressure-relaxation",s.pressure_relaxation);
    s.momentum_linear.rel_tol=number("--momentum-tolerance",s.momentum_linear.rel_tol);s.pressure_linear.rel_tol=number("--pressure-tolerance",s.pressure_linear.rel_tol);
    s.momentum_linear.abs_tol=number("--absolute-tolerance",s.momentum_linear.abs_tol);s.pressure_linear.abs_tol=number("--absolute-tolerance",s.pressure_linear.abs_tol);
    if(s.sst) {s.sst->k_relaxation=s.sst->omega_relaxation=number("--turbulence-relaxation",s.sst->k_relaxation);s.sst->linear.abs_tol=number("--absolute-tolerance",s.sst->linear.abs_tol);s.sst->linear.rel_tol=number("--turbulence-tolerance",s.sst->linear.rel_tol);}
    if(s.sa) {s.sa->relaxation=number("--turbulence-relaxation",s.sa->relaxation);s.sa->linear.abs_tol=number("--absolute-tolerance",s.sa->linear.abs_tol);s.sa->linear.rel_tol=number("--turbulence-tolerance",s.sa->linear.rel_tol);}
    // Elbow follows its case's flow/scheme (Coupled) unless --scheme overrides.
    if(const auto scheme=str("--scheme","");!scheme.empty()) {
      KUTTA_REQUIRE(scheme=="simplec" || scheme=="coupled","invalid scheme");
      s.coupling=scheme=="coupled"?PressureVelocityCoupling::Coupled:PressureVelocityCoupling::Simplec;
    }
    s.coupled_linear.rel_tol=number("--coupled-tolerance",s.coupled_linear.rel_tol);
    s.coupled_linear.max_iter=integer("--coupled-iterations",s.coupled_linear.max_iter);
    SteadySolver solver(*dev,std::move(mesh),std::move(bc),s);const auto n=solver.mesh().n_cells;
    if(initialization) setup::initialize(solver,*initialization);
    else if(channel) bench::initialize_channel(solver);
    else if(model=="sa") solver.initialize_sa_fields(std::vector<Vec3d>(n,{1,0,0}),std::vector<double>(n,0),std::vector<double>(n,3*s.viscosity));
    else solver.initialize_fields(std::vector<Vec3d>(n,{1,0,0}),std::vector<double>(n,0),std::vector<double>(n,2.25e-7),std::vector<double>(n,25e-6/s.viscosity));
    using Clock=std::chrono::steady_clock;dev->synchronize();const auto start=Clock::now();
    std::cout<<std::setprecision(12)<<"CASE source="<<KUTTA_BENCH_COMMIT<<" name="<<name<<" model="<<model<<" backend="<<backend<<" cells="<<n<<" dt="<<s.pseudo_time_step<<" momentum_relaxation="<<s.momentum_relaxation<<" pressure_relaxation="<<s.pressure_relaxation
      <<" scheme="<<(s.coupling==PressureVelocityCoupling::Coupled?"coupled":"simplec")
      <<" coupled_bytes_per_cell="<<double(solver.coupled_solver_bytes())/n<<" continuity_target="<<continuity_target<<std::endl;
    std::deque<std::pair<int,Sample>> window;bool converged=false;
    while(solver.report().iteration<iterations) {
      try {solver.iterate(std::min(interval,iterations-solver.report().iteration));}
      catch(const Error& e) {
        std::cout<<"FAILURE completed_iterations="<<solver.report().iteration
          <<" seconds="<<std::chrono::duration<double>(Clock::now()-start).count()<<" reason="<<e.what()<<std::endl;
        return 1;
      }
      const auto r=solver.monitor();
      Sample value;
      if(channel || elbow) {std::vector<std::string> walls;for(const auto& p:solver.mesh().patches) if(p.kind==PatchKind::Wall) walls.push_back(p.name);value.drag=solver.forces(walls).force().x;}
      else value=model=="sa"?sa::sample(*dev,solver,c):sample(*dev,solver,c);
      while(window.size()>1 && window[1].first<=r.iteration-100) window.pop_front();
      const bool has_window=!window.empty() && window.front().first<=r.iteration-100;
      const double delta=has_window?change(window.front().second,value,c):std::numeric_limits<double>::infinity();
      const double step_change=window.empty()?std::numeric_limits<double>::infinity():change(window.back().second,value,c);
      const double d=std::max({r.u.scaled,r.v.scaled,r.w.scaled,r.k.scaled,r.omega.scaled,r.nu_tilde.scaled,r.continuity.scaled});
      // V6's recorded raw mass target and steady transport criterion differ
      // from the TMR gates. WFGC remains an explicitly allowed diagnostic.
      converged=elbow?(std::max({r.u.scaled,r.v.scaled,r.w.scaled,r.k.scaled,r.omega.scaled})<1e-6 && r.continuity.raw<continuity_target):
        (has_window && d<=1e-5 && delta<=.001);
      std::cout<<"STEP it="<<r.iteration<<" seconds="<<std::chrono::duration<double>(Clock::now()-start).count()<<" defect="<<d<<" continuity="<<r.continuity.scaled<<" raw="<<r.continuity.raw<<" scale="<<r.continuity.fluent_normalization<<" corrected="<<r.corrected_mass_imbalance<<" u="<<r.u.scaled<<" v="<<r.v.scaled<<" k="<<r.k.scaled<<" omega="<<r.omega.scaled<<" sa="<<r.nu_tilde.scaled<<" change="<<delta<<" step_change="<<step_change<<" drag="<<value.drag<<" cf0="<<value.cf[0]<<" cf1="<<value.cf[1]<<" cf2="<<value.cf[2]<<" inner_u="<<r.momentum[0].iterations<<" inner_sa="<<r.turbulence[0].iterations<<" inner_p="<<r.pressure.iterations<<" p_ratio="<<r.pressure.final_residual/std::max(r.pressure.initial_residual,1e-30)<<" converged="<<converged<<std::endl;
      if(converged) break;
      if(std::chrono::duration<double>(Clock::now()-start).count()>=max_seconds) {std::cout<<"BUDGET max_seconds="<<max_seconds<<std::endl;break;}
      window.push_back({r.iteration,value});
    }
    return converged?0:2;
  } catch(const std::exception& e) {std::cerr<<e.what()<<std::endl;return 1;}
}
