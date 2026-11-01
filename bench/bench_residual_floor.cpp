// Where a Fluent case's continuity residual stops falling, and why
// (validation/results/residual_floor/REPORT.md). Prints:
//  - FLUENT_*: the saved residual rows (raw kg/s and Fluent's own scale) and
//    whether the saved cell fields are float32 values, i.e. how Fluent ran;
//  - KUTTA_HISTORY: Kutta's raw continuity, its first-five scale, the
//    corrected imbalance summed in float on the device and recomputed in
//    double from the stored face fluxes, sum |F|, the per-sample field change
//    and the last step's x-momentum and pressure Krylov iterations;
//  - KUTTA_FLOOR: the residuals re-evaluated at the final state and at that state
//    moved by one float32 ulp (velocity, pressure, both), and the closure of
//    the stored float area vectors.
// It never changes the solver; it only reads it. The legacy cylinder .dat must be
// the CR-restored copy (bench_restore_stripped_cr). A trailing float32 argument
// turns SolverSettings::paired_state off, reproducing the float32 floor. The
// probes re-seed hi alone, so with the pair they measure hi's own residual.
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/setup/fluent.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

using namespace kutta;
using namespace kutta::detail;
namespace {
void fluent_reference(const fluent::FluentCase& c,const std::filesystem::path& data) {
  const auto d=fluent::read_data(data,c);
  for(const auto* name:{"SV_U","SV_V","SV_W","SV_P"}) if(const auto f=d.cell_fields.find(name);f!=d.cell_fields.end()) {
    std::size_t total=0,exact=0;
    for(const auto& s:f->second.sections) for(const double v:s.values) {++total;exact+=double(float(v))==v;}
    std::cout<<"FLUENT_PRECISION cells field="<<name<<" values="<<total<<" float32_exact="<<exact<<'\n';
  }
  for(const auto* name:{"continuity","x-velocity"}) if(const auto r=d.residuals.find(name);r!=d.residuals.end() && !r->second.rows.empty()) {
    const auto& h=r->second;
    for(const std::size_t i:{std::size_t(0),h.rows.size()-1})
      std::cout<<"FLUENT_RESIDUAL "<<name<<" iteration="<<h.iterations[i]<<" raw="<<h.rows[i][0]<<" scale="<<h.rows[i][1]<<" scaled="<<h.rows[i][0]/h.rows[i][1]<<'\n';
  }
}
struct State {std::vector<float> u,v,w,p,k,second;};
State download(Device& dev,const fv::SteadySolver& s,const fv::SolverSettings& settings) {
  State out{dev.download(s.velocity().x),dev.download(s.velocity().y),dev.download(s.velocity().z),dev.download(s.pressure()),{},{}};
  if(settings.sst || settings.realizable_ke) out.k=dev.download(s.kinetic_energy());
  if(settings.sst) out.second=dev.download(s.specific_dissipation());
  if(settings.realizable_ke) out.second=dev.download(s.dissipation());
  return out;
}
// Re-seed the solver at a state (solver cell order) and evaluate its residuals.
fv::IterationReport evaluate(fv::SteadySolver& s,const fv::SolverSettings& settings,const State& x) {
  const auto& perm=s.permutation().cells_new_to_old;const auto n=x.u.size();
  std::vector<Vec3d> velocity(n);std::vector<double> p(n),k(x.k.size()),second(x.second.size());
  for(std::size_t i=0;i<n;++i) {
    velocity[perm[i]]={x.u[i],x.v[i],x.w[i]};p[perm[i]]=x.p[i];
    if(!k.empty()) {k[perm[i]]=x.k[i];second[perm[i]]=x.second[i];}
  }
  if(settings.realizable_ke) s.initialize_ke_fields(velocity,p,k,second);else s.initialize_fields(velocity,p,k,second);
  return s.monitor();
}
}

int main(int argc,char** argv) {
  try {
    KUTTA_REQUIRE(argc==6 || (argc==7 && std::string_view(argv[6])=="float32"),
      "Usage: bench_residual_floor case data|- cpu|metal approximations iterations [float32]");
    const auto c=fluent::read_case(argv[1]);const auto mesh=read_mesh(argv[1]);
    std::cout<<std::setprecision(8);
    if(std::string_view(argv[2])!="-") fluent_reference(c,argv[2]);
    auto problem=setup::fluent_setup(c,mesh,setup::parse_approximations(argv[4]));
    problem.settings.paired_state=argc==6;
    const auto& settings=problem.settings;
    KUTTA_REQUIRE(std::string_view(argv[3])=="cpu" || std::string_view(argv[3])=="metal","expected cpu|metal");
    auto dev=Device::create(std::string_view(argv[3])=="cpu"?Backend::Cpu:Backend::Metal);
    const int iterations=std::stoi(argv[5]);KUTTA_REQUIRE(iterations>0,"iterations must be positive");
    fv::SteadySolver solver(*dev,mesh,problem.boundaries,settings);
    setup::initialize(solver,setup::initial_fields(c,mesh));
    const auto& m=solver.mesh();const auto n=std::size_t(m.n_cells);
    std::cout<<"KUTTA_INIT continuity_raw="<<solver.monitor().continuity.raw<<'\n';
    State previous=download(*dev,solver,settings);
    double flux_sum=0;
    // The monitor's own schedule: iterations 1-5 and every monitor interval.
    for(int it=1;it<=iterations;++it) {
      solver.iterate(1);
      if(it>5 && it%settings.monitor_interval!=0 && it!=iterations) continue;
      const auto& r=solver.report();
      const auto mass=dev->download(solver.mass_flux());
      std::vector<double> net(n,0);flux_sum=0;
      for(Index f=0;f<m.n_faces();++f) {
        net[m.owner[f]]+=mass[f];flux_sum+=std::abs(double(mass[f]));
        if(f<m.n_internal_faces()) {net[m.neighbour[f]]-=mass[f];flux_sum+=std::abs(double(mass[f]));}
      }
      double corrected=0;for(const double x:net) corrected+=std::abs(x);
      const auto now=download(*dev,solver,settings);
      double du=0,dp=0;
      for(std::size_t i=0;i<n;++i) {
        du=std::max({du,double(std::abs(now.u[i]-previous.u[i])),double(std::abs(now.v[i]-previous.v[i])),double(std::abs(now.w[i]-previous.w[i]))});
        dp=std::max(dp,double(std::abs(now.p[i]-previous.p[i])));
      }
      previous=now;
      if(it<=5 || it%(5*settings.monitor_interval)==0 || it==iterations)
        std::cout<<"KUTTA_HISTORY iteration="<<it<<" continuity_raw="<<r.continuity.raw<<" scale="<<r.continuity.fluent_normalization
          <<" scaled="<<r.continuity.scaled<<" corrected_f32="<<r.corrected_mass_imbalance<<" corrected_f64="<<corrected
          <<" sum_abs_flux="<<flux_sum<<" x_velocity_raw="<<r.u.raw<<" x_velocity_scaled="<<r.u.scaled
          <<" max_velocity_change="<<du<<" max_pressure_change="<<dp
          <<" u_krylov="<<r.momentum[0].iterations<<" p_krylov="<<r.pressure.iterations<<'\n';
    }
    const auto final_report=solver.report();
    const auto x=download(*dev,solver,settings);
    {
      // A uniform 1 m/s flow's imbalance from the float area vectors' closure.
      const auto& dm=solver.device_mesh();const auto ax=dev->download(dm.area.x),ay=dev->download(dm.area.y),az=dev->download(dm.area.z);
      std::vector<Vec3d> closure(n);
      for(Index f=0;f<m.n_faces();++f) {
        const Vec3d a{ax[f],ay[f],az[f]};closure[m.owner[f]]=closure[m.owner[f]]+a;
        if(f<m.n_internal_faces()) closure[m.neighbour[f]]=closure[m.neighbour[f]]-a;
      }
      double l1=0;for(const auto& s:closure) l1+=std::abs(s.x)+std::abs(s.y)+std::abs(s.z);
      std::cout<<"KUTTA_FLOOR area_closure_l1="<<l1<<" m^2 (rho*U*this bounds a uniform flow's imbalance)\n";
    }
    const double eps=std::numeric_limits<float>::epsilon();
    std::mt19937 rng(7);
    auto nudge=[&](std::vector<float>& values) {for(auto& v:values) v=std::nextafter(v,(rng()&1)?INFINITY:-INFINITY);};
    for(const auto* kind:{"none","velocity","pressure","both","both","both"}) {
      auto y=x;const std::string_view k=kind;
      if(k=="velocity" || k=="both") {nudge(y.u);nudge(y.v);nudge(y.w);}
      if(k=="pressure" || k=="both") nudge(y.p);
      const auto r=evaluate(solver,settings,y);
      std::cout<<"KUTTA_FLOOR perturbation_ulp="<<kind<<" continuity_raw="<<r.continuity.raw<<" x_velocity_raw="<<r.u.raw<<'\n';
    }
    const double scale=final_report.continuity.fluent_normalization;
    std::cout<<"KUTTA_FINAL iteration="<<final_report.iteration<<" continuity_raw="<<final_report.continuity.raw<<" scale="<<scale
      <<" scaled="<<final_report.continuity.scaled<<" raw_over_eps_sum_abs_flux="<<final_report.continuity.raw/(eps*flux_sum)
      <<" scale_over_sum_abs_flux="<<scale/flux_sum<<'\n';
    return 0;
  } catch(const std::exception& e) {std::cerr<<"bench_residual_floor: "<<e.what()<<'\n';return 1;}
}
