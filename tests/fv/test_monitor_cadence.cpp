#include "kernels.hpp"
#include "../backend_helpers.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"

using namespace kutta;
using namespace kutta::fv;

// Monitoring only measures. With the automatic pseudo time step, the step used
// to be refreshed in monitor(), so runs at different monitor intervals diverged.
// The differentiable limiter's thresholds are refreshed with the flux, never by
// monitor(); SST covers them under SIMPLEC and Coupled.
TEST_CASE("Monitor interval does not change the solution for any turbulence model") {
  enum class Model {sst,sa,realizable_ke};
  struct Case {Model model;Limiter limiter=Limiter::Standard;PressureVelocityCoupling coupling=PressureVelocityCoupling::Simplec;bool paired=true;};
  constexpr auto coupled=PressureVelocityCoupling::Coupled;
  const auto mesh=make_box_hex(6,4,1,{0,0,0},{1,.4,.1});
  for(auto [model,limiter,coupling,paired]:{Case{Model::sst},Case{Model::sa},Case{Model::realizable_ke},
      Case{Model::sst,Limiter::Differentiable},Case{Model::sst,Limiter::Standard,coupled},
      Case{Model::sst,Limiter::Differentiable,coupled},Case{Model::sst,Limiter::Differentiable,coupled,false}}) {
    const TurbulenceCondition inflow=model==Model::sst?TurbulenceCondition{TurbulenceValues{.01,4}}:
      model==Model::sa?TurbulenceCondition{ModifiedViscosity{.003}}:TurbulenceCondition{KeValues{.01,.01}};
    std::vector<PatchCondition> bc;
    for(const auto& p:mesh.patches) {
      if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0},inflow}});
      else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0,Vec3d{1,0,0},BackflowPressure::Static,inflow}});
      else if(p.name=="ymin" || p.name=="ymax") bc.push_back({p.name,Wall{}});
      else bc.push_back({p.name,Symmetry{}});
    }
    KUTTA_FOR_EACH_BACKEND(dev) {
      const auto run=[&](int interval) {
        SolverSettings s;s.viscosity=.001;s.monitor_interval=interval;s.limiter=limiter;s.coupling=coupling;s.paired_state=paired;
        if(model==Model::sst) s.sst=SstSettings{};
        else if(model==Model::sa) s.sa=SaSettings{};
        else s.realizable_ke=RealizableKeSettings{};
        REQUIRE(s.pseudo_time_step==0);  // automatic stepping is what monitor() used to disturb
        SteadySolver solver(*dev,mesh,bc,s);
        const std::vector<Vec3d> u(mesh.n_cells,{.5,0,0});const std::vector<double> p(mesh.n_cells,0);
        if(model==Model::sst) solver.initialize_fields(u,p,std::vector<double>(mesh.n_cells,.01),std::vector<double>(mesh.n_cells,4));
        else if(model==Model::sa) solver.initialize_sa_fields(u,p,std::vector<double>(mesh.n_cells,.003));
        else solver.initialize_ke_fields(u,p,std::vector<double>(mesh.n_cells,.01),std::vector<double>(mesh.n_cells,.01));
        solver.iterate(12);
        const auto& v=solver.velocity();
        std::vector<std::vector<float>> fields{dev->download(v.x),dev->download(v.y),dev->download(v.z),
          dev->download(solver.pressure()),dev->download(solver.mass_flux())};
        if(model==Model::sst) {fields.push_back(dev->download(solver.kinetic_energy()));fields.push_back(dev->download(solver.specific_dissipation()));}
        else if(model==Model::sa) fields.push_back(dev->download(solver.modified_viscosity()));
        else {fields.push_back(dev->download(solver.kinetic_energy()));fields.push_back(dev->download(solver.dissipation()));}
        return fields;
      };
      // Name the case and field, so a platform-only failure says where it is.
      INFO("model ",int(model)," limiter ",int(limiter)," coupling ",int(coupling)," paired ",paired," backend ",dev->name());
      const auto every=run(1),fifth=run(5);
      REQUIRE(every.size()==fifth.size());
      for(std::size_t f=0;f<every.size();++f) {
        std::size_t differ=0;
        for(std::size_t i=0;i<every[f].size();++i) differ+=every[f][i]!=fifth[f][i];
        CAPTURE(f);CHECK(differ==0);
      }
    }
  }
}
