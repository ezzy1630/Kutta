// Bounded offline history of the baseline SST plate. Public field downloads
// are diagnostic only; this instrumentation never enters SteadySolver.
#include "../../../validation/cases/tmr.hpp"

int main(int argc,char** argv) {
  using namespace kutta;using namespace fv;using namespace validation::tmr;
  try {
    KUTTA_REQUIRE(argc==4,"usage: plate_trace cpu|metal tmr-directory output.csv");
    const std::string backend=argv[1];
    KUTTA_REQUIRE(backend=="cpu" || backend=="metal","unsupported diagnostic backend");
    auto dev=Device::create(backend=="cpu"?Backend::Cpu:Backend::Metal);
    auto mesh=read_plot3d(grid(argv[2],Case::FlatPlate,0),span);
    SolverSettings settings;settings.viscosity=1./5e6;settings.sst=SstSettings{};
    settings.monitor_interval=100;
    settings.sst->k_convection=settings.sst->omega_convection=ConvectionScheme::FirstOrderUpwind;
    std::vector<PatchCondition> bc;
    const TurbulenceValues turbulence{2.25e-7,25e-6/settings.viscosity};
    for(const auto& p:mesh.patches) {
      if(p.kind==PatchKind::Wall) bc.push_back({p.name,Wall{}});
      else if(p.source_type=="subsonic_inflow_pt") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0},turbulence}});
      else if(p.source_type=="back_pressure" || p.source_type=="farfield_riem")
        bc.push_back({p.name,PressureOutlet{0,Vec3d{1,0,0},BackflowPressure::Static,turbulence}});
      else {KUTTA_REQUIRE(p.kind==PatchKind::Symmetry,"unexpected patch");bc.push_back({p.name,Symmetry{}});}
    }
    SteadySolver solver(*dev,std::move(mesh),std::move(bc),settings);
    const int n=solver.mesh().n_cells,cell=795;
    const auto point=solver.geometry().cell_centroid.at(cell);
    KUTTA_REQUIRE(n==816 && std::abs(point.x+.26638596013494364)<1e-12 &&
      std::abs(point.y-.000123022928060775)<1e-12,"captured cell identity mismatch");
    solver.initialize_fields(std::vector<Vec3d>(n,{1,0,0}),std::vector<double>(n,0),
      std::vector<double>(n,turbulence.k),std::vector<double>(n,turbulence.omega));
    // Match the frozen V5 monitor cadence, including its explicit monitor call.
    for(int i=0;i<29;++i) {solver.iterate(100);solver.monitor();}
    auto previous_u=dev->download(solver.velocity().x),previous_v=dev->download(solver.velocity().y);
    auto previous_p=dev->download(solver.pressure());
    std::ofstream out(argv[3]);KUTTA_REQUIRE(out,"cannot open trace output");
    out<<std::setprecision(17)<<"iteration,cell_u,cell_v,cell_p,delta_u,delta_v,delta_p,max_delta_u,max_delta_v,max_delta_p,pressure_rhs_l2,continuity_raw_if_monitored\n";
    for(int iteration=2901;iteration<=3000;++iteration) {
      solver.iterate(1);
      const auto u=dev->download(solver.velocity().x),v=dev->download(solver.velocity().y),p=dev->download(solver.pressure());
      double du=0,dv=0,dp=0;
      for(int c=0;c<n;++c) {du=std::max(du,std::abs(double(u[c])-previous_u[c]));
        dv=std::max(dv,std::abs(double(v[c])-previous_v[c]));dp=std::max(dp,std::abs(double(p[c])-previous_p[c]));}
      const auto& r=solver.report();
      out<<iteration<<','<<u[cell]<<','<<v[cell]<<','<<p[cell]<<','<<double(u[cell])-previous_u[cell]<<','
        <<double(v[cell])-previous_v[cell]<<','<<double(p[cell])-previous_p[cell]<<','<<du<<','<<dv<<','<<dp<<','<<r.pressure.initial_residual<<',';
      if(r.residual_iteration==iteration) out<<r.continuity.raw;
      out<<'\n';previous_u=u;previous_v=v;previous_p=p;
    }
    const auto r=solver.monitor();
    std::cout<<std::setprecision(17)<<"backend="<<backend<<" iterations="<<r.iteration<<" continuity="<<r.continuity.scaled
      <<" raw="<<r.continuity.raw<<" u="<<r.u.scaled<<" v="<<r.v.scaled<<'\n';
    return 0;
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
