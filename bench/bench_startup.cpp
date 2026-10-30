// Deliberately synchronizing diagnostic of the real startup path, never a gate.
#include <iostream>
#include <charconv>
#include "kutta/core/parse.hpp"
#include "../tests/linalg/instrumented_device.hpp"
#include "../validation/cases/tmr.hpp"
#include "kernels.hpp"
#include "kutta/mesh/detail.hpp"

int main(int argc,char** argv) {
  using namespace kutta;using namespace kutta::detail;
  try {
    KUTTA_REQUIRE(argc>=5 && argc<=7,"usage: bench_startup <tmr-root> <cpu|metal> <level> <steps> [uniform|hybrid|laminar] [dt]");
    const auto integer=[](const char* text) {
      int value=0;const auto [end,error]=std::from_chars(text,text+std::char_traits<char>::length(text),value);
      KUTTA_REQUIRE(error==std::errc{} && *end=='\0',"invalid integer");return value;
    };
    const int level=integer(argv[3]),steps=integer(argv[4]);
    KUTTA_REQUIRE(level>=0 && level<=2 && steps>0,"invalid level or step count");
    const std::string backend=argv[2];
    KUTTA_REQUIRE(backend=="cpu" || backend=="metal","unknown backend");
    test_linalg::InstrumentedDevice dev(backend=="cpu"?Backend::Cpu:Backend::Metal);
    auto mesh=read_plot3d(validation::tmr::grid(argv[1],validation::tmr::Case::Bump,level),validation::tmr::span);
    renumber_cells_rcm(mesh);auto geom=compute_geometry(mesh);
    double skew=0,min_height=std::numeric_limits<double>::infinity();
    for(Index f=0;f<mesh.n_faces();++f) {
      const auto r=geom.face_centroid[f]-geom.cell_centroid[mesh.owner[f]],area=geom.face_area[f];
      const double height=std::abs(dot(r,area))/norm(area);min_height=std::min(min_height,height);
      if(f>=mesh.n_internal_faces()) skew=std::max(skew,norm(r)/height);
    }
    std::cout<<std::setprecision(12)<<"GEOMETRY max_boundary_offset_over_normal_distance="<<skew<<" min_normal_distance="<<min_height<<'\n';
    int step=0;
    auto range=[&](const char* label,BufferHandle& buffer) {
      const auto v=dev.read<float>(buffer);
      const auto lo=std::min_element(v.begin(),v.end()),hi=std::max_element(v.begin(),v.end());
      double square=0;for(auto x:v) square+=double(x)*x;
      std::cout<<' '<<label<<"_min="<<*lo<<' '<<label<<"_max="<<*hi<<' '<<label<<"_l2="<<std::sqrt(square);
    };
    dev.after_launch=[&](KernelId id,std::size_t,std::span<const KernelArg> a) {
      if(id==KernelId::fv_pressure_faces) {
        const auto owner=dev.read<Index>(*a[0].buffer),neighbour=dev.read<Index>(*a[1].buffer);
        const auto px=dev.read<float>(*a[2].buffer),py=dev.read<float>(*a[3].buffer),pz=dev.read<float>(*a[4].buffer);
        const auto nx=dev.read<float>(*a[5].buffer),ny=dev.read<float>(*a[6].buffer),nz=dev.read<float>(*a[7].buffer);
        const auto gx=dev.read<float>(*a[9].buffer),gy=dev.read<float>(*a[10].buffer),gz=dev.read<float>(*a[11].buffer);
        const auto p=dev.read<float>(*a[8].buffer),jump=dev.read<float>(*a[14].buffer);
        double original=0,retained=0,ratio=0,remote_original=0,remote_retained=0;
        for(Index f=0;f<a[15].get<Index>();++f) {
          const auto c=owner[f],n=neighbour[f];
          const Vec3d rp{px[f],py[f],pz[f]},rn{nx[f],ny[f],nz[f]},gp{gx[c],gy[c],gz[c]},gn{gx[n],gy[n],gz[n]};
          const auto d=rp-rn,skew=(rp+rn)*.5;
          const double shared=p[c]-p[n]+dot((gp+gn)*.5,d);
          original=std::max(original,std::abs(dot(gp-gn,skew)));
          retained=std::max(retained,std::abs(jump[f]-shared));ratio=std::max(ratio,norm(skew)/norm(d));
          if(norm(skew)>2*norm(d)) {
            remote_original=std::max(remote_original,std::abs(dot(gp-gn,skew)));
            remote_retained=std::max(remote_retained,std::abs(jump[f]-shared));
          }
        }
        std::cout<<"PRESSURE_SKEW step="<<step<<" original_term="<<original<<" retained_term="<<retained
          <<" remote_original="<<remote_original<<" remote_retained="<<remote_retained<<" max_offset_over_link="<<ratio<<'\n';
      }
      if(id==KernelId::fv_floor) {
        std::cout<<"STARTUP step="<<step<<" kernel=floor value="<<a[1].get<float>();
        range("field",*a[0].buffer);std::cout<<'\n';
      }
      if(id==KernelId::fv_momentum_faces) {
        std::cout<<"STARTUP step="<<step<<" kernel=face_velocity";
        range("cell_v",*a[10].buffer);range("face_v",*a[26].buffer);std::cout<<'\n';
      }
      if(id==KernelId::fv_rhie_chow) {
        std::cout<<"STARTUP step="<<step<<" kernel=flux";
        range("physical",*a[8].buffer);range("inverse",*a[9].buffer);
        range("flux",*a[14].buffer);std::cout<<'\n';
      }
      if(id==KernelId::fv_pressure_assemble) {
        std::cout<<"STARTUP step="<<step<<" kernel=pressure";
        range("rhs",*a[10].buffer);std::cout<<'\n';
      }
      if(id==KernelId::fv_correct_velocity) {
        std::cout<<"STARTUP step="<<step++<<" kernel=correct";
        range("correction",*a[6].buffer);range("u",*a[7].buffer);range("v",*a[8].buffer);
        range("p",*a[10].buffer);std::cout<<'\n';
        const auto v=dev.read<float>(*a[8].buffer),gy=dev.read<float>(*a[4].buffer),inv=dev.read<float>(*a[1].buffer);
        const auto j=std::max_element(v.begin(),v.end(),[](float x,float y){return std::abs(x)<std::abs(y);})-v.begin();
        std::cout<<"WORST cell="<<j<<" x="<<geom.cell_centroid[j].x<<" y="<<geom.cell_centroid[j].y<<" v="<<v[j]<<" gy="<<gy[j]<<" d="<<geom.cell_volume[j]*inv[j]<<'\n';
      }
    };
    const auto input=read_plot3d(validation::tmr::grid(argv[1],validation::tmr::Case::Bump,level),validation::tmr::span);
    const std::string mode=argc>5?argv[5]:"uniform";
    KUTTA_REQUIRE(mode=="uniform" || mode=="hybrid" || mode=="laminar","invalid initialization/model diagnostic");
    std::vector<fv::PatchCondition> bc;fv::TurbulenceValues turbulence{2.25e-7,75};
    for(const auto& p:input.patches) {
      if(p.kind==PatchKind::Wall) bc.push_back({p.name,fv::Wall{}});
      else if(p.source_type=="subsonic_inflow_pt") bc.push_back({p.name,fv::VelocityInlet{Vec3d{1,0,0},turbulence}});
      else if(p.source_type=="back_pressure") bc.push_back({p.name,fv::PressureOutlet{0,Vec3d{1,0,0},fv::BackflowPressure::Static,turbulence}});
      else bc.push_back({p.name,fv::Symmetry{}});
    }
    fv::SolverSettings settings;settings.viscosity=1./3e6;settings.pseudo_time_step=.01;
    if(argc>6) {const auto dt=parse_double(argv[6]);KUTTA_REQUIRE(dt && *dt>0,"invalid pseudo time step");settings.pseudo_time_step=*dt;}
    settings.momentum_relaxation=.5;settings.pressure_relaxation=.3;settings.monitor_interval=100;
    if(mode!="laminar") {
      settings.sst=fv::SstSettings{};
      settings.sst->k_convection=fv::ConvectionScheme::FirstOrderUpwind;
      settings.sst->omega_convection=fv::ConvectionScheme::FirstOrderUpwind;
    }
    fv::SteadySolver solver(dev,input,bc,settings);
    std::vector<Vec3d> u(input.n_cells,{1,0,0});std::vector<double> p(input.n_cells,0),k(input.n_cells,turbulence.k),w(input.n_cells,turbulence.omega);
    if(mode=="hybrid") solver.initialize_hybrid({});
    else if(mode=="laminar") solver.initialize_fields(u,p);
    else solver.initialize_fields(u,p,k,w);
    for(int i=0;i<steps;++i) {
      solver.iterate(1);const auto& r=solver.report();
      std::cout<<"LINEAR iteration="<<r.iteration;
      for(const auto& m:r.momentum) std::cout<<" momentum="<<m.converged<<','<<m.initial_residual<<','<<m.final_residual;
      std::cout<<" pressure="<<r.pressure.converged<<','<<r.pressure.initial_residual<<','<<r.pressure.final_residual<<'\n';
    }
    std::cout<<"END completed="<<solver.report().iteration<<'\n';
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
