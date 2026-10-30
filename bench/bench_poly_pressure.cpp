#include <iostream>
#include <iomanip>
#include "kutta/core/parse.hpp"
#include "../tests/fv/physics_meshes.hpp"
#include "kutta/fv/steady_solver.hpp"
using namespace kutta;
using namespace kutta::fv;
using namespace kutta::detail;
int main(int argc,char** argv) {
  const auto rings=argc>1?parse_double(argv[1]):std::optional<double>(8);
  const auto rotation=argc>2?parse_double(argv[2]):std::optional<double>(0);
  const auto steps=argc>3?parse_double(argv[3]):std::optional<double>(1200);
  const auto viscosity=argc>4?parse_double(argv[4]):std::optional<double>(1);
  if(argc>5 || !rings || !rotation || !steps || !viscosity ||
      *rings<4 || *rings>128 || std::fmod(*rings,2)!=0 ||
      (*rotation!=0 && *rotation!=1) || *steps<1 || *steps>100000 ||
      std::floor(*steps)!=*steps || !std::isfinite(*viscosity) || *viscosity<=0) {
    std::cerr<<"Usage: bench_poly_pressure [even rings 4..128] [solid 0|1] [steps 1..100000] [positive viscosity]\n";
    return 2;
  }
  const int nr=static_cast<int>(*rings),iterations=static_cast<int>(*steps);
  const bool solid=*rotation==1;
  auto pressure=[&](Vec3d x) {const auto r=std::hypot(x.x,x.y);return solid?.5*r*r:r*r/18.-8./9*std::log(r)-8./(9*r*r);};
  auto exact=[&](Vec3d x) {const double a=solid?1.:-1./3+4./(3*(x.x*x.x+x.y*x.y));return Vec3d{-a*x.y,a*x.x,0};};
  for(auto backend:Device::available()) {
    auto dev=Device::create(backend);auto mesh=test::physics_annulus(nr,8*nr,true);std::vector<PatchCondition> bc;
    for(const auto& patch:mesh.patches) {
      if(patch.name=="inner" || (solid && patch.name=="outer")) bc.push_back({patch.name,Wall{{},Rotation{{},{0,0,1}}}});
      else if(patch.name=="outer") bc.push_back({patch.name,Wall{}});
      else bc.push_back({patch.name,Symmetry{}});
    }
    SolverSettings opt;opt.viscosity=*viscosity;opt.limit=false;opt.pseudo_time_step=.1;
    opt.momentum_linear={1e-5,1e-11,200};opt.pressure_linear={1e-5,1e-12,300};
    SteadySolver s(*dev,mesh,bc,opt);s.initialize(exact,pressure);s.iterate(iterations);
    const auto r=s.monitor();const auto p=dev->download(s.pressure());const auto u=dev->download(s.velocity().x),v=dev->download(s.velocity().y);
    double mean=0,volume=0,ep=0,eu=0;
    std::vector<double> means(nr),weights(nr),lo(nr,1e30),hi(nr,-1e30);
    for(Index c=0;c<mesh.n_cells;++c) {const auto x=s.geometry().cell_centroid[c];const auto vol=s.geometry().cell_volume[c];mean+=vol*(p[c]-pressure(x));volume+=vol;}
    mean/=volume;
    for(Index c=0;c<mesh.n_cells;++c) {
      const auto x=s.geometry().cell_centroid[c];const auto vol=s.geometry().cell_volume[c];const auto d=p[c]-pressure(x)-mean;
      ep+=vol*d*d;const auto vel=exact(x);eu+=vol*(std::pow(u[c]-vel.x,2)+std::pow(v[c]-vel.y,2));
      const auto ring=std::clamp(int((std::hypot(x.x,x.y)-1)*nr),0,nr-1);means[ring]+=vol*d;weights[ring]+=vol;lo[ring]=std::min(lo[ring],d);hi[ring]=std::max(hi[ring],d);
    }
    std::cout<<std::setprecision(10)<<"POLY_SOLVE backend="<<backend_name(backend)<<" nr="<<nr<<" cells="<<mesh.n_cells<<" solid="<<solid<<" viscosity="<<*viscosity<<" steps="<<iterations
      <<" pressure_L2="<<std::sqrt(ep/volume)/std::abs(pressure({2,0,0})-pressure({1,0,0}))<<" velocity_L2="<<std::sqrt(eu/volume)
      <<" momentum="<<r.u.scaled<<" continuity_raw="<<r.continuity.raw<<" pressure_rings=";
    for(int i=0;i<nr;++i) std::cout<<means[i]/weights[i]<<":"<<lo[i]<<":"<<hi[i]<<",";std::cout<<std::endl;
  }
}
