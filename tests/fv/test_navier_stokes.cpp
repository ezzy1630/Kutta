#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>

#include "../backend_helpers.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

using namespace kutta;
using namespace kutta::fv;

namespace {
double f(double x) {return x*x*(1-x)*(1-x);}
double df(double x) {return 2*x-6*x*x+4*x*x*x;}
double ddf(double x) {return 2-12*x+12*x*x;}
double dddf(double x) {return -12+24*x;}
Vec3d exact_velocity(Vec3d p) {return {100*f(p.x)*df(p.y),-100*df(p.x)*f(p.y),0};}
double exact_pressure(Vec3d p) {return std::cos(std::numbers::pi*p.x)*std::cos(std::numbers::pi*p.y);}
Vec3d source(Vec3d p) {
  const auto q=exact_velocity(p);const double pi=std::numbers::pi,mu=.1;
  const double ux=100*df(p.x)*df(p.y),uy=100*f(p.x)*ddf(p.y);
  const double vx=-100*ddf(p.x)*f(p.y),vy=-100*df(p.x)*df(p.y);
  const double lapu=100*(ddf(p.x)*df(p.y)+f(p.x)*dddf(p.y));
  const double lapv=-100*(dddf(p.x)*f(p.y)+df(p.x)*ddf(p.y));
  return {q.x*ux+q.y*uy-mu*lapu-pi*std::sin(pi*p.x)*std::cos(pi*p.y),
    q.x*vx+q.y*vy-mu*lapv-pi*std::cos(pi*p.x)*std::sin(pi*p.y),0};
}
std::pair<double,double> manufactured(Device& dev,int n,bool tet) {
  auto mesh=tet?make_box_tet(n,n,n,{0,0,0},{1,1,1}):make_box_hex(n,n,1,{0,0,0},{1,1,1./n});
  std::vector<PatchCondition> bc;
  for(auto& patch:mesh.patches)
    if(patch.name=="zmin" || patch.name=="zmax") bc.push_back({patch.name,Symmetry{}});
    else bc.push_back({patch.name,Wall{}});
  SolverSettings settings;settings.viscosity=.1;settings.pseudo_time_scale=3;settings.limit=false;
  settings.momentum_linear={1e-3,1e-10,200};settings.pressure_linear={1e-4,1e-12,300};
  SteadySolver solver(dev,std::move(mesh),bc,settings);
  solver.initialize(exact_velocity,exact_pressure);solver.set_source(source);
  double total_area=0;
  for(auto area:solver.geometry().face_area) total_area+=kutta::detail::norm(area);
  int it=0;
  for(;it<2000;it+=20) {
    solver.iterate(20);const auto r=solver.report();
    if(it%200==0) std::cout<<"FV_V2_PROGRESS backend="<<backend_name(dev.backend())<<" tet="<<tet<<" n="<<n
      <<" iter="<<r.iteration<<" ru="<<r.u.scaled<<" rv="<<r.v.scaled<<" rc="<<r.continuity.scaled<<" rc_raw="<<r.continuity.raw<<std::endl;
    const double mass_roundoff=32*std::numeric_limits<float>::epsilon()*settings.density*r.maximum_velocity*total_area;
    if(std::max({r.u.scaled,r.v.scaled,r.w.scaled})<1e-7 &&
       (r.continuity.scaled<1e-5 || r.continuity.raw<mass_roundoff)) break;
  }
  CHECK(it<2000);
  auto u=dev.download(solver.velocity().x),v=dev.download(solver.velocity().y),w=dev.download(solver.velocity().z),p=dev.download(solver.pressure());
  const auto& g=solver.geometry();
  double eu=0,ep=0,vol=0,mean=0;
  for(Index c=0;c<solver.mesh().n_cells;++c) {mean+=g.cell_volume[c]*(p[c]-exact_pressure(g.cell_centroid[c]));vol+=g.cell_volume[c];}
  mean/=vol;
  for(Index c=0;c<solver.mesh().n_cells;++c) {
    auto exact=exact_velocity(g.cell_centroid[c]);
    eu+=g.cell_volume[c]*(std::pow(u[c]-exact.x,2)+std::pow(v[c]-exact.y,2)+std::pow(w[c]-exact.z,2));
    ep+=g.cell_volume[c]*std::pow(p[c]-exact_pressure(g.cell_centroid[c])-mean,2);
  }
  eu=std::sqrt(eu/vol);ep=std::sqrt(ep/vol);
  std::cout<<"FV_V2 backend="<<backend_name(dev.backend())<<" tet="<<tet<<" n="<<n
    <<" iterations="<<solver.report().iteration<<" u_L2="<<eu<<" p_L2="<<ep
    <<" ru="<<solver.report().u.scaled<<" rv="<<solver.report().v.scaled<<" rw="<<solver.report().w.scaled
    <<" rc_raw="<<solver.report().continuity.raw<<" corrected_mass="<<solver.report().corrected_mass_imbalance<<std::endl;
  return {eu,ep};
}
}

void check_refinement(Device& dev,bool tet) {
  double old_u=0,old_p=0;
  for(int n:tet?std::vector<int>{4,8,16}:std::vector<int>{8,16,32}) {
    auto [eu,ep]=manufactured(dev,n,tet);
    if(old_u>0) {
      double ou=std::log2(old_u/eu),op=std::log2(old_p/ep);
      std::cout<<"FV_V2_ORDER backend="<<backend_name(dev.backend())<<" tet="<<tet<<" n="<<n
        <<" u="<<ou<<" p="<<op<<std::endl;
      if(n==(tet?16:32)) {CHECK(ou>=1.9);CHECK(op>=1.8);}
    }
    old_u=eu;old_p=ep;
  }
}
TEST_CASE("FV V2 manufactured Navier Stokes hex refinement") {
  KUTTA_FOR_EACH_BACKEND(dev) check_refinement(*dev,false);
}
TEST_CASE("FV V2 manufactured Navier Stokes tet refinement") {
  KUTTA_FOR_EACH_BACKEND(dev) check_refinement(*dev,true);
}
