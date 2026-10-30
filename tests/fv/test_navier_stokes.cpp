#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>

#include "../backend_helpers.hpp"
#include "pressure_operator.hpp"
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
struct Errors { double velocity,pressure,wall_pressure,wall_viscous,wall_total; };
Errors manufactured(Device& dev,int n,bool tet,int max_iterations=2000) {
  using namespace kutta::detail;
  auto mesh=tet?make_box_tet(n,n,n,{0,0,0},{1,1,1}):make_box_hex(n,n,1,{0,0,0},{1,1,1./n});
  // Select x=0, 0<=y<=1/2 as real wall patches, keeping every face in its
  // original position. Whole-wall pressure cancels for this MMS.
  const auto initial_geometry=compute_geometry(mesh);
  std::vector<Patch> split;std::vector<std::string> force_patches;
  for(const auto& patch:mesh.patches) {
    if(patch.name!="xmin") {split.push_back(patch);continue;}
    for(Index face=patch.start;face<patch.start+patch.size;) {
      const bool selected=initial_geometry.face_centroid[face].y<.5;
      auto part=patch;part.start=face;part.name="xmin_"+std::to_string(face);
      do {++face;} while(face<patch.start+patch.size &&
        (initial_geometry.face_centroid[face].y<.5)==selected);
      part.size=face-part.start;split.push_back(part);
      if(selected) force_patches.push_back(part.name);
    }
  }
  mesh.patches=std::move(split);
  std::vector<PatchCondition> bc;
  for(auto& patch:mesh.patches)
    if(patch.name=="zmin" || patch.name=="zmax") bc.push_back({patch.name,Symmetry{}});
    else bc.push_back({patch.name,Wall{}});
  SolverSettings settings;settings.viscosity=.1;settings.limit=false;
  settings.momentum_linear={1e-3,1e-10,200};settings.pressure_linear={1e-4,1e-12,300};
  SteadySolver solver(dev,std::move(mesh),bc,settings);
  solver.initialize(exact_velocity,exact_pressure);solver.set_source(source);
  double total_area=0;
  for(auto area:solver.geometry().face_area) total_area+=kutta::detail::norm(area);
  int it=0;
  for(;it<max_iterations;it+=20) {
    solver.iterate(20);const auto r=solver.report();
    if(it%200==0) std::cout<<"FV_V2_PROGRESS backend="<<backend_name(dev.backend())<<" tet="<<tet<<" n="<<n
      <<" iter="<<r.iteration<<" ru="<<r.u.scaled<<" rv="<<r.v.scaled<<" rc="<<r.continuity.scaled<<" rc_raw="<<r.continuity.raw<<std::endl;
    const double mass_roundoff=32*std::numeric_limits<float>::epsilon()*settings.density*r.maximum_velocity*total_area;
    if(std::max({r.u.scaled,r.v.scaled,r.w.scaled})<1e-7 &&
       (r.continuity.scaled<1e-5 || r.continuity.raw<mass_roundoff)) break;
  }
  CHECK(it<max_iterations);
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
  std::vector<bool> boundary(solver.mesh().n_cells),wall(solver.mesh().n_cells);
  for(const auto& patch:solver.mesh().patches) for(Index face=patch.start;face<patch.start+patch.size;++face) {
    boundary[solver.mesh().owner[face]]=true;
    if(patch.kind==PatchKind::Wall) wall[solver.mesh().owner[face]]=true;
  }
  double ei=0,vi=0,ew=0,vw=0,emax=0;Index worst=0;
  for(Index c=0;c<solver.mesh().n_cells;++c) {
    const double error=p[c]-exact_pressure(g.cell_centroid[c])-mean;
    if(!boundary[c]) {ei+=g.cell_volume[c]*error*error;vi+=g.cell_volume[c];}
    if(!wall[c]) {ew+=g.cell_volume[c]*error*error;vw+=g.cell_volume[c];}
    if(std::abs(error)>emax) {emax=std::abs(error);worst=c;}
  }
  const auto point=g.cell_centroid[worst];
  std::cout<<"FV_V2_BOUNDARY backend="<<backend_name(dev.backend())<<" tet="<<tet<<" n="<<n
    <<" interior_p_L2="<<(vi>0?std::sqrt(ei/vi):-1)<<" excluding_walls_p_L2="<<std::sqrt(ew/vw)
    <<" p_Linf="<<emax<<" at="<<point.x<<","<<point.y<<","<<point.z
    <<" boundary="<<boundary[worst]<<" wall="<<wall[worst]<<std::endl;
  eu=std::sqrt(eu/vol);ep=std::sqrt(ep/vol);
  std::cout<<"FV_V2 backend="<<backend_name(dev.backend())<<" tet="<<tet<<" n="<<n
    <<" iterations="<<solver.report().iteration<<" u_L2="<<eu<<" p_L2="<<ep
    <<" ru="<<solver.report().u.scaled<<" rv="<<solver.report().v.scaled<<" rw="<<solver.report().w.scaled
    <<" rc_raw="<<solver.report().continuity.raw<<" corrected_mass="<<solver.report().corrected_mass_imbalance<<std::endl;
  // Analytic wall-on-solid traction is p*n - mu*(grad u + grad u^T)*n.
  // At xmin: pressure=(-cos(pi*y),0,0), viscous=(0,-20*f(y),0).
  // Integrals over y=[0,1/2]: -1/pi and -20*(1/60)=-1/3.
  // These closed-form antiderivatives avoid a reference quadrature error.
  const double span=tet?1.:1./n;
  const auto force=solver.forces(force_patches,mean);
  const Vec3d exact_p{-span/std::numbers::pi,0,0},exact_v{0,-span/3,0};
  const double efp=kutta::detail::norm(force.pressure-exact_p)/span;
  const double efv=kutta::detail::norm(force.viscous-exact_v)/span;
  const double eft=kutta::detail::norm(force.force()-exact_p-exact_v)/span;
  std::cout<<"FV_V2_FORCE backend="<<backend_name(dev.backend())<<" tet="<<tet<<" n="<<n
    <<" pressure_x="<<force.pressure.x/span<<" viscous_y="<<force.viscous.y/span
    <<" pressure_error="<<efp<<" viscous_error="<<efv<<" total_error="<<eft<<std::endl;
  return {eu,ep,efp,efv,eft};
}
}

void check_refinement(Device& dev,bool tet) {
  Errors old{};double old_operator=0;
  for(int n:tet?std::vector<int>{4,8,16}:std::vector<int>{8,16,32}) {
    const auto e=manufactured(dev,n,tet);
    const double momentum=test::pressure_operator_error(dev,n,tet).momentum;
    if(old.velocity>0) {
      const double ou=test::order(old.velocity,e.velocity),op=test::order(old.pressure,e.pressure);
      const double operator_order=test::order(old_operator,momentum);
      const double ofp=test::order(old.wall_pressure,e.wall_pressure);
      const double ofv=test::order(old.wall_viscous,e.wall_viscous);
      const double oft=test::order(old.wall_total,e.wall_total);
      std::cout<<"FV_V2_ORDER backend="<<backend_name(dev.backend())<<" tet="<<tet<<" n="<<n
        <<" u="<<ou<<" p="<<op<<" momentum_operator="<<operator_order
        <<" wall_pressure="<<ofp<<" wall_viscous="<<ofv<<" wall_total="<<oft<<std::endl;
      if(n==(tet?16:32)) {
        CHECK(ou>=1.9);
        if(tet) {CHECK(operator_order>=.9);CHECK(op>=operator_order-.1);}
        else CHECK(op>=1.8);
        CHECK(ofp>=1.8);CHECK(ofv>=1.8);CHECK(oft>=1.8);
      }
    }
    old=e;old_operator=momentum;
  }
}

TEST_CASE("FV V2 manufactured Navier Stokes hex refinement" * doctest::test_suite("slow")) {
  KUTTA_FOR_EACH_BACKEND(dev) check_refinement(*dev,false);
}
TEST_CASE("FV V2 manufactured Navier Stokes tet refinement" * doctest::test_suite("slow")) {
  KUTTA_FOR_EACH_BACKEND(dev) check_refinement(*dev,true);
}

TEST_CASE("FV laminar tet manufactured solution converges without SST") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    // The 384-cell V2 case converges in 160 steps with the verified laminar
    // operator. Omitting extrapolated pressure observations stalls this solve.
    const auto error=manufactured(*dev,4,true,200);
    CHECK(error.velocity<.18);
    CHECK(error.pressure<.28);
  }
}
