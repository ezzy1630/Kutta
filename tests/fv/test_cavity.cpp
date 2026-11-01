#include "../fluent_examples.hpp"
#include <array>
#include <iomanip>
#include "../backend_helpers.hpp"
#include "../../validation/cases/cavity.hpp"

using namespace kutta;
namespace {
const std::filesystem::path data=std::filesystem::path(KUTTA_SOURCE_DIR)/"validation/data";
// Use actual cell widths (65,129,257 are not exact factors of two).
double refinement_order(double e0,double e1,int n0,int n1) {
  return std::log(e0/e1)/std::log(double(n1)/n0);
}
double extrapolate(double coarse,double fine,double ratio,double order) {
  return fine+(fine-coarse)/(std::pow(ratio,order)-1);
}
}

TEST_CASE("FV cavity primary vortex diagnostic recovers analytic streamfunction") {
  std::array<double,2> psi_error{},omega_error{};
  for(int level=0;level<2;++level) {
    const int n=level?33:17;
    const auto mesh=make_box_hex(n,n,1,{0,0,0},{1,1,1});
    const auto g=compute_geometry(mesh);
    const auto f=[](double x){return x*x*(1-x)*(1-x);};
    const auto df=[](double x){return 2*x-6*x*x+4*x*x*x;};
    std::vector<float> u,v;
    for(auto p:g.cell_centroid) {u.push_back(-f(p.x)*df(p.y));v.push_back(df(p.x)*f(p.y));}
    const auto vortex=validation::primary_vortex(g,u,v,n);
    CHECK(vortex.x==doctest::Approx(.5).epsilon(1e-6));
    CHECK(vortex.y==doctest::Approx(.5).epsilon(1e-6));
    psi_error[level]=std::abs(vortex.psi-1./256);
    omega_error[level]=std::abs(vortex.omega-.125);
  }
  CHECK(refinement_order(psi_error[0],psi_error[1],17,33)>=1.9);
  CHECK(refinement_order(omega_error[0],omega_error[1],17,33)>=1.9);
}

TEST_CASE("FV V3 Re1000 cavity refinement against Erturk" * doctest::test_suite("slow")) {
  KUTTA_FOR_EACH_BACKEND(dev) {
    std::array<validation::CavityResult,3> results;
    std::array<double,3> eu{},ev{};
    const std::array<int,3> grids{65,129,257};
    for(int level=0;level<3;++level) {
      const int n=grids[level];
      auto& r=results[level];
      r=validation::cavity(*dev,make_box_hex(n,n,1,{0,0,0},{1,1,1}),1000,10000,n,2e-6);
      REQUIRE(r.converged);
      eu[level]=validation::reference_error(r.u,data/"erturk2005_cavity/erturk2005_u_vertical_centerline_Re1000.csv",0,1);
      ev[level]=validation::reference_error(r.v,data/"erturk2005_cavity/erturk2005_v_horizontal_centerline_Re1000.csv",0,1);
      validation::compare_cavity(r,data,1000,dev->backend(),std::to_string(n));
      std::cout<<std::setprecision(12)<<"FV_CAVITY_REFINEMENT backend="<<backend_name(dev->backend())
        <<" n="<<n<<" u_error="<<eu[level]<<" v_error="<<ev[level]
        <<" u_center="<<validation::interpolate(r.u,.5)<<" v_center="<<validation::interpolate(r.v,.5)
        <<" psi="<<r.vortex.psi<<" omega="<<r.vortex.omega<<" x="<<r.vortex.x<<" y="<<r.vortex.y;
      if(level) std::cout<<" u_order="<<refinement_order(eu[level-1],eu[level],grids[level-1],n)
        <<" v_order="<<refinement_order(ev[level-1],ev[level],grids[level-1],n);
      std::cout<<std::endl;
    }
    CHECK(eu[2]<eu[1]);CHECK(ev[2]<ev[1]);
    CHECK(refinement_order(eu[1],eu[2],129,257)>=1.6);
    CHECK(refinement_order(ev[1],ev[2],129,257)>=1.6);
    // Second-order Richardson estimate on the finest pair; the measured
    // three-grid profile slopes above check that this assumption is credible.
    const auto& a=results[1];const auto& b=results[2];
    const auto rich=[](double x,double y){return extrapolate(x,y,257./129,2);};
    std::cout<<"FV_CAVITY_RICHARDSON backend="<<backend_name(dev->backend())
      <<" u_center="<<rich(validation::interpolate(a.u,.5),validation::interpolate(b.u,.5))
      <<" v_center="<<rich(validation::interpolate(a.v,.5),validation::interpolate(b.v,.5))
      <<" psi="<<rich(a.vortex.psi,b.vortex.psi)<<" omega="<<rich(a.vortex.omega,b.vortex.omega)
      <<" x="<<rich(a.vortex.x,b.vortex.x)<<" y="<<rich(a.vortex.y,b.vortex.y)<<std::endl;
    // Preserve the source rows (including blank centre coordinates) verbatim
    // in the report: extrapolated and spectral references are diagnostics.
    std::ifstream reference(data/"erturk2005_cavity/erturk2005_primary_vortex_Re1000.csv");
    REQUIRE(reference.good());
    for(std::string line;std::getline(reference,line);) std::cout<<"FV_CAVITY_PRIMARY_REFERENCE "<<line<<'\n';
  }
}

TEST_CASE("FV V3 Re100 and Re400 cavity reported deviations" * doctest::test_suite("slow")) {
  KUTTA_FOR_EACH_BACKEND(dev) for(int re:{100,400}) {
    auto r=validation::cavity(*dev,make_box_hex(129,129,1,{0,0,0},{1,1,1}),re);
    CHECK(r.converged);
    validation::compare_cavity(r,data,re,dev->backend(),"129x129");
  }
}

TEST_CASE("FV V4 Fluent cavity matches its actual grid" * doctest::test_suite("slow")) {
  const auto file=fluent_examples()/"lid_driven_cavity/lid_driven_cavity.msh";
  if(!std::filesystem::exists(file)) {MESSAGE("SKIP: external Fluent cavity fixture absent");return;}
  auto mesh=read_fluent_msh(file);
  std::map<long long,double> x,y;
  for(auto p:mesh.points) {x.emplace(std::llround(p.x*1e10),p.x);y.emplace(std::llround(p.y*1e10),p.y);}
  std::vector<double> xs,ys;
  for(auto [key,value]:x) xs.push_back(value);
  for(auto [key,value]:y) ys.push_back(value);
  const int nx=int(xs.size())-1,ny=int(ys.size())-1;
  REQUIRE(nx*ny==mesh.n_cells);
  auto matched=make_box_hex(nx,ny,1,{0,0,0},{1,1,1});
  for(auto& p:matched.points) {p.x=xs[std::llround(p.x*nx)];p.y=ys[std::llround(p.y*ny)];}
  KUTTA_FOR_EACH_BACKEND(dev) for(int re:{100,400,1000}) {
    auto imported=validation::cavity(*dev,mesh,re);
    auto generated=validation::cavity(*dev,matched,re);
    CHECK(imported.converged);CHECK(generated.converged);
    validation::compare_cavity(imported,data,re,dev->backend(),"Fluent");
    for(const auto& [coordinate,value]:imported.u)
      CHECK(std::abs(value-validation::interpolate(generated.u,coordinate))<2e-6);
    for(const auto& [coordinate,value]:imported.v)
      CHECK(std::abs(value-validation::interpolate(generated.v,coordinate))<2e-6);
  }
}
