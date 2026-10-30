#include "../backend_helpers.hpp"
#include "support.hpp"
#include "kutta/fv/turbulence/sa.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kernels.hpp"

using namespace kutta;
using namespace sst_test;
namespace {
bool near(double a,double b) {return std::abs(a-b)<=3e-5*std::max(std::abs(b),1e-20);}
std::vector<std::vector<float>> snapshot(Device& d,const SaFields& f) {
  std::vector<std::vector<float>> result;
  for(auto* b:{&f.mu_t,&f.source,&f.sink,&f.diffusivity,&f.face_diffusivity,&f.fv1,&f.fv2,&f.modified_rate,&f.fw}) result.push_back(d.download(*b));
  return result;
}
}
TEST_CASE("SA transport coefficients match independent double equations and backend oracle") {
  std::vector<std::vector<float>> cpu;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=box(*dev,8);State s(*dev,m);
    const std::vector<float> q{0,1e-9f,1e-4f,.001f,.003f,.01f,1,100};
    dev->upload<float>(s.k,q);dev->fill(s.gu.y,.00002f);dev->fill(s.gk.x,.1f);
    auto strain=filled(*dev,8,20),rotation=filled(*dev,8,20);
    for(auto form:{SaFormulation::Fluent,SaFormulation::TmrStandard}) {
      SaOptions o;o.formulation=form;SaModel model(*dev,m,o);
      model.update({s.rho,s.mu,s.k,s.gk,strain,rotation});const auto& f=model.fields();
      const auto fv=dev->download(f.fv1),source=dev->download(f.source),sink=dev->download(f.sink);
      const auto rate=dev->download(f.modified_rate),fw=dev->download(f.fw),mt=dev->download(f.mu_t);
      const auto distance=dev->download(m.wall_distance);
      for(int i=0;i<8;++i) {
        const double v=q[i],chi=v/double(float(.001)),cv=o.cv1,kp=o.kappa;
        const double f1=std::pow(chi,3)/(std::pow(chi,3)+std::pow(cv,3));
        double st=20+v/(kp*kp*distance[i]*distance[i])*(1-chi/(1+chi*f1));
        if(form==SaFormulation::TmrStandard) st=std::max(st,6.);
        double r=v==0?0:v/(st*kp*kp*distance[i]*distance[i]);
        if(form==SaFormulation::TmrStandard) r=std::min(r,10.);
        double g=r+o.cw2*(std::pow(r,6)-r),wall=g*std::pow((1+std::pow(o.cw3,6))/(std::pow(g,6)+std::pow(o.cw3,6)),1./6);
        double ft=form==SaFormulation::TmrStandard?1.2*std::exp(-.5*chi*chi):0;
        double growth=o.cb1*(1-ft)*st-(o.cb1/(kp*kp)+(1+o.cb2)/o.sigma)*wall*v/(distance[i]*distance[i])+o.cb1*ft/(kp*kp)*v/(distance[i]*distance[i]);
        CHECK(near(fv[i],f1));CHECK(near(rate[i],st));CHECK(near(fw[i],wall));
        const auto local_source=[&](double q) {
          const double chi=q/double(float(.001)),f1=std::pow(chi,3)/(std::pow(chi,3)+std::pow(cv,3));
          double rate=20+q/(kp*kp*distance[i]*distance[i])*(1-chi/(1+chi*f1));
          if(form==SaFormulation::TmrStandard) rate=std::max(rate,6.);
          double r=q==0?0:q/(rate*kp*kp*distance[i]*distance[i]);
          r=std::clamp(r,-1e4,form==SaFormulation::TmrStandard?10.:1e4);
          const double g=r+o.cw2*(std::pow(r,6)-r);
          const double fw=g*std::pow((1+std::pow(o.cw3,6))/(std::pow(g,6)+std::pow(o.cw3,6)),1./6);
          const double ft=form==SaFormulation::TmrStandard?1.2*std::exp(-.5*chi*chi):0;
          return (o.cb1*(1-ft)*rate-((o.cb1/(kp*kp)+(1+o.cb2)/o.sigma)*fw-o.cb1*ft/(kp*kp))*q/(distance[i]*distance[i]))*q;
        };
        const double h=v>0?v*1e-5:1e-11;
        const double derivative=v>0?(local_source(v+h)-local_source(v-h))/(2*h):(local_source(h)-local_source(0))/h;
        const double implicit=std::max({-growth,-derivative,0.});
        CHECK(near(source[i],(growth+implicit)*v+o.cb2/o.sigma*double(float(.1))*double(float(.1))));
        CHECK(near(sink[i],implicit));CHECK(near(mt[i],std::min(v*f1,double(o.max_viscosity_ratio)*double(float(.001)))));
      }
      auto first=snapshot(*dev,f);model.update({s.rho,s.mu,s.k,s.gk,strain,rotation});CHECK(first==snapshot(*dev,f));
      if(form==SaFormulation::Fluent) {
        if(dev->backend()==Backend::Cpu) cpu=first;
        else {REQUIRE(cpu.size()==first.size());for(std::size_t a=0;a<first.size();++a) for(std::size_t i=0;i<first[a].size();++i) CHECK(near(first[a][i],cpu[a][i]));}
      }
    }
    SaOptions alternate;alternate.production=SaProduction::StrainVorticity;
    SaModel model(*dev,m,alternate);dev->fill(strain,0);dev->fill(rotation,20);
    model.update({s.rho,s.mu,s.k,s.gk,strain,rotation});CHECK(dev->download(model.fields().modified_rate)[0]==-20);
    alternate.wall=SaWallTreatment::FluentAutomatic;CHECK_THROWS(SaModel(*dev,m,alternate));
  }
}
TEST_CASE("SA boundary methods invert damped viscosity and distinguish outlet backflow") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=box(*dev);const int n=m.n_boundary_faces();
    std::vector<Index> kinds{0,1,2,2,3,4},spec(n,1);
    std::vector<float> flow{0,-1,1,-1,0,-1};
    auto kind=dev->upload<Index>(kinds),method=dev->upload<Index>(spec);
    auto mass=dev->upload<float>(flow),value=filled(*dev,n,10),q=filled(*dev,1,.001f);
    auto b=dev->alloc<float>(n);auto transport=dev->alloc<Index>(n);
    kern::sa_boundary(*dev,n,m.owner,kind,method,value,value,value,value,value,mass,q,b,transport,0,.001f,7.1f);
    auto bc=dev->download(b);CHECK(bc[0]==0);CHECK(bc[2]==.001f);CHECK(bc[4]==.001f);
    for(int j:{1,3,5}) {double chi=bc[j]/.001;CHECK(near(chi*std::pow(chi,3)/(std::pow(chi,3)+std::pow(7.1,3)),10));}
    CHECK(dev->download(transport)==std::vector<Index>{0,0,1,0,1,0});
    dev->fill(method,0);dev->fill(value,.007f);
    kern::sa_boundary(*dev,n,m.owner,kind,method,value,value,value,value,value,mass,q,b,transport,0,.001f,7.1f);
    CHECK(dev->download(b)[1]==.007f);
    dev->fill(value,0);dev->fill(method,1);
    kern::sa_boundary(*dev,n,m.owner,kind,method,value,value,value,value,value,mass,q,b,transport,0,.001f,7.1f);
    CHECK(dev->download(b)[1]==0);
    dev->fill(method,2);dev->fill(value,.05f);
    auto length=filled(*dev,n,.07f),u=filled(*dev,n,2),zero=filled(*dev,n,0);
    kern::sa_boundary(*dev,n,m.owner,kind,method,value,length,u,zero,zero,mass,q,b,transport,0,.001f,7.1f);
    CHECK(near(dev->download(b)[1],std::pow(.09,.25)*std::sqrt(1.5)*.05*.07*2));
  }
}
TEST_CASE("SA real transport path preserves laminar Couette and rejects conflicting model data") {
  using namespace kutta::fv;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=make_box_hex(2,12,1,{0,0,0},{1,1,.1});std::vector<PatchCondition> bc;
    for(const auto& p:m.patches) {
      if(p.name=="ymax") bc.push_back({p.name,Wall{{1,0,0}}});
      else if(p.name=="ymin") bc.push_back({p.name,Wall{}});
      else if(p.name=="xmin" || p.name=="xmax") bc.push_back({p.name,VelocityInlet{std::function<Vec3d(Vec3d)>([](Vec3d x){return Vec3d{x.y,0,0};}),ModifiedViscosity{0}}});
      else bc.push_back({p.name,Symmetry{}});
    }
    SolverSettings o;o.sa=SaSettings{};o.viscosity=.01;
    SteadySolver solver(*dev,m,bc,o);
    std::vector<Vec3d> u;std::vector<double> p(m.n_cells,0),q(m.n_cells,0);
    for(auto x:compute_geometry(m).cell_centroid) u.push_back({x.y,0,0});
    solver.initialize_sa_fields(u,p,q);solver.iterate(12);auto r=solver.monitor();
    CHECK(r.nu_tilde.raw==0);CHECK(dev->max_abs(solver.modified_viscosity(),m.n_cells)==0);
    auto ux=dev->download(solver.velocity().x);
    for(Index c=0;c<m.n_cells;++c) CHECK(std::abs(ux[c]-solver.geometry().cell_centroid[c].y)<2e-5);
    CHECK(near(solver.forces(std::array<std::string,1>{"ymax"}).viscous.x,-.001));
    const auto stress=dev->download(solver.sa_wall().shear.x);
    for(const auto& patch:solver.mesh().patches) if(patch.name=="ymax")
      for(Index f=patch.start;f<patch.start+patch.size;++f) CHECK(near(stress[f-solver.mesh().n_internal_faces()],.01));
    o.sst=SstSettings{};CHECK_THROWS(SteadySolver(*dev,m,bc,o));
  }
}
