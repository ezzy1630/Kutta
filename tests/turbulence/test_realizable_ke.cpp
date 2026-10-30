#include "../backend_helpers.hpp"
#include "support.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/fv/turbulence/realizable_ke.hpp"
#include "kernels.hpp"
using namespace kutta;
using namespace sst_test;
namespace {
bool near(double a,double b,double tolerance=3e-5) {return std::abs(a-b)<=tolerance*std::max(std::abs(b),1e-20);}
std::vector<std::vector<float>> snapshot(Device& d,const RealizableKeFields& f) {
  std::vector<std::vector<float>> result;
  for(const auto* b:{&f.mu_t,&f.cmu,&f.strain,&f.k_source,&f.k_sink,&f.epsilon_source,&f.epsilon_sink,&f.diffusivity_k,&f.diffusivity_epsilon,
      &f.face_diffusivity_k,&f.face_diffusivity_epsilon,&f.wall_coefficient,&f.wall_y_star,&f.wall_tau,&f.wall_epsilon,&f.wall_production,
      &f.wall_shear.x,&f.wall_shear.y,&f.wall_shear.z,&f.fixed_epsilon}) result.push_back(d.download(*b));
  return result;
}
}
TEST_CASE("Realizable KE variable Cmu and epsilon equations match double oracle and all backends") {
  std::vector<std::vector<float>> cpu;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=box(*dev,8);for(auto& p:m.patches) p.kind=PatchKind::Symmetry;State s(*dev,m);
    const std::vector<float> ux{0,0,1,-2,3,0,1e-12f,1e3f},uy{0,5,0,1,2,-10,0,0},vx{0,0,0,-1,1,10,0,0};
    dev->upload<float>(s.gu.x,ux);dev->upload<float>(s.gu.y,uy);dev->upload<float>(s.gv.x,vx);
    dev->fill(s.gv.y,0);dev->fill(s.gw.z,0);
    RealizableKeModel model(*dev,m);model.update({s.velocity,{s.gu,s.gv,s.gw},s.rho,s.mu,s.k,s.omega,s.wall_velocity});
    const auto fields=snapshot(*dev,model.fields());const auto& f=model.fields();
    const auto cmu=dev->download(f.cmu),mt=dev->download(f.mu_t),ks=dev->download(f.k_source),es=dev->download(f.epsilon_source),ed=dev->download(f.epsilon_sink);
    for(int i=0;i<8;++i) {
      double a=ux[i],b=.5*(uy[i]+vx[i]),ss=a*a+2*b*b,rate=std::sqrt(2*ss);
      double W=ss>0?(a*a*a+3*a*b*b)/std::pow(ss,1.5):0;
      double As=std::sqrt(6.)*std::cos(std::acos(std::clamp(std::sqrt(6.)*W,-1.,1.))/3);
      double us=std::sqrt(ss+.5*std::pow(uy[i]-vx[i],2)),c=1/(4.04+As*.5*us),visc=std::min(c*.5,100000*double(float(.001)));
      double eta=rate*.5,C1=std::max(.43,eta/(eta+5));
      CHECK(near(cmu[i],c));CHECK(near(mt[i],visc));CHECK(near(ks[i],(visc+double(float(.001)))*rate*rate));
      CHECK(near(es[i],C1*rate*2));CHECK(near(ed[i],double(float(1.9))*2/(1+std::sqrt(double(float(.001))*2))));
    }
    model.update({s.velocity,{s.gu,s.gv,s.gw},s.rho,s.mu,s.k,s.omega,s.wall_velocity});CHECK(fields==snapshot(*dev,model.fields()));
    if(dev->backend()==Backend::Cpu) cpu=fields;
    else {REQUIRE(cpu.size()==fields.size());for(std::size_t j=0;j<fields.size();++j) for(std::size_t i=0;i<fields[j].size();++i) CHECK(near(fields[j][i],cpu[j][i]));}
    RealizableKeOptions o;o.wall=KeWallTreatment::Enhanced;CHECK_THROWS(RealizableKeModel(*dev,m,o));
  }
}
TEST_CASE("Realizable KE scalable walls clamp every y-star formula and retain moving-wall tangency") {
  std::vector<std::vector<float>> cpu;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=box(*dev);State s(*dev,m);dev->fill(s.velocity.x,15);dev->fill(s.wall_velocity.x,2);
    const std::vector<float> distance{.0001f,.001f,.01f,.1f,1,10};dev->upload<float>(m.wall_normal_distance,distance);
    RealizableKeModel model(*dev,m);model.update({s.velocity,{s.gu,s.gv,s.gw},s.rho,s.mu,s.k,s.omega,s.wall_velocity});const auto& f=model.fields();
    auto coefficient=dev->download(f.wall_coefficient),epsilon=dev->download(f.wall_epsilon),ystar=dev->download(f.wall_y_star);
    auto tx=dev->download(f.wall_shear.x),ty=dev->download(f.wall_shear.y),tz=dev->download(f.wall_shear.z);
    const double us=std::pow(.09,.25),mu=double(float(.001));
    for(int b=0;b<6;++b) {
      const double y=us*distance[b]/mu,effective=std::max(y,11.225),c=us*.4187/std::log(9.793*effective);
      CHECK(near(ystar[b],y));CHECK(near(coefficient[b],c));CHECK(near(epsilon[b],std::pow(us,4)/(.4187*effective*mu)));
      const auto a=dev->download(m.area.x)[b],d=dev->download(m.area.y)[b],e=dev->download(m.area.z)[b];
      CHECK(std::abs(a*tx[b]+d*ty[b]+e*tz[b])<1e-6);
    }
    CHECK(coefficient[0]==coefficient[1]);
    const auto values=snapshot(*dev,f);
    if(dev->backend()==Backend::Cpu) cpu=values;
    else {REQUIRE(values.size()==cpu.size());for(std::size_t j=0;j<values.size();++j) for(std::size_t i=0;i<values[j].size();++i) CHECK(near(values[j][i],cpu[j][i]));}
    dev->fill(s.velocity.x,2);model.update({s.velocity,{s.gu,s.gv,s.gw},s.rho,s.mu,s.k,s.omega,s.wall_velocity});
    CHECK(dev->max_abs(f.wall_tau,6)==0);CHECK(dev->download(f.wall_coefficient)[0]>0);
    dev->fill(s.k,0);model.update({s.velocity,{s.gu,s.gv,s.gw},s.rho,s.mu,s.k,s.omega,s.wall_velocity});
    // A fixed wall row must be attainable after the transport floor is applied.
    CHECK(dev->download(f.fixed_epsilon)[0]==1e-20f);
    for(float value:dev->download(f.wall_epsilon)) CHECK(value==1e-20f);
  }
}
TEST_CASE("Realizable KE explicit EWT algebra and Kader asymptotes are portable") {
  std::vector<std::vector<float>> cpu;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=box(*dev,8);State s(*dev,m);const std::vector<float> distance{1e-7f,1e-5f,.01f,.1f,.2f,1,2,10};
    dev->upload<float>(m.wall_distance,distance);auto outer=filled(*dev,8,.5f);auto f=allocate_ewt_fields(*dev,8);
    evaluate_ewt(*dev,m,s.rho,s.mu,s.k,outer,{20,.09f},f);
    const auto viscosity=dev->download(f.viscosity),epsilon=dev->download(f.epsilon),weight=dev->download(f.outer_weight);
    for(int i=0;i<8;++i) {
      double y=distance[i],re=y/double(float(.001)),cl=.4187*std::pow(double(float(.09)),-.75);
      double x=(re-200)/(20/std::atanh(.98)),tail=std::exp(-2*std::abs(x)),blend=x>=0?1/(1+tail):tail/(1+tail);
      double mt=blend*.5+(1-blend)*double(float(.09))*y*cl*(-std::expm1(-re/70));
      double eps=1/(y*cl*(-std::expm1(-re/(2*cl))));
      CHECK(near(weight[i],blend));CHECK(near(viscosity[i],mt));CHECK(near(epsilon[i],eps));
    }
    auto yp=dev->upload<float>(std::vector<float>{0,.001f,1,10,30,100,300,1000});auto up=dev->alloc<float>(8),du=dev->alloc<float>(8);
    kern::rke_ewt_law(*dev,8,yp,up,du);auto u=dev->download(up),derivative=dev->download(du);
    CHECK(u[0]==0);CHECK(derivative[0]==1);CHECK(near(u[1],.001));CHECK(near(u[7],std::log(9793)/.4187));
    for(int i=0;i<8;++i) {CHECK(std::isfinite(u[i]));CHECK(derivative[i]>0);}
    const std::vector<std::vector<float>> values{viscosity,epsilon,weight,u,derivative};
    if(dev->backend()==Backend::Cpu) cpu=values;
    else for(std::size_t j=0;j<values.size();++j) for(std::size_t i=0;i<values[j].size();++i) {
      INFO("field "<<j<<" index "<<i<<" GPU "<<values[j][i]<<" CPU "<<cpu[j][i]);CHECK(near(values[j][i],cpu[j][i]));
    }
  }
}
TEST_CASE("Realizable KE inlet boundary conversions use epsilon and hydraulic diameter") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=box(*dev);const int n=6;auto kind=dev->alloc<Index>(n);dev->fill(kind,1);
    auto methods=dev->upload<Index>(std::vector<Index>{0,1,2,3,0,0});
    auto a=filled(*dev,n,.05f),b=filled(*dev,n,1),mass=filled(*dev,n,-1),u=filled(*dev,n,2),v=filled(*dev,n,0);
    auto k=filled(*dev,1,.2f),epsilon=filled(*dev,1,.3f),bk=dev->alloc<float>(n),be=dev->alloc<float>(n);auto transport=dev->alloc<Index>(n);
    kern::rke_boundary(*dev,n,m.owner,kind,methods,a,b,mass,u,v,v,k,epsilon,bk,be,transport,0,1,.001f);
    auto kinetic=dev->download(bk),eps=dev->download(be);double ke=1.5*std::pow(2*double(float(.05)),2);
    CHECK(near(kinetic[0],double(float(.05))));CHECK(near(eps[0],1));
    CHECK(near(kinetic[1],ke));CHECK(near(eps[1],std::pow(.09,.75)*std::pow(ke,1.5)));
    CHECK(near(eps[2],.09*ke*ke/double(float(.001))));CHECK(near(eps[3],eps[1]/.07));
  }
}
TEST_CASE("Realizable KE coupled steady path transports positive decaying turbulence on CPU and Metal") {
  using namespace kutta::fv;std::vector<float> cpu;double cpu_force=0;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=make_box_hex(12,2,1,{0,0,0},{1,.2,.1});std::vector<PatchCondition> bc;
    for(const auto& p:m.patches) {
      if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0},KeValues{.01,.01}}});
      else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0,Vec3d{1,0,0},BackflowPressure::Static,KeValues{.01,.01}}});
      else bc.push_back({p.name,Symmetry{}});
    }
    SolverSettings o;o.realizable_ke=RealizableKeSettings{};o.viscosity=.001;
    SteadySolver solver(*dev,m,bc,o);
    solver.initialize_ke_fields(std::vector<Vec3d>(m.n_cells,{1,0,0}),std::vector<double>(m.n_cells,0),std::vector<double>(m.n_cells,.01),std::vector<double>(m.n_cells,.01));
    solver.iterate(30);auto r=solver.monitor();CHECK(std::isfinite(r.k.scaled));CHECK(std::isfinite(r.epsilon.scaled));
    auto k=dev->download(solver.kinetic_energy());for(float x:k) {CHECK(x>0);CHECK(x<.01);}
    for(float x:dev->download(solver.dissipation())) {CHECK(std::isfinite(x));CHECK(x>0);}
    if(dev->backend()==Backend::Cpu) cpu=k;else for(std::size_t i=0;i<k.size();++i) CHECK(near(k[i],cpu[i],1e-4));
    for(auto& condition:bc) if(condition.patch=="ymin" || condition.patch=="ymax") condition.condition=Wall{};
    SteadySolver wall_solver(*dev,m,bc,o);
    wall_solver.initialize_ke_fields(std::vector<Vec3d>(m.n_cells,{1,0,0}),std::vector<double>(m.n_cells,0),std::vector<double>(m.n_cells,.01),std::vector<double>(m.n_cells,.01));
    wall_solver.iterate(20);const auto force=wall_solver.forces(std::array<std::string,2>{"ymin","ymax"}).viscous.x;
    CHECK(force>0);CHECK(std::isfinite(wall_solver.monitor().epsilon.scaled));
    if(dev->backend()==Backend::Cpu) cpu_force=force;else CHECK(near(force,cpu_force,1e-4));
  }
}
