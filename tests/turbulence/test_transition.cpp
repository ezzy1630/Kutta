// Fluent's Intermittency Transition Model (src/kutta/fv/turbulence/transition.md):
// kernels against independent double-precision formulas, the SST coupling,
// the wall-normal direction and the Fluent setup adapter.
#include <algorithm>
#include <cmath>
#include "../backend_helpers.hpp"
#include "../fluent_examples.hpp"
#include "support.hpp"
#include "kernels.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/setup/convergence.hpp"
#include "kutta/setup/fluent.hpp"
#include "kutta/setup/settings.hpp"

using namespace kutta;
using namespace sst_test;
namespace {
bool near(double a,double b,double tolerance=2e-5,double floor=1e-30) {
  return std::abs(a-b)<=tolerance*std::max({std::abs(a),std::abs(b),floor});
}
// Fluent TG 4.8.2 defaults (CPG in Fluent's rescaled lambda units).
struct Constants {double ctu1=100,ctu2=1000,ctu3=1,cpg1=1,cpg2=-.5,cpg3=0,lim1=1.5,lim2=3;};
double fpg(double lambda,const Constants& c) {
  const double f=lambda>=0?std::min(1+c.cpg1*lambda,c.lim1):std::min(1+c.cpg2*lambda+c.cpg3*std::min(lambda+1,0.),c.lim2);
  return std::max(f,0.);
}
struct GammaTerms {double production,destruction,critical,onset;};
GammaTerms gamma_terms(double rho,double mu,double k,double w,double d,double dvdy,double s,double vorticity,double g,const Constants& c) {
  const double tu=std::min(100*std::sqrt(2*k/3)/(w*d),100.);
  const double lambda=std::clamp(-.1111*dvdy*d*d/(mu/rho)+.1875,-10.,10.);
  const double critical=c.ctu1+c.ctu2*std::exp(-c.ctu3*tu*fpg(lambda,c));
  const double rev=rho*d*d*s/mu,rt=rho*k/(mu*w);
  const double onset=std::max(std::min(rev/(2.2*critical),2.)-std::max(1-std::pow(rt/3.5,3),0.),0.);
  const double turb=std::exp(-std::pow(rt/2,4));
  return {100*rho*s*g*(1-g)*onset,.06*rho*vorticity*g*turb*(50*g-1),critical,onset};
}
// P_k^lim (Menter et al. 2015; NASA TMR SST-2003-Menter-gamma-2015).
double production_limit(double g,double rho,double mu,double mut,double d,double s,double vorticity,double ck,double csep) {
  const double on=std::min(std::max(rho*d*d*s/mu/(2.2*1100)-1,0.),3.);
  return 5*ck*std::max(g-.2,0.)*(1-g)*on*std::max(3*csep*mu-mut,0.)*s*vorticity;
}
}

TEST_CASE("Intermittency sources match independent double formulas") {
  // Cells span laminar (R_T small), onset (Re_V large), freestream (Omega~0),
  // favourable/adverse/clipped lambda and gamma at both ends.
  const std::vector<double> rho{1,1,1.2,1,1,1,1,1,1,1,1,1},mu{1e-5,1e-5,2e-5,1e-5,1e-5,1e-3,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5};
  const std::vector<double> k{1e-6,1e-4,1e-3,3e-3,1e-8,1e-2,1e-5,1e-4,1e-4,1e-4,2e-3,1e-4};
  const std::vector<double> w{50,10,400,20,1,2,100,30,30,30,1e4,30};
  const std::vector<double> d{1e-4,2e-3,5e-4,1e-2,3e-3,.1,1e-3,2e-3,2e-3,2e-3,1e-5,2e-3};
  const std::vector<double> dvdy{0,-20,5,0,300,0,-1e4,1e3,-50,50,0,8};
  const std::vector<double> s{2000,300,800,1,500,10,1e4,300,300,300,5e4,300};
  const std::vector<double> o{2000,300,800,0,500,10,1e4,300,300,300,5e4,300};
  const std::vector<double> g{.02,.5,1,1,1e-7,.9,.3,.6,.6,.6,.05,.6};
  const std::size_t n=k.size();
  for(const auto& c:{Constants{},Constants{150,800,1.3,1.2,2,-.7,1.7,2.5}}) {
    std::vector<std::vector<float>> backends;
    KUTTA_FOR_EACH_BACKEND(dev) {
      auto R=scalar(*dev,rho),M=scalar(*dev,mu),K=scalar(*dev,k),W=scalar(*dev,w),D=scalar(*dev,d),V=scalar(*dev,dvdy);
      auto S=scalar(*dev,s),O=scalar(*dev,o),G=scalar(*dev,g);
      auto source=dev->alloc<float>(n),sink=dev->alloc<float>(n),critical=dev->alloc<float>(n),onset=dev->alloc<float>(n);
      kern::transition_sources(*dev,n,R,M,K,W,D,V,S,O,G,source,sink,critical,onset,float(c.ctu1),float(c.ctu2),float(c.ctu3),
        float(c.cpg1),float(c.cpg2),float(c.cpg3),float(c.lim1),float(c.lim2),1e-14f,1e-20f);
      const auto a=dev->download(source),b=dev->download(sink),rc=dev->download(critical),on=dev->download(onset);
      for(std::size_t i=0;i<n;++i) {
        INFO("cell "<<i);
        // Float inputs are what the kernel sees.
        const auto f=[](double x){return double(float(x));};
        const auto t=gamma_terms(f(rho[i]),f(mu[i]),f(k[i]),f(w[i]),f(d[i]),f(dvdy[i]),f(s[i]),f(o[i]),f(g[i]),c);
        CHECK(near(rc[i],t.critical));CHECK(near(on[i],t.onset,2e-5,1e-6));
        CHECK(a[i]>=0);CHECK(b[i]>=0);
        // The linearized pair reproduces P - E at the current gamma, to float
        // rounding of its two (possibly cancelling) terms.
        const double scale=std::max({t.production,std::abs(t.destruction),double(a[i]),1e-30});
        CHECK(std::abs(a[i]-b[i]*f(g[i])-(t.production-t.destruction))<=3e-5*scale);
      }
      std::vector<float> all;for(const auto* v:{&a,&b,&rc,&on}) all.insert(all.end(),v->begin(),v->end());
      backends.push_back(std::move(all));
    }
    for(std::size_t j=1;j<backends.size();++j)
      for(std::size_t i=0;i<backends[0].size();++i) CHECK(near(backends[j][i],backends[0][i],1e-5,1e-12));
  }
  // F_PG floor and branches, and the lambda bound, in Fluent's units.
  Constants negative;negative.cpg2=2;
  CHECK(fpg(-3,negative)==0);CHECK(fpg(.2,Constants{})==1.2);CHECK(fpg(5,Constants{})==1.5);CHECK(fpg(-1,Constants{})==1.5);
  KUTTA_FOR_EACH_BACKEND(dev) {
    std::vector<double> lam_dvdy{-1e6,1e6,0};
    auto one=scalar(*dev,{1,1,1}),V=scalar(*dev,lam_dvdy),D=scalar(*dev,{1,1,1}),mu=scalar(*dev,{1,1,1});
    auto K=scalar(*dev,{1.5e-6,1.5e-6,1.5e-6}),W=scalar(*dev,{1,1,1}),S=scalar(*dev,{0,0,0}),G=scalar(*dev,{1,1,1});
    auto source=dev->alloc<float>(3),sink=dev->alloc<float>(3),critical=dev->alloc<float>(3),onset=dev->alloc<float>(3);
    kern::transition_sources(*dev,3,one,mu,K,W,D,V,S,S,G,source,sink,critical,onset,100,1000,1,2,-.5f,0,1.5f,3,1e-14f,1e-20f);
    const auto rc=dev->download(critical);
    // Tu_L = .1: lambda clipped to +10 (F_PG = min(21, 1.5)) and -10 (F_PG = min(6, 3)).
    CHECK(near(rc[0],100+1000*std::exp(-.1*1.5),1e-5));
    CHECK(near(rc[1],100+1000*std::exp(-.1*3),1e-5));
    CHECK(near(rc[2],100+1000*std::exp(-.1*std::min(1+2*.1875,1.5)),1e-5));
  }
}

TEST_CASE("Intermittency wall-normal derivative, boundary values and diffusivity") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    // n . grad(n . U) for a general tensor and unit normals.
    const std::vector<Vec3d> normal{{0,1,0},{1,0,0},{.6,.8,0},{0,.28,.96}};
    const double G[3][3]{{.3,-40,2},{7,-.5,1.5},{-3,9,.2}};
    std::vector<double> c[9];
    for(std::size_t i=0;i<normal.size();++i) for(int r=0;r<3;++r) for(int q=0;q<3;++q) c[3*r+q].push_back(G[r][q]*(1+i));
    auto n=vector(*dev,normal);std::vector<Buffer<float>> g;for(auto& v:c) g.push_back(scalar(*dev,v));
    auto out=dev->alloc<float>(normal.size());
    kern::transition_normal_derivative(*dev,normal.size(),n.x,n.y,n.z,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],out);
    const auto dvdy=dev->download(out);
    for(std::size_t i=0;i<normal.size();++i) {
      const double v[3]{normal[i].x,normal[i].y,normal[i].z};double expected=0;
      for(int r=0;r<3;++r) for(int q=0;q<3;++q) expected+=v[r]*G[r][q]*(1+i)*v[q];
      CHECK(near(dvdy[i],expected,2e-6,1));
    }
    // Inflow faces (inlet; outlet or pressure inlet with backflow) are fixed.
    auto m=box(*dev);const int nb=m.n_boundary_faces();
    std::vector<Index> kinds{0,1,2,2,3,4};std::vector<float> flow{0,-1,1,-1,0,-1};
    auto kind=dev->upload<Index>(kinds);auto mass=dev->upload<float>(flow),inflow=filled(*dev,nb,.7f),gamma=filled(*dev,1,.2f);
    auto b=dev->alloc<float>(nb);auto transport=dev->alloc<Index>(nb);
    kern::transition_boundary(*dev,nb,m.owner,kind,inflow,mass,gamma,b,transport,0);
    CHECK(dev->download(b)==std::vector<float>{.2f,.7f,.2f,.7f,.2f,.7f});
    CHECK(dev->download(transport)==std::vector<Index>{1,0,1,0,1,0});
    auto big=box(*dev,4);std::vector<Index> wall(big.n_boundary_faces(),0);wall[0]=1;
    auto wall_mask=dev->upload<Index>(wall);auto mu=scalar(*dev,{1,2,3,4}),mut=scalar(*dev,{10,20,30,40});
    auto face=dev->alloc<float>(big.n_faces);
    kern::transition_face_diffusion(*dev,big.n_faces,big.owner,big.neighbour,big.weight,wall_mask,mu,mut,big.n_internal_faces,face,1/1.3f);
    const auto fd=dev->download(face);const auto weight=dev->download(big.weight);
    const auto owner=dev->download(big.owner),neighbour=dev->download(big.neighbour);
    for(Index f=0;f<big.n_faces;++f) {
      const auto cell=[&](int i){return (1.+i)+10*(1.+i)/1.3;};
      if(f<big.n_internal_faces) CHECK(near(fd[f],weight[f]*cell(owner[f])+(1-weight[f])*cell(neighbour[f]),1e-6));
      else if(f==big.n_internal_faces) CHECK(fd[f]==0);
      else CHECK(near(fd[f],cell(owner[f]),1e-6));
    }
    auto clamp=scalar(*dev,{-1,0,.5,1,2});kern::transition_clamp(*dev,5,clamp,1e-7f,1);
    CHECK(dev->download(clamp)==std::vector<float>{1e-7f,1e-7f,.5f,1,1});
  }
}

TEST_CASE("SST intermittency coupling scales k terms, adds P_k^lim and protects F1") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=box(*dev,8);State s(*dev,m);
    dev->upload<float>(s.k,std::vector<float>{1e-7f,1e-5f,1e-4f,1e-3f,1e-2f,1e-6f,1e-4f,1});
    dev->upload<float>(s.mu,std::vector<float>{1e-5f,1e-5f,2e-5f,1e-5f,1e-3f,1e-6f,1e-5f,1e-3f});
    dev->upload<float>(s.omega,std::vector<float>{5,40,3,20,2,80,1,9});
    dev->fill(s.gu.y,300);dev->fill(s.gv.x,-20);dev->fill(s.gk.x,.03f);dev->fill(s.go.y,-4);
    SstModel plain(*dev,m);plain.update(s.inputs());
    const auto reference=snapshot(*dev,plain.fields());
    SstOptions options;options.intermittency=SstIntermittency{1.3f,.7f};SstModel coupled(*dev,m,options);
    const std::vector<float> gammas{.01f,.15f,.21f,.5f,.9f,1,.3f,.6f};
    auto gamma=dev->upload<float>(gammas);auto in=s.inputs();in.intermittency=&gamma;
    CHECK_THROWS_WITH(coupled.update(s.inputs()),doctest::Contains("intermittency input"));
    CHECK_THROWS_WITH(plain.update(in),doctest::Contains("intermittency input"));
    coupled.update(in);const auto& f=coupled.fields();
    const auto k=dev->download(s.k),mu=dev->download(s.mu),d=dev->download(m.wall_distance);
    const auto f1=dev->download(f.f1),dk=dev->download(f.diffusivity_k),dw=dev->download(f.diffusivity_omega);
    const auto ks=dev->download(f.k_source),kk=dev->download(f.k_sink),limit=dev->download(f.production_limit);
    const auto mut=dev->download(f.mu_t),strain=dev->download(f.strain),vorticity=dev->download(f.vorticity);
    const auto ws=dev->download(f.omega_source),wk=dev->download(f.omega_sink);
    // Fields not touched by the coupling are bitwise SST's.
    CHECK(mut==reference[0]);CHECK(dev->download(f.production_k)==reference[26]);CHECK(strain==reference[19]);
    int protected_cells=0;
    for(int i=0;i<8;++i) {
      INFO("cell "<<i);
      const double f3=std::exp(-std::pow(double(d[i])*std::sqrt(double(k[i]))/mu[i]/120,8));
      const double blend=std::max(double(reference[21][i]),f3);
      CHECK(near(f1[i],blend,1e-5));
      protected_cells+=f3>reference[21][i];
      CHECK(near(dk[i],mu[i]+mut[i]*(blend/1.176+(1-blend)),1e-5));
      CHECK(near(dw[i],mu[i]+mut[i]*(blend/2+(1-blend)/1.168),1e-5));
      const double plim=production_limit(gammas[i],1,mu[i],mut[i],d[i],strain[i],vorticity[i],1.3,.7);
      CHECK(near(limit[i],plim,3e-5,1e-30));
      CHECK(near(ks[i],gammas[i]*double(reference[1][i])+plim,3e-5));
      CHECK(near(kk[i],std::clamp(double(gammas[i]),.1,1.)*reference[2][i],1e-6));
      if(f1[i]==reference[21][i]) {CHECK(ws[i]==reference[3][i]);CHECK(wk[i]==reference[4][i]);}
    }
    CHECK(protected_cells>0);CHECK(protected_cells<8);
    CHECK(std::any_of(limit.begin(),limit.end(),[](float v){return v>0;}));
    CHECK_FALSE(plain.fields().production_limit.handle());
  }
}

TEST_CASE("Wall direction points from the nearest wall point to each centroid") {
  auto mesh=make_box_hex(4,5,3,{0,0,0},{1,2,1});
  // Walls only at ymin and xmax; the others are symmetry planes.
  for(auto& p:mesh.patches) if(p.name!="ymin" && p.name!="xmax") p.kind=PatchKind::Symmetry;
  const auto g=compute_geometry(mesh);
  KUTTA_FOR_EACH_BACKEND(dev) {
    const auto n=wall_direction(*dev,mesh,g);
    const auto x=dev->download(n.x),y=dev->download(n.y),z=dev->download(n.z);
    for(Index c=0;c<mesh.n_cells;++c) {
      const auto p=g.cell_centroid[c];
      const Vec3d expected=p.y<1-p.x?Vec3d{0,1,0}:Vec3d{-1,0,0};
      INFO("cell "<<c<<" at "<<p.x<<","<<p.y);
      CHECK(near(x[c],expected.x,1e-6,1));CHECK(near(y[c],expected.y,1e-6,1));CHECK(near(z[c],expected.z,1e-6,1));
    }
  }
}

#if KUTTA_HAS_HDF5
namespace {
void set(fluent::SExpr& alist,const std::string& settings) {
  for(auto replacement:fluent::parse_sexprs(settings)) {
    const auto key=replacement.elements.front().atom_text();
    auto it=std::find_if(alist.elements.begin(),alist.elements.end(),[&](const auto& e){return !e.elements.empty() && e.elements.front().atom_text()==key;});
    if(it==alist.elements.end()) alist.elements.push_back(std::move(replacement));else *it=std::move(replacement);
  }
}
// The elbow case with case-config flags, rpvars and every inlet/outlet zone edited.
fluent::FluentCase edited(const fluent::FluentCase& original,const std::string& config,const std::string& variables="",const std::string& zones="") {
  auto sections=original.raw_sections;auto& rp=sections.at("Rampant Variables")[0].elements[1];
  for(auto& e:rp.elements) if(!e.elements.empty() && e.elements[0].atom_text()==std::string_view("case-config")) set(e.elements[1],config);
  set(rp,variables);
  if(!zones.empty()) for(auto& thread:sections.at("Thread Variables")) {
    const auto type=thread.elements[1].elements[1].atom_text();  // (id type name ...)
    if(type==std::string_view("velocity-inlet") || type==std::string_view("pressure-outlet")) set(thread.elements[2],zones);
  }
  auto c=fluent::extract_case(std::move(sections),original.solver_version,original.origin);c.mesh=original.mesh;return c;
}
}
TEST_CASE("Fluent intermittency transition settings are read, checked and solved") {
  const auto path=fluent_examples()/"mixing_elbow/mixing_elbow.cas.h5";
  if(!std::filesystem::exists(path)) {MESSAGE("SKIP: external Fluent elbow case absent");return;}
  const auto original=fluent::read_case(path);const auto mesh=read_fluent_cff(path);
  const auto allowed=setup::parse_approximations("wfgc");
  CHECK(original.viscous.transition==fluent::TransitionModel::none);
  CHECK_FALSE(setup::fluent_setup(original,mesh,allowed).settings.sst->transition);
  const std::string on="(rp-trans-one-eqn? . #t) (rp-transition-model . gamma-transport-eqn)";
  auto c=edited(original,on);
  REQUIRE(c.viscous.transition==fluent::TransitionModel::intermittency);
  const auto problem=setup::fluent_setup(c,mesh,allowed);
  REQUIRE(problem.settings.sst->transition);
  const auto& t=*problem.settings.sst->transition;
  // The fixture stores Fluent's defaults, equal to the Theory Guide's list.
  CHECK(t.model.ctu1==100);CHECK(t.model.ctu2==1000);CHECK(t.model.ctu3==1);CHECK(t.model.cpg1==1);CHECK(t.model.cpg2==-.5f);
  CHECK(t.model.cpg3==0);CHECK(t.model.cpg1_lim==1.5f);CHECK(t.model.cpg2_lim==3);CHECK(t.model.csep==1);CHECK(t.model.ck==1);
  CHECK(t.model.sigma==1);CHECK(t.model.minimum==1e-7f);CHECK(t.model.maximum==1);
  CHECK(t.convection==fv::ConvectionScheme::SecondOrderUpwind);CHECK(t.relaxation==.75);  // intermit/pseudo-relax
  CHECK(std::any_of(problem.convergence_notes.begin(),problem.convergence_notes.end(),[](const auto& s){return s.find("dudn-method")!=s.npos;}));
  for(const auto& b:problem.boundaries) if(const auto* in=std::get_if<fv::VelocityInlet>(&b.condition)) CHECK(in->intermittency==1);
  CHECK(setup::solved_residuals(c.viscous.model,true).back()=="intermit");
  {
    auto copy=c;setup::RunControls run;const auto inventory=setup::settings_inventory(copy,run);
    const auto find=[&](const std::string& p){return std::find_if(inventory.begin(),inventory.end(),[&](const auto& s){return s.path==p;});};
    REQUIRE(find("models/transition")!=inventory.end());CHECK(find("models/transition")->value==Json(std::string("gamma-transport-eqn")));
    CHECK(find("methods/intermittency-scheme")!=inventory.end());CHECK(find("controls/intermit-relaxation")!=inventory.end());
    CHECK(std::any_of(run.convergence.residuals.begin(),run.convergence.residuals.end(),[](const auto& r){return r.equation=="intermit";}));
  }
  // User-tuned constants and inflow intermittency are read, not replaced.
  const auto tuned=setup::fluent_setup(edited(original,on,"(trans-one-eqn-ctu1 . 150.) (trans-one-eqn-csep . 1.5) (intermit-sig . 1.2)","(intermit (constant . .4) (profile \"\" \"\"))"),mesh,allowed);
  CHECK(tuned.settings.sst->transition->model.ctu1==150);CHECK(tuned.settings.sst->transition->model.csep==1.5f);
  CHECK(tuned.settings.sst->transition->model.sigma==1.2f);
  for(const auto& b:tuned.boundaries) {
    if(const auto* in=std::get_if<fv::VelocityInlet>(&b.condition)) CHECK(in->intermittency==.4);
    if(const auto* out=std::get_if<fv::PressureOutlet>(&b.condition)) CHECK(out->intermittency==.4);
  }
  // Unsupported physics and inconsistent selectors are named.
  const auto rejects=[&](const fluent::FluentCase& bad,const std::string& path) {
    try {(void)setup::fluent_setup(bad,mesh,allowed);FAIL("accepted "<<path);}
    catch(const setup::RejectedSettings& e) {
      CHECK(std::any_of(e.issues().begin(),e.issues().end(),[&](const auto& i){return i.path==path;}));
    }
  };
  rejects(edited(original,on,"(trans-one-eqn-crossflow? . #t)"),"rp/trans-one-eqn-crossflow?");
  rejects(edited(original,on,"(trans-one-eqn-dudn-method . 2)"),"rp/trans-one-eqn-dudn-method");
  rejects(edited(original,on,"(trans-intermit-smoothing? . #t)"),"rp/trans-intermit-smoothing?");
  rejects(edited(original,on,"(intermit-sig . 0)"),"rp/intermit-sig");
  rejects(edited(original,on,"","(intermit (constant . 1.5) (profile \"\" \"\"))"),"zone/"+std::to_string(std::find_if(original.boundaries.begin(),original.boundaries.end(),[](const auto& z){return z.type==fluent::ZoneType::velocity_inlet;})->id)+"/intermit");
  rejects(edited(original,"(rp-trans-sst? . #t)"),"case-config/rp-trans-sst?");
  rejects(edited(original,"(rp-transition-model . gamma-algebraic)"),"case-config/rp-transition-model");
  rejects(edited(original,"(rp-trans-one-eqn? . #t)"),"case-config/rp-transition-model");  // selector says none
  // Standard initialization takes intermit/default; hybrid keeps gamma = 1.
  auto standard=edited(original,on,"(hyb-init? . #f) (intermit/default . .8)");
  const auto fields=std::get<setup::CellFields>(setup::initial_fields(standard,mesh));
  REQUIRE(fields.intermittency.size()==std::size_t(mesh.n_cells));CHECK(fields.intermittency[0]==.8);
  // A short coupled solve stays finite and agrees across backends.
  std::vector<float> oracle;
  KUTTA_FOR_EACH_BACKEND(dev) {
    fv::SteadySolver solver(*dev,mesh,problem.boundaries,problem.settings);
    setup::initialize(solver,setup::initial_fields(c,mesh));
    solver.iterate(4);const auto r=solver.monitor();
    CHECK(std::isfinite(r.intermittency.fluent_scaled));CHECK(r.intermittency.raw>0);
    const auto gamma=dev->download(solver.intermittency());
    CHECK(std::all_of(gamma.begin(),gamma.end(),[](float v){return v>=1e-7f && v<=1;}));
    CHECK(*std::min_element(gamma.begin(),gamma.end())<1);  // laminar destruction acts at the walls
    std::vector<float> result=gamma;const auto k=dev->download(solver.kinetic_energy());result.insert(result.end(),k.begin(),k.end());
    if(oracle.empty()) oracle=result;
    else for(std::size_t i=0;i<result.size();++i) CHECK(result[i]==doctest::Approx(oracle[i]).epsilon(2e-3).scale(1e-3));
  }
}
#endif
