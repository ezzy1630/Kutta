#include <cmath>
#include <limits>
#include "kernels.hpp"

#include "../backend_helpers.hpp"
#include "support.hpp"

using namespace sst_test;
namespace {
double first(Device& d, const Buffer<float>& b) { return d.download(b).front(); }
bool near(double a, double b, double tol = 3e-6) {
  return std::isfinite(a) && std::isfinite(b) && std::abs(a-b) <= tol*std::max(std::abs(b),1e-25);
}
double gamma1994(double f1) {
  return f1*(.075/.09 - .41*.41/(2*.3)) + (1-f1)*(.0828/.09 - .41*.41/(1.168*.3));
}
}  // namespace

TEST_CASE("fused SST constitutive update equals staged kernels on every backend") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    constexpr int n=97;
    std::vector<float> positive(n),signed_values(n),distance(n);
    for(int i=0;i<n;++i) {
      positive[i]=std::pow(10.f,float(i%17-8));
      signed_values[i]=float(i%13-6)*.017f;
      distance[i]=std::pow(10.f,float(i%11-5));
    }
    distance.back()=std::numeric_limits<float>::infinity();
    auto rho=filled(*dev,n,1.2f),mu=filled(*dev,n,1.8e-5f);
    auto k=dev->upload<float>(positive),omega=dev->upload<float>(positive);
    auto d=dev->upload<float>(distance),g=dev->upload<float>(signed_values);
    auto rate=filled(*dev,n,7.3f),vorticity=filled(*dev,n,11.1f);
    std::array<Buffer<float>,6> staged,fused;
    for(auto* out:{&staged,&fused}) for(auto& q:*out) q=dev->alloc<float>(n);
    for(int low_re:{0,1}) for(int use_vorticity:{0,1}) for(float cap:{2.f,1e5f}) {
      kern::sst_blending(*dev,n,rho,mu,k,omega,d,g,g,g,g,g,g,
        staged[0],staged[1],staged[2],1e-10f,1e-14f,1e-20f);
      kern::sst_viscosity(*dev,n,rho,mu,k,omega,rate,vorticity,staged[0],staged[1],
        staged[3],staged[4],staged[5],low_re,use_vorticity,cap,1e-14f,1e-20f);
      kern::sst_blending_viscosity(*dev,n,rho,mu,k,omega,d,g,g,g,g,g,g,rate,vorticity,
        fused[0],fused[1],fused[2],fused[3],fused[4],fused[5],
        1e-10f,1e-14f,1e-20f,low_re,use_vorticity,cap);
      for(std::size_t j=0;j<staged.size();++j) CHECK(dev->download(staged[j])==dev->download(fused[j]));
    }
  }
}

TEST_CASE("SST constants, F1/F2 limits, signed cross diffusion and implicit sinks") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = box(*dev); State s(*dev,m);
    SstOptions o; o.wall_production = false;
    SstModel model(*dev,m,o);
    dev->fill(m.wall_distance,1e-5f);
    model.update(s.inputs()); const auto& f = model.fields();
    CHECK(first(*dev,f.f1) == 1); CHECK(first(*dev,f.f2) == 1);
    CHECK(near(first(*dev,f.mu_t),.5));
    CHECK(near(first(*dev,f.diffusivity_k),.001+.5/1.176));
    CHECK(near(first(*dev,f.diffusivity_omega),.001+.5/2));
    CHECK(near(first(*dev,f.k_sink),.18));
    CHECK(near(first(*dev,f.omega_sink),.30));
    dev->fill(m.wall_distance,1e6f);
    model.update(s.inputs());
    CHECK(first(*dev,f.f1) < 1e-18); CHECK(first(*dev,f.f2) < 1e-8);
    CHECK(near(first(*dev,f.diffusivity_k),.501));
    CHECK(near(first(*dev,f.diffusivity_omega),.001+.5/1.168));
    CHECK(near(first(*dev,f.omega_sink),.3312));
    // A negative gradient dot product must not be replaced by CDkw's F1 floor.
    dev->fill(m.wall_distance,2); dev->fill(s.k,.1f);
    dev->fill(s.gk.x,1); dev->fill(s.go.x,-1);
    model.update(s.inputs());
    const auto cd = -1/1.168;
    const auto a = std::sqrt(.1)/(.09*2*2);
    const auto expected_f1 = std::tanh(std::pow(a,4));
    CHECK(near(first(*dev,f.f1),expected_f1));
    CHECK(near(first(*dev,f.cdkw),cd));
    CHECK(first(*dev,f.k_source) == 0);
    CHECK(near(first(*dev,f.omega_source),(.075*expected_f1+.0828*(1-expected_f1))*4));
    CHECK(near(first(*dev,f.omega_source)-2*first(*dev,f.omega_sink),
      -(.075*expected_f1+.0828*(1-expected_f1))*4 + (1-expected_f1)*cd));
    dev->fill(s.go.x,1); model.update(s.inputs());
    CHECK(first(*dev,f.cdkw) > 0);
    CHECK(near(first(*dev,f.omega_source),(1-first(*dev,f.f1))/1.168+(.075*first(*dev,f.f1)+.0828*(1-first(*dev,f.f1)))*4));
  }
}

TEST_CASE("SST production switches, viscosity limiter and gamma choices") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = box(*dev); State s(*dev,m);
    dev->fill(s.gu.x,5000); dev->fill(s.gv.y,-5000);
    SstOptions o; o.wall_production = false;
    for (bool limiter : {false,true}) for (bool limit_w : {false,true}) {
      o.production_limiter = limiter; o.limit_omega_production = limit_w;
      SstModel model(*dev,m,o); model.update(s.inputs()); const auto& f = model.fields();
      const auto mt = first(*dev,f.mu_t), f1 = first(*dev,f.f1);
      const auto raw = mt*1e8;
      CHECK(near(first(*dev,f.strain),10000)); CHECK(first(*dev,f.vorticity) == 0);
      CHECK(near(first(*dev,f.k_source),limiter ? 1.8 : raw));
      CHECK(near(first(*dev,f.omega_source),gamma1994(f1)*(limiter && limit_w ? 1.8 : raw)/mt+(.075*f1+.0828*(1-f1))*4));
      o.gamma_set = SstGammaSet::Menter2003;
      SstModel alt(*dev,m,o); alt.update(s.inputs());
      CHECK(near(first(*dev,alt.fields().omega_source),
        (f1*5/9+(1-f1)*.44)*(limiter && limit_w ? 1.8 : raw)/mt+(.075*f1+.0828*(1-f1))*4));
      o.gamma_set = SstGammaSet::Computed1994;
    }
    o.production_limiter = false; o.kato_launder = true;
    SstModel kato(*dev,m,o); kato.update(s.inputs()); CHECK(first(*dev,kato.fields().k_source) == 0);
    o.kato_launder = false; o.viscosity_limiter = SstViscosityLimiter::Vorticity;
    SstModel vort_limiter(*dev,m,o); vort_limiter.update(s.inputs());
    CHECK(near(first(*dev,vort_limiter.fields().mu_t),.5));
    // Pure rotation has S=0, Omega=20. Production is selectable independently
    // of the viscosity limiter, so the two switches cannot be conflated.
    dev->fill(s.gu.x,0); dev->fill(s.gv.y,0); dev->fill(s.gu.y,-10); dev->fill(s.gv.x,10);
    o.vorticity_production = true;
    SstModel vort_prod(*dev,m,o); vort_prod.update(s.inputs());
    CHECK(first(*dev,vort_prod.fields().strain) == 0);
    CHECK(near(first(*dev,vort_prod.fields().vorticity),20));
    CHECK(near(first(*dev,vort_prod.fields().k_source),400*first(*dev,vort_prod.fields().mu_t)));
    o.vorticity_production = false; o.max_viscosity_ratio = 2;
    dev->fill(s.gu.y,0); dev->fill(s.gv.x,0);
    SstModel capped(*dev,m,o); capped.update(s.inputs());
    CHECK(near(first(*dev,capped.fields().mu_t),.002));
  }
}

TEST_CASE("SST optional low-Re coefficients and explicit floors") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = box(*dev); State s(*dev,m);
    dev->fill(s.k,1e-8f); dev->fill(s.mu,1); dev->fill(s.omega,1);
    dev->fill(s.gu.y,1); // Re_T=1e-8; alpha*_0=.024, beta*_0=.024
    SstOptions o; o.low_re = true; o.production_limiter = false; o.wall_production = false;
    SstModel model(*dev,m,o); model.update(s.inputs()); const auto& f = model.fields();
    CHECK(near(first(*dev,f.mu_t),.024e-8));
    CHECK(near(first(*dev,f.k_sink),.024));
    CHECK(near(first(*dev,f.omega_source),gamma1994(first(*dev,f.f1))/9+(.075*first(*dev,f.f1)+.0828*(1-first(*dev,f.f1)))));
    // Test floors through the actual update path, without repairing input arrays.
    dev->fill(s.k,-1); dev->fill(s.omega,0);
    o.low_re = false; o.k_floor = 1e-6f; o.omega_floor = 1e-4f;
    SstModel floors(*dev,m,o); floors.update(s.inputs());
    CHECK(first(*dev,s.k) == -1); CHECK(first(*dev,s.omega) == 0);
    for (const auto& array : snapshot(*dev,floors.fields()))
      for (float x : array) CHECK(std::isfinite(x));
    CHECK(near(first(*dev,floors.fields().k_sink),9e-6));
  }
}

TEST_CASE("SST wall sublayer, log layer, moving walls and zero relative velocity") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = box(*dev); State s(*dev,m);
    dev->fill(s.k,0); dev->fill(s.velocity.x,1e-6f);
    SstModel model(*dev,m); model.update(s.inputs()); const auto& f = model.fields();
    // Face y-min has tangential U along x; x-normal faces have no shear.
    const auto tau = dev->download(f.wall_tau);
    CHECK(tau[0] == 0); CHECK(tau[1] == 0);
    CHECK(near(tau[2],.001*1e-6/.5));
    const auto omega = first(*dev,f.fixed_omega);
    CHECK(near(omega,2*.001/(.072*.5*.5),1e-4));
    CHECK(dev->download(f.fixed_omega_mask) == std::vector<Index>{1});
    CHECK(dev->download(f.selected_wall_face) == std::vector<Index>{0});
    dev->fill(s.mu,1e-5f); dev->fill(s.k,1); dev->fill(s.velocity.x,10);
    model.update(s.inputs());
    const auto us = std::sqrt(1e-5*10/.5+.3), yp = us*.5/1e-5;
    const auto log_up = std::log(9.793*yp)/.4187;
    const auto wall_tau = dev->download(f.wall_tau)[2];
    CHECK(near(wall_tau,us*10/log_up));
    CHECK(near(dev->download(f.wall_y_plus)[2],yp));
    CHECK(near(dev->download(f.wall_y_plus_utau)[2],std::sqrt(wall_tau)*.5/1e-5));
    CHECK(near(dev->download(f.wall_shear.x)[2],-wall_tau));
    const double wl = 2e-5/(.072*.25), wt = us/(.3*.4187*.5);
    const double wall_w = std::pow(std::pow(wl,1.3)+std::pow(wt,1.3),1/1.3);
    // Nearest/tied face 0 is normal to U, so inspect the selected y-face by
    // making it closer than the other five faces.
    auto dy = dev->download(m.wall_normal_distance); dy[2] = .49f;
    dev->upload<float>(m.wall_normal_distance,dy); model.update(s.inputs());
    CHECK(dev->download(f.selected_wall_face) == std::vector<Index>{2});
    CHECK(wall_w > wt); // finite viscous contribution to the log branch
    const auto p_lam = 1.0/2 * std::pow(dev->download(f.wall_tau)[2]/1e-5,2);
    const auto p_log = std::pow(dev->download(f.wall_tau)[2],2)/1e-5 /
                       (.4187*dev->download(f.wall_y_plus)[2]);
    CHECK(near(first(*dev,f.wall_production),1/(1/p_lam+1/p_log)));
    CHECK(near(first(*dev,f.k_source),std::min(first(*dev,f.wall_production),1.8)));
    CHECK(near(first(*dev,f.fixed_omega),std::pow(
      std::pow(2e-5/(.072*.49*.49),1.3)+
      std::pow(std::sqrt(1e-5*10/.49+.3)/(.3*.4187*.49),1.3),1/1.3)));
    // Translating the wall with the fluid removes all relative shear.
    dev->fill(s.wall_velocity.x,10); model.update(s.inputs());
    for (float v : dev->download(f.wall_tau)) CHECK(v == 0);
    CHECK(first(*dev,f.wall_production) == 0);
  }
}

TEST_CASE("SST wall law and omega correlation blend continuously") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = box(*dev); State s(*dev,m);
    // k=0 gives y+ = sqrt(U*d/nu); scan through the log clamp and layer blend.
    dev->fill(s.k,0);
    auto distances = dev->download(m.wall_normal_distance); distances[2] = .4f;
    dev->upload<float>(m.wall_normal_distance,distances);
    SstModel model(*dev,m);
    for (double yp : {.2,1.,4.,8.,11.,30.,300.}) {
      const auto speed = yp*yp*.001/.5;
      dev->fill(s.velocity.x,float(speed*(1-1e-5))); model.update(s.inputs());
      auto before = dev->download(model.fields().wall_tau)[2];
      auto wb = first(*dev,model.fields().fixed_omega);
      dev->fill(s.velocity.x,float(speed*(1+1e-5))); model.update(s.inputs());
      auto after = dev->download(model.fields().wall_tau)[2];
      auto wa = first(*dev,model.fields().fixed_omega);
      CHECK(after >= before); CHECK(after/before < 1.0001);
      CHECK(std::abs(wa/wb-1) < 1e-4);
    }
  }
}

TEST_CASE("SST F1 and F2 blend continuously across their transition") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = box(*dev); State s(*dev,m);
    SstModel model(*dev,m);
    for (float distance : {.001f,.01f,.1f,1.f,2.f,4.f,8.f,16.f,100.f}) {
      dev->fill(m.wall_distance,distance*(1-1e-5f)); model.update(s.inputs());
      const auto a1 = first(*dev,model.fields().f1), a2 = first(*dev,model.fields().f2);
      dev->fill(m.wall_distance,distance*(1+1e-5f)); model.update(s.inputs());
      const auto b1 = first(*dev,model.fields().f1), b2 = first(*dev,model.fields().f2);
      CHECK(a1 >= b1); CHECK(a2 >= b2);
      CHECK(a1-b1 < 1e-4); CHECK(a2-b2 < 1e-4);
    }
  }
}

TEST_CASE("SST interpolation, no-wall constraints and boundary masking") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = box(*dev,2); State s(*dev,m);
    for (auto& p : m.patches) p.kind = PatchKind::Symmetry;
    auto vis = std::vector<float>{.001f,.002f}; dev->upload<float>(s.mu,vis);
    auto ke = std::vector<float>{1,2}; dev->upload<float>(s.k,ke);
    SstModel model(*dev,m); model.update(s.inputs()); const auto& f = model.fields();
    CHECK(dev->download(f.fixed_omega_mask) == std::vector<Index>{0,0});
    CHECK(dev->download(f.selected_wall_face) == std::vector<Index>{-1,-1});
    for (float x : dev->download(f.fixed_omega)) CHECK(x == 0);
    for (float x : dev->download(f.wall_tau)) CHECK(x == 0);
    const auto gc = dev->download(f.diffusivity_k), gf = dev->download(f.face_diffusivity_k);
    CHECK(near(gf[0],.5*(gc[0]+gc[1])));
    for (Index b = m.n_internal_faces; b < m.n_faces; ++b)
      CHECK(near(gf[b],gc[dev->download(m.owner)[b]]));
  }
}

TEST_CASE("SST every kernel matches CPU and repeats deterministically") {
  auto cpu = Device::create(Backend::Cpu);
  for (bool low_re : {false,true}) {
    auto run = [&](Device& d) {
      auto m = box(d,7); State s(d,m);
      s.velocity = constant_vector(d,7,{.7,-.4,1.2});
      s.gu = constant_vector(d,7,{.2,-.7,.8});
      s.gv = constant_vector(d,7,{1.2,-.3,.9});
      s.gw = constant_vector(d,7,{-.6,.1,.1});
      s.gk = constant_vector(d,7,{.003,-.002,.001});
      s.go = constant_vector(d,7,{-.6,.3,.2});
      s.wall_velocity = constant_vector(d,m.n_boundary_faces(),{.1,.02,.03});
      std::vector<float> k{1e-7f,.2f,.6f,1.2f,2,10,100}; d.upload<float>(s.k,k);
      SstOptions o; o.low_re = low_re;
      SstModel model(d,m,o); model.update(s.inputs());
      const auto result = snapshot(d,model.fields());
      const auto mask = d.download(model.fields().fixed_omega_mask);
      CHECK(mask == std::vector<Index>(7,1));
      model.update(s.inputs()); CHECK(result == snapshot(d,model.fields()));
      return result;
    };
    const auto reference = run(*cpu);
    KUTTA_FOR_EACH_BACKEND(dev) {
      const auto result = run(*dev);
      REQUIRE(result.size() == reference.size());
      for (std::size_t a = 0; a < result.size(); ++a) {
        INFO("field " << a);
        REQUIRE(result[a].size() == reference[a].size());
        for (std::size_t i = 0; i < result[a].size(); ++i) {
          INFO("entry " << i); CHECK(near(result[a][i],reference[a][i],1e-5));
        }
      }
    }
  }
}

TEST_CASE("SST unsupported and invalid options fail explicitly") {
  auto d = Device::create(Backend::Cpu); auto m = box(*d);
  SstOptions o; o.wall_omega = SstWallOmega::Tabulated;
  CHECK_THROWS_WITH_AS(SstModel(*d,m,o),doctest::Contains("tabulated"),Error);
  o = {}; o.curvature_correction = true;
  CHECK_THROWS_WITH_AS(SstModel(*d,m,o),doctest::Contains("curvature"),Error);
  o = {}; o.compressibility = true; CHECK_THROWS_AS(SstModel(*d,m,o),Error);
  o = {}; o.kato_launder = true; o.vorticity_production = true; CHECK_THROWS_AS(SstModel(*d,m,o),Error);
  o = {}; o.cdkw_floor = 0; CHECK_THROWS_AS(SstModel(*d,m,o),Error);
  o = {}; o.cross_diffusion_clip = 0; CHECK_THROWS_AS(SstModel(*d,m,o),Error);
  o = {}; o.production_clip = std::numeric_limits<float>::quiet_NaN(); CHECK_THROWS_AS(SstModel(*d,m,o),Error);
  State s(*d,m); SstModel model(*d,m);
  s.k = filled(*d,2,1); CHECK_THROWS_AS(model.update(s.inputs()),Error);
}

TEST_CASE("SST floor cross-diffusion overflow probe remains finite") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = box(*dev); State s(*dev,m);
    // No wall constraint: represents a layer-2 transport cell.
    for (auto& patch : m.patches) patch.kind = PatchKind::Symmetry;
    SstOptions o;
    CHECK(o.k_floor == 1e-14f); CHECK(o.omega_floor == 1e-20f);
    for (bool clip : {false,true}) for (bool low_re : {false,true}) {
      o.limit_cross_diffusion = clip; o.low_re = low_re;
      SstModel model(*dev,m,o);
      for (float ke : {0.f,o.k_floor,1.f}) for (float w : {0.f,o.omega_floor,2.f})
        for (float sign : {-1.f,1.f}) {
          INFO("clip=" << clip << " low_re=" << low_re << " k=" << ke << " omega=" << w << " sign=" << sign);
          dev->fill(s.k,ke); dev->fill(s.omega,w);
          dev->fill(s.gk.x,1e5f); dev->fill(s.go.x,sign*1e14f);
          model.update(s.inputs());
          for (const auto& array : snapshot(*dev,model.fields()))
            for (float value : array) CHECK(std::isfinite(value));
          if (clip) {
            const double we = std::max(w,o.omega_floor);
            const double beta = .075*first(*dev,model.fields().f1)+.0828*(1-first(*dev,model.fields().f1));
            const double bound = o.cross_diffusion_clip*beta*we*we;
            const double cross = std::min((1-first(*dev,model.fields().f1))*std::abs(first(*dev,model.fields().cdkw)),bound);
            if (sign > 0) CHECK(near(first(*dev,model.fields().omega_source)/bound,(cross+beta*we*we)/bound));
            else CHECK(near(first(*dev,model.fields().omega_sink),2*beta*we+cross/we));
          }
        }
    }
    // Exercise nonzero production at floors with every limiter setting too.
    dev->fill(s.gu.y,1e4f); dev->fill(s.velocity.x,10);
    for (bool walls : {false,true}) {
      for (auto& patch : m.patches) patch.kind = walls ? PatchKind::Wall : PatchKind::Symmetry;
      for (bool low_re : {false,true}) for (bool limiter : {false,true})
        for (bool limit_w : {false,true}) {
          o = {}; o.low_re = low_re; o.production_limiter = limiter; o.limit_omega_production = limit_w;
          SstModel producing(*dev,m,o);
          dev->fill(s.k,0); dev->fill(s.omega,0); producing.update(s.inputs());
          const auto zeros = snapshot(*dev,producing.fields());
          CHECK(first(*dev,s.k) == 0); CHECK(first(*dev,s.omega) == 0);
          dev->fill(s.k,o.k_floor); dev->fill(s.omega,o.omega_floor); producing.update(s.inputs());
          CHECK(snapshot(*dev,producing.fields()) == zeros);
          for (float ke : {0.f,o.k_floor}) for (float w : {0.f,o.omega_floor,2.f}) {
            dev->fill(s.k,ke); dev->fill(s.omega,w);
            producing.update(s.inputs());
            for (const auto& array : snapshot(*dev,producing.fields()))
              for (float value : array) CHECK(std::isfinite(value));
          }
        }
    }
    // Wall production and fixed omega must also tolerate zero/floored inputs.
    for (auto& patch : m.patches) patch.kind = PatchKind::Wall;
    for (float ke : {0.f,o.k_floor}) for (float w : {0.f,o.omega_floor}) {
      dev->fill(s.k,ke); dev->fill(s.omega,w); dev->fill(s.velocity.x,10);
      SstModel wall(*dev,m); wall.update(s.inputs());
      for (const auto& array : snapshot(*dev,wall.fields()))
        for (float value : array) CHECK(std::isfinite(value));
    }
  }
}

TEST_CASE("SST signed cross-diffusion clip scales with omega destruction") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = box(*dev); State s(*dev,m);
    for (auto& p : m.patches) p.kind = PatchKind::Symmetry;
    dev->fill(m.wall_distance,2); dev->fill(s.k,.1f);
    dev->fill(s.gk.x,1); dev->fill(s.go.x,1);
    SstOptions o; o.cross_diffusion_clip = .25f;
    for (bool clip : {false,true}) for (float sign : {-1.f,1.f}) {
      o.limit_cross_diffusion = clip; dev->fill(s.go.x,sign);
      SstModel model(*dev,m,o); model.update(s.inputs()); const auto& f = model.fields();
      const double f1 = first(*dev,f.f1), beta = .075*f1+.0828*(1-f1);
      const double raw = (1-f1)*sign/1.168;
      const double cross = clip ? std::clamp(raw,-o.cross_diffusion_clip*beta*4,o.cross_diffusion_clip*beta*4) : raw;
      CHECK(near(first(*dev,f.cdkw),sign/1.168)); // F1 diagnostic is not physically clipped
      CHECK(near(first(*dev,f.omega_source),std::max(cross,0.)+beta*4));
      CHECK(near(first(*dev,f.omega_sink),beta*4+std::max(-cross,0.)/2));
    }
  }
}

TEST_CASE("SST F1 cross branch and transport match an independent double oracle") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = box(*dev); State s(*dev,m);
    for (auto& p : m.patches) p.kind = PatchKind::Symmetry;
    // Menter 1994 multiplier sigma_w2=.856; Fluent stores its rounded divisor
    // 1.168. Use the exact Fluent reciprocal for tight transport comparisons,
    // and separately bound the small rounding difference against Menter.
    const double sigma_w2 = 1/1.168, r = 1.3, mu = .002, ke = .1, w = 2, y = .5, rate = 3;
    dev->fill(s.rho,float(r)); dev->fill(s.mu,float(mu)); dev->fill(s.k,float(ke));
    dev->fill(s.omega,float(w)); dev->fill(m.wall_distance,float(y)); dev->fill(s.gu.y,float(rate));
    dev->fill(s.gk.x,1);
    SstOptions o; o.wall_production = false;
    for (double dot : {2.,4.,8.,-4.}) {
      dev->fill(s.go.x,float(dot));
      SstModel model(*dev,m,o); model.update(s.inputs()); const auto& f = model.fields();
      const auto cd = 2*r*sigma_w2*dot/w;
      const auto a = std::max(std::sqrt(ke)/(.09*w*y),500*mu/(r*y*y*w));
      const auto c = 4*r*sigma_w2*ke/(std::max(cd,double(o.cdkw_floor))*y*y);
      if (dot > 0) REQUIRE(c < a); // exercise the previously untested F1 branch
      const auto f1 = std::tanh(std::pow(std::min(a,c),4));
      const auto f2 = std::tanh(std::pow(std::max(2*std::sqrt(ke)/(.09*w*y),500*mu/(r*y*y*w)),2));
      const auto mt = r*.31*ke/std::max(.31*w,rate*f2);
      const auto pk = std::min(mt*rate*rate,10*r*.09*ke*w);
      const auto beta = .075*f1+.0828*(1-f1), cross = (1-f1)*cd;
      CHECK(near(first(*dev,f.f1),f1)); CHECK(near(first(*dev,f.f2),f2));
      CHECK(near(first(*dev,f.cdkw),cd)); CHECK(near(first(*dev,f.cdkw),2*r*.856*dot/w,3e-4));
      CHECK(near(first(*dev,f.mu_t),mt));
      CHECK(near(first(*dev,f.k_source),pk)); CHECK(near(first(*dev,f.k_sink),r*.09*w));
      CHECK(near(first(*dev,f.omega_source),gamma1994(f1)*r*pk/mt+std::max(cross,0.)+r*beta*w*w));
      CHECK(near(first(*dev,f.omega_sink),2*r*beta*w+std::max(-cross,0.)/w));
    }
  }
}

TEST_CASE("SST implicit wall shear, reconstruction omega and structural zero k flux") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = box(*dev); State s(*dev,m);
    m.patches.front().kind = PatchKind::Symmetry;
    SstOptions o;
    SstModel model(*dev,m,o); const auto& f = model.fields();
    for (float ke : {0.f,1.f}) for (float speed : {0.f,1e-12f,1e-6f,10.f}) {
      dev->fill(s.k,ke); dev->fill(s.mu,1e-5f);
      // Relative normal motion must not enter the tangential implicit law.
      s.velocity = constant_vector(*dev,1,{speed+2,3,0});
      s.wall_velocity = constant_vector(*dev,m.n_boundary_faces(),{2,2,0});
      model.update(s.inputs());
      auto coeff = dev->download(f.wall_shear_coefficient), tau = dev->download(f.wall_tau);
      auto face_w = dev->download(f.wall_face_omega), face_k = dev->download(f.face_diffusivity_k);
      auto fw_diff = dev->download(f.face_diffusivity_omega), cell_k = dev->download(f.diffusivity_k);
      auto sx = dev->download(f.wall_shear.x), sy = dev->download(f.wall_shear.y), sz = dev->download(f.wall_shear.z);
      auto ax = dev->download(m.area.x), ay = dev->download(m.area.y), az = dev->download(m.area.z);
      const auto actual_u = dev->download(s.velocity.x)[0]-2.;
      for (const auto& patch : m.patches) for (Index face = patch.start; face < patch.start+patch.size; ++face) {
        const auto b = face-m.n_internal_faces;
        if (patch.kind != PatchKind::Wall) {
          CHECK(coeff[b] == 0); CHECK(face_w[b] == 0);
          CHECK(near(face_k[face],cell_k[0]));
          continue;
        }
        CHECK(face_k[face] == 0); CHECK(near(fw_diff[face],1e-5));
        const Vec3d normal{ax[face],ay[face],az[face]}, relative{actual_u,1,0};
        const auto tangent = relative-normal*(dot(relative,normal)/dot(normal,normal));
        const auto u = norm(tangent), d = .5, nu = 1e-5;
        const auto us = std::sqrt(nu*u/d+.3*std::max(double(ke),double(o.k_floor)));
        const auto yp = us*d/nu, log_up = std::log(9.793*std::max(yp,.2))/.4187;
        const auto expected = std::pow(std::pow(nu/d,4)+std::pow(us/log_up,4),.25);
        CHECK(near(coeff[b],expected)); CHECK(std::isfinite(coeff[b])); CHECK(coeff[b] > 0);
        CHECK(near(tau[b],expected*u));
        const Vec3d implicit_force = tangent*(-coeff[b]*norm(normal));
        const Vec3d explicit_force{sx[b]*norm(normal),sy[b]*norm(normal),sz[b]*norm(normal)};
        CHECK(norm(implicit_force-explicit_force) <= 1e-6*std::max(norm(implicit_force),1e-25));
        const auto expected_face_w = std::pow(std::pow(6*nu/(.072*d*d),1.3)+
          std::pow(2*us/(.3*.4187*d),1.3),1/1.3);
        CHECK(near(face_w[b],expected_face_w));
      }
    }
    // Exact zero tangential speed in both the viscous and k-driven log regimes.
    for (float ke : {0.f,1.f}) {
      dev->fill(s.k,ke); s.velocity = constant_vector(*dev,1,{2,2,0});
      model.update(s.inputs()); auto zero = dev->download(f.wall_shear_coefficient);
      dev->fill(s.velocity.z,1e-6f); model.update(s.inputs());
      auto small = dev->download(f.wall_shear_coefficient);
      CHECK(near(zero[2],small[2],1e-4)); CHECK(zero[2] > 0);
      if (ke == 0) CHECK(near(zero[2],1e-5/.5));
    }
  }
}
