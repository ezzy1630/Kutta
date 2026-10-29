#include <cmath>
#include <limits>

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
    CHECK(near(first(*dev,f.omega_sink),.15));
    dev->fill(m.wall_distance,1e6f);
    model.update(s.inputs());
    CHECK(first(*dev,f.f1) < 1e-18); CHECK(first(*dev,f.f2) < 1e-8);
    CHECK(near(first(*dev,f.diffusivity_k),.501));
    CHECK(near(first(*dev,f.diffusivity_omega),.001+.5/1.168));
    CHECK(near(first(*dev,f.omega_sink),.1656));
    // A negative gradient dot product must not be replaced by CDkw's F1 floor.
    dev->fill(m.wall_distance,2); dev->fill(s.k,.1f);
    dev->fill(s.gk.x,1); dev->fill(s.go.x,-1);
    model.update(s.inputs());
    const auto cd = -1/1.168;
    const auto a = std::sqrt(.1)/(.09*2*2);
    const auto expected_f1 = std::tanh(std::pow(a,4));
    CHECK(near(first(*dev,f.f1),expected_f1));
    CHECK(near(first(*dev,f.cdkw),cd));
    CHECK(first(*dev,f.k_source) == 0); CHECK(first(*dev,f.omega_source) == 0);
    CHECK(near(-2*first(*dev,f.omega_sink),
      -(.075*expected_f1+.0828*(1-expected_f1))*4 + (1-expected_f1)*cd));
    dev->fill(s.go.x,1); model.update(s.inputs());
    CHECK(first(*dev,f.cdkw) > 0);
    CHECK(near(first(*dev,f.omega_source),(1-first(*dev,f.f1))/1.168));
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
      CHECK(near(first(*dev,f.omega_source),gamma1994(f1)*(limiter && limit_w ? 1.8 : raw)/mt));
      o.gamma_set = SstGammaSet::Menter2003;
      SstModel alt(*dev,m,o); alt.update(s.inputs());
      CHECK(near(first(*dev,alt.fields().omega_source),
        (f1*5/9+(1-f1)*.44)*(limiter && limit_w ? 1.8 : raw)/mt));
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
    CHECK(near(first(*dev,f.omega_source),gamma1994(first(*dev,f.f1))/9));
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
  o = {}; o.production_clip = std::numeric_limits<float>::quiet_NaN(); CHECK_THROWS_AS(SstModel(*d,m,o),Error);
  State s(*d,m); SstModel model(*d,m);
  s.k = filled(*d,2,1); CHECK_THROWS_AS(model.update(s.inputs()),Error);
}
