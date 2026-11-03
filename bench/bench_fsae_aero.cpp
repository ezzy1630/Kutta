// Runs one FSAE aero validation solve (validation/cases/fsae_aero.hpp) and
// prints a RESULT line; tools/fsae_aero/run.py drives the sweeps.
#include <iostream>
#include <map>
#include <set>
#include <sstream>

#include "../validation/cases/fsae_aero.hpp"

int main(int argc, char** argv) {
  using namespace kutta;
  using namespace validation::fsae_aero;
  try {
    std::map<std::string, std::string> args;
    for (int i = 1; i < argc; ++i) {
      std::string key = argv[i];
      KUTTA_REQUIRE(key.starts_with("--"), "expected --option, got " + key);
      const bool is_flag = key == "--moving-ground" || key == "--coupled" || key == "--downforce" || key == "--cells" || key == "--lift-z" || key == "--hybrid" || key == "--transition" || key == "--laminar" || key == "--transition-harmonic" || key == "--kato-launder";
      KUTTA_REQUIRE(is_flag || i + 1 < argc, "missing value for " + key);
      args[key] = is_flag ? "1" : argv[++i];
    }
    std::set<std::string> used;
    const auto text = [&](const std::string& k, std::string d = {}) {
      if(!args.contains(k)) return d;
      used.insert(k);return args.at(k);
    };
    const auto flag = [&](const std::string& k) {used.insert(k);return args.contains(k);};
    const auto number = [&](const std::string& k, double d) {
      if (!args.contains(k)) return d;
      used.insert(k);
      std::istringstream in(args.at(k));
      in.imbue(std::locale::classic());
      double v = 0;
      in >> v;
      KUTTA_REQUIRE(in && in.eof(), "invalid number for " + k);
      return v;
    };
    Options o;
    o.mesh = text("--mesh");
    const auto limiter = text("--limiter", "standard");
    KUTTA_REQUIRE(limiter == "standard" || limiter == "differentiable", "--limiter is standard or differentiable");
    o.limiter = limiter == "standard" ? fv::Limiter::Standard : fv::Limiter::Differentiable;
    o.limiter_coefficient = number("--limiter-coefficient", o.limiter_coefficient);
    const auto pressure_scheme = text("--pressure-scheme", "second-order");
    KUTTA_REQUIRE(pressure_scheme == "standard" || pressure_scheme == "second-order", "--pressure-scheme is standard or second-order");
    o.pressure_interpolation = pressure_scheme == "standard" ? fv::PressureInterpolation::Standard : fv::PressureInterpolation::SecondOrder;
    const auto boundary_gradient = text("--pressure-boundary-gradient", o.pressure_boundary_gradient == fv::PressureBoundaryGradient::Automatic ? "automatic"
      : o.pressure_boundary_gradient == fv::PressureBoundaryGradient::NormalFoot ? "normal-foot" : "extrapolated");
    KUTTA_REQUIRE(boundary_gradient == "automatic" || boundary_gradient == "extrapolated" || boundary_gradient == "normal-foot",
                  "--pressure-boundary-gradient is automatic, extrapolated or normal-foot");
    o.pressure_boundary_gradient = boundary_gradient == "normal-foot" ? fv::PressureBoundaryGradient::NormalFoot
      : boundary_gradient == "extrapolated" ? fv::PressureBoundaryGradient::Extrapolated : fv::PressureBoundaryGradient::Automatic;
    const auto laminar_reconstruction = text("--laminar-reconstruction", o.laminar_reconstruction == fv::LaminarReconstruction::Automatic ? "automatic"
      : o.laminar_reconstruction == fv::LaminarReconstruction::Cubic ? "cubic" : "projected");
    KUTTA_REQUIRE(laminar_reconstruction == "automatic" || laminar_reconstruction == "projected" || laminar_reconstruction == "cubic",
                  "--laminar-reconstruction is automatic, projected or cubic");
    o.laminar_reconstruction = laminar_reconstruction == "cubic" ? fv::LaminarReconstruction::Cubic
      : laminar_reconstruction == "projected" ? fv::LaminarReconstruction::Projected : fv::LaminarReconstruction::Automatic;
    const auto transport = text("--sst-transport", "fluent");
    KUTTA_REQUIRE(transport == "fluent" || transport == "cfl3d", "--sst-transport is fluent or cfl3d");
    if (transport == "cfl3d") {
      // Original SST-V transport (CFL3D Appendix H and NASA TMR). Wall
      // treatment stays Fluent's; this is not a complete CFL3D wall model.
      o.sst_model.viscosity_limiter = SstViscosityLimiter::Vorticity;
      o.sst_model.vorticity_production = true;
      o.sst_model.production_clip = 20;
      o.sst_model.limit_omega_production = false;
      o.sst_model.cdkw_floor = 1e-20f;
      o.sst_model.limit_cross_diffusion = false;
      o.sst_model.max_viscosity_ratio = 1e10f;
    }
    const auto viscosity = text("--sst-viscosity", o.sst_model.viscosity_limiter == SstViscosityLimiter::Strain ? "strain" : "vorticity");
    const auto production = text("--sst-production", o.sst_model.vorticity_production ? "vorticity" : "strain");
    KUTTA_REQUIRE(viscosity == "strain" || viscosity == "vorticity", "--sst-viscosity is strain or vorticity");
    KUTTA_REQUIRE(production == "strain" || production == "vorticity", "--sst-production is strain or vorticity");
    o.sst_model.viscosity_limiter = viscosity == "strain" ? SstViscosityLimiter::Strain : SstViscosityLimiter::Vorticity;
    o.sst_model.vorticity_production = production == "vorticity";
    const auto turbulence_convection = text("--turbulence-convection", "second-order");
    KUTTA_REQUIRE(turbulence_convection == "first-order" || turbulence_convection == "second-order", "--turbulence-convection is first-order or second-order");
    o.turbulence_convection = turbulence_convection == "first-order" ? fv::ConvectionScheme::FirstOrderUpwind : fv::ConvectionScheme::SecondOrderUpwind;
    o.poor_cell_quality = number("--poor-cell-quality", o.poor_cell_quality);
    o.output = text("--out");
    o.reynolds = number("--re", 0);
    o.length = number("--length", 1);
    o.reference_area = number("--area", 1);
    o.moving_ground = flag("--moving-ground");
    const auto ground = text("--ground", "wall");
    KUTTA_REQUIRE(ground == "wall" || ground == "prescribed-zero", "--ground is wall or prescribed-zero");
    o.ground_condition = ground == "wall" ? GroundCondition::Wall : GroundCondition::PrescribedZeroTurbulence;
    o.coupled = flag("--coupled");
    o.hybrid = flag("--hybrid");
    o.transition = flag("--transition");
    o.intermittency_relaxation = number("--urf-gamma", o.intermittency_relaxation);
    o.transition_harmonic = flag("--transition-harmonic");
    o.kato_launder = flag("--kato-launder");
    o.omega_production_molecular_ratio = number("--omega-production-molecular-ratio", 0);
    o.laminar = flag("--laminar");
    KUTTA_REQUIRE(o.laminar || o.laminar_reconstruction != fv::LaminarReconstruction::Cubic,
                  "--laminar-reconstruction cubic requires --laminar");
    o.write_cells = flag("--cells");
    if (flag("--downforce")) o.lift_direction = {0, -1, 0};
    if (flag("--lift-z")) o.lift_direction = {0, 0, 1};
    o.intensity = number("--tu", o.intensity);
    o.viscosity_ratio = number("--ratio", o.viscosity_ratio);
    o.max_iterations = int(number("--iterations", o.max_iterations));
    o.check_every = int(number("--check", o.check_every));
    o.window = int(number("--window", o.window));
    o.max_seconds = number("--max-seconds", o.max_seconds);
    o.settle = number("--settle", o.settle);
    o.pseudo_time_step = number("--dt", o.pseudo_time_step);
    o.pseudo_time_scale = number("--time-scale", o.pseudo_time_scale);
    o.momentum_relaxation = number("--urf-u", o.momentum_relaxation);
    o.pressure_relaxation = number("--urf-p", o.pressure_relaxation);
    if (args.contains("--walls")) {
      o.force_walls.clear();
      std::istringstream in(text("--walls"));
      for (std::string w; std::getline(in, w, ',');) o.force_walls.push_back(w);
    }
    const auto backend = text("--backend", "metal");
    KUTTA_REQUIRE(backend == "cpu" || backend == "metal" || backend == "cuda", "unknown backend");
    for(const auto& argument:args) KUTTA_REQUIRE(used.contains(argument.first), "unknown option: " + argument.first);
    auto dev = Device::create(backend == "cpu" ? Backend::Cpu : backend == "metal" ? Backend::Metal : Backend::Cuda);
    const auto r = run(*dev, o);
    std::cout << "RESULT {\"source\":\"" << KUTTA_BENCH_COMMIT << "\",\"device\":\"" << dev->name() << "\",\"backend\":\"" << backend
              << "\",\"pressure_scheme\":\"" << pressure_scheme << "\",\"pressure_boundary_gradient\":\"" << boundary_gradient
              << "\",\"laminar_reconstruction\":\"" << laminar_reconstruction << "\",\"result\":";
    print(std::cout, r);
    std::cout << "}\n";
    return r.failure.empty() ? 0 : 2;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
