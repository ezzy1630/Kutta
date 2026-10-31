// Runs one FSAE aero validation solve (validation/cases/fsae_aero.hpp) and
// prints a RESULT line; tools/fsae_aero/run.py drives the sweeps.
#include <iostream>
#include <map>
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
      const bool flag = key == "--moving-ground" || key == "--coupled" || key == "--downforce" || key == "--cells" || key == "--lift-z" || key == "--hybrid";
      KUTTA_REQUIRE(flag || i + 1 < argc, "missing value for " + key);
      args[key] = flag ? "1" : argv[++i];
    }
    const auto text = [&](const std::string& k, std::string d = {}) { return args.contains(k) ? args.at(k) : d; };
    const auto number = [&](const std::string& k, double d) {
      if (!args.contains(k)) return d;
      std::istringstream in(args.at(k));
      in.imbue(std::locale::classic());
      double v = 0;
      in >> v;
      KUTTA_REQUIRE(in && in.eof(), "invalid number for " + k);
      return v;
    };
    Options o;
    o.mesh = text("--mesh");
    o.output = text("--out");
    o.reynolds = number("--re", 0);
    o.length = number("--length", 1);
    o.reference_area = number("--area", 1);
    o.moving_ground = args.contains("--moving-ground");
    o.coupled = args.contains("--coupled");
    o.hybrid = args.contains("--hybrid");
    o.write_cells = args.contains("--cells");
    if (args.contains("--downforce")) o.lift_direction = {0, -1, 0};
    if (args.contains("--lift-z")) o.lift_direction = {0, 0, 1};
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
      std::istringstream in(args.at("--walls"));
      for (std::string w; std::getline(in, w, ',');) o.force_walls.push_back(w);
    }
    const auto backend = text("--backend", "metal");
    KUTTA_REQUIRE(backend == "cpu" || backend == "metal" || backend == "cuda", "unknown backend");
    auto dev = Device::create(backend == "cpu" ? Backend::Cpu : backend == "metal" ? Backend::Metal : Backend::Cuda);
    const auto r = run(*dev, o);
    std::cout << "RESULT {\"source\":\"" << KUTTA_BENCH_COMMIT << "\",\"device\":\"" << dev->name() << "\",\"backend\":\"" << backend
              << "\",\"result\":";
    print(std::cout, r);
    std::cout << "}\n";
    return r.failure.empty() ? 0 : 2;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
