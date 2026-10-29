#include <charconv>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>

#include "kernels.hpp"
#include "kutta/linalg/amg.hpp"
#include "../tests/linalg/matrices.hpp"

namespace {
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
int positive(std::string_view text) {
  int value = 0;
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
  KUTTA_REQUIRE(error == std::errc{} && end == text.data() + text.size() && value > 0,
                "benchmark: expected a positive integer");
  return value;
}
std::string json_string(std::string_view input) {
  std::string result = "\"";
  for (char c : input) {
    if (c == '\\' || c == '"') result += '\\';
    if (static_cast<unsigned char>(c) < 32) result += ' ';
    else result += c;
  }
  return result + '"';
}
void run(kutta::Device& dev, int size, double anisotropy, const kutta::AmgOptions& options,
         int max_iter) {
  using namespace kutta;
  using namespace kutta::test_linalg;
  KUTTA_REQUIRE(size >= 2 && std::int64_t(size) * size <=
                    std::numeric_limits<Index>::max() / size / 7,
                "benchmark: cube side is outside the 32-bit CSR range");
  auto a = poisson(size, size, size, anisotropy);
  auto matrix = upload(dev, a);
  auto original_values = dev.upload<float>(std::vector<float>(a.val.begin(), a.val.end()));
  std::vector<float> b(a.pattern.n_rows);
  // Deterministic forcing, independent of A and of the selected solver. Unlike
  // an A*random-x RHS, this retains the low modes AMG is intended to remove.
  for (Index i = 0; i < a.pattern.n_rows; ++i)
    b[i] = float(int(hash(i + 53) % 2001) - 1000) / 1024.0f;
  auto db = dev.upload<float>(b);
  auto x = dev.alloc<float>(b.size());
  for (auto kind : {SolverKind::PcgJacobi, SolverKind::AmgPcg}) {
    dev.synchronize();
    auto start = Clock::now();
    auto solver = kind == SolverKind::AmgPcg ? make_amg_solver(kind, dev, a, options) :
                                             make_solver(kind, dev, a);
    dev.synchronize();
    double setup = seconds(start);
    // The same pattern with scaled values models a new outer iteration.
    // Restore before solving so the forcing and matrix are identical.
    kern::scale_copy(dev, matrix.nnz, matrix.val, original_values, 1.25f);
    dev.synchronize(); start = Clock::now();
    solver->refresh(dev, matrix); dev.synchronize();
    double refresh = seconds(start);
    dev.copy(matrix.val, original_values); solver->refresh(dev, matrix);
    dev.fill(x, 0); dev.synchronize(); start = Clock::now();
    auto stats = solver->solve(dev, matrix, x, db, {1e-6, 1e-12, max_iter});
    dev.synchronize(); double solve = seconds(start);
    const char* name = kind == SolverKind::PcgJacobi ? "PcgJacobi" : "AmgPcg";
    std::cout << std::setprecision(9) << "backend=" << backend_name(dev.backend())
              << " device=" << dev.name() << " n=" << size << " anisotropy=" << anisotropy
              << " solver=" << name << " setup_s=" << setup << " refresh_s=" << refresh
              << " iterations=" << stats.iterations << " solve_s=" << solve
              << " relative_residual=" << stats.final_residual / stats.initial_residual
              << " converged=" << stats.converged << " " << solver->describe() << '\n';
    std::cout << "{\"benchmark\":\"linalg\",\"backend\":" << json_string(backend_name(dev.backend()))
              << ",\"device\":" << json_string(dev.name()) << ",\"n\":" << size
              << ",\"rows\":" << a.pattern.n_rows << ",\"nnz\":" << a.pattern.nnz()
              << ",\"anisotropy\":" << anisotropy << ",\"solver\":" << json_string(name)
              << ",\"description\":" << json_string(solver->describe())
              << ",\"rel_tol\":1e-6,\"setup_s\":" << setup << ",\"refresh_s\":" << refresh
              << ",\"iterations\":" << stats.iterations << ",\"solve_s\":" << solve
              << ",\"initial_residual\":" << stats.initial_residual
              << ",\"final_residual\":" << stats.final_residual
              << ",\"converged\":" << (stats.converged ? "true" : "false") << "}" << std::endl;
  }
}
}  // namespace

int main(int argc, char** argv) {
  try {
    int small = 64, large = 128, anisotropic = 64, max_iter = 2000;
    kutta::AmgOptions options;
    for (int i = 1; i < argc; ++i) {
      std::string_view arg = argv[i];
      KUTTA_REQUIRE(i + 1 < argc, "benchmark: option needs a value");
      std::string_view value = argv[++i];
      if (arg == "--small") small = positive(value);
      else if (arg == "--large") large = positive(value);
      else if (arg == "--anisotropic") anisotropic = positive(value);
      else if (arg == "--max-iter") max_iter = positive(value);
      else if (arg == "--passes") options.aggregation_passes = positive(value);
      else if (arg == "--cycle") {
        KUTTA_REQUIRE(value == "v" || value == "k", "benchmark: cycle must be v or k");
        options.cycle = value == "v" ? kutta::AmgCycle::V : kutta::AmgCycle::K;
      } else throw kutta::Error("benchmark: unknown option");
    }
    for (auto backend : kutta::Device::available()) {
      auto dev = kutta::Device::create(backend);
      run(*dev, small, 1, options, max_iter);
      if (large != small) run(*dev, large, 1, options, max_iter);
      run(*dev, anisotropic, 1000, options, max_iter);
    }
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
