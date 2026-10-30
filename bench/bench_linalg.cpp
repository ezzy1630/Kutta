#include <charconv>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>

#include "kernels.hpp"
#include "kutta/core/parse.hpp"
#include "kutta/linalg/amg.hpp"
#include "../tests/linalg/matrices.hpp"
#include "../tests/linalg/instrumented_device.hpp"

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
void run(kutta::test_linalg::InstrumentedDevice& dev, int size, double anisotropy, const kutta::AmgOptions& options,
         int max_iter, double wind = 0, int ordering = 0, double tolerance = 1e-6) {
  using namespace kutta;
  using namespace kutta::test_linalg;
  KUTTA_REQUIRE(size >= 2 && std::int64_t(size) * size <=
                    std::numeric_limits<Index>::max() / size / 7,
                "benchmark: cube side is outside the 32-bit CSR range");
  auto a = poisson(size, size, size, anisotropy, false, false, wind);
  std::vector<Index> order(a.pattern.n_rows);
  std::iota(order.begin(), order.end(), 0);
  if (ordering != 0) {
    if (ordering == 1) std::reverse(order.begin(), order.end());
    else std::sort(order.begin(), order.end(), [](Index i, Index j) {
      return std::pair{hash(i), i} < std::pair{hash(j), j};
    });
    a = permute(a, order);
  }
  auto matrix = upload(dev, a);
  auto original_values = dev.upload<float>(std::vector<float>(a.val.begin(), a.val.end()));
  std::vector<float> b(a.pattern.n_rows);
  // Deterministic forcing, independent of A and of the selected solver. Unlike
  // an A*random-x RHS, this retains the low modes AMG is intended to remove.
  for (Index i = 0; i < a.pattern.n_rows; ++i)
    b[i] = float(int(hash(order[i] + 53) % 2001) - 1000) / 1024.0f;
  auto db = dev.upload<float>(b);
  auto x = dev.alloc<float>(b.size());
  for (auto kind : {wind == 0 ? SolverKind::PcgJacobi : SolverKind::BicgstabJacobi,
                    wind == 0 ? SolverKind::AmgPcg : SolverKind::AmgBicgstab}) {
    const auto outside_bytes = dev.live_bytes;
    dev.peak_bytes = outside_bytes;
    dev.synchronize();
    auto start = Clock::now();
    auto solver = (kind == SolverKind::AmgPcg || kind == SolverKind::AmgBicgstab) ? make_amg_solver(kind, dev, a, options) :
                                             make_solver(kind, dev, a);
    dev.synchronize();
    double setup = seconds(start);
    const auto solver_bytes = dev.live_bytes - outside_bytes;
    const auto setup_peak_bytes = dev.peak_bytes - outside_bytes;
    // The same pattern with scaled values models a new outer iteration.
    // Restore before solving so the forcing and matrix are identical.
    kern::scale_copy(dev, matrix.nnz, matrix.val, original_values, 1.25f);
    dev.synchronize(); start = Clock::now();
    solver->refresh(dev, matrix); dev.synchronize();
    double refresh = seconds(start);
    dev.copy(matrix.val, original_values); solver->refresh(dev, matrix);
    dev.fill(x, 0); dev.synchronize(); start = Clock::now();
    dev.reductions = dev.downloads = dev.device_reductions = dev.uploads = 0;
    auto stats = solver->solve(dev, matrix, x, db, {tolerance, 1e-12, max_iter});
    dev.synchronize(); double solve = seconds(start);
    const char* name = kind == SolverKind::PcgJacobi ? "PcgJacobi" :
        kind == SolverKind::AmgPcg ? "AmgPcg" :
        kind == SolverKind::BicgstabJacobi ? "BicgstabJacobi" : "AmgBicgstab";
    std::cout << std::setprecision(9) << "backend=" << backend_name(dev.backend())
              << " device=" << dev.name() << " n=" << size << " anisotropy=" << anisotropy
              << " solver=" << name << " setup_s=" << setup << " refresh_s=" << refresh
              << " iterations=" << stats.iterations << " solve_s=" << solve
              << " relative_residual=" << stats.final_residual / stats.initial_residual
              << " reductions=" << dev.reductions << " scalar_polls=" << dev.downloads
              << " uploads=" << dev.uploads << " device_reductions=" << dev.device_reductions << " bytes_per_row=" << double(solver_bytes) / a.pattern.n_rows
              << " setup_peak_bytes_per_row=" << double(setup_peak_bytes) / a.pattern.n_rows
              << " wind=" << wind << " ordering=" << ordering << " converged=" << stats.converged << " " << solver->describe() << '\n';
    std::cout << "{\"benchmark\":\"linalg\",\"backend\":" << json_string(backend_name(dev.backend()))
              << ",\"device\":" << json_string(dev.name()) << ",\"n\":" << size
              << ",\"rows\":" << a.pattern.n_rows << ",\"nnz\":" << a.pattern.nnz()
              << ",\"wind\":" << wind << ",\"ordering\":" << ordering
              << ",\"reductions\":" << dev.reductions
              << ",\"scalar_polls\":" << dev.downloads << ",\"uploads\":" << dev.uploads
              << ",\"device_reductions\":" << dev.device_reductions
              << ",\"solver_bytes\":" << solver_bytes << ",\"setup_peak_bytes\":" << setup_peak_bytes
              << ",\"anisotropy\":" << anisotropy << ",\"solver\":" << json_string(name)
              << ",\"description\":" << json_string(solver->describe())
              << ",\"source_commit\":" << json_string(KUTTA_BENCH_COMMIT)
              << ",\"rel_tol\":" << tolerance << ",\"setup_s\":" << setup << ",\"refresh_s\":" << refresh
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
    int wind = 0, ordering = 0;
    double tolerance = 1e-6;
    std::string_view selected_backend = "all";
    for (int i = 1; i < argc; ++i) {
      std::string_view arg = argv[i];
      KUTTA_REQUIRE(i + 1 < argc, "benchmark: option needs a value");
      std::string_view value = argv[++i];
      if (arg == "--small") small = positive(value);
      else if (arg == "--large") large = positive(value);
      else if (arg == "--anisotropic") anisotropic = positive(value);
      else if (arg == "--max-iter") max_iter = positive(value);
      else if (arg == "--backend") selected_backend = value;
      else if (arg == "--tolerance") {
        const auto parsed = kutta::parse_double(value);
        KUTTA_REQUIRE(parsed && std::isfinite(*parsed) && *parsed > 0, "benchmark: invalid tolerance");
        tolerance = *parsed;
      }
      else if (arg == "--wind") wind = positive(value);
      else if (arg == "--ordering") ordering = positive(value);
      else if (arg == "--sweeps") options.pre_sweeps = options.post_sweeps = positive(value);
      else if (arg == "--passes") options.aggregation_passes = positive(value);
      else if (arg == "--cycle") {
        KUTTA_REQUIRE(value == "v" || value == "k", "benchmark: cycle must be v or k");
        options.cycle = value == "v" ? kutta::AmgCycle::V : kutta::AmgCycle::K;
      } else throw kutta::Error("benchmark: unknown option");
    }
    for (auto backend : kutta::Device::available()) {
      if (selected_backend != "all" && selected_backend != kutta::backend_name(backend)) continue;
      kutta::test_linalg::InstrumentedDevice dev(backend);
      run(dev, small, 1, options, max_iter, wind, ordering, tolerance);
      if (large != small) run(dev, large, 1, options, max_iter, wind, ordering, tolerance);
      if (wind == 0) run(dev, anisotropic, 1000, options, max_iter, 0, 0, tolerance);

    }
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
