#include <bit>
#include <cstdint>

#include "../backend_helpers.hpp"
#include "instrumented_device.hpp"
#include "kutta/linalg/amg.hpp"
#include "matrices.hpp"

using namespace kutta;
using namespace kutta::test_linalg;

TEST_CASE("A predicted exit removes a large Metal readback without changing solve bits") {
  // Above the unconditional small-system speculation limit. A well-conditioned
  // system exits in its first batch, so the previous solve predicts that batch.
  const auto host = poisson(131073, 1);
  const auto forcing = rhs(host, random_solution(host.pattern.n_rows));
  KUTTA_FOR_EACH_BACKEND(selected) {
    InstrumentedDevice dev(selected->backend());
    auto matrix = upload(dev, host);
    auto b = dev.upload<float>(forcing), x = dev.alloc<float>(b.size()), y = dev.alloc<float>(b.size());
    for (auto kind : {SolverKind::AmgPcg, SolverKind::AmgBicgstab}) {
      INFO(int(kind));
      auto warm = make_solver(kind, dev, host), fresh = make_solver(kind, dev, host);
      const SolveControl control{1e-2, 1e-12, 16};
      dev.fill(x, 0);
      REQUIRE(warm->solve(dev, matrix, x, b, control).converged);
      dev.fill(x, 0);dev.fill(y, 0);
      const auto start = dev.host_waits;
      const auto baseline = fresh->solve(dev, matrix, x, b, control);
      const auto baseline_waits = dev.host_waits - start;
      const auto middle = dev.host_waits;
      const auto candidate = warm->solve(dev, matrix, y, b, control);
      CHECK(dev.host_waits - middle + (dev.backend() == Backend::Metal ? 1 : 0) == baseline_waits);
      CHECK(baseline.iterations == candidate.iterations);
      CHECK(baseline.converged == candidate.converged);
      CHECK(std::bit_cast<std::uint64_t>(baseline.initial_residual) ==
            std::bit_cast<std::uint64_t>(candidate.initial_residual));
      CHECK(std::bit_cast<std::uint64_t>(baseline.final_residual) ==
            std::bit_cast<std::uint64_t>(candidate.final_residual));
      const auto a = dev.download(x), c = dev.download(y);
      for (std::size_t i = 0; i < a.size(); ++i)
        CHECK(std::bit_cast<std::uint32_t>(a[i]) == std::bit_cast<std::uint32_t>(c[i]));
    }
  }
}
