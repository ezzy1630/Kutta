#include <cmath>
#include <vector>

#include "../backend_helpers.hpp"
#include "kernels.hpp"

using namespace kutta;

TEST_CASE("blas kernels and reductions match the reference on every backend") {
  const std::size_t n = 100'003;  // not a multiple of any block size
  std::vector<float> x(n), y(n);
  for (std::size_t i = 0; i < n; ++i) {
    x[i] = std::sin(0.001f * float(i));
    y[i] = 1.0f + 0.5f * std::cos(0.002f * float(i));
  }
  std::vector<float> expect(n);
  for (std::size_t i = 0; i < n; ++i) expect[i] = y[i] + 2.5f * x[i];
  double dot_ref = 0;
  for (std::size_t i = 0; i < n; ++i) dot_ref += double(x[i]) * double(expect[i]);

  KUTTA_FOR_EACH_BACKEND(dev) {
    auto bx = dev->upload<float>(x);
    auto by = dev->upload<float>(y);
    kern::axpy(*dev, n, by, bx, 2.5f);
    const auto got = dev->download(by);
    for (std::size_t i = 0; i < n; ++i) REQUIRE(got[i] == doctest::Approx(expect[i]).epsilon(1e-6));
    CHECK(dev->dot(bx, by, n) == doctest::Approx(dot_ref).epsilon(1e-9));
    CHECK(dev->max_abs(bx, n) == doctest::Approx(1.0).epsilon(1e-6));

    auto bz = dev->alloc<float>(n);
    dev->fill(bz, 3.0f);
    CHECK(dev->sum(bz, n) == doctest::Approx(3.0 * double(n)).epsilon(1e-12));
  }
}
