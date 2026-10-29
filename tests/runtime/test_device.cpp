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

TEST_CASE("device scalar reductions preserve low terms and ordered scratch reuse") {
  // Cancellation and products that cannot be represented by one float.
  const std::size_t n = 100003;
  std::vector<float> a(n), b(n), c(n), d(n);
  double sum = 0, dot = 0, dot2 = 0;
  for (std::size_t i = 0; i < n; ++i) {
    a[i] = i % 3 == 0 ? 1024.0f : (i % 3 == 1 ? -1024.0f : 0.000123f);
    b[i] = 1.00000011920928955078125f;
    c[i] = float(int(i % 31) - 15) / 37;
    d[i] = float(int(i % 17) - 8) / 19;
    sum += a[i]; dot += double(a[i]) * b[i]; dot2 += double(c[i]) * d[i];
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto x = dev->upload<float>(a), y = dev->upload<float>(b);
    auto u = dev->upload<float>(c), v = dev->upload<float>(d);
    auto out = dev->alloc<float>(10);
    dev->fill(out, 42);
    dev->sum_into(x, n, out, 1);
    dev->dot_into(x, y, n, out, 3);
    dev->dot2_into(x, y, u, v, n, out, 5);
    const auto got = dev->download(out);
    auto pair = [&](int i) { return double(got[i]) + got[i + 1]; };
    CHECK(got[0] == 42); CHECK(got[9] == 42);
    CHECK(pair(1) == doctest::Approx(sum).epsilon(1e-10));
    CHECK(pair(3) == doctest::Approx(dot).epsilon(1e-10));
    CHECK(pair(5) == doctest::Approx(dot).epsilon(1e-10));
    CHECK(pair(7) == doctest::Approx(dot2).epsilon(1e-10));
    dev->dot2_into(x, y, u, v, n, out, 5);
    CHECK(dev->download(out) == got); // fixed reduction order
    dev->dot2_into(x, y, u, v, 0, out);
    const auto empty = dev->download(out);
    for (int i = 0; i < 4; ++i) CHECK(empty[i] == 0);
    CHECK_THROWS_AS(dev->dot_into(x, y, n + 1, out), Error);
    CHECK_THROWS_AS(dev->dot2_into(x, y, u, v, n, out, 7), Error);
    CHECK_THROWS_AS(dev->sum_into(x, n, x), Error);
    auto bad = dev->upload<float>(std::vector<float>{std::numeric_limits<float>::infinity()});
    dev->sum_into(bad, 1, out);
    CHECK_FALSE(std::isfinite(dev->download(out)[0]));
  }
}
