#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

#include "../backend_helpers.hpp"
#include "kernels.hpp"

using namespace kutta;

namespace {
void close_vectors(const std::vector<float>& actual, const std::vector<float>& expected) {
  REQUIRE(actual.size() == expected.size());
  for (std::size_t i = 0; i < actual.size(); ++i)
    REQUIRE(std::abs(double(actual[i]) - expected[i]) <= 2e-6 * std::max(1.0, std::abs(double(expected[i]))));
}
void close_reduction(double actual, double expected) {
  CHECK(std::abs(actual - expected) <= 1e-9 * std::max(1.0, std::abs(expected)));
}
}

TEST_CASE("compiled GPU backends are exercised, not silently skipped") {
#if KUTTA_HAS_METAL
  auto device = Device::create(Backend::Metal);
  CHECK(device->backend() == Backend::Metal);
#endif
#if KUTTA_HAS_CUDA
  // A too-old driver makes CUDA drop out of available(); fail loudly instead.
  auto cuda = Device::create(Backend::Cuda);
  CHECK(cuda->backend() == Backend::Cuda);
#endif
  for (auto backend : Device::available()) {
    auto device = Device::create(backend);
    std::cout << "runtime test backend: " << backend_name(backend) << " (" << device->name() << ")\n";
  }
}

TEST_CASE("all generated BLAS kernels match the CPU oracle including tail grids") {
  auto cpu = Device::create(Backend::Cpu);
  for (std::size_t n : {std::size_t{1}, std::size_t{31}, std::size_t{257}, std::size_t{100003}}) {
    std::vector<float> x(n), y(n);
    for (std::size_t i = 0; i < n; ++i) {
      x[i] = std::sin(float(i) * .17f);
      y[i] = std::cos(float(i) * .11f);
    }
    auto cx = cpu->upload<float>(x), cy = cpu->upload<float>(y), cz = cpu->alloc<float>(n);
    KUTTA_FOR_EACH_BACKEND(dev) {
      auto dx = dev->upload<float>(x), dy = dev->upload<float>(y), dz = dev->alloc<float>(n);
      auto reset = [&] { cpu->upload<float>(cy, y); dev->upload<float>(dy, y); };
      reset();
      kern::axpy(*cpu, n, cy, cx, .37f); kern::axpy(*dev, n, dy, dx, .37f);
      close_vectors(dev->download(dy), cpu->download(cy));
      reset();
      kern::xpby(*cpu, n, cy, cx, -.91f); kern::xpby(*dev, n, dy, dx, -.91f);
      close_vectors(dev->download(dy), cpu->download(cy));
      reset();
      kern::axpbyz(*cpu, n, cz, cx, cy, .37f, -.91f);
      kern::axpbyz(*dev, n, dz, dx, dy, .37f, -.91f);
      close_vectors(dev->download(dz), cpu->download(cz));
      kern::scale_copy(*cpu, n, cz, cx, -.91f); kern::scale_copy(*dev, n, dz, dx, -.91f);
      close_vectors(dev->download(dz), cpu->download(cz));
      kern::mul(*cpu, n, cz, cx, cy); kern::mul(*dev, n, dz, dx, dy);
      close_vectors(dev->download(dz), cpu->download(cz));
    }
  }
}

TEST_CASE("transfers fills copies and queued launches preserve call order and lifetime") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    constexpr std::size_t n = 513;
    std::vector<float> input(n, 2.0f);
    auto a = dev->upload<float>(input), b = dev->alloc<float>(n);
    CHECK(dev->download(a) == input);
    dev->fill(b, -3.25f);
    CHECK(dev->download(b) == std::vector<float>(n, -3.25f));
    dev->copy(b, a);
    kern::axpy(*dev, n, b, a, 2.0f);
    CHECK(dev->download(b) == std::vector<float>(n, 6.0f));
    for (int i = 0; i < 600; ++i) {
      kern::xpby(*dev, n, b, a, 0.0f);
      kern::axpy(*dev, n, b, a, 2.0f);
      if (i % 137 == 0) CHECK(dev->download(b) == std::vector<float>(n, 6.0f));
    }
    dev->copy(a, b);
    dev->fill(b, 9.0f);
    CHECK(dev->download(a) == std::vector<float>(n, 6.0f));
    CHECK(dev->download(b) == std::vector<float>(n, 9.0f));
    kern::scale_copy(*dev, n, b, a, 2.0f);
    std::fill(input.begin(), input.end(), 4.0f);
    dev->upload<float>(a, input);
    std::fill(input.begin(), input.end(), -99.0f);
    CHECK(dev->download(b) == std::vector<float>(n, 12.0f));
    CHECK(dev->download(a) == std::vector<float>(n, 4.0f));
    // Pending commands retain their resources after the host handle is destroyed.
    {
      auto temporary = dev->alloc<float>(n);
      dev->fill(temporary, 7.0f);
      dev->copy(b, temporary);
    }
    CHECK(dev->sum(b, n) == 7.0 * n);
    dev->copy(b, b);
    CHECK(dev->download(b) == std::vector<float>(n, 7.0f));

    std::vector<std::int32_t> integers(n, -123456789);
    auto ia = dev->upload<std::int32_t>(integers), ib = dev->alloc<std::int32_t>(n);
    dev->copy(ib, ia);
    CHECK(dev->download(ib) == integers);
    dev->fill(ia, std::int32_t{-7654321});
    CHECK(dev->download(ia) == std::vector<std::int32_t>(n, -7654321));
    std::vector<std::uint32_t> words{0, 1, 0xffffffffU, 0x80000000U};
    auto wa = dev->upload<std::uint32_t>(words), wb = dev->alloc<std::uint32_t>(words.size());
    dev->copy(wb, wa);
    CHECK(dev->download(wb) == words);
  }
}

TEST_CASE("zero sized operations are no-ops and reductions validate bounds") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto a = dev->alloc<float>(0), b = dev->alloc<float>(0);
    dev->upload<float>(a, {});
    dev->fill(a, 2.0f); dev->copy(b, a);
    kern::axpy(*dev, 0, b, a, 1.0f);
    kern::xpby(*dev, 0, b, a, 1.0f);
    kern::axpbyz(*dev, 0, b, a, a, 1.0f, 1.0f);
    kern::scale_copy(*dev, 0, b, a, 1.0f);
    kern::mul(*dev, 0, b, a, a);
    CHECK(dev->download(b).empty());
    CHECK(dev->sum(a, 0) == 0);
    CHECK(dev->dot(a, b, 0) == 0);
    CHECK(dev->max_abs(a, 0) == 0);
    CHECK_THROWS_AS(dev->sum(a, 1), Error);
    CHECK_THROWS_AS(dev->dot(a, b, 1), Error);
    CHECK_THROWS_AS(dev->max_abs(a, 1), Error);
    CHECK_THROWS_AS(dev->launch(static_cast<KernelId>(65535), 0, {}), Error);
    dev->synchronize();
  }
}

TEST_CASE("compensated reductions retain mixed magnitudes cancellation and product residuals") {
  constexpr std::size_t n = 1'000'003;
  std::vector<float> x(n), y(n);
  const std::array<float, 8> pattern{16384.0f, .125f, -16384.0f, .25f, 1.0001f, -1.0001f, .00001f, .75f};
  double sum = 0, dot = 0, maximum = 0;
  for (std::size_t i = 0; i < n; ++i) {
    x[i] = pattern[i % pattern.size()];
    y[i] = i % 8 < 4 ? 1.0f : 1.0003f;
    sum += double(x[i]); dot += double(x[i]) * double(y[i]);
    maximum = std::max(maximum, std::abs(double(x[i])));
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto a = dev->upload<float>(x), b = dev->upload<float>(y);
    close_reduction(dev->sum(a, n), sum);
    close_reduction(dev->dot(a, b, n), dot);
    CHECK(dev->max_abs(a, n) == maximum);
    close_reduction(dev->sum(a, 17), [&] { double s = 0; for (int i = 0; i < 17; ++i) s += x[i]; return s; }());
  }
}

TEST_CASE("reductions keep low terms when high products and sums cancel") {
  constexpr std::size_t n = 1'000'002;
  std::vector<float> x(n), y(n);
  for (std::size_t i = 0; i < n; ++i) {
    x[i] = i % 2 == 0 ? 1.000244140625f : -1.0f;
    y[i] = i % 2 == 0 ? 1.000244140625f : 1.00048828125f;
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto a = dev->upload<float>(x), b = dev->upload<float>(y);
    // Each pair cancels its rounded product, leaving exactly 2^-24.
    const double expected = double(n / 2) * std::ldexp(1.0, -24);
    CHECK(std::abs(dev->dot(a, b, n) - expected) <= expected * 1e-9);
    for (std::size_t i = 0; i < n; ++i)
      x[i] = i % 3 == 0 ? 33554432.0f : (i % 3 == 1 ? 1.0f : -33554432.0f);
    dev->upload<float>(a, x);
    CHECK(dev->sum(a, n) == double(n / 3));
    // Restore the shared host fixture for the next backend.
    for (std::size_t i = 0; i < n; ++i) x[i] = i % 2 == 0 ? 1.000244140625f : -1.0f;
  }
}

TEST_CASE("reductions reach one hundred million elements at double reference accuracy") {
  constexpr std::size_t n = 100'000'000;
  std::vector<float> x(n);
  double sum = 0, dot = 0;
  for (std::size_t i = 0; i < n; ++i) {
    x[i] = 1.0f + float(i % 127) / 128.0f;
    sum += double(x[i]); dot += double(x[i]) * double(x[i]);
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto a = dev->upload<float>(x);
    close_reduction(dev->sum(a, n), sum);
    close_reduction(dev->dot(a, a, n), dot);
    CHECK(dev->max_abs(a, n) == 1.984375);
  }
}

TEST_CASE("identical queued computations and reductions are bitwise deterministic") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    std::array<std::vector<float>, 2> outputs;
    std::array<double, 2> sums{}, dots{}, maxima{};
    for (int run = 0; run < 2; ++run) {
      auto x = dev->alloc<float>(100003), y = dev->alloc<float>(100003);
      dev->fill(x, .1234567f); dev->fill(y, -.7654321f);
      for (int step = 0; step < 300; ++step) {
        kern::xpby(*dev, x.size(), y, x, .75f);
        kern::axpy(*dev, x.size(), x, y, -.03125f);
      }
      outputs[run] = dev->download(x);
      sums[run] = dev->sum(x, x.size());
      dots[run] = dev->dot(x, y, x.size());
      maxima[run] = dev->max_abs(y, y.size());
    }
    CHECK(std::memcmp(outputs[0].data(), outputs[1].data(), outputs[0].size() * sizeof(float)) == 0);
    CHECK(std::memcmp(&sums[0], &sums[1], sizeof(double)) == 0);
    CHECK(std::memcmp(&dots[0], &dots[1], sizeof(double)) == 0);
    CHECK(std::memcmp(&maxima[0], &maxima[1], sizeof(double)) == 0);
  }
}

TEST_CASE("reductions propagate NaN and Inf instead of hiding them") {
  const float nan = std::nanf("");
  const float inf = INFINITY;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto with_nan = dev->upload<float>(std::vector<float>{1, nan, 1, 1});
    CHECK(std::isnan(dev->max_abs(with_nan, 4)));
    CHECK(std::isnan(dev->sum(with_nan, 4)));
    auto with_inf = dev->upload<float>(std::vector<float>{1, inf, 1});
    CHECK(dev->sum(with_inf, 3) == inf);
    CHECK(dev->max_abs(with_inf, 3) == inf);
  }
}

TEST_CASE("math wrappers are accurate on every backend") {
  // Guards against fast-math variants (e.g. Metal's default fast tanh near 0):
  // each backend must agree with the CPU oracle running the same float kernel.
  const std::size_t n = 4096;
  std::vector<float> x(n);
  for (std::size_t i = 0; i < n; ++i) x[i] = -4.0f + 8.0f * float(i) / float(n - 1);
  auto cpu = Device::create(Backend::Cpu);
  auto cpu_x = cpu->upload<float>(x);
  auto cpu_out = cpu->alloc<float>(n);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto bx = dev->upload<float>(x);
    auto out = dev->alloc<float>(n);
    for (int fn = 0; fn < 5; ++fn) {
      kern::test_math(*cpu, n, cpu_out, cpu_x, fn);
      kern::test_math(*dev, n, out, bx, fn);
      const auto ref = cpu->download(cpu_out);
      const auto got = dev->download(out);
      for (std::size_t i = 0; i < n; ++i) {
        INFO("fn " << fn << " x " << x[i]);
        REQUIRE(std::fabs(got[i] - ref[i]) <= 4.0f * std::numeric_limits<float>::epsilon() * std::fabs(ref[i]));
      }
    }
  }
}
