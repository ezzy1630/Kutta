#include <array>
#include <bit>
#include <cstdint>
#include <limits>

#include "../backend_helpers.hpp"
#include "kernels.hpp"

using namespace kutta;

TEST_CASE("AMG zero-start sweep preserves the full CSR sweep's bits") {
  // Seven finite entries per row span product-splitting and underflow ranges.
  // RHS values include signed zero, subnormals, infinities and a NaN payload.
  const std::array<float, 12> rhs_values{0.f, -0.f, 1.f, -1.f, 1e30f, -1e30f,
      1e-30f, -1e-30f, std::numeric_limits<float>::denorm_min(),
      std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(),
      std::bit_cast<float>(std::uint32_t{0x7fc01234})};
  const std::array<float, 8> matrix_values{0.f, -0.f, 1.f, -1.f, 1e30f, -1e30f, 1e-30f, -1e-30f};
  for (Index n : {513, 2053}) {
    std::vector<Index> rows(n + 1), cols(7 * n);
    std::vector<float> values(7 * n), b(n), inverse(n);
    for (Index i = 0; i < n; ++i) {
      rows[i] = 7 * i;
      b[i] = rhs_values[i % rhs_values.size()];
      inverse[i] = matrix_values[(i / rhs_values.size()) % matrix_values.size()];
      for (Index j = 0; j < 7; ++j) {
        cols[7 * i + j] = (i + j * 31) % n;
        values[7 * i + j] = matrix_values[(i + j) % matrix_values.size()];
      }
    }
    rows[n] = 7 * n;
    KUTTA_FOR_EACH_BACKEND(dev) {
      auto dr = dev->upload<Index>(rows), dc = dev->upload<Index>(cols);
      auto da = dev->upload<float>(values), db = dev->upload<float>(b), di = dev->upload<float>(inverse);
      auto x = dev->alloc<float>(n), expected = dev->alloc<float>(n), actual = dev->alloc<float>(n);
      auto status = dev->alloc<Index>(3);
      for (Index phase : {0, 1, 2}) {
        dev->fill(x, 0); dev->fill(status, phase);
        dev->fill(expected, 42); dev->fill(actual, 42);
        kern::amg_smooth(*dev, n, dr, dc, da, x, db, di, expected, status);
        kern::amg_smooth_zero(*dev, n, dr, dc, da, x, db, di, actual, status);
        const auto a = dev->download(expected), c = dev->download(actual);
        for (Index i = 0; i < n; ++i)
          CHECK(std::bit_cast<std::uint32_t>(a[i]) == std::bit_cast<std::uint32_t>(c[i]));
      }
    }
  }
}
