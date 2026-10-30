#include <cmath>
#include <limits>

#include "../backend_helpers.hpp"
#include "kernels.hpp"
#include "kutta/linalg/csr.hpp"
#include "matrices.hpp"

using namespace kutta;

TEST_CASE("linalg kernels match explicit host calculations on every backend") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto r = dev->upload<float>(std::vector<float>{2, -4, 6});
    auto x = dev->upload<float>(std::vector<float>{1, 2, 3});
    auto inverse = dev->upload<float>(std::vector<float>{0.5f, 0.25f, 2});
    auto out = dev->alloc<float>(3);
    kern::krylov_multiply(*dev, 3, r, inverse, out);
    CHECK(dev->download(out) == std::vector<float>{1, -1, 12});
    kern::krylov_smooth(*dev, 3, r, inverse, x, 0.5f);
    CHECK(dev->download(x) == std::vector<float>{1.5f, 1.5f, 9});
    kern::krylov_shift(*dev, 3, x, 1.0f);
    CHECK(dev->download(x) == std::vector<float>{0.5f, 0.5f, 8});
    kern::krylov_direction(*dev, 3, r, x, 2);
    CHECK(dev->download(x) == std::vector<float>{3, -3, 22});
    kern::krylov_step(*dev, 3, x, out, r, inverse, 2);
    CHECK(dev->download(x) == std::vector<float>{7, -11, 34});
    CHECK(dev->download(out) == std::vector<float>{0, -1.5f, 8});
    kern::krylov_bicg_direction(*dev, 3, r, inverse, x, 0.5f, 2);
    CHECK(dev->download(x) == std::vector<float>{5, -9.75f, 21});
    kern::krylov_difference(*dev, 3, r, inverse, out, 2);
    CHECK(dev->download(out) == std::vector<float>{1, -4.5f, 2});
    kern::krylov_bicg_step(*dev, 3, x, r, inverse, 2, 4);
    CHECK(dev->download(x) == std::vector<float>{11, -16.75f, 41});

    auto rows = dev->upload<Index>(std::vector<Index>{0, 2, 5, 7});
    auto diag = dev->upload<Index>(std::vector<Index>{0, 3, 6});
    auto values = dev->upload<float>(std::vector<float>{2, -1, -1, 4, -2, -2, 8});
    auto l1 = dev->alloc<float>(3), invalid = dev->alloc<float>(3);
    kern::krylov_diagonal(*dev, 3, rows, diag, values, out, invalid, 0);
    CHECK(dev->download(out) == std::vector<float>{0.5f, 0.25f, 0.125f});
    kern::krylov_diagonal(*dev, 3, rows, diag, values, l1, invalid, 1);
    CHECK(dev->download(l1) == dev->download(out));
    CHECK(dev->max_abs(invalid, 3) == 0);
    dev->fill(values, 0);
    kern::krylov_diagonal(*dev, 3, rows, diag, values, out, invalid, 0);
    CHECK(dev->download(invalid) == std::vector<float>{1, 1, 1});
  }
}

TEST_CASE("AMG restriction prolongation and Galerkin gather match host sums") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto offsets = dev->upload<Index>(std::vector<Index>{0, 2, 5});
    auto members = dev->upload<Index>(std::vector<Index>{0, 2, 1, 3, 4});
    auto map = dev->upload<Index>(std::vector<Index>{0, 1, 0, 1, 1});
    auto fine = dev->upload<float>(std::vector<float>{1, 2, 3, 4, 5});
    auto coarse = dev->alloc<float>(2);
    kern::amg_restrict(*dev, 2, offsets, members, fine, coarse);
    CHECK(dev->download(coarse) == std::vector<float>{4, 11});
    kern::amg_prolong(*dev, 5, map, coarse, fine);
    CHECK(dev->download(fine) == std::vector<float>{5, 13, 7, 15, 16});
    kern::amg_galerkin(*dev, 2, offsets, members, fine, coarse);
    CHECK(dev->download(coarse) == std::vector<float>{12, 44});
  }
}

TEST_CASE("compensated CSR residual resolves product and cancellation error") {
  CsrMatrixHost a;
  a.pattern.n_rows = 3;
  a.pattern.row_ptr = {0, 3, 6, 9};
  a.pattern.col = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  a.pattern.diag = {0, 4, 8};
  a.val = {6, -1, -1, -1, 6, -1, -1, -1, 6};
  std::vector<float> xx = {std::nextafter(1.0f, 2.0f), -1.0f, 1.0f};
  std::vector<float> bb = {6, -8, 6};
  std::vector<float> reference(3);
  for (Index i = 0; i < 3; ++i) {
    double value = bb[i];
    for (Index j = a.pattern.row_ptr[i]; j < a.pattern.row_ptr[i + 1]; ++j)
      value -= a.val[j] * double(xx[a.pattern.col[j]]);
    reference[i] = static_cast<float>(value);
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto matrix = upload(*dev, a);
    auto x = dev->upload<float>(xx), b = dev->upload<float>(bb), r = dev->alloc<float>(3);
    residual(*dev, matrix, x, b, r);
    auto result = dev->download(r);
    for (std::size_t i = 0; i < result.size(); ++i)
      CHECK(result[i] == doctest::Approx(reference[i]).epsilon(1e-6));
  }
}

TEST_CASE("fused coarse Jacobi preserves the ordinary sweep arithmetic") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    for (Index count : {1, 17, 24}) {
      auto a = kutta::test_linalg::poisson(count, 1);
      auto matrix = upload(*dev, a);
      auto forcing = kutta::test_linalg::random_solution(count);
      auto b = dev->upload<float>(forcing), x = dev->alloc<float>(count);
      auto y = dev->alloc<float>(count), r = dev->alloc<float>(count);
      std::vector<float> inverse(count);
      for (Index row = 0; row < count; ++row) inverse[row] = 1.0f / float(a.val[a.pattern.diag[row]]);
      auto di = dev->upload<float>(inverse);
      auto status = dev->alloc<Index>(1); dev->fill(status, 0);
      for (int sweeps : {1, 5, 32}) {
        dev->fill(x, 0); dev->fill(y, 0);
        for (int sweep = 0; sweep < sweeps; ++sweep) {
          residual(*dev, matrix, y, b, r);
          kern::krylov_smooth(*dev, count, r, di, y, 2.0f / 3.0f);
        }
        kern::amg_coarse_solve(*dev, 1, matrix.row_ptr, matrix.col, matrix.val,
                               x, b, di, status, count, sweeps);
        CHECK(dev->download(x) == dev->download(y));
        dev->fill(x, 0);
        for (int sweep = 0; sweep < sweeps; ++sweep) {
          kern::amg_smooth(*dev, count, matrix.row_ptr, matrix.col, matrix.val,
                            x, b, di, r, status);
          std::swap(x, r);
        }
        CHECK(dev->download(x) == dev->download(y));
      }
      dev->fill(status, 2); dev->fill(x, 42);
      kern::amg_coarse_solve(*dev, 1, matrix.row_ptr, matrix.col, matrix.val,
                             x, b, di, status, count, 32);
      CHECK(dev->download(x) == std::vector<float>(count, 42));
    }
  }
}

TEST_CASE("device coefficient guards hold vectors on zero or nonfinite rho and omega") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto status = dev->alloc<Index>(3);
    auto scalars = dev->alloc<float>(24);
    auto x = dev->upload<float>(std::vector<float>{2, -3});
    auto p = dev->upload<float>(std::vector<float>{1, 1});
    auto r = dev->upload<float>(std::vector<float>{4, 5});
    for (int slot : {0, 12}) {
      for (float bad : {0.0f, std::numeric_limits<float>::quiet_NaN(),
                         std::numeric_limits<float>::infinity()}) {
        std::vector<float> values(24, 0);
        values[0] = values[8] = values[12] = values[14] = 1;
        values[slot] = bad;
        dev->upload<float>(scalars, values);
        dev->fill(status, 0);
        kern::solver_bicg_beta(*dev, 1, scalars, status);
        kern::solver_step(*dev, 2, x, r, p, p, scalars, status);
        CHECK(dev->download(status)[0] == 2);
        CHECK(dev->download(x) == std::vector<float>{2, -3});
        CHECK(dev->download(r) == std::vector<float>{4, 5});
      }
    }
  }
}

TEST_CASE("AMG Galerkin preserves weak terms across strong cancellation") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto offsets=dev->upload<Index>(std::vector<Index>{0,3});
    auto entries=dev->upload<Index>(std::vector<Index>{0,1,2});
    auto fine=dev->upload<float>(std::vector<float>{1e8f,1,-1e8f});auto coarse=dev->alloc<float>(1);
    kern::amg_galerkin(*dev,1,offsets,entries,fine,coarse);
    CHECK(dev->download(coarse)[0]==1);
    auto rows=dev->upload<Index>(std::vector<Index>{0,1});
    auto cols=dev->upload<Index>(std::vector<Index>{0});
    auto members=dev->upload<Index>(std::vector<Index>{0});
    auto fine_rows=dev->upload<Index>(std::vector<Index>{0,3});
    auto fine_cols=dev->upload<Index>(std::vector<Index>{0,0,0});
    auto aggregate=dev->upload<Index>(std::vector<Index>{0});
    kern::amg_galerkin_rows(*dev,1,rows,cols,rows,members,fine_rows,fine_cols,aggregate,fine,coarse);
    CHECK(dev->download(coarse)[0]==1);
  }
}
