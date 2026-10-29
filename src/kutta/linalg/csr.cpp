#include "kutta/linalg/csr.hpp"

#include <cmath>
#include <limits>

#include "kernels.hpp"
#include "kutta/linalg/detail.hpp"

namespace kutta {
namespace {
void validate(const CsrPattern& p) {
  KUTTA_REQUIRE(p.n_rows >= 0, "CSR: negative row count");
  KUTTA_REQUIRE(p.col.size() <= std::size_t(std::numeric_limits<Index>::max()),
                "CSR: exceeds 32-bit entry limit");
  KUTTA_REQUIRE(p.row_ptr.size() == std::size_t(p.n_rows) + 1 &&
                    p.diag.size() == std::size_t(p.n_rows), "CSR: invalid array sizes");
  KUTTA_REQUIRE(p.row_ptr.front() == 0 && p.row_ptr.back() == p.nnz(),
                "CSR: invalid row offsets");
  for (Index i = 0; i < p.n_rows; ++i) {
    KUTTA_REQUIRE(p.row_ptr[i] >= 0 && p.row_ptr[i] <= p.row_ptr[i + 1] &&
                      p.row_ptr[i + 1] <= p.nnz(), "CSR: invalid row range");
    Index previous = -1;
    for (Index j = p.row_ptr[i]; j < p.row_ptr[i + 1]; ++j) {
      KUTTA_REQUIRE(p.col[j] > previous && p.col[j] < p.n_rows,
                    "CSR: columns must be in range, unique and sorted");
      previous = p.col[j];
    }
    KUTTA_REQUIRE(p.diag[i] >= p.row_ptr[i] && p.diag[i] < p.row_ptr[i + 1] &&
                      p.col[p.diag[i]] == i, "CSR: missing or incorrect diagonal");
  }
}

void check_vectors(const CsrMatrix& a, const Buffer<float>& x, const Buffer<float>& y) {
  KUTTA_REQUIRE(x.size() == std::size_t(a.n_rows) && y.size() == x.size(),
                "CSR: vector size mismatch");
  KUTTA_REQUIRE(x.handle() != y.handle(), "CSR: input/output aliasing is unsupported");
}
}  // namespace

CsrMatrix upload(Device& dev, const CsrPattern& p) {
  validate(p);
  CsrMatrix a;
  a.n_rows = p.n_rows;
  a.nnz = p.nnz();
  a.row_ptr = dev.upload<Index>(p.row_ptr);
  a.col = p.col.empty() ? dev.alloc<Index>(0) : dev.upload<Index>(p.col);
  a.diag = p.diag.empty() ? dev.alloc<Index>(0) : dev.upload<Index>(p.diag);
  a.val = dev.alloc<float>(p.col.size());
  dev.fill(a.val, 0.0f);
  return a;
}

void linalg_detail::validate_host(const CsrMatrixHost& m) {
  validate(m.pattern);
  KUTTA_REQUIRE(m.val.size() == m.pattern.col.size(), "CSR: value count mismatch");
  for (double v : m.val)
    KUTTA_REQUIRE(std::isfinite(v) && std::abs(v) <= std::numeric_limits<float>::max(),
                  "CSR: nonfinite or out-of-range value");
}

CsrMatrix upload(Device& dev, const CsrMatrixHost& m) {
  KUTTA_REQUIRE(m.val.size() == m.pattern.col.size(), "CSR: value count mismatch");
  std::vector<float> values;
  values.reserve(m.val.size());
  for (double v : m.val) {
    KUTTA_REQUIRE(std::isfinite(v) && std::abs(v) <= std::numeric_limits<float>::max(),
                  "CSR: nonfinite or out-of-range value");
    values.push_back(static_cast<float>(v));
  }
  auto a = upload(dev, m.pattern);
  if (!values.empty()) dev.upload<float>(a.val, values);
  return a;
}

void spmv(Device& dev, const CsrMatrix& a, const Buffer<float>& x, Buffer<float>& y) {
  check_vectors(a, x, y);
  kern::csr_spmv(dev, a.n_rows, a.row_ptr, a.col, a.val, x, y);
}

void residual(Device& dev, const CsrMatrix& a, const Buffer<float>& x,
              const Buffer<float>& b, Buffer<float>& r) {
  check_vectors(a, x, r);
  KUTTA_REQUIRE(b.size() == r.size() && b.handle() != r.handle(),
                "CSR: residual RHS size or alias mismatch");
  kern::csr_residual(dev, a.n_rows, a.row_ptr, a.col, a.val, x, b, r);
}
}  // namespace kutta
