#pragma once

// Sparse matrices. The finite-volume operators all share one sparsity pattern
// per mesh (cell adjacency); it is built once on the host, uploaded once, and
// only the values change every outer iteration.

#include <vector>

#include "kutta/core/types.hpp"
#include "kutta/runtime/device.hpp"

namespace kutta {

// Host CSR pattern. Columns sorted ascending within each row; the diagonal is
// always present.
struct CsrPattern {
  Index n_rows = 0;
  std::vector<Index> row_ptr{0};  // size n_rows+1
  std::vector<Index> col;         // size nnz
  std::vector<Index> diag;        // size n_rows: index into col of (i, i)

  Index nnz() const { return static_cast<Index>(col.size()); }
};

struct CsrMatrixHost {
  CsrPattern pattern;
  std::vector<double> val;  // size nnz
};

// Device CSR matrix. `val` is the only part expected to change after setup.
struct CsrMatrix {
  Index n_rows = 0;
  Index nnz = 0;
  Buffer<std::int32_t> row_ptr;
  Buffer<std::int32_t> col;
  Buffer<std::int32_t> diag;
  Buffer<float> val;
};

CsrMatrix upload(Device& dev, const CsrPattern& pattern);
CsrMatrix upload(Device& dev, const CsrMatrixHost& m);

// y = A x
void spmv(Device& dev, const CsrMatrix& a, const Buffer<float>& x, Buffer<float>& y);

// r = b - A x
void residual(Device& dev, const CsrMatrix& a, const Buffer<float>& x, const Buffer<float>& b,
              Buffer<float>& r);

}  // namespace kutta
