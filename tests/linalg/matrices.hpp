#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>

#include "kutta/linalg/csr.hpp"

namespace kutta::test_linalg {
inline std::uint32_t hash(std::uint32_t x) {
  x ^= x >> 16; x *= 0x7feb352dU; x ^= x >> 15; x *= 0x846ca68bU; return x ^ (x >> 16);
}

inline CsrMatrixHost poisson(Index nx, Index ny, Index nz = 1, double anisotropy = 1,
                            bool random = false, bool neumann = false, double wind = 0) {
  CsrMatrixHost a;
  auto& p = a.pattern;
  p.n_rows = nx * ny * nz;
  p.diag.resize(p.n_rows);
  const Index dims[] = {nx, ny, nz};
  const Index strides[] = {1, nx, nx * ny};
  const int dimension = nz == 1 ? 2 : 3;
  for (Index i = 0; i < p.n_rows; ++i) {
    double diagonal = 0;
    std::vector<std::pair<Index, double>> row;
    for (int d = 0; d < dimension; ++d) {
      Index coord = (i / strides[d]) % dims[d];
      for (int sign : {-1, 1}) {
        Index j = i + sign * strides[d];
        bool inside = coord + sign >= 0 && coord + sign < dims[d];
        double weight = d == 0 ? anisotropy : 1;
        if (random && inside) {
          auto key = hash(std::uint32_t(std::min(i, j)) ^
                          (hash(std::uint32_t(std::max(i, j))) + 0x9e3779b9U));
          weight *= 0.1 + 9.9 * double(key % 1000) / 999;
        }
        if (inside) row.emplace_back(j, -weight);
        if (inside || !neumann) diagonal += weight;
      }
    }
    if (wind != 0) {
      diagonal += wind;
      if (i % nx != 0) {
        auto left = std::find_if(row.begin(), row.end(), [i](auto e) { return e.first == i - 1; });
        left->second -= wind;
      }
    }
    row.emplace_back(i, diagonal);
    std::sort(row.begin(), row.end());
    for (auto [j, value] : row) {
      if (i == j) p.diag[i] = static_cast<Index>(p.col.size());
      p.col.push_back(j);
      // Compare against exactly the operator stored by the device.
      a.val.push_back(static_cast<float>(value));
    }
    p.row_ptr.push_back(static_cast<Index>(p.col.size()));
  }
  return a;
}

inline std::vector<double> multiply(const CsrMatrixHost& a, const std::vector<double>& x) {
  std::vector<double> b(a.pattern.n_rows);
  for (Index i = 0; i < a.pattern.n_rows; ++i)
    for (Index j = a.pattern.row_ptr[i]; j < a.pattern.row_ptr[i + 1]; ++j)
      b[i] += a.val[j] * x[a.pattern.col[j]];
  return b;
}
inline std::vector<float> random_solution(Index n) {
  std::vector<float> x(n);
  for (Index i = 0; i < n; ++i) x[i] = float(int(hash(i + 19) % 2001) - 1000) / 1024.0f;
  return x;
}
inline std::vector<float> rhs(const CsrMatrixHost& a, const std::vector<float>& x) {
  auto b = multiply(a, std::vector<double>(x.begin(), x.end()));
  return std::vector<float>(b.begin(), b.end());
}
inline double dot(const std::vector<double>& a, const std::vector<double>& b) {
  return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}
inline void project(std::vector<double>& x) {
  const double mean = std::accumulate(x.begin(), x.end(), 0.0) / double(x.size());
  for (auto& v : x) v -= mean;
}

// Independent host double-precision CG reference, used only for small SPD tests.
inline std::vector<double> reference_cg(const CsrMatrixHost& a, const std::vector<float>& b,
                                      bool neumann = false) {
  std::vector<double> x(b.size()), r(b.begin(), b.end()), p = r;
  if (neumann) { project(r); p = r; }
  double rr = dot(r, r), start = rr;
  for (int iter = 0; iter < 10000 && rr > 1e-26 * start; ++iter) {
    auto q = multiply(a, p);
    double alpha = rr / dot(p, q);
    for (std::size_t i = 0; i < x.size(); ++i) { x[i] += alpha * p[i]; r[i] -= alpha * q[i]; }
    if (neumann) project(r);
    double next = dot(r, r);
    for (std::size_t i = 0; i < x.size(); ++i) p[i] = r[i] + (next / rr) * p[i];
    rr = next;
  }
  if (neumann) project(x);
  return x;
}

// Partial-pivoted double dense solve, deliberately independent of Krylov/AMG.
inline std::vector<double> reference_dense(const CsrMatrixHost& a, const std::vector<float>& rhs) {
  const auto n = a.pattern.n_rows;
  std::vector<double> matrix(std::size_t(n) * n), b(rhs.begin(), rhs.end());
  for (Index i = 0; i < n; ++i)
    for (Index j = a.pattern.row_ptr[i]; j < a.pattern.row_ptr[i + 1]; ++j)
      matrix[std::size_t(i) * n + a.pattern.col[j]] = a.val[j];
  for (Index k = 0; k < n; ++k) {
    Index pivot = k;
    for (Index i = k + 1; i < n; ++i)
      if (std::abs(matrix[std::size_t(i) * n + k]) >
          std::abs(matrix[std::size_t(pivot) * n + k])) pivot = i;
    for (Index j = k; j < n; ++j) std::swap(matrix[std::size_t(k) * n + j], matrix[std::size_t(pivot) * n + j]);
    std::swap(b[k], b[pivot]);
    for (Index i = k + 1; i < n; ++i) {
      double factor = matrix[std::size_t(i) * n + k] / matrix[std::size_t(k) * n + k];
      for (Index j = k + 1; j < n; ++j) matrix[std::size_t(i) * n + j] -= factor * matrix[std::size_t(k) * n + j];
      b[i] -= factor * b[k];
    }
  }
  std::vector<double> x(n);
  for (Index i = n; i-- > 0;) {
    double v = b[i];
    for (Index j = i + 1; j < n; ++j) v -= matrix[std::size_t(i) * n + j] * x[j];
    x[i] = v / matrix[std::size_t(i) * n + i];
  }
  return x;
}
}  // namespace kutta::test_linalg
