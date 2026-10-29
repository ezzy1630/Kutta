#pragma once

#include <cmath>
#include <memory>
#include <limits>

#include "kernels.hpp"
#include "kutta/linalg/solver.hpp"

namespace kutta::linalg_detail {
inline double norm(Device& dev, const Buffer<float>& x) {
  return std::sqrt(dev.dot(x, x, x.size()));
}
inline void project(Device& dev, Buffer<float>& x) {
  if (!x.empty()) kern::krylov_shift(dev, x.size(), x,
      static_cast<float>(dev.sum(x, x.size()) / double(x.size())));
}

inline bool coefficient(double x) {
  return std::isfinite(x) && std::abs(x) <= std::numeric_limits<float>::max() &&
         (x == 0 || std::abs(x) >= std::numeric_limits<float>::denorm_min());
}

// Structure and nullspace are fixed for a solver's lifetime. A refresh may
// change coefficients, but must retain that pattern and the nullspace choice.
void validate_host(const CsrMatrixHost& a);
std::vector<float> diagonal_values(const CsrMatrixHost& a, bool l1, bool nullspace);
void check_nullspace(const CsrMatrixHost& a, NullSpace choice);

class Preconditioner {
 public:
  virtual ~Preconditioner() = default;
  virtual void refresh(Device& dev, const CsrMatrix& a) = 0;
  virtual void apply(Device& dev, const CsrMatrix& a, const Buffer<float>& r,
                     Buffer<float>& z) = 0;
  virtual std::string describe() const = 0;
  virtual bool variable() const { return false; }
};

std::unique_ptr<LinearSolver> krylov_solver(SolverKind kind, Device& dev,
    const CsrMatrixHost& a, std::unique_ptr<Preconditioner> preconditioner,
    NullSpace nullspace);
}  // namespace kutta::linalg_detail
