#pragma once

#include <cmath>
#include <memory>

#include "kernels.hpp"
#include "kutta/linalg/solver.hpp"

namespace kutta::linalg_detail {
inline double norm(Device& dev, const Buffer<float>& x) {
  return std::sqrt(dev.dot(x, x, x.size()));
}
inline void project(Device& dev, Buffer<float>& x, Buffer<float>& scalar) {
  if (x.empty()) return;
  dev.sum_into(x, x.size(), scalar);
  const float count = static_cast<float>(x.size());
  kern::solver_project(dev, x.size(), x, scalar, count,
                       static_cast<float>(double(x.size()) - double(count)));
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
                     Buffer<float>& z, const Buffer<Index>& status) = 0;
  virtual std::string describe() const = 0;
  virtual bool variable() const { return false; }
};

std::unique_ptr<LinearSolver> krylov_solver(SolverKind kind, Device& dev,
    const CsrMatrixHost& a, std::unique_ptr<Preconditioner> preconditioner,
    NullSpace nullspace, std::shared_ptr<KrylovWorkspace> workspace = {});
}  // namespace kutta::linalg_detail
