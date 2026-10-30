#pragma once

// Iterative linear solvers. A solver is set up once per (pattern, solver kind)
// and then reused every outer iteration with new matrix values; refresh()
// recomputes value-dependent data (AMG coarse operators, Jacobi diagonals) on
// the device without rebuilding structure.
//
// Residual norms are reported Fluent/AMG-style: ||b - A x||_2, accumulated in
// double (see Device reductions), relative to the initial residual.

#include <memory>
#include <string>

#include "kutta/linalg/csr.hpp"

namespace kutta {

struct SolveControl {
  double rel_tol = 1e-2;  // outer SIMPLE loops only need ~1-2 orders per solve
  double abs_tol = 1e-12;
  int max_iter = 200;
  int residual_check_interval = 8;  // updates between host convergence checks
};

struct SolveStats {
  int iterations = 0;
  double initial_residual = 0;
  double final_residual = 0;
  bool converged = false;
};

enum class SolverKind : std::uint8_t {
  Jacobi,          // damped Jacobi sweeps (smoother / cheap scalar transport)
  PcgJacobi,       // SPD only
  BicgstabJacobi,  // nonsymmetric (momentum, turbulence)
  AmgPcg,          // SPD pressure: aggregation AMG preconditioned CG
  AmgBicgstab,     // nonsymmetric with AMG preconditioner
};

enum class NullSpace : std::uint8_t { None, Constant };
enum class OperatorStructure : std::uint8_t { General, SymmetricMMatrix };

// Krylov work vectors (residual, directions, best iterate...) are rebuilt by
// every solve. Solvers of one size whose calls never run concurrently may share
// one workspace instead of each holding about nine vectors. Preconditioner
// state (Jacobi diagonal, AMG hierarchy) is never shared.
class KrylovWorkspace {
 public:
  enum Slot { R, Z, P, Q, Shadow, S, T, Zs, Best, Scratch, Count };
  explicit KrylovWorkspace(Index rows) : rows_(rows) {}
  Index rows() const { return rows_; }
  // Allocated on first use; every user must pass the same Device.
  Buffer<float>& vector(Device& dev, Slot slot) {
    KUTTA_REQUIRE(!device_ || device_ == &dev, "Krylov workspace shared across devices");
    device_ = &dev;
    auto& v = vectors_[slot];
    if (!v.handle()) v = dev.alloc<float>(rows_);  // a handle even when empty
    return v;
  }
 private:
  Index rows_;
  Device* device_ = nullptr;
  Buffer<float> vectors_[Count];
};

struct SolverOptions {
  NullSpace nullspace = NullSpace::None;  // the equation, not rounded row sums, decides
  // Caller guarantees symmetry, nonpositive off-diagonals and stored row
  // diagonal dominance, including on refresh. AMG preserves this construction.
  OperatorStructure structure = OperatorStructure::General;
  std::shared_ptr<KrylovWorkspace> workspace;  // private when empty
};

class LinearSolver {
 public:
  virtual ~LinearSolver() = default;
  virtual void refresh(Device& dev, const CsrMatrix& a) = 0;
  virtual SolveStats solve(Device& dev, const CsrMatrix& a, Buffer<float>& x,
                           const Buffer<float>& b, const SolveControl& control) = 0;
  virtual std::string describe() const = 0;  // e.g. "AMG-PCG, 6 levels, op complexity 1.31"
};

// Builds structure (and for AMG the aggregation hierarchy) from the host matrix.
// Host values initialize diagonals/coarse operators and choose AMG aggregates.
std::unique_ptr<LinearSolver> make_solver(SolverKind kind, Device& dev, const CsrMatrixHost& a);

std::unique_ptr<LinearSolver> make_solver(SolverKind kind, Device& dev, const CsrMatrixHost& a,
                                           const SolverOptions& options);

}  // namespace kutta
