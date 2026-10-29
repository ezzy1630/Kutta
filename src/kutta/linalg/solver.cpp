#include "kutta/linalg/detail.hpp"

#include <algorithm>
#include <limits>

namespace kutta::linalg_detail {
void check_nullspace(const CsrMatrixHost& a, NullSpace choice) {
  KUTTA_REQUIRE(choice == NullSpace::None || choice == NullSpace::Constant,
                "solver: invalid nullspace option");
#ifndef NDEBUG
  if (choice != NullSpace::Constant) return;
  // Diagnostic only: float assembly cannot preserve exact row sums. Never use
  // this tolerance to decide whether an invertible equation should be projected.
  for (Index i = 0; i < a.pattern.n_rows; ++i) {
    double sum = 0, magnitude = 0;
    for (Index j = a.pattern.row_ptr[i]; j < a.pattern.row_ptr[i + 1]; ++j) {
      const double v = static_cast<float>(a.val[j]);
      sum += v; magnitude += std::abs(v);
      const auto col = a.pattern.col[j];
      const auto first = a.pattern.col.begin() + a.pattern.row_ptr[col];
      const auto last = a.pattern.col.begin() + a.pattern.row_ptr[col + 1];
      const auto opposite = std::lower_bound(first, last, i);
      KUTTA_REQUIRE(opposite != last && *opposite == i &&
          static_cast<float>(a.val[opposite - a.pattern.col.begin()]) == v,
          "constant nullspace: expected symmetric pressure operator");
    }
    KUTTA_REQUIRE(std::abs(sum) <= 8 * std::numeric_limits<float>::epsilon() * magnitude,
                  "constant nullspace: inconsistent pressure row sum");
  }
#else
  (void)a;
#endif
}

std::vector<float> diagonal_values(const CsrMatrixHost& a, bool l1, bool nullspace) {
  std::vector<float> inverse(a.pattern.n_rows);
  for (Index i = 0; i < a.pattern.n_rows; ++i) {
    const float d = static_cast<float>(a.val[a.pattern.diag[i]]);
    KUTTA_REQUIRE(!l1 || d >= 0, "AMG: requires positive diagonals");
    KUTTA_REQUIRE(d != 0 || (nullspace && a.pattern.n_rows == 1),
                  "solver: invalid values or zero diagonal");
    float off = 0;
    if (l1) for (Index j = a.pattern.row_ptr[i]; j < a.pattern.row_ptr[i + 1]; ++j)
      if (j != a.pattern.diag[i]) off += std::abs(static_cast<float>(a.val[j]));
    const float scale = l1 ? std::max(std::abs(d), 0.5f * (std::abs(d) + off)) : d;
    inverse[i] = scale != 0 ? 1.0f / scale : 0;
    KUTTA_REQUIRE(std::isfinite(scale) && std::isfinite(inverse[i]),
                  "solver: invalid values or diagonal inverse");
  }
  return inverse;
}

namespace {
class Jacobi final : public Preconditioner {
 public:
  Jacobi(Device& dev, const CsrMatrixHost& a)
      : inverse_(dev.upload<float>(diagonal_values(a, false, false))) {}
  void refresh(Device& dev, const CsrMatrix& a) override {
    auto invalid = dev.alloc<float>(a.n_rows);
    kern::krylov_diagonal(dev, a.n_rows, a.row_ptr, a.diag, a.val, inverse_, invalid, 0);
    KUTTA_REQUIRE(dev.max_abs(invalid, a.n_rows) == 0, "Jacobi: invalid values or zero diagonal");
  }
  void apply(Device& dev, const CsrMatrix&, const Buffer<float>& r, Buffer<float>& z) override {
    kern::krylov_multiply(dev, r.size(), r, inverse_, z);
  }
  std::string describe() const override { return "Jacobi"; }
 private:
  Buffer<float> inverse_;
};

class Krylov final : public LinearSolver {
 public:
  Krylov(SolverKind kind, Device& dev, const CsrMatrixHost& a,
         std::unique_ptr<Preconditioner> p, NullSpace nullspace)
      : kind_(kind), device_(&dev), n_(a.pattern.n_rows), nnz_(a.pattern.nnz()),
        nullspace_(nullspace == NullSpace::Constant), preconditioner_(std::move(p)),
        r_(dev.alloc<float>(n_)), z_(dev.alloc<float>(n_)), p_(dev.alloc<float>(kind == SolverKind::Jacobi ? 0 : n_)),
        q_(dev.alloc<float>(kind == SolverKind::Jacobi ? 0 : n_)),
        shadow_(dev.alloc<float>(bicg(kind) ? n_ : 0)), s_(dev.alloc<float>(bicg(kind) ? n_ : 0)),
        t_(dev.alloc<float>(bicg(kind) ? n_ : 0)), zs_(dev.alloc<float>(bicg(kind) ? n_ : 0)),
        best_(dev.alloc<float>(n_)) {}

  void refresh(Device& dev, const CsrMatrix& a) override {
    check_matrix(dev, a);
    preconditioner_->refresh(dev, a);
  }

  SolveStats solve(Device& dev, const CsrMatrix& a, Buffer<float>& x,
                   const Buffer<float>& b, const SolveControl& c) override {
    check_matrix(dev, a);
    KUTTA_REQUIRE(x.size() == std::size_t(n_) && b.size() == x.size() &&
                      x.handle() != b.handle(), "solve: vector size or alias mismatch");
    KUTTA_REQUIRE(c.max_iter >= 0 && std::isfinite(c.rel_tol) && c.rel_tol >= 0 &&
                      std::isfinite(c.abs_tol) && c.abs_tol >= 0 && c.residual_check_interval > 0, "solve: invalid control");
    KUTTA_REQUIRE(std::isfinite(norm(dev, b)) && std::isfinite(norm(dev, x)),
                  "solve: nonfinite input vector");
    if (nullspace_ && n_ != 0) {
      KUTTA_REQUIRE(std::abs(dev.sum(b, n_)) <=
                        1e-6 * std::sqrt(double(n_)) * norm(dev, b),
                    "solve: incompatible RHS for constant Neumann nullspace");
      project(dev, x);
    }
    SolveStats stats;
    residual(dev, a, x, b, r_);
    stats.initial_residual = norm(dev, r_);
    stats.final_residual = stats.initial_residual;
    const double target = std::max(c.abs_tol, c.rel_tol * stats.initial_residual);
    if (stats.final_residual <= target) { stats.converged = true; return stats; }
    if (nullspace_) project(dev, r_);
    best_residual_ = stats.initial_residual;
    dev.copy(best_, x);
    if (kind_ == SolverKind::Jacobi) {
      for (int i = 0; i < c.max_iter; ++i) {
        apply(dev, a, r_, z_);
        kern::axpy(dev, n_, x, z_, 2.0f / 3.0f);
        ++stats.iterations;
        residual(dev, a, x, b, r_);
        if (stats.iterations % c.residual_check_interval == 0 || stats.iterations == c.max_iter) {
          if (verified(dev, a, x, b, target, stats) || diverged(stats)) break;
        }
      }
    } else if (kind_ == SolverKind::PcgJacobi || kind_ == SolverKind::AmgPcg) {
      pcg(dev, a, x, b, c, target, stats);
    } else {
      bicgstab(dev, a, x, b, c, target, stats);
    }
    // All exits, including breakdown/max_iter, report the actual residual.
    verified(dev, a, x, b, target, stats);
    if (!std::isfinite(stats.final_residual) || stats.final_residual > best_residual_) {
      dev.copy(x, best_);
      verified(dev, a, x, b, target, stats);
    }
    return stats;
  }

  std::string describe() const override {
    const char* name = kind_ == SolverKind::Jacobi ? "damped Jacobi" :
        (kind_ == SolverKind::PcgJacobi || kind_ == SolverKind::AmgPcg ?
          (preconditioner_->variable() ? "flexible PCG" : "PCG") : "BiCGStab");
    return std::string(name) + ", " + preconditioner_->describe() +
        (nullspace_ ? ", constant-nullspace projection" : "");
  }

 private:
  static bool bicg(SolverKind kind) {
    return kind == SolverKind::BicgstabJacobi || kind == SolverKind::AmgBicgstab;
  }
  void check_matrix(Device& dev, const CsrMatrix& a) const {
    KUTTA_REQUIRE(&dev == device_ && a.n_rows == n_ && a.nnz == nnz_ &&
                      a.val.size() == std::size_t(nnz_), "solver: device or matrix mismatch");
  }
  void apply(Device& dev, const CsrMatrix& a, const Buffer<float>& r, Buffer<float>& z) {
    preconditioner_->apply(dev, a, r, z);
    if (nullspace_) project(dev, z);
  }
  bool verified(Device& dev, const CsrMatrix& a, Buffer<float>& x, const Buffer<float>& b,
                double target, SolveStats& stats) {
    if (nullspace_) project(dev, x);
    residual(dev, a, x, b, r_);
    stats.final_residual = norm(dev, r_);
    stats.converged = std::isfinite(stats.final_residual) && stats.final_residual <= target;
    if (stats.final_residual < best_residual_) {
      best_residual_ = stats.final_residual;
      dev.copy(best_, x);
    }
    if (nullspace_) project(dev, r_);
    return stats.converged;
  }
  bool diverged(const SolveStats& stats) const {
    return !std::isfinite(stats.final_residual) ||
           stats.final_residual > 1e3 * stats.initial_residual;
  }
  void pcg(Device& dev, const CsrMatrix& a, Buffer<float>& x, const Buffer<float>& b,
           const SolveControl& c, double target, SolveStats& stats) {
    apply(dev, a, r_, z_);
    dev.copy(p_, z_);
    double rho = dev.dot(r_, z_, n_);
    int restarts = 0;
    double last_verified = stats.initial_residual;
    while (stats.iterations < c.max_iter) {
      spmv(dev, a, p_, q_);
      const double pq = dev.dot(p_, q_, n_);
      if (!(rho > 0 && pq > 0) || (!coefficient(rho / pq) || rho / pq == 0)) break;
      const double alpha = rho / pq;
      kern::krylov_step(dev, n_, x, r_, p_, q_, static_cast<float>(alpha));
      ++stats.iterations;
      const double rn = norm(dev, r_);
      if (rn <= target || !std::isfinite(rn) || rn > 1e3 * stats.initial_residual ||
          stats.iterations % 50 == 0) {
        if (verified(dev, a, x, b, target, stats) || diverged(stats)) break;
        // Reliable residual replacement restarts the recurrence. Twelve
        // stagnant replacements indicate a float precision floor.
        restarts = stats.final_residual >= 0.9 * last_verified ? restarts + 1 : 0;
        last_verified = stats.final_residual;
        if (restarts >= 12) break;
        apply(dev, a, r_, z_);
        rho = dev.dot(r_, z_, n_);
        dev.copy(p_, z_);
        continue;
      }
      apply(dev, a, r_, z_);
      const double next = dev.dot(r_, z_, n_);
      // One-term flexible CG for the nonlinear K-cycle. Fixed V/Jacobi
      // preconditioners use the ordinary Fletcher-Reeves recurrence.
      const double beta = preconditioner_->variable() ?
          -alpha * dev.dot(z_, q_, n_) / rho : next / rho;
      if (!(next > 0) || !std::isfinite(beta) || std::abs(beta) >
          std::numeric_limits<float>::max()) break;
      kern::krylov_direction(dev, n_, z_, p_, static_cast<float>(beta));
      rho = next;
    }
  }
  void bicgstab(Device& dev, const CsrMatrix& a, Buffer<float>& x, const Buffer<float>& b,
                const SolveControl& c, double target, SolveStats& stats) {
    dev.copy(shadow_, r_);
    dev.fill(p_, 0.0f);
    dev.fill(q_, 0.0f);
    double old_rho = 1, alpha = 1, omega = 1;
    int restarts = 0;
    double last_verified = stats.initial_residual;
    while (stats.iterations < c.max_iter) {
      const double rho = dev.dot(shadow_, r_, n_);
      if (rho == 0 || !std::isfinite(rho) || omega == 0 || !coefficient(omega)) break;
      const double beta = (rho / old_rho) * (alpha / omega);
      if (!coefficient(beta)) break;
      kern::krylov_bicg_direction(dev, n_, r_, q_, p_, static_cast<float>(beta),
                                static_cast<float>(omega));
      apply(dev, a, p_, z_);
      spmv(dev, a, z_, q_);
      const double rv = dev.dot(shadow_, q_, n_);
      if (rv == 0 || !std::isfinite(rv) || !coefficient(rho / rv)) break;
      alpha = rho / rv;
      kern::krylov_difference(dev, n_, r_, q_, s_, static_cast<float>(alpha));
      ++stats.iterations;
      const double sn = norm(dev, s_);
      if (sn <= target) {
        kern::krylov_bicg_step(dev, n_, x, z_, z_, static_cast<float>(alpha), 0.0f);
      } else {
        apply(dev, a, s_, zs_);
        spmv(dev, a, zs_, t_);
        const double tt = dev.dot(t_, t_, n_);
        if (!(tt > 0) || !std::isfinite(tt)) break;
        omega = dev.dot(t_, s_, n_) / tt;
        if (omega == 0 || !coefficient(omega)) break;
        kern::krylov_bicg_step(dev, n_, x, z_, zs_, static_cast<float>(alpha),
                             static_cast<float>(omega));
        kern::krylov_difference(dev, n_, s_, t_, r_, static_cast<float>(omega));
      }
      const double rn = sn <= target ? sn : norm(dev, r_);
      if (rn <= target || !std::isfinite(rn) || rn > 1e3 * stats.initial_residual ||
          stats.iterations % 50 == 0) {
        if (verified(dev, a, x, b, target, stats) || diverged(stats)) break;
        restarts = stats.final_residual >= 0.9 * last_verified ? restarts + 1 : 0;
        last_verified = stats.final_residual;
        if (restarts >= 12) break;
        dev.copy(shadow_, r_);
        dev.fill(p_, 0.0f);
        dev.fill(q_, 0.0f);
        old_rho = alpha = omega = 1;
      } else old_rho = rho;
    }
  }

  SolverKind kind_;
  Device* device_;
  Index n_, nnz_;
  bool nullspace_;
  std::unique_ptr<Preconditioner> preconditioner_;
  Buffer<float> r_, z_, p_, q_, shadow_, s_, t_, zs_, best_;
  double best_residual_ = 0;
};
}  // namespace

std::unique_ptr<LinearSolver> krylov_solver(SolverKind kind, Device& dev,
    const CsrMatrixHost& a, std::unique_ptr<Preconditioner> p, NullSpace nullspace) {
  if (!p) p = std::make_unique<Jacobi>(dev, a);
  auto solver = std::make_unique<Krylov>(kind, dev, a, std::move(p), nullspace);
  return solver;
}
}  // namespace kutta::linalg_detail
