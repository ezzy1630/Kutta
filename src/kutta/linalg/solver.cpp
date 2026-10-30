#include "kutta/linalg/detail.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <optional>

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
  Jacobi(Device& dev, const CsrMatrixHost& a, KrylovWorkspace& workspace)
      : inverse_(dev.upload<float>(diagonal_values(a, false, false))),
        invalid_(workspace.vector(dev, KrylovWorkspace::Scratch)) {}
  void refresh(Device& dev, const CsrMatrix& a) override {
    kern::krylov_diagonal(dev, a.n_rows, a.row_ptr, a.diag, a.val, inverse_, invalid_, 0);
    KUTTA_REQUIRE(dev.max_abs(invalid_, a.n_rows) == 0, "Jacobi: invalid values or zero diagonal");
  }
  void apply(Device& dev, const CsrMatrix&, const Buffer<float>& r, Buffer<float>& z,
             const Buffer<Index>& status) override {
    kern::krylov_multiply_guarded(dev, r.size(), r, inverse_, z, status);
  }
  std::string describe() const override { return "Jacobi"; }
 private:
  Buffer<float> inverse_;
  Buffer<float>& invalid_;  // refresh-only scratch
};

class Krylov final : public LinearSolver {
 public:
  Krylov(SolverKind kind, Device& dev, const CsrMatrixHost& a,
         std::unique_ptr<Preconditioner> p, NullSpace nullspace,
         std::shared_ptr<KrylovWorkspace> workspace)
      : kind_(kind), device_(&dev), n_(a.pattern.n_rows), nnz_(a.pattern.nnz()),
        nullspace_(nullspace == NullSpace::Constant), workspace_(std::move(workspace)),
        preconditioner_(std::move(p)), unused_(dev.alloc<float>(0)),
        r_(vector(KrylovWorkspace::R, true)), z_(vector(KrylovWorkspace::Z, true)),
        p_(vector(KrylovWorkspace::P, kind != SolverKind::Jacobi)),
        q_(vector(KrylovWorkspace::Q, kind != SolverKind::Jacobi)),
        shadow_(vector(KrylovWorkspace::Shadow, bicg(kind))), s_(vector(KrylovWorkspace::S, bicg(kind))),
        t_(vector(KrylovWorkspace::T, bicg(kind))), zs_(vector(KrylovWorkspace::Zs, bicg(kind))),
        best_(vector(KrylovWorkspace::Best, true)), scalars_(dev.alloc<float>(24)),
        status_(dev.alloc<Index>(3)), projection_(dev.alloc<float>(nullspace_ ? 2 : 0)),
        speculate_(kind != SolverKind::Jacobi && n_ <= kSpeculativeRows) {}

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
                      std::isfinite(c.abs_tol) && c.abs_tol >= 0 && c.residual_check_interval > 0,
                  "solve: invalid control");
    using Request = Device::HostReduction;
    speculated_.reset();
    // Entry checks share one host wait. Without a nullspace x is not changed
    // before the residual, so r_ (scratch until here) joins the same batch.
    if (!nullspace_) residual(dev, a, x, b, r_);
    const auto entry = nullspace_ ?
        dev.reduce_many({Request::dot(b, b, b.size()), Request::dot(x, x, x.size()),
                         Request::sum(b, n_)}) :
        dev.reduce_many({Request::dot(b, b, b.size()), Request::dot(x, x, x.size()),
                         Request::dot(r_, r_, r_.size())});
    const double bnorm = std::sqrt(entry[0]);
    KUTTA_REQUIRE(std::isfinite(bnorm) && std::isfinite(std::sqrt(entry[1])),
                  "solve: nonfinite input vector");
    SolveStats stats;
    if (nullspace_) {
      if (n_ != 0) {
        KUTTA_REQUIRE(std::abs(entry[2]) <= 1e-6 * std::sqrt(double(n_)) * bnorm,
                      "solve: incompatible RHS for constant Neumann nullspace");
        project(dev, x, projection_);
      }
      residual(dev, a, x, b, r_);
      stats.initial_residual = norm(dev, r_);
    } else stats.initial_residual = std::sqrt(entry[2]);
    stats.final_residual = stats.initial_residual;
    verified_current_ = true;
    const double target = std::max(c.abs_tol, c.rel_tol * stats.initial_residual);
    if (!std::isfinite(stats.initial_residual)) return stats;
    if (stats.final_residual <= target) { stats.converged = true; return stats; }
    if (nullspace_) project(dev, r_, projection_);
    best_residual_ = stats.initial_residual;
    dev.copy(best_, x);
    if (kind_ != SolverKind::Jacobi) {
      // Thresholds use the same double host calculation as the public norm API.
      // Upload only at the solve boundary; iteration state is reset by kernels.
      std::array<float, 24> scalars{};
      auto pair = [&](int offset, double value) {
        scalars[offset] = static_cast<float>(value);
        scalars[offset + 1] = std::isfinite(scalars[offset]) ?
            static_cast<float>(value - double(scalars[offset])) : 0;
      };
      pair(20, target * target);
      pair(22, 1e6 * stats.initial_residual * stats.initial_residual);
      kern::solver_limits(dev, scalars.size(), scalars_, scalars[20], scalars[21],
                          scalars[22], scalars[23]);
    }
    reset(dev, 0);
    if (kind_ == SolverKind::Jacobi) {
      for (int i = 0; i < c.max_iter; ++i) {
        apply(dev, a, r_, z_);
        kern::axpy(dev, n_, x, z_, 2.0f / 3.0f);
        verified_current_ = false;
        ++stats.iterations;
        if (stats.iterations % c.residual_check_interval == 0 || stats.iterations == c.max_iter) {
          if (verified(dev, a, x, b, target, stats) || diverged(stats)) break;
        } else residual(dev, a, x, b, r_);
      }
    } else if (kind_ == SolverKind::PcgJacobi || kind_ == SolverKind::AmgPcg) {
      pcg(dev, a, x, b, c, target, stats);
    } else {
      bicgstab(dev, a, x, b, c, target, stats);
    }
    // All exits, including breakdown/max_iter, report the actual residual.
    // Candidate exits have just checked this exact x. Constant projection can
    // round x again, so retain its historical final verification sequence.
    if (!verified_current_ || nullspace_) verified(dev, a, x, b, target, stats);
    if (!std::isfinite(stats.final_residual) || stats.final_residual > best_residual_) {
      dev.copy(x, best_);
      verified(dev, a, x, b, target, stats);
    }
    previous_iterations_ = stats.converged ? stats.iterations : 0;
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
  Buffer<float>& vector(KrylovWorkspace::Slot slot, bool used) {
    return used ? workspace_->vector(*device_, slot) : unused_;
  }
  static bool bicg(SolverKind kind) {
    return kind == SolverKind::BicgstabJacobi || kind == SolverKind::AmgBicgstab;
  }
  void check_matrix(Device& dev, const CsrMatrix& a) const {
    KUTTA_REQUIRE(&dev == device_ && a.n_rows == n_ && a.nnz == nnz_ &&
                      a.val.size() == std::size_t(nnz_), "solver: device or matrix mismatch");
  }
  void apply(Device& dev, const CsrMatrix& a, const Buffer<float>& r, Buffer<float>& z) {
    preconditioner_->apply(dev, a, r, z, status_);
    if (nullspace_) project(dev, z, projection_);
  }
  bool verified(Device& dev, const CsrMatrix& a, Buffer<float>& x, const Buffer<float>& b,
                double target, SolveStats& stats) {
    if (speculated_) {
      // The same kernels already ran on copies of this exact x; adopt them.
      if (nullspace_) dev.copy(x, speculative_x());
      dev.copy(r_, speculative_r());
      stats.final_residual = *speculated_;
      speculated_.reset();
    } else {
      if (nullspace_) project(dev, x, projection_);
      residual(dev, a, x, b, r_);
      stats.final_residual = norm(dev, r_);
    }
    verified_current_ = true;
    stats.converged = std::isfinite(stats.final_residual) && stats.final_residual <= target;
    if (stats.final_residual < best_residual_) {
      best_residual_ = stats.final_residual;
      dev.copy(best_, x);
    }
    if (nullspace_) project(dev, r_, projection_);
    return stats.converged;
  }
  bool diverged(const SolveStats& stats) const {
    return !std::isfinite(stats.final_residual) ||
           stats.final_residual > 1e3 * stats.initial_residual;
  }
  // Reads the batch status. On small systems the reliable residual of the
  // current x is computed alongside, so a candidate exit costs one host wait
  // instead of two. Copies keep x and r_ untouched until verified() adopts it.
  std::vector<Index> poll(Device& dev, const CsrMatrix& a, const Buffer<float>& x,
                          const Buffer<float>& b) {
    if (speculate_) {
      const Buffer<float>* source = &x;
      if (nullspace_) {
        dev.copy(speculative_x(), x);
        project(dev, speculative_x(), projection_);
        source = &speculative_x();
      }
      auto& r = speculative_r();
      residual(dev, a, *source, b, r);
      speculated_ = std::sqrt(dev.reduce_many({Device::HostReduction::dot(r, r, r.size())})[0]);
    }
    return dev.download(status_);
  }
  // Work vectors that no later step reads before rewriting once a batch has
  // been queued: BiCGStab recomputes s and t (and PCG q and z) from the
  // carried r, p, v and x. Batch reductions may still read them, but only
  // guarded kernels consume those scalars. Holding the speculative residual
  // there keeps Krylov storage unchanged.
  Buffer<float>& speculative_r() { return bicg(kind_) ? t_ : q_; }
  Buffer<float>& speculative_x() { return bicg(kind_) ? s_ : z_; }
  void reset(Device& dev, int iterations) {
    kern::solver_reset(dev, 1, scalars_, status_, iterations);
  }
  // A candidate pauses device updates until the full b-A*x is checked. Batches
  // may contain no-op updates after that point; they never inflate iterations.
  bool replacement(Device& dev, const CsrMatrix& a, Buffer<float>& x,
                   const Buffer<float>& b, double target, SolveStats& stats,
                   int& restarts, double& last_verified) {
    if (verified(dev, a, x, b, target, stats) || diverged(stats)) return false;
    restarts = stats.final_residual >= 0.9 * last_verified ? restarts + 1 : 0;
    last_verified = stats.final_residual;
    return restarts < 12;
  }
  void pcg(Device& dev, const CsrMatrix& a, Buffer<float>& x, const Buffer<float>& b,
           const SolveControl& c, double target, SolveStats& stats) {
    auto restart = [&] {
      reset(dev, stats.iterations);
      apply(dev, a, r_, z_);
      dev.copy(p_, z_);
      dev.dot_into(r_, z_, n_, scalars_, 0);
    };
    restart();
    int restarts = 0;
    double last_verified = stats.initial_residual;
    int check_interval = batch_interval(c);
    while (stats.iterations < c.max_iter) {
      const int batch = std::min(check_interval, c.max_iter - stats.iterations);
      verified_current_ = false;
      speculated_.reset();
      for (int j = 0; j < batch; ++j) {
        kern::csr_spmv_guarded(dev, n_, a.row_ptr, a.col, a.val, p_, q_, status_);
        dev.dot_into(p_, q_, n_, scalars_, 2);
        kern::solver_pcg_alpha(dev, 1, scalars_, status_);
        kern::solver_step(dev, n_, x, r_, p_, q_, scalars_, status_);
        dev.dot_into(r_, r_, n_, scalars_, 18);
        kern::solver_check(dev, 1, scalars_, status_, 1);
        if (stats.iterations + j + 1 == c.max_iter) continue;
        apply(dev, a, r_, z_);
        if (preconditioner_->variable()) dev.dot2_into(r_, z_, z_, q_, n_, scalars_, 4);
        else dev.dot_into(r_, z_, n_, scalars_, 4);
        kern::solver_pcg_beta(dev, 1, scalars_, status_, preconditioner_->variable() ? 1 : 0);
        kern::solver_direction(dev, n_, z_, p_, scalars_, status_);
      }
      const auto state = poll(dev, a, x, b);
      stats.iterations = state[1];
      if (state[0] == 2) break;
      if (state[0] == 1) {
        if (!replacement(dev, a, x, b, target, stats, restarts, last_verified) ||
            stats.iterations == c.max_iter) break;
        // Near a float residual floor, repeated candidates can stop after one
        // update. Avoid queueing a full batch behind each reliable replacement.
        if (stats.final_residual <= 10 * target) check_interval = 1;
        restart();
      }
    }
  }
  void bicgstab(Device& dev, const CsrMatrix& a, Buffer<float>& x, const Buffer<float>& b,
                const SolveControl& c, double target, SolveStats& stats) {
    auto restart = [&] {
      reset(dev, stats.iterations);
      dev.copy(shadow_, r_);
      dev.fill(p_, 0.0f);
      dev.fill(q_, 0.0f);
    };
    restart();
    int restarts = 0;
    double last_verified = stats.initial_residual;
    int check_interval = batch_interval(c);
    while (stats.iterations < c.max_iter) {
      const int batch = std::min(check_interval, c.max_iter - stats.iterations);
      verified_current_ = false;
      speculated_.reset();
      for (int j = 0; j < batch; ++j) {
        dev.dot_into(shadow_, r_, n_, scalars_, 0);
        kern::solver_bicg_beta(dev, 1, scalars_, status_);
        kern::solver_bicg_direction(dev, n_, r_, q_, p_, scalars_, status_);
        apply(dev, a, p_, z_);
        kern::csr_spmv_guarded(dev, n_, a.row_ptr, a.col, a.val, z_, q_, status_);
        dev.dot_into(shadow_, q_, n_, scalars_, 2);
        kern::solver_bicg_alpha(dev, 1, scalars_, status_);
        kern::solver_bicg_difference(dev, n_, r_, q_, s_, scalars_, status_);
        dev.dot_into(s_, s_, n_, scalars_, 16);
        kern::solver_bicg_s(dev, 1, scalars_, status_);
        apply(dev, a, s_, zs_);
        kern::csr_spmv_guarded(dev, n_, a.row_ptr, a.col, a.val, zs_, t_, status_);
        dev.dot2_into(t_, t_, t_, s_, n_, scalars_, 2);
        kern::solver_bicg_omega(dev, 1, scalars_, status_);
        kern::solver_bicg_step(dev, n_, x, r_, z_, zs_, s_, t_, scalars_, status_);
        dev.dot_into(r_, r_, n_, scalars_, 18);
        kern::solver_check(dev, 1, scalars_, status_, 0);
      }
      const auto state = poll(dev, a, x, b);
      stats.iterations = state[1];
      if (state[0] == 2) break;
      if (state[0] == 1) {
        if (!replacement(dev, a, x, b, target, stats, restarts, last_verified) ||
            stats.iterations == c.max_iter) break;
        // Near a float residual floor, repeated candidates can stop after one
        // update. Avoid queueing a full batch behind each reliable replacement.
        if (stats.final_residual <= 10 * target) check_interval = 1;
        restart();
      }
    }
  }

  int batch_interval(const SolveControl& c) const {
    // The caller's interval bounds updates between host polls. A reused solver
    // often converges in fewer updates: smaller batches avoid frozen AMG cycles
    // and reductions without changing device convergence tests or arithmetic.
    return previous_iterations_ > 0 && previous_iterations_ < c.residual_check_interval ?
        previous_iterations_ + 1 : c.residual_check_interval;
  }

  SolverKind kind_;
  Device* device_;
  Index n_, nnz_;
  bool nullspace_;
  // Every solve rebuilds these before reading them, so they may be shared.
  // Declared first: a Jacobi preconditioner refers to its scratch slot.
  std::shared_ptr<KrylovWorkspace> workspace_;
  std::unique_ptr<Preconditioner> preconditioner_;
  Buffer<float> unused_;
  Buffer<float> &r_, &z_, &p_, &q_, &shadow_, &s_, &t_, &zs_, &best_;
  Buffer<float> scalars_;
  Buffer<Index> status_;
  Buffer<float> projection_;
  // Polling with a speculative residual costs one extra residual and dot per
  // batch; worthwhile only where a host wait costs more than that work.
  static constexpr Index kSpeculativeRows = 1 << 17;
  bool speculate_;
  std::optional<double> speculated_;  // norm of speculative_r() for the current x
  double best_residual_ = 0;
  int previous_iterations_ = 0;
  bool verified_current_ = false;
};
}  // namespace

std::unique_ptr<LinearSolver> krylov_solver(SolverKind kind, Device& dev,
    const CsrMatrixHost& a, std::unique_ptr<Preconditioner> p, NullSpace nullspace,
    std::shared_ptr<KrylovWorkspace> workspace) {
  if (!workspace) workspace = std::make_shared<KrylovWorkspace>(a.pattern.n_rows);
  KUTTA_REQUIRE(workspace->rows() == a.pattern.n_rows, "Krylov workspace size mismatch");
  if (!p) p = std::make_unique<Jacobi>(dev, a, *workspace);
  auto solver = std::make_unique<Krylov>(kind, dev, a, std::move(p), nullspace, std::move(workspace));
  return solver;
}
}  // namespace kutta::linalg_detail
