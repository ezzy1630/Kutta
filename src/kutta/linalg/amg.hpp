#pragma once

#include "kutta/linalg/solver.hpp"

namespace kutta {
enum class AmgCycle { V, K };

struct AmgOptions {
  AmgCycle cycle = AmgCycle::K;
  int pre_sweeps = 2;
  int post_sweeps = 2;
  int coarse_sweeps = 32;
  Index coarse_rows = 24;
  int max_levels = 24;
  int aggregation_passes = 3;
  double strength_threshold = 0.25;
};

// Same lifecycle as make_solver; only AmgPcg/AmgBicgstab are accepted.
// K uses two recursive inner Krylov corrections (flexible CG for SPD,
// GCR for nonsymmetric). V is a fixed linear symmetric preconditioner
// when pre_sweeps == post_sweeps and the matrix is SPD.
std::unique_ptr<LinearSolver> make_amg_solver(SolverKind kind, Device& dev,
    const CsrMatrixHost& a, const AmgOptions& options);
}  // namespace kutta
