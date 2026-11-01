#pragma once

// Fluent-style pressure-based Coupled solve (Theory Guide 23.4.3.2): momentum
// and continuity in one block system per outer iteration, with implicit
// pressure gradient and implicit Rhie-Chow pressure dissipation. Turbulence
// stays segregated. The system is solved in delta form against the steady
// SIMPLEC residual, so both schemes share one discrete fixed point.
//
// Storage: the fine level is matrix-free (momentum off-diagonals come from the
// face mass flux and viscosity, pressure coupling from face areas and the
// pressure least-squares deltas); per cell it stores diagonals and a 3x3
// gradient-sensitivity tensor, plus one Rhie-Chow coefficient per face. Coarse
// levels are piecewise-constant aggregation Galerkin operators in a compact
// structured block CSR (8 floats per entry, 3 per row).

#include <memory>
#include <span>
#include <string>
#include <vector>

#include "kutta/fv/operators.hpp"
#include "kutta/linalg/solver.hpp"

namespace kutta::fv {

struct CoupledOptions {
  // The coarsest level (aggregation often stalls near 100 rows) gets parallel
  // block-Jacobi sweeps: a single-thread Gauss-Seidel there cost about half of
  // the Metal block solve for the same Krylov iteration counts.
  int pre_sweeps = 1, post_sweeps = 2, coarsest_sweeps = 8;
  float smoother_damping = .7f;
  // Fluent's flow Courant number for the Coupled scheme (default 200; Theory
  // Guide eq. 18.4-24 in Fluent 12): implicit relaxation a_P (1 + 1/CFL) of
  // the momentum Jacobian. It changes only the convergence path.
  double courant = 200;
  Index coarsest_rows = 48;
  int aggregation_passes = 3;
};

// Per outer iteration inputs. Residuals are b - A u of the assembled momentum
// rows at the current state; mass is the current Rhie-Chow face mass flux.
struct CoupledResiduals {
  const Buffer<float> &u, &v, &w;
  const Buffer<float>& mass;           // kg/s per face, out of owner
  const Buffer<float>& viscosity;      // face momentum diffusivity
  const Buffer<float>& physical_diagonal;
  // Cell weights of the face velocity interpolation (u_f = (a_P u_P +
  // a_N u_N)/(a_P + a_N)); the continuity-velocity block follows them.
  const Buffer<float>& face_weights;
  double density = 1, velocity_scale = 1;
  // Signed net outflow per cell to use instead of summing mass (the
  // float-pair continuity); null sums mass.
  const Buffer<float>* net_outflow = nullptr;
};

// The low words of the float-pair state (SolverSettings::paired_state) and
// kpair_add's slacks.
struct CoupledPairs {
  DeviceVec3& velocity_lo;
  Buffer<float>& pressure_lo;
  float velocity_slack = 0, pressure_slack = 0;
};

class CoupledSolver {
 public:
  // boundary_kind follows fv_flow.kern.h (2 and 4 are pressure-specified).
  // reference >= 0 anchors a closed domain's pressure level in that cell.
  // pressure_lsq is the residual's pressure reconstruction; its own-cell
  // sensitivity enters the Jacobian, except in cells flagged in
  // low_order_pressure (empty: none), whose face pressures the residual takes
  // from the cell value (SolverSettings::poor_cell_quality). All references
  // must outlive the solver.
  CoupledSolver(Device& dev, const DeviceMesh& mesh, const Stencil& stencil,
                const PolyMesh& host, const MeshGeometry& geometry,
                const Buffer<Index>& boundary_kind, const BoundaryLeastSquares& pressure_lsq,
                Index reference, std::span<const Index> low_order_pressure = {}, CoupledOptions options = {});
  ~CoupledSolver();
  CoupledSolver(const CoupledSolver&) = delete;
  CoupledSolver& operator=(const CoupledSolver&) = delete;

  // Captures one component's assembled momentum diagonal (pseudo time
  // included, no implicit relaxation): the physical (unrelaxed, no pseudo
  // time) part is floored at the off-diagonal sum and relaxed by the Courant
  // number. Off-diagonals are shared and never copied.
  void capture_diagonal(int component, const CsrMatrix& momentum, const Buffer<float>& physical);
  // Assembles the Jacobian and RHS, then solves for the update.
  SolveStats solve(const CoupledResiduals& residuals, const SolveControl& control);
  // u += alpha_u du, p += alpha_p dp (Fluent's explicit relaxation), then
  // re-gauges p at the reference cell of a closed domain.
  void update(DeviceVec3& velocity, Buffer<float>& pressure, double alpha_u, double alpha_p);
  // The same on the float pairs hi + lo.
  void update(DeviceVec3& velocity, Buffer<float>& pressure, const CoupledPairs& pairs, double alpha_u, double alpha_p);
  // Bytes of device memory this solver owns (excludes shared mesh/flow data).
  std::size_t bytes() const;
  std::string describe() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace kutta::fv
