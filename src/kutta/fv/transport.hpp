#pragma once

#include "kutta/fv/operators.hpp"
#include "kutta/linalg/solver.hpp"

namespace kutta::fv {

enum class ConvectionScheme { FirstOrderUpwind, SecondOrderUpwind };
enum class ScalarBoundary : Index { FixedValue=0, ZeroGradient=1 };

struct TransportSettings {
  ConvectionScheme convection=ConvectionScheme::SecondOrderUpwind;
  bool limit=true;
  float density=1;
  float pseudo_time_step=0;  // zero disables pseudo time
  float relaxation=1;
  SolveControl linear{1e-4,1e-12,200};
};

// Solver policy is kept here so the parallel linalg lane can change it once.
inline constexpr SolverKind transport_solver_kind=SolverKind::BicgstabJacobi;
inline constexpr SolverKind pressure_solver_kind=SolverKind::AmgPcg;

class TransportEquation {
  std::shared_ptr<CsrMatrix> matrix_storage_;
 public:
  TransportEquation(Device& dev, const DeviceMesh& mesh, const Stencil& stencil,
                    const PolyMesh& host, const MeshGeometry& geometry,
                    std::shared_ptr<CsrMatrix> shared_matrix={});
  void assemble(const Buffer<float>& field, const Buffer<float>& face_values,
                const Buffer<Index>& boundary_kind, const Buffer<float>& mass_flux,
                const Buffer<float>& diffusivity, const Buffer<float>& source,
                const TransportSettings& settings);
  // grad must already contain the LSQ reconstruction of field/face_values.
  void assemble_reconstructed(const Buffer<float>& field,const Buffer<float>& face_values,
                const Buffer<Index>& boundary_kind,const Buffer<float>& mass_flux,
                const Buffer<float>& diffusivity,const Buffer<float>& source,
                const TransportSettings& settings);
  SolveStats solve(Buffer<float>& field, const SolveControl& control);
  CsrMatrix& matrix;
  Buffer<float> rhs, diagonal, consistent_inverse;
  DeviceVec3 grad;
  Buffer<float> limiter, face_correction;
 private:
  Device& dev_;
  const DeviceMesh& mesh_;
  const Stencil& stencil_;
  std::unique_ptr<LinearSolver> solver_;
};

}  // namespace kutta::fv
