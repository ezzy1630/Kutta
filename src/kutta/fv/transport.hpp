#pragma once

#include <optional>

#include "kutta/fv/operators.hpp"
#include "kutta/fv/face_projection.hpp"
#include "kutta/linalg/solver.hpp"

namespace kutta::fv {

enum class ConvectionScheme { FirstOrderUpwind, SecondOrderUpwind };
enum class ScalarBoundary : Index { FixedValue=0, ZeroGradient=1 };
enum class TransportPreconditioner { Jacobi, Amg };

struct TransportSettings {
  ConvectionScheme convection=ConvectionScheme::SecondOrderUpwind;
  bool limit=true;
  float density=1;
  float pseudo_time_step=0;  // zero disables pseudo time
  float relaxation=1;
  // The low words of a float-pair field (fv_paired_limiter); null for float.
  const Buffer<float>* limiter_low=nullptr;
  std::optional<PairedGradient> limiter_paired_gradient;
  // fv_limiter's eps2: the differentiable limiter's squared threshold
  // (SolverSettings::limiter); negative is Barth-Jespersen.
  float limiter_epsilon2=-1;
  SolveControl linear{1e-4,1e-12,200};
};

struct TransportTerms {
  const Buffer<float>* sink=nullptr;      // per volume, implicit sink*phi
  const Buffer<float>* diagonal=nullptr;  // integrated per-cell contribution
  const Buffer<float>* rhs=nullptr;       // integrated per-cell contribution
};

// Solver policy is kept here so the parallel linalg lane can change it once.
inline constexpr SolverKind transport_solver_kind=SolverKind::BicgstabJacobi;
inline constexpr SolverKind pressure_solver_kind=SolverKind::AmgPcg;

// Equation assemblies and physical residuals must execute sequentially on
// one Device. RHS, diagonals, gradients and limiter fields remain per equation.
struct TransportWorkspace {
  TransportWorkspace(Device& dev,Index cells,Index faces)
      :device(&dev),zero(dev.alloc<float>(cells)),face_correction(dev.alloc<float>(faces)) {dev.fill(zero,0);}
  Device* device;
  Buffer<float> zero,face_correction;
};

class TransportEquation {
  std::shared_ptr<CsrMatrix> matrix_storage_;
  std::shared_ptr<TransportWorkspace> scratch_;
 public:
  TransportEquation(Device& dev, const DeviceMesh& mesh, const Stencil& stencil,
                    const PolyMesh& host, const MeshGeometry& geometry,
                    std::shared_ptr<CsrMatrix> shared_matrix={},
                    std::shared_ptr<const FaceReconstruction> diffusion_reconstruction={},
                    std::shared_ptr<KrylovWorkspace> workspace={},
                    std::shared_ptr<TransportWorkspace> scratch={});
  void assemble(const Buffer<float>& field, const Buffer<float>& face_values,
                const Buffer<Index>& boundary_kind, const Buffer<float>& mass_flux,
                const Buffer<float>& diffusivity, const Buffer<float>& source,
                const TransportSettings& settings,TransportTerms terms={});
  // grad must already contain the LSQ reconstruction of field/face_values.
  void assemble_reconstructed(const Buffer<float>& field,const Buffer<float>& face_values,
                const Buffer<Index>& boundary_kind,const Buffer<float>& mass_flux,
                const Buffer<float>& diffusivity,const Buffer<float>& source,
                const TransportSettings& settings,TransportTerms terms={});
  // Null restores the compact secondary-gradient diffusion.
  void set_reconstruction(std::shared_ptr<const FaceReconstruction> reconstruction);
  // Rebuild at the next solve from its assembled operator, then refresh values
  // on later solves. Both choices use nonsymmetric BiCGStab.
  void set_preconditioner(TransportPreconditioner preconditioner);
  // Call before assembling a paired momentum equation. No iteration-time
  // allocation; ordinary scalar solves retain their existing float path.
  void enable_compensated_assembly();
  // Solve a scalar increment from its compensated physical defect.
  void enable_compensated_solve();
  void physical_residual(const Buffer<float>& field,Buffer<float>& result,const Buffer<Index>* excluded=nullptr);
  // An inner equation already below this physical steady defect need not
  // chase float roundoff. A changed source is checked again at every solve.
  void set_steady_tolerance(double tolerance);
  const Buffer<float>& rhs_low() const {return rhs_low_.empty()?zero_:rhs_low_;}
  const Buffer<float>& diagonal_low() const {return diagonal_low_.empty()?zero_:diagonal_low_;}
  SolveStats solve(Buffer<float>& field, const SolveControl& control,const Buffer<Index>* excluded=nullptr);
  CsrMatrix& matrix;
  Buffer<float> rhs, diagonal, consistent_inverse;
  DeviceVec3 grad;
  Buffer<float> limiter;
  Buffer<float>& face_correction;
 private:
  Device& dev_;
  const DeviceMesh& mesh_;
  const Stencil& stencil_;
  std::unique_ptr<LinearSolver> solver_;
  std::shared_ptr<KrylovWorkspace> workspace_;
  std::optional<TransportPreconditioner> pending_preconditioner_;
  Buffer<float>& zero_;
  Buffer<float> rhs_low_,diagonal_low_,scalar_delta_,scalar_residual_;
  struct SteadyStopping {double tolerance;Buffer<float> residual,absolute,scale;};
  std::optional<SteadyStopping> steady_stopping_;
  std::shared_ptr<const FaceReconstruction> diffusion_reconstruction_;
  Buffer<float> diffusion_correction_;
};

}  // namespace kutta::fv
