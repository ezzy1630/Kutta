#pragma once

#include <functional>
#include <array>
#include <optional>
#include <span>
#include <variant>
#include <utility>

#include "kutta/fv/transport.hpp"
#include "kutta/fv/zone_physics.hpp"
#include "kutta/fv/interior_jump.hpp"
#include "kutta/fv/turbulence/sst.hpp"
#include "kutta/fv/turbulence/sa.hpp"
#include "kutta/fv/turbulence/realizable_ke.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta::fv {

using VelocityProfile=std::variant<Vec3d,std::function<Vec3d(Vec3d)>>;
struct Wall {
  Vec3d velocity{};
  std::optional<Rotation> rotation;
  bool relative_to_zone=false;
};
struct TurbulenceValues { double k=0,omega=0; };
struct IntensityLength { double intensity=0,length=0; };
struct IntensityViscosityRatio { double intensity=0,ratio=0; };
struct ModifiedViscosity { double value=0; };
struct SaViscosityRatio { double ratio=0; };
struct KeValues {double k=0,epsilon=0;};
struct IntensityDiameter {double intensity=0,diameter=0;};
using TurbulenceCondition=std::variant<TurbulenceValues,IntensityLength,IntensityViscosityRatio,ModifiedViscosity,SaViscosityRatio,KeValues,IntensityDiameter>;
struct VelocityInlet {
  VelocityProfile velocity=Vec3d{};
  std::optional<TurbulenceCondition> turbulence;
};
struct NormalBackflow {};
using BackflowProfile=std::variant<NormalBackflow,Vec3d,std::function<Vec3d(Vec3d)>>;
enum class BackflowPressure { Static, Total };
struct PressureOutlet {
  double pressure=0; BackflowProfile backflow_velocity=NormalBackflow{};
  BackflowPressure backflow_pressure=BackflowPressure::Static;
  std::optional<TurbulenceCondition> turbulence;
};
struct Symmetry {};
struct PressureInlet {
  double total_pressure=0;
  // Absent means inward face normal. Otherwise global Cartesian direction.
  std::optional<Vec3d> direction;
  std::optional<TurbulenceCondition> turbulence;
};
struct MassFlowInlet {
  double mass_flow=0; // kg/s, positive into the domain, for the WHOLE patch
  std::optional<Vec3d> direction;
  std::optional<TurbulenceCondition> turbulence;
};
using BoundaryCondition=std::variant<Wall,VelocityInlet,PressureOutlet,Symmetry,PressureInlet,MassFlowInlet>;
struct PatchCondition { std::string patch; BoundaryCondition condition; };

// Initialization controls affect no steady equation. The adapter validates
// stored Fluent settings before constructing this value.
struct HybridInitialization {
  int iterations=10;
  std::array<double,2> explicit_relaxation{1,1},implicit_relaxation{1,1},tolerance{1e-6,1e-6},initial_scalar{0,0};
  std::array<bool,2> solve{true,true};
  bool averaged_turbulence=true,linear_potential=false,constant_speed=false;
  double turbulence_intensity=.05,turbulence_viscosity_ratio=10;
  // Only populated when the case requests specified initial inlet pressures.
  std::vector<std::pair<std::string,double>> inlet_pressure;
};

struct SstSettings {
  SstOptions model;
  ConvectionScheme k_convection=ConvectionScheme::SecondOrderUpwind;
  ConvectionScheme omega_convection=ConvectionScheme::SecondOrderUpwind;
  double k_relaxation=.75,omega_relaxation=.75;
  SolveControl linear{1e-3,1e-12,100};
};
struct SaSettings {
  SaOptions model;
  ConvectionScheme convection=ConvectionScheme::SecondOrderUpwind;
  double relaxation=.75;
  SolveControl linear{1e-3,1e-12,100};
};
struct RealizableKeSettings {
  RealizableKeOptions model;
  ConvectionScheme k_convection=ConvectionScheme::SecondOrderUpwind,epsilon_convection=ConvectionScheme::SecondOrderUpwind;
  double k_relaxation=.75,epsilon_relaxation=.75;
  SolveControl linear{1e-3,1e-12,100};
};

// Fluent's flow/scheme. SIMPLE replaces SIMPLEC's consistent inverse with the
// relaxed diagonal's; Coupled solves momentum and continuity as one
// pseudo-transient block system (fv/coupled_solver.hpp). Rhie-Chow always uses
// the unrelaxed physical diagonal, so all three share one discrete fixed point.
enum class PressureVelocityCoupling { Simplec, Simple, Coupled };
// Fluent pressure-velocity options (Theory Guide, pressure-based solver).
// Unlike the coupling, these face interpolations change the fixed point.
// Face velocity in the face mass flux, TG 23-48 (distance) and 23-49 (momentum).
enum class FaceVelocityWeights { Momentum, Distance };
// Face pressure in the momentum pressure force, TG 23-43 (second order) and
// 23-42 (standard, momentum-coefficient weighted).
enum class PressureInterpolation { SecondOrder, Standard };
// Fluent's cell gradient method (Theory Guide 23.3.3): least-squares cell
// based, or Green-Gauss cell based (operators.hpp, green_gauss_cell). It
// changes the fixed point. Green-Gauss serves laminar and k-epsilon flows.
enum class GradientScheme { LeastSquares, GreenGaussCell };

struct SolverSettings {
  double density=1, viscosity=.01;
  PressureVelocityCoupling coupling=PressureVelocityCoupling::Simplec;
  ConvectionScheme convection=ConvectionScheme::SecondOrderUpwind;
  bool limit=true;
  // SIMPLEC: implicit momentum and pressure-correction relaxation. Coupled:
  // Fluent's pseudo-time explicit relaxation of the block update.
  double momentum_relaxation=.7, pressure_relaxation=1;
  // False is Fluent's classic steady formulation: no pseudo-time term, only
  // the implicit under-relaxation factors.
  bool pseudo_transient=true;
  FaceVelocityWeights face_velocity=FaceVelocityWeights::Momentum;
  PressureInterpolation pressure_interpolation=PressureInterpolation::SecondOrder;
  GradientScheme gradient=GradientScheme::LeastSquares;
  double pseudo_time_scale=1;
  // Zero selects the conservative global automatic time step, refreshed every
  // iteration from the current flow. Its velocity scale is the mean cell
  // speed for Coupled (Fluent TG 23-103) and the maximum for SIMPLE/SIMPLEC.
  double pseudo_time_step=0;
  SolveControl momentum_linear{1e-2,1e-10,100};
  SolveControl pressure_linear{1e-3,1e-12,300};
  // Coupled block solve: one order of residual reduction per outer
  // iteration, Fluent's default AMG termination criterion.
  SolveControl coupled_linear{.1,1e-12,20};
  int monitor_interval=20;
  // Closed (pure-Neumann) domains only: gauge pressure is zero in the cell
  // containing this point, as at Fluent's reference pressure location.
  Vec3d reference_point{};
  std::optional<SstSettings> sst;
  // BEGIN FSAE physics: Fluent source IDs survive the owned RCM permutation.
  std::vector<CellZonePhysics> cell_zone_physics;
  std::vector<InteriorCondition> interior_conditions;
  // END FSAE physics
  // Experimental resolved-wall SA; mutually exclusive with SST. The Fluent
  // case adapter remains closed until model and automatic-wall gates pass.
  std::optional<SaSettings> sa;
  std::optional<RealizableKeSettings> realizable_ke;
};

struct ResidualValue {
  double raw=0,scaled=0; // legacy steady scaling, retained for stopping criteria
  double fluent_scaled=0,fluent_normalization=0;
};
struct IterationReport {
  int iteration=0;
  std::optional<int> residual_iteration;
  ResidualValue u,v,w,continuity,k,omega,nu_tilde,epsilon;
  double corrected_mass_imbalance=0;
  double maximum_velocity=0, pseudo_time_step=0;
  SolveStats momentum[3],pressure,turbulence[2];
};
struct ReferenceValues { double density=1,velocity=1,area=1,pressure=0; };
struct ForceReport {
  Vec3d pressure,viscous;
  Vec3d area_vector;
  double reference_pressure=0;
  Vec3d force() const;
  double coefficient(Vec3d direction,const ReferenceValues& reference) const;
};

// Owns an RCM-renumbered copy of the mesh. Device must outlive the solver.
// No host field transfers inside iterate(), except linear solvers and monitors.
class SteadySolver {
 public:
  SteadySolver(Device& dev,PolyMesh mesh,std::vector<PatchCondition> boundaries,
               SolverSettings settings={});
  ~SteadySolver();
  SteadySolver(const SteadySolver&)=delete;
  SteadySolver& operator=(const SteadySolver&)=delete;

  void initialize(const std::function<Vec3d(Vec3d)>& velocity,
                  const std::function<double(Vec3d)>& pressure);
  // Arrays follow the ORIGINAL input mesh cell order, before owned RCM.
  void initialize_fields(std::span<const Vec3d> velocity,std::span<const double> pressure,
      std::span<const double> k={},std::span<const double> omega={});
  void initialize_hybrid(const HybridInitialization& settings);
  void initialize_sa_fields(std::span<const Vec3d> velocity,std::span<const double> pressure,
                            std::span<const double> modified_viscosity);
  void initialize_ke_fields(std::span<const Vec3d> velocity,std::span<const double> pressure,
                           std::span<const double> k,std::span<const double> epsilon);
  // Force per unit volume, not acceleration. Called once at case setup.
  void set_source(const std::function<Vec3d(Vec3d)>& source);
  void iterate(int count);
  IterationReport monitor();
  ForceReport forces(std::span<const std::string> wall_patches,
                     double reference_pressure=0);
  const PolyMesh& mesh() const;
  const MeshGeometry& geometry() const;
  const DeviceMesh& device_mesh() const;
  const DeviceVec3& velocity() const;
  const Buffer<float>& pressure() const;
  const Buffer<float>& mass_flux() const;
  const Buffer<float>& kinetic_energy() const;
  const Buffer<float>& specific_dissipation() const;
  const SstFields& turbulence() const;
  // Wall shear stress on the fluid per boundary face (zero off walls) for the
  // active closure; force on the wall is -shear*|A|. Laminar refreshes it.
  const DeviceVec3& wall_shear() const;
  const SaFields& sa_turbulence() const;
  const SaWallFields& sa_wall() const;
  const RealizableKeFields& realizable_turbulence() const;
  const Buffer<float>& dissipation() const;
  const Buffer<float>& modified_viscosity() const;
  const MeshPermutation& permutation() const;
  const Buffer<float>& boundary_pressure() const;
  const IterationReport& report() const;
  // Device bytes the Coupled block solver owns beyond SIMPLEC; zero for SIMPLEC.
  std::size_t coupled_solver_bytes() const;
 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
}  // namespace kutta::fv
