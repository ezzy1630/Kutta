#pragma once

#include <functional>
#include <optional>
#include <span>
#include <variant>

#include "kutta/fv/transport.hpp"
#include "kutta/fv/turbulence/sst.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta::fv {

using VelocityProfile=std::variant<Vec3d,std::function<Vec3d(Vec3d)>>;
struct Wall { Vec3d velocity{}; };
struct TurbulenceValues { double k=0,omega=0; };
struct IntensityLength { double intensity=0,length=0; };
struct IntensityViscosityRatio { double intensity=0,ratio=0; };
using TurbulenceCondition=std::variant<TurbulenceValues,IntensityLength,IntensityViscosityRatio>;
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
using BoundaryCondition=std::variant<Wall,VelocityInlet,PressureOutlet,Symmetry>;
struct PatchCondition { std::string patch; BoundaryCondition condition; };

struct SstSettings {
  SstOptions model;
  ConvectionScheme k_convection=ConvectionScheme::SecondOrderUpwind;
  ConvectionScheme omega_convection=ConvectionScheme::SecondOrderUpwind;
  double k_relaxation=.75,omega_relaxation=.75;
  SolveControl linear{1e-3,1e-12,100};
};

struct SolverSettings {
  double density=1, viscosity=.01;
  ConvectionScheme convection=ConvectionScheme::SecondOrderUpwind;
  bool limit=true;
  double momentum_relaxation=.7, pressure_relaxation=1;
  double pseudo_time_scale=1;
  // Zero selects the conservative global automatic time scale, updated only
  // at monitor intervals (where speed reductions are already permitted).
  double pseudo_time_step=0;
  SolveControl momentum_linear{1e-2,1e-10,100};
  SolveControl pressure_linear{1e-3,1e-12,300};
  int monitor_interval=20;
  // Closed (pure-Neumann) domains only: gauge pressure is zero in the cell
  // containing this point, as at Fluent's reference pressure location.
  Vec3d reference_point{};
  std::optional<SstSettings> sst;
};

struct ResidualValue {
  double raw=0,scaled=0; // legacy steady scaling, retained for stopping criteria
  double fluent_scaled=0,fluent_normalization=0;
};
struct IterationReport {
  int iteration=0;
  std::optional<int> residual_iteration;
  ResidualValue u,v,w,continuity,k,omega;
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
  const MeshPermutation& permutation() const;
  const Buffer<float>& boundary_pressure() const;
  const IterationReport& report() const;
 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
}  // namespace kutta::fv
