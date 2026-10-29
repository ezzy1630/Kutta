#pragma once

#include <functional>
#include <span>
#include <variant>

#include "kutta/fv/transport.hpp"

namespace kutta::fv {

using VelocityProfile=std::variant<Vec3d,std::function<Vec3d(Vec3d)>>;
struct Wall { Vec3d velocity{}; };
struct VelocityInlet { VelocityProfile velocity=Vec3d{}; };
struct PressureOutlet { double pressure=0; VelocityProfile backflow_velocity=Vec3d{}; };
struct Symmetry {};
using BoundaryCondition=std::variant<Wall,VelocityInlet,PressureOutlet,Symmetry>;
struct PatchCondition { std::string patch; BoundaryCondition condition; };

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
  Index reference_cell=0;  // symmetric pin for closed, pure-Neumann pressure
};

struct ResidualValue { double raw=0,scaled=0; };
struct IterationReport {
  int iteration=0;
  ResidualValue u,v,w,continuity;
  double corrected_mass_imbalance=0;
  double maximum_velocity=0, pseudo_time_step=0;
  SolveStats momentum[3],pressure;
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
  const IterationReport& report() const;
 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
}  // namespace kutta::fv
