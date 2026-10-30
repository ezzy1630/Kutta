#pragma once

#include <array>
#include <optional>
#include <span>

#include "kutta/fv/device_mesh.hpp"

namespace kutta::fv {
struct Rotation { Vec3d origin{},angular_velocity{}; }; // radians/second, global axes
struct RigidFrame {
  Vec3d translation{};
  Rotation rotation;
  Vec3d velocity(Vec3d point) const;
};
// Fluent's Cartesian principal directions: direction 2 is orthogonalized
// against direction 1, then direction 3 is their right-handed cross product.
struct PorousResistance {
  Vec3d viscous{},inertial{}; // 1/alpha [1/m^2], C2 [1/m]
  Vec3d direction1{1,0,0},direction2{0,1,0};
  double porosity=1;
  bool relative_velocity=true;
};
struct CellZonePhysics {
  int source_zone_id=-1;
  std::optional<RigidFrame> frame;
  std::optional<PorousResistance> porous;
};

std::array<Vec3d,3> porous_basis(const PorousResistance& resistance);
// Resolves already-absolute parent motion at the child's origin.
RigidFrame compose_frames(const RigidFrame& parent,const RigidFrame& child);
std::vector<RigidFrame> cell_frames(const PolyMesh& mesh,std::span<const CellZonePhysics> zones);

// All field-dependent work stays in generated Device kernels. Host geometry
// is used only once to construct tensors and rigid-frame face fluxes.
class ZonePhysics {
 public:
  ZonePhysics(Device& dev,const PolyMesh& mesh,const MeshGeometry& geometry,
              std::span<const CellZonePhysics> zones,double density,double viscosity);
  void relative_flux(const Buffer<Index>& boundary_kind,Index internal,Buffer<float>& mass) const;
  void momentum(int component,const DeviceVec3& velocity,const Buffer<float>& volume,
                Buffer<float>& diagonal,Buffer<float>& rhs) const;
  const Buffer<float>& frame_normal_flux() const { return frame_flux_; }
  bool has_frame_motion() const { return moving_; }
 private:
  Device& dev_;
  Index n_,nf_,zones_=0;
  float rho_,mu_;
  bool active_=false,moving_=false,relative_motion_=false;
  Buffer<Index> zone_index_;
  Buffer<float> coefficients_,frame_flux_,relative_velocity_;
};
} // namespace kutta::fv
