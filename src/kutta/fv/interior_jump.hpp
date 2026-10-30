#pragma once
#include <variant>
#include <optional>
#include <array>
#include "kutta/fv/operators.hpp"

namespace kutta::fv {
struct FanCurve {
  std::vector<double> coefficients; // c0 + c1*|v_normal| + ... [Pa]
  bool reverse=false;
  std::optional<std::array<double,2>> velocity_range;
};
struct PorousJump {
  double permeability=1,thickness=0,inertial=0; // m^2, m, 1/m
};
struct InteriorCondition {
  int source_zone_id=-1;
  std::variant<FanCurve,PorousJump> condition;
};
class InteriorJumps {
 public:
  InteriorJumps(Device& dev,const PolyMesh& mesh,const MeshGeometry& geometry,
                std::span<const InteriorCondition> conditions,double density,double viscosity);
  void update(const Buffer<float>& mass);
  void gradient(const DeviceMesh& mesh,const BoundaryLeastSquares& lsq,DeviceVec3& gradient) const;
  void momentum(const DeviceMesh& mesh,DeviceVec3& source) const;
  void flux(const DeviceMesh& mesh,const Stencil& stencil,const Buffer<float>& diagonal,
            Buffer<float>& mass,Buffer<float>& coefficient) const;
 private:
  Device& dev_;
  Index n_,nf_,active_faces_=0,active_cells_=0;
  bool active_=false;
  float rho_,mu_;
  Buffer<Index> kinds_,counts_,faces_,cells_,offsets_,incidence_,entries_;
  Buffer<float> parameters_,area_,jumps_,slopes_,previous_;
};
} // namespace kutta::fv
