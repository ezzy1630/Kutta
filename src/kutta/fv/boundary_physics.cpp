#include "kutta/fv/boundary_physics.hpp"
#include <cmath>
#include "kernels.hpp"
#include "kutta/mesh/detail.hpp"

namespace kutta::fv {
using namespace detail;
Vec3d wall_velocity(const Wall& wall,Vec3d point,const RigidFrame& frame) {
  auto value=wall.velocity;
  if(wall.rotation) value=value+cross(wall.rotation->angular_velocity,point-wall.rotation->origin);
  if(wall.relative_to_zone) value=value+frame.velocity(point);
  KUTTA_REQUIRE(finite(value),"wall: nonfinite velocity");return value;
}
Vec3d mass_flow_velocity(const MassFlowInlet& inlet,Vec3d area,double patch_area,
    double projected_area,double rho) {
  KUTTA_REQUIRE(std::isfinite(inlet.mass_flow) && inlet.mass_flow>=0 && patch_area>0 && rho>0,
    "mass-flow-inlet: nonnegative mass flow and positive area/density required");
  if(inlet.direction) {
    KUTTA_REQUIRE(finite(*inlet.direction) && norm(*inlet.direction)>0 && projected_area>0,
      "mass-flow-inlet: nonzero inward direction required");
    const auto d=*inlet.direction/norm(*inlet.direction);
    KUTTA_REQUIRE(dot(d,area)<0,"mass-flow-inlet: direction must enter every face");
    return d*(inlet.mass_flow/(rho*projected_area));
  }
  return area*(-inlet.mass_flow/(rho*patch_area*norm(area)));
}
void pressure_inlet_velocity(Device& dev,const DeviceMesh& mesh,const Buffer<Index>& kind,
    const DeviceVec3& direction,const Buffer<float>& mass,DeviceVec3& boundary,double rho) {
  kern::fv_pressure_inlet_velocity(dev,mesh.n_boundary_faces(),kind,
    mesh.area.x,mesh.area.y,mesh.area.z,direction.x,direction.y,direction.z,mass,
    boundary.x,boundary.y,boundary.z,mesh.n_internal_faces,float(rho));
}
} // namespace kutta::fv
