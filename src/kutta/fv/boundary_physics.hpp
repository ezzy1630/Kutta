#pragma once
#include "kutta/fv/steady_solver.hpp"

namespace kutta::fv {
Vec3d wall_velocity(const Wall& wall,Vec3d point,const RigidFrame& frame);
Vec3d mass_flow_velocity(const MassFlowInlet& inlet,Vec3d area,double patch_area,
                        double projected_area,double density);
// Overrides pressure-inlet advected momentum only. Pressure/flux remain coupled
// through the pressure boundary and pressure-correction operator.
void pressure_inlet_velocity(Device& dev,const DeviceMesh& mesh,const Buffer<Index>& kind,
    const DeviceVec3& direction,const Buffer<float>& mass,DeviceVec3& boundary,double density);
} // namespace kutta::fv
