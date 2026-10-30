#pragma once
#include "kutta/setup/fluent.hpp"

namespace kutta::setup {
std::vector<fv::CellZonePhysics> fluent_zone_physics(const fluent::FluentCase& input);
fv::Wall fluent_wall(const fluent::BoundaryZone& wall);
} // namespace kutta::setup
