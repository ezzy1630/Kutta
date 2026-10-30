#pragma once
#include <string>
#include <vector>
#include "kutta/core/types.hpp"
namespace kutta {
struct OrientedFace {
  Index face=-1;
  int direction=1; // source normal = direction * Kutta owner-outward normal
};
struct InternalFaceZone {
  std::string name,source_type;
  int source_zone_id=-1;
  std::vector<OrientedFace> faces;
};
} // namespace kutta
