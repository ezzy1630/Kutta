#pragma once

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

namespace kutta::detail {
struct Plot3dBlock {
  std::array<int, 3> dims;
  Index point_offset;
  bool extruded;
};
void apply_neutral_map(const std::filesystem::path& path,
                       const std::vector<Plot3dBlock>& blocks,
                       std::vector<Vec3d>& points, std::vector<Face>& faces,
                       std::vector<Patch>& patches);
}  // namespace kutta::detail
