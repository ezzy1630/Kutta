#pragma once
#include "kutta/core/json_value.hpp"
#include "view.hpp"

namespace kutta::view {
struct Geometry {
  std::vector<float> positions;
  std::vector<std::uint32_t> indices;
  std::vector<Field> fields;
  // Nonempty selects KUTL: concatenated polyline indices and L+1 offsets.
  std::vector<std::uint32_t> line_offsets;
};
struct FieldRange {
  std::string name, unit;
  std::uint32_t components = 1;
  std::optional<double> min, max;
};
struct Artifact {
  Json::Object descriptor;
  std::vector<FieldRange> fields;
};
Artifact write_geometry(const std::filesystem::path &, const Geometry &);
std::vector<std::uint32_t> triangulate(std::span<const Vec3d>);
void append_deltas(std::vector<Field> &);
} // namespace kutta::view
