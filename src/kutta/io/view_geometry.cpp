#include "view_geometry.hpp"
#include "binary.hpp"
#include <cmath>
#include <limits>
#include <set>

namespace kutta::view {
Artifact write_geometry(const std::filesystem::path &path, const Geometry &g) {
  const auto vertices = g.positions.size() / 3;
  const bool lines = !g.line_offsets.empty();
  KUTTA_REQUIRE(g.positions.size() % 3 == 0 && vertices <= UINT32_MAX &&
                    g.indices.size() <= UINT32_MAX,
                "geometry exceeds uint32 capacity");
  for (float p : g.positions)
    KUTTA_REQUIRE(std::isfinite(p), "nonfinite float32 geometry position");
  for (auto i : g.indices)
    KUTTA_REQUIRE(i < vertices, "geometry index out of bounds");
  if (lines) {
    KUTTA_REQUIRE(g.line_offsets.front() == 0 &&
                      g.line_offsets.back() == g.indices.size(),
                  "invalid polyline offsets");
    for (std::size_t i = 1; i < g.line_offsets.size(); ++i)
      KUTTA_REQUIRE(g.line_offsets[i] >= g.line_offsets[i - 1],
                    "unordered polyline offsets");
  } else
    KUTTA_REQUIRE(g.indices.size() % 3 == 0,
                  "triangle index count must be a multiple of three");
  binary::Writer w(path);
  for (auto n :
       {lines ? 0x4c54554bu : 0x5654554bu, 1u, 32u, std::uint32_t(vertices),
        std::uint32_t(g.indices.size()), std::uint32_t(g.fields.size()),
        lines ? std::uint32_t(g.line_offsets.size() - 1) : 0u, 0u})
    w.scalar(n);
  w.array<float>(g.positions);
  w.array<std::uint32_t>(g.indices);
  if (lines)
    w.array<std::uint32_t>(g.line_offsets);
  Artifact out;
  Json::Array descriptors;
  std::uint64_t offset = 32 + g.positions.size() * 4 + g.indices.size() * 4 +
                         (lines ? g.line_offsets.size() * 4 : 0);
  std::set<std::string> names;
  for (const auto &f : g.fields) {
    KUTTA_REQUIRE(names.insert(f.name).second &&
                      (f.components == 1 || f.components == 3) &&
                      f.values.size() == vertices * f.components,
                  "invalid vertex field: " + f.name);
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (std::size_t i = 0; i < f.values.size(); i += f.components) {
      for (std::uint32_t c = 0; c < f.components; ++c)
        KUTTA_REQUIRE(!std::isinf(f.values[i + c]),
                      "infinite output field: " + f.name);
      const auto v = f.components == 1 ? double(f.values[i])
                                       : std::hypot(double(f.values[i]),
                                                    double(f.values[i + 1]),
                                                    double(f.values[i + 2]));
      if (std::isfinite(v)) {
        lo = std::min(lo, v);
        hi = std::max(hi, v);
      }
    }
    const auto min =
                   std::isfinite(lo) ? std::optional<double>(lo) : std::nullopt,
               max =
                   std::isfinite(hi) ? std::optional<double>(hi) : std::nullopt;
    descriptors.emplace_back(
        Json::Object{{"name", f.name},
                     {"unit", f.unit},
                     {"components", int(f.components)},
                     {"location", "vertex"},
                     {"offset_bytes", double(offset)},
                     {"min", min ? Json(*min) : Json(nullptr)},
                     {"max", max ? Json(*max) : Json(nullptr)}});
    out.fields.push_back({f.name, f.unit, f.components, min, max});
    w.array<float>(f.values);
    offset += f.values.size() * 4;
  }
  w.close();
  out.descriptor = {{"file", path.filename().string()},
                    {"vertices", double(vertices)},
                    {"indices", double(g.indices.size())},
                    {"fields", std::move(descriptors)},
                    {"primitive", lines ? "polylines" : "triangles"}};
  if (lines)
    out.descriptor.emplace("line_count", double(g.line_offsets.size() - 1));
  return out;
}
} // namespace kutta::view
