#include "kutta/mesh/detail.hpp"
#include "post.hpp"
#include <algorithm>
#include <cmath>
#include <set>

namespace kutta::post {
using namespace detail;
namespace {
struct Vertex {
  Vec3d point;
  std::uint64_t id;
  double scalar = 0;
};
double scalar(const Sample &s, int component) {
  if (s.components == 1) {
    KUTTA_REQUIRE(component <= 0,
                  "scalar iso field has no requested vector component");
    return s.value[0];
  }
  return component < 0 ? std::hypot(s.value[0], s.value[1], s.value[2])
                       : s.value[component];
}
} // namespace
view::Geometry isosurface(const Volume &volume, const Iso &spec) {
  const auto description = volume.info(spec.field);
  KUTTA_REQUIRE(!description.boundary_only,
                "iso field is boundary-only: " + spec.field);
  const auto &m = volume.snapshot().mesh;
  view::Geometry out;
  std::vector<std::string> names = spec.fields;
  if (std::find(names.begin(), names.end(), spec.field) == names.end())
    names.push_back(spec.field);
  for (const auto &name : names) {
    const auto f = volume.info(name);
    KUTTA_REQUIRE(!f.boundary_only,
                  "iso output field is boundary-only: " + name);
    out.fields.push_back({name, f.unit, f.components, {}});
  }
  const auto iso_column = std::size_t(
      std::find(names.begin(), names.end(), spec.field) - names.begin());
  std::vector<double> scalar_nodes(m.points.size());
  std::vector<bool> scalar_valid(m.points.size());
  struct Cache {
    std::vector<float> values;
    std::vector<bool> valid;
  };
  std::vector<Cache> caches(names.size());
  std::vector<Sample> center_values(names.size());
  const auto value = [&](std::size_t field, const Vertex &v) {
    if (v.id >= m.points.size())
      return center_values[field];
    auto &cache = caches[field];
    const auto count = out.fields[field].components;
    if (cache.valid.empty()) {
      cache.valid.resize(m.points.size());
      cache.values.resize(m.points.size() * count);
    }
    const auto n = std::size_t(v.id);
    if (!cache.valid[n]) {
      const auto sample = volume.node(names[field], Index(n));
      for (std::uint32_t k = 0; k < count; ++k)
        cache.values[n * count + k] = float(sample.value[k]);
      cache.valid[n] = true;
    }
    Sample result;
    result.components = count;
    for (std::uint32_t k = 0; k < count; ++k)
      result.value[k] = cache.values[n * count + k];
    return result;
  };
  const auto node = [&](Index n) {
    if (!scalar_valid[n]) {
      scalar_nodes[n] = scalar(volume.node(spec.field, n), spec.component);
      scalar_valid[n] = true;
    }
    return Vertex{m.points[n], std::uint64_t(n), scalar_nodes[n]};
  };
  std::map<std::pair<std::uint64_t, std::uint64_t>, std::uint32_t> edges;
  std::set<std::array<std::uint32_t, 3>> triangles;
  const auto intersection = [&](const Vertex &a, const Vertex &b) {
    const double difference = b.scalar - a.scalar;
    double t = difference == 0 ? 0 : (spec.value - a.scalar) / difference;
    t = std::clamp(t, 0.0, 1.0);
    auto key = std::minmax(a.id, b.id);
    std::pair<std::uint64_t, std::uint64_t> identity{key.first, key.second};
    if (t <= 1e-12)
      identity = {a.id, a.id};
    else if (t >= 1 - 1e-12)
      identity = {b.id, b.id};
    const auto found = edges.find(identity);
    if (found != edges.end())
      return found->second;
    KUTTA_REQUIRE(out.positions.size() / 3 < UINT32_MAX,
                  "iso surface exceeds uint32 capacity");
    const auto id = std::uint32_t(out.positions.size() / 3);
    const auto p = a.point + (b.point - a.point) * t;
    out.positions.insert(out.positions.end(),
                         {float(p.x), float(p.y), float(p.z)});
    for (std::size_t f = 0; f < names.size(); ++f) {
      const auto first = value(f, a), second = value(f, b);
      for (std::uint32_t k = 0; k < out.fields[f].components; ++k)
        out.fields[f].values.push_back(
            float((1 - t) * first.value[k] + t * second.value[k]));
    }
    edges.emplace(identity, id);
    return id;
  };
  constexpr int edge_pairs[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                    {1, 2}, {1, 3}, {2, 3}};
  for (Index cell = 0; cell < m.n_cells; ++cell) {
    Vertex center{volume.decomposition_center(cell),
                  std::uint64_t(m.points.size()) + std::uint64_t(cell), 0};
    for (std::size_t f = 0; f < names.size(); ++f)
      center_values[f] = volume.sample(names[f], cell, center.point);
    center.scalar = scalar(center_values[iso_column], spec.component);
    for (Index j = volume.cell_faces().offsets[cell];
         j < volume.cell_faces().offsets[cell + 1]; ++j) {
      const auto face = volume.cell_faces().faces[j];
      const auto ns = nodes(m, face);
      const auto face_triangles = volume.face_triangles(face);
      for (std::size_t t = 0; t < face_triangles.size(); t += 3) {
        Vertex v[4] = {center, node(ns[face_triangles[t]]),
                       node(ns[face_triangles[t + 1]]),
                       node(ns[face_triangles[t + 2]])};
        const auto a = v[1].point - v[0].point, b = v[2].point - v[0].point,
                   c = v[3].point - v[0].point;
        const double determinant = dot(a, cross(b, c));
        KUTTA_REQUIRE(
            determinant * volume.cell_faces().signs[j] >=
                -1e-12 * norm(a) * norm(b) * norm(c),
            "iso decomposition is inverted at cell " + std::to_string(cell) +
                ", face " + std::to_string(face) + ", relative determinant " +
                std::to_string(determinant * volume.cell_faces().signs[j] /
                               (norm(a) * norm(b) * norm(c))));
        if (std::abs(determinant) <= 1e-18 * norm(a) * norm(b) * norm(c))
          continue;
        bool finite = true;
        double lo = INFINITY, hi = -INFINITY;
        for (const auto &x : v) {
          finite &= std::isfinite(x.scalar);
          lo = std::min(lo, x.scalar);
          hi = std::max(hi, x.scalar);
        }
        if (!finite || spec.value < lo || spec.value > hi || lo == hi)
          continue;
        const auto normal = (cross(b, c) * (v[1].scalar - v[0].scalar) +
                             cross(c, a) * (v[2].scalar - v[0].scalar) +
                             cross(a, b) * (v[3].scalar - v[0].scalar)) /
                            determinant;
        if (norm(normal) == 0)
          continue;
        std::vector<std::uint32_t> polygon;
        for (const auto &edge : edge_pairs) {
          const auto &p = v[edge[0]];
          const auto &q = v[edge[1]];
          if ((p.scalar <= spec.value && q.scalar >= spec.value) ||
              (p.scalar >= spec.value && q.scalar <= spec.value)) {
            if (p.scalar == q.scalar) {
              for (const auto *x : {&p, &q}) {
                const auto id = intersection(*x, *x);
                if (std::find(polygon.begin(), polygon.end(), id) ==
                    polygon.end())
                  polygon.push_back(id);
              }
            } else {
              const auto id = intersection(p, q);
              if (std::find(polygon.begin(), polygon.end(), id) ==
                  polygon.end())
                polygon.push_back(id);
            }
          }
        }
        if (polygon.size() < 3)
          continue;
        const auto position = [&](auto id) {
          const auto i = std::size_t(id) * 3;
          return Vec3d{out.positions[i], out.positions[i + 1],
                       out.positions[i + 2]};
        };
        Vec3d mean;
        for (auto id : polygon)
          mean = mean + position(id);
        mean = mean / double(polygon.size());
        const auto n = normal / norm(normal);
        const auto u =
            cross(n, std::abs(n.x) < .8 ? Vec3d{1, 0, 0} : Vec3d{0, 1, 0});
        const auto w = cross(n, u);
        std::sort(polygon.begin(), polygon.end(), [&](auto x, auto y) {
          const auto p = position(x) - mean, q = position(y) - mean;
          return std::atan2(dot(p, w), dot(p, u)) <
                 std::atan2(dot(q, w), dot(q, u));
        });
        for (std::size_t i = 1; i + 1 < polygon.size(); ++i) {
          std::array triangle{polygon[0], polygon[i], polygon[i + 1]};
          if (norm(cross(position(triangle[1]) - position(triangle[0]),
                         position(triangle[2]) - position(triangle[0]))) <=
              volume.tolerance() * volume.tolerance())
            continue;
          auto identity = triangle;
          std::sort(identity.begin(), identity.end());
          if (triangles.insert(identity).second)
            out.indices.insert(out.indices.end(), triangle.begin(),
                               triangle.end());
        }
      }
    }
  }
  return out;
}
} // namespace kutta::post
