#include "kutta/mesh/detail.hpp"
#include "post.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace kutta::post {
using namespace detail;
namespace {
constexpr double nan = std::numeric_limits<double>::quiet_NaN();
const view::Field *find(std::span<const view::Field> fields,
                        std::string_view name) {
  const auto it = std::find_if(fields.begin(), fields.end(),
                               [&](const auto &f) { return f.name == name; });
  return it == fields.end() ? nullptr : &*it;
}

double component(Vec3d p, int j) { return j == 0 ? p.x : j == 1 ? p.y : p.z; }
Vec3d solve(double a[3][3], Vec3d rhs) {
  // Jacobi pseudoinverse handles genuinely rank-deficient one-layer meshes;
  // no arbitrary regularization biases a resolved linear direction.
  double v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  for (int sweep = 0; sweep < 12; ++sweep)
    for (int p = 0; p < 2; ++p)
      for (int q = p + 1; q < 3; ++q) {
        if (std::abs(a[p][q]) <=
            1e-16 * (std::abs(a[p][p]) + std::abs(a[q][q])))
          continue;
        const double angle = .5 * std::atan2(2 * a[p][q], a[q][q] - a[p][p]),
                     c = std::cos(angle), s = std::sin(angle);
        const double app = c * c * a[p][p] - 2 * c * s * a[p][q] +
                           s * s * a[q][q],
                     aqq = s * s * a[p][p] + 2 * c * s * a[p][q] +
                           c * c * a[q][q];
        for (int k = 0; k < 3; ++k)
          if (k != p && k != q) {
            const double x = a[k][p], y = a[k][q];
            a[k][p] = a[p][k] = c * x - s * y;
            a[k][q] = a[q][k] = s * x + c * y;
          }
        a[p][p] = app;
        a[q][q] = aqq;
        a[p][q] = a[q][p] = 0;
        for (int k = 0; k < 3; ++k) {
          const double x = v[k][p], y = v[k][q];
          v[k][p] = c * x - s * y;
          v[k][q] = s * x + c * y;
        }
      }
  const auto largest = std::max({a[0][0], a[1][1], a[2][2]});
  Vec3d out;
  for (int j = 0; j < 3; ++j)
    if (a[j][j] > largest * 1e-12 && a[j][j] > 0) {
      const Vec3d axis{v[0][j], v[1][j], v[2][j]};
      out = out + axis * (dot(axis, rhs) / a[j][j]);
    }
  return out;
}
struct Bounds {
  Vec3d lo{INFINITY, INFINITY, INFINITY}, hi{-INFINITY, -INFINITY, -INFINITY};
  void add(Vec3d p) {
    lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
    hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
  }
  void add(const Bounds &b) {
    add(b.lo);
    add(b.hi);
  }
  bool contains(Vec3d p, double e) const {
    return p.x >= lo.x - e && p.x <= hi.x + e && p.y >= lo.y - e &&
           p.y <= hi.y + e && p.z >= lo.z - e && p.z <= hi.z + e;
  }
};
struct Node {
  Bounds box;
  Index left = -1, right = -1, begin = 0, end = 0;
};
struct Spatial {
  const PolyMesh &mesh;
  MeshGeometry geometry;
  kutta::CellFaces faces;
  double scale = 1, eps = 1e-10;
  mutable std::vector<Bounds> boxes;
  mutable std::vector<Node> tree;
  mutable std::vector<Index> order;
  mutable std::vector<std::int8_t> convex;
  explicit Spatial(const PolyMesh &m)
      : mesh(m), geometry(compute_geometry(m)), faces(build_cell_faces(m)),
        convex(m.n_cells, -1) {
    Bounds b;
    for (auto p : m.points)
      b.add(p);
    scale = norm(b.hi - b.lo);
    eps = std::max(scale * 1e-10, std::max(norm(b.lo), norm(b.hi)) * 8 *
                                      std::numeric_limits<double>::epsilon());
  }
  Index build(Index begin, Index end) const {
    const auto id = index(tree.size());
    tree.push_back({});
    Node n;
    n.begin = begin;
    n.end = end;
    for (Index i = begin; i < end; ++i)
      n.box.add(boxes[order[i]]);
    if (end - begin > 8) {
      const auto size = n.box.hi - n.box.lo;
      const int axis = size.x >= size.y && size.x >= size.z ? 0
                       : size.y >= size.z                   ? 1
                                                            : 2;
      const auto mid = begin + (end - begin) / 2;
      std::nth_element(order.begin() + begin, order.begin() + mid,
                       order.begin() + end, [&](Index a, Index b) {
                         return component(geometry.cell_centroid[a], axis) <
                                component(geometry.cell_centroid[b], axis);
                       });
      n.left = build(begin, mid);
      n.right = build(mid, end);
    }
    tree[id] = n;
    return id;
  }
  void index_cells() const {
    if (!tree.empty())
      return;
    boxes.resize(mesh.n_cells);
    order.resize(mesh.n_cells);
    std::iota(order.begin(), order.end(), 0);
    for (Index c = 0; c < mesh.n_cells; ++c)
      for (Index j = faces.offsets[c]; j < faces.offsets[c + 1]; ++j)
        for (auto n : nodes(mesh, faces.faces[j]))
          boxes[c].add(mesh.points[n]);
    tree.reserve(std::size_t(mesh.n_cells) / 2 + 1);
    build(0, mesh.n_cells);
  }
};
std::optional<double> hit_triangle(Vec3d from, Vec3d delta, Vec3d a, Vec3d b,
                                   Vec3d c, double eps) {
  const auto e1 = b - a, e2 = c - a, h = cross(delta, e2);
  const auto det = dot(e1, h);
  if (std::abs(det) <= 1e-14 * norm(delta) * norm(e1) * norm(e2))
    return {};
  const auto s = from - a;
  const double u = dot(s, h) / det;
  const auto q = cross(s, e1);
  const double v = dot(delta, q) / det, t = dot(e2, q) / det;
  if (u < -eps || v < -eps || u + v > 1 + eps || t < -eps || t > 1 + eps)
    return {};
  return std::clamp(t, 0.0, 1.0);
}
std::string_view base_name(std::string_view name) {
  const auto at = name.find('@');
  return name.substr(0, at);
}
} // namespace
struct Volume::Impl {
  const view::Snapshot &snapshot;
  std::shared_ptr<Spatial> spatial;
  mutable std::map<std::string, std::vector<float>> gradients;
  mutable std::vector<Index> point_offsets, point_cells;
  explicit Impl(const view::Snapshot &s, std::shared_ptr<Spatial> shared = {})
      : snapshot(s), spatial(shared ? std::move(shared)
                                    : std::make_shared<Spatial>(s.mesh)) {}
  const std::vector<float> &gradient(const view::Field &field) const {
    const auto existing = gradients.find(field.name);
    if (existing != gradients.end())
      return existing->second;
    const auto &m = spatial->mesh;
    const auto &g = spatial->geometry;
    const auto &cf = spatial->faces;
    const auto *boundary = find(snapshot.boundary, field.name);
    std::vector<float> out(std::size_t(m.n_cells) * field.components * 3);
    for (Index cell = 0; cell < m.n_cells; ++cell) {
      double a[3][3]{};
      Vec3d rhs[3]{};
      for (Index j = cf.offsets[cell]; j < cf.offsets[cell + 1]; ++j) {
        const auto f = cf.faces[j];
        Index other = -1;
        Vec3d delta;
        std::size_t sample;
        if (f < m.n_internal_faces()) {
          other = m.owner[f] == cell ? m.neighbour[f] : m.owner[f];
          delta = g.cell_centroid[other] - g.cell_centroid[cell];
          sample = std::size_t(other) * field.components;
        } else {
          if (!boundary)
            continue;
          delta = g.face_centroid[f] - g.cell_centroid[cell];
          sample = std::size_t(f - m.n_internal_faces()) * field.components;
        }
        const auto &values = other >= 0 ? field.values : boundary->values;
        bool valid = true;
        for (std::uint32_t k = 0; k < field.components; ++k)
          valid &= std::isfinite(values[sample + k]) &&
                   std::isfinite(
                       field.values[std::size_t(cell) * field.components + k]);
        if (!valid)
          continue;
        const double length2 = dot(delta, delta);
        if (length2 == 0)
          continue;
        const double weight = 1 / length2;
        for (int x = 0; x < 3; ++x)
          for (int y = 0; y < 3; ++y)
            a[x][y] += component(delta, x) * component(delta, y) * weight;
        for (std::uint32_t k = 0; k < field.components; ++k)
          rhs[k] =
              rhs[k] +
              delta *
                  ((values[sample + k] -
                    field.values[std::size_t(cell) * field.components + k]) *
                   weight);
      }
      for (std::uint32_t k = 0; k < field.components; ++k) {
        double copy[3][3];
        for (int row = 0; row < 3; ++row)
          for (int col = 0; col < 3; ++col)
            copy[row][col] = a[row][col];
        const auto v = solve(copy, rhs[k]);
        const auto i = (std::size_t(cell) * field.components + k) * 3;
        out[i] = float(v.x);
        out[i + 1] = float(v.y);
        out[i + 2] = float(v.z);
      }
    }
    return gradients.emplace(field.name, std::move(out)).first->second;
  }
  Sample raw(const view::Field &f, Index cell, Vec3d point) const {
    Sample out;
    out.components = f.components;
    const auto &grad = gradient(f);
    const auto delta = point - spatial->geometry.cell_centroid[cell];
    for (std::uint32_t k = 0; k < f.components; ++k) {
      const auto i = (std::size_t(cell) * f.components + k) * 3;
      out.value[k] = f.values[std::size_t(cell) * f.components + k] +
                     dot({grad[i], grad[i + 1], grad[i + 2]}, delta);
    }
    return out;
  }
};
Volume::Volume(const view::Snapshot &s, const Volume *same_geometry)
    : impl_(std::make_unique<Impl>(
          s, same_geometry ? same_geometry->impl_->spatial : nullptr)) {}
Volume::~Volume() = default;
const view::Snapshot &Volume::snapshot() const { return impl_->snapshot; }
const MeshGeometry &Volume::geometry() const {
  return impl_->spatial->geometry;
}
const kutta::CellFaces &Volume::cell_faces() const {
  return impl_->spatial->faces;
}
double Volume::scale() const { return impl_->spatial->scale; }
double Volume::tolerance() const { return impl_->spatial->eps; }
std::vector<std::uint32_t> Volume::face_triangles(Index face) const {
  const auto &m = snapshot().mesh;
  std::vector<Vec3d> points;
  for (auto n : nodes(m, face))
    points.push_back(m.points[n]);
  auto result = view::triangulate(points);
  const auto quality = [&](const auto &triangles) {
    double worst = INFINITY;
    for (std::size_t t = 0; t < triangles.size(); t += 3) {
      const auto a = points[triangles[t]],
                 normal = cross(points[triangles[t + 1]] - a,
                                points[triangles[t + 2]] - a);
      const auto size = norm(normal);
      if (size == 0)
        continue;
      worst = std::min(
          worst,
          dot(a - geometry().cell_centroid[m.owner[face]], normal) / size);
      if (face < m.n_internal_faces())
        worst =
            std::min(worst, dot(geometry().cell_centroid[m.neighbour[face]] - a,
                                normal) /
                                size);
    }
    return worst;
  };
  double best = quality(result);
  if (best >= 0)
    return result;
  // A warped polygon has no unique diagonal. Select one decomposition for
  // both neighbours, preferring positive tetrahedra about their centres.
  // Rotated ear clipping also preserves concave projected polygons.
  for (std::size_t shift = 1; shift < points.size(); ++shift) {
    std::vector<Vec3d> rotated;
    for (std::size_t i = 0; i < points.size(); ++i)
      rotated.push_back(points[(i + shift) % points.size()]);
    auto triangles = view::triangulate(rotated);
    for (auto &i : triangles)
      i = std::uint32_t((i + shift) % points.size());
    const auto score = quality(triangles);
    if (score > best) {
      best = score;
      result = std::move(triangles);
    }
    if (best >= 0)
      break;
  }
  return result;
}
Vec3d Volume::decomposition_center(Index cell) const {
  struct Plane {
    Vec3d point, normal;
  };
  std::vector<Plane> planes;
  const auto &m = snapshot().mesh;
  const auto &cf = cell_faces();
  auto center = geometry().cell_centroid[cell];
  bool inside = true;
  for (Index j = cf.offsets[cell]; j < cf.offsets[cell + 1]; ++j) {
    const auto face = cf.faces[j];
    const auto ns = nodes(m, face);
    const auto triangles = face_triangles(face);
    for (std::size_t t = 0; t < triangles.size(); t += 3) {
      const auto a = m.points[ns[triangles[t]]];
      auto n = cross(m.points[ns[triangles[t + 1]]] - a,
                     m.points[ns[triangles[t + 2]]] - a) *
               double(cf.signs[j]);
      if (norm(n) == 0)
        continue;
      n = n / norm(n);
      planes.push_back({a, n});
      inside &= dot(center - a, n) < 0;
    }
  }
  if (inside)
    return center;
  // A Fluent face-pyramid centroid can lie outside a warped triangle's
  // halfspace. Find a point in the triangulated cell kernel, without changing
  // any solver geometry or allowing overlapping, inverted tetrahedra.
  const double margin = std::cbrt(geometry().cell_volume[cell]) * 1e-10;
  for (int iteration = 0; iteration < 2048; ++iteration) {
    double worst = 0;
    for (const auto &p : planes) {
      const double distance = dot(center - p.point, p.normal) + margin;
      worst = std::max(worst, distance);
      if (distance > 0)
        center = center - p.normal * distance;
    }
    if (worst <= margin * .01)
      return center;
  }
  throw Error("iso decomposition has no positive cell kernel at cell " +
              std::to_string(cell));
}
FieldInfo Volume::info(std::string_view name) const {
  const auto base = base_name(name);
  const auto suffix = name.substr(base.size());
  KUTTA_REQUIRE(suffix.empty() || suffix == "@fluent" || suffix == "@delta",
                "unknown field suffix: " + std::string(name));
  if (suffix == "@delta") {
    const auto reference = info(std::string(base) + "@fluent");
    const auto primary = info(base);
    KUTTA_REQUIRE(primary.components == reference.components,
                  "delta field shape mismatch");
    return {std::string(name), primary.unit, primary.components,
            primary.boundary_only};
  }
  if (const auto *f = find(snapshot().cells, name))
    return {f->name, f->unit, f->components, false};
  if (const auto *f = find(snapshot().boundary, name))
    return {f->name, f->unit, f->components, true};
  const auto velocity = std::string("velocity") + std::string(suffix),
             pressure = std::string("pressure") + std::string(suffix);
  if (base == "velocity_magnitude" || base == "velocity_x" ||
      base == "velocity_y" || base == "velocity_z" || base == "vorticity" ||
      base == "vorticity_magnitude" || base == "q_criterion" ||
      base == "helicity") {
    KUTTA_REQUIRE(find(snapshot().cells, velocity),
                  "velocity is unavailable for " + std::string(name));
    return {std::string(name),
            base == "q_criterion"                                  ? "1/s^2"
            : base == "vorticity" || base == "vorticity_magnitude" ? "1/s"
            : base == "helicity"                                   ? "m/s^2"
                                                                   : "m/s",
            base == "vorticity" ? 3u : 1u, false};
  }
  if (base == "pressure_coefficient" || base == "total_pressure" ||
      base == "total_pressure_coefficient") {
    KUTTA_REQUIRE(find(snapshot().cells, pressure),
                  "pressure is unavailable for " + std::string(name));
    const auto &r = snapshot().properties;
    if (base != "pressure_coefficient")
      KUTTA_REQUIRE(r.density && *r.density > 0 &&
                        find(snapshot().cells, velocity),
                    "total pressure requires density and velocity");
    if (base != "total_pressure")
      KUTTA_REQUIRE(r.reference_density && *r.reference_density > 0 &&
                        r.velocity && *r.velocity > 0 && r.pressure,
                    "coefficient requires saved reference density, speed and "
                    "pressure (or post reference overrides)");
    return {std::string(name), base == "total_pressure" ? "Pa" : "1", 1, false};
  }
  if (base == "wall_shear_magnitude") {
    KUTTA_REQUIRE(find(snapshot().boundary,
                       std::string("wall_shear") + std::string(suffix)),
                  "wall shear is unavailable");
    return {std::string(name), "Pa", 1, true};
  }
  if (base == "skin_friction" || base == "skin_friction_vector") {
    const auto &r = snapshot().properties;
    KUTTA_REQUIRE(r.reference_density && *r.reference_density > 0 &&
                      r.velocity && *r.velocity > 0,
                  "skin friction requires reference density and speed");
    KUTTA_REQUIRE(find(snapshot().boundary,
                       std::string("wall_shear") + std::string(suffix)),
                  "wall shear is unavailable for " + std::string(name));
    return {std::string(name), "1", base == "skin_friction_vector" ? 3u : 1u,
            true};
  }
  throw Error("unknown or unavailable post-processing field: " +
              std::string(name));
}
std::vector<FieldInfo> Volume::available() const {
  std::map<std::string, FieldInfo> fields;
  for (const auto *source : {&snapshot().cells, &snapshot().boundary})
    for (const auto &f : *source)
      fields.emplace(f.name, FieldInfo{f.name, f.unit, f.components,
                                       source == &snapshot().boundary &&
                                           !find(snapshot().cells, f.name)});
  for (const auto *source : {&snapshot().cells, &snapshot().boundary})
    for (const auto &field : *source)
      if (field.name.ends_with("@fluent"))
        try {
          auto info_delta =
              info(field.name.substr(0, field.name.size() - 7) + "@delta");
          fields.insert_or_assign(info_delta.name, info_delta);
        } catch (const Error &) {
        }
  for (auto suffix : {"", "@fluent", "@delta"})
    for (auto name :
         {"velocity_magnitude", "velocity_x", "velocity_y", "velocity_z",
          "vorticity", "vorticity_magnitude", "q_criterion", "helicity",
          "wall_shear_magnitude", "pressure_coefficient", "total_pressure",
          "total_pressure_coefficient", "skin_friction",
          "skin_friction_vector"})
      try {
        auto f = info(std::string(name) + suffix);
        fields.insert_or_assign(f.name, f);
      } catch (const Error &) {
      }
  std::vector<FieldInfo> out;
  for (auto &[name, f] : fields) {
    (void)name;
    out.push_back(std::move(f));
  }
  return out;
}
std::array<Vec3d, 3> Volume::velocity_gradient(Index cell,
                                               std::string_view suffix) const {
  const auto *f =
      find(snapshot().cells, std::string("velocity") + std::string(suffix));
  KUTTA_REQUIRE(f && f->components == 3, "velocity field unavailable");
  const auto &g = impl_->gradient(*f);
  std::array<Vec3d, 3> out;
  for (int i = 0; i < 3; ++i) {
    const auto at = (std::size_t(cell) * 3 + i) * 3;
    out[i] = {g[at], g[at + 1], g[at + 2]};
  }
  return out;
}
Sample Volume::sample(std::string_view name, Index cell, Vec3d point) const {
  KUTTA_REQUIRE(cell >= 0 && cell < impl_->spatial->mesh.n_cells,
                "sample cell out of range");
  const auto description = info(name);
  const auto base = base_name(name);
  const auto suffix = name.substr(base.size());
  if (suffix == "@delta") {
    auto out = sample(base, cell, point);
    const auto ref = sample(std::string(base) + "@fluent", cell, point);
    for (std::uint32_t i = 0; i < out.components; ++i)
      out.value[i] -= ref.value[i];
    return out;
  }
  const auto velocity_name = std::string("velocity") + std::string(suffix);
  // Derived quantities are evaluated from reconstructed primitives, rather
  // than interpolating a nonlinear value precomputed at a cell centre.
  if (base == "velocity_magnitude" || base == "velocity_x" ||
      base == "velocity_y" || base == "velocity_z") {
    const auto u = sample(velocity_name, cell, point);
    return {{{base == "velocity_magnitude"
                  ? std::hypot(u.value[0], u.value[1], u.value[2])
                  : u.value[base.back() == 'x'   ? 0
                            : base.back() == 'y' ? 1
                                                 : 2],
              0, 0}},
            1};
  }
  if (base == "vorticity" || base == "vorticity_magnitude" ||
      base == "q_criterion" || base == "helicity") {
    const auto g = velocity_gradient(cell, suffix);
    const Vec3d curl{g[2].y - g[1].z, g[0].z - g[2].x, g[1].x - g[0].y};
    if (base == "vorticity")
      return {{{curl.x, curl.y, curl.z}}, 3};
    if (base == "vorticity_magnitude")
      return {{{norm(curl), 0, 0}}, 1};
    if (base == "helicity") {
      const auto u = sample(velocity_name, cell, point);
      return {{{dot({u.value[0], u.value[1], u.value[2]}, curl), 0, 0}}, 1};
    }
    double strain = 0, rotation = 0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        const double a = component(g[i], j), b = component(g[j], i);
        strain += .25 * (a + b) * (a + b);
        rotation += .25 * (a - b) * (a - b);
      }
    return {{{.5 * (rotation - strain), 0, 0}}, 1};
  }
  if (base == "pressure_coefficient" || base == "total_pressure" ||
      base == "total_pressure_coefficient") {
    const auto &r = snapshot().properties;
    if (base == "pressure_coefficient" &&
        (!r.reference_density || !r.velocity || !r.pressure)) {
      const auto *f = find(snapshot().cells, name);
      KUTTA_REQUIRE(f, "coefficient reference unavailable");
      return impl_->raw(*f, cell, point);
    }
    double p =
        sample(std::string("pressure") + std::string(suffix), cell, point)
            .value[0];
    if (base != "pressure_coefficient") {
      const auto u = sample(velocity_name, cell, point);
      p += .5 * *r.density *
           (u.value[0] * u.value[0] + u.value[1] * u.value[1] +
            u.value[2] * u.value[2]);
    }
    if (base != "total_pressure")
      p = (p - *r.pressure) /
          (.5 * *r.reference_density * *r.velocity * *r.velocity);
    return {{{p, 0, 0}}, 1};
  }
  const auto *f = find(snapshot().cells, name);
  KUTTA_REQUIRE(f && !description.boundary_only,
                "field is boundary-only: " + std::string(name));
  return impl_->raw(*f, cell, point);
}
Sample Volume::boundary(std::string_view name, Index face, Vec3d point) const {
  const auto &m = impl_->spatial->mesh;
  KUTTA_REQUIRE(face >= m.n_internal_faces() && face < m.n_faces(),
                "expected boundary face");
  const auto base = base_name(name), suffix = name.substr(base.size());
  if (suffix == "@delta") {
    auto out = boundary(base, face, point);
    const auto ref = boundary(std::string(base) + "@fluent", face, point);
    for (std::uint32_t i = 0; i < out.components; ++i)
      out.value[i] -= ref.value[i];
    return out;
  }

  if (base == "wall_shear_magnitude") {
    const auto shear =
        boundary(std::string("wall_shear") + std::string(suffix), face, point);
    return {
        {{std::hypot(shear.value[0], shear.value[1], shear.value[2]), 0, 0}},
        1};
  }
  if (base == "skin_friction" || base == "skin_friction_vector") {
    (void)info(name);
    const auto shear =
        boundary(std::string("wall_shear") + std::string(suffix), face, point);
    const auto &r = snapshot().properties;
    const double q = .5 * *r.reference_density * *r.velocity * *r.velocity;
    if (base == "skin_friction")
      return {{{std::hypot(shear.value[0], shear.value[1], shear.value[2]) / q,
                0, 0}},
              1};
    auto out = shear;
    for (auto &v : out.value)
      v /= q;
    return out;
  }
  if (base == "velocity_magnitude" || base == "velocity_x" ||
      base == "velocity_y" || base == "velocity_z") {
    const auto u =
        boundary(std::string("velocity") + std::string(suffix), face, point);
    return {{{base == "velocity_magnitude"
                  ? std::hypot(u.value[0], u.value[1], u.value[2])
                  : u.value[base.back() == 'x'   ? 0
                            : base.back() == 'y' ? 1
                                                 : 2],
              0, 0}},
            1};
  }
  if (base == "pressure_coefficient" || base == "total_pressure" ||
      base == "total_pressure_coefficient") {
    const auto &r = snapshot().properties;
    if (base == "pressure_coefficient" &&
        (!r.reference_density || !r.velocity || !r.pressure)) {
      const auto *f = find(snapshot().boundary, name);
      KUTTA_REQUIRE(f, "coefficient reference unavailable");
      return {{{f->values[face - m.n_internal_faces()], 0, 0}}, 1};
    }
    (void)info(name);
    double p =
        boundary(std::string("pressure") + std::string(suffix), face, point)
            .value[0];
    if (base != "pressure_coefficient") {
      const auto u =
          boundary(std::string("velocity") + std::string(suffix), face, point);
      p += .5 * *r.density *
           (u.value[0] * u.value[0] + u.value[1] * u.value[1] +
            u.value[2] * u.value[2]);
    }
    if (base != "total_pressure")
      p = (p - *r.pressure) /
          (.5 * *r.reference_density * *r.velocity * *r.velocity);
    return {{{p, 0, 0}}, 1};
  }
  if (const auto *f = find(snapshot().boundary, name)) {
    Sample out;
    out.components = f->components;
    for (std::uint32_t i = 0; i < f->components; ++i)
      out.value[i] =
          f->values[std::size_t(face - m.n_internal_faces()) * f->components +
                    i];
    return out;
  }
  if (base == "helicity") {
    const auto u = boundary(std::string("velocity") + std::string(suffix), face,
                            point),
               omega = sample(std::string("vorticity") + std::string(suffix),
                              m.owner[face], point);
    return {{{u.value[0] * omega.value[0] + u.value[1] * omega.value[1] +
                  u.value[2] * omega.value[2],
              0, 0}},
            1};
  }
  return sample(name, m.owner[face], point);
}
Sample Volume::node(std::string_view name, Index node) const {
  const auto &m = impl_->spatial->mesh;
  if (impl_->point_offsets.empty()) {
    impl_->point_offsets.assign(m.points.size() + 1, 0);
    std::vector<Index> ids;
    for (Index c = 0; c < m.n_cells; ++c) {
      ids.clear();
      for (Index j = cell_faces().offsets[c]; j < cell_faces().offsets[c + 1];
           ++j)
        for (auto n : nodes(m, cell_faces().faces[j]))
          ids.push_back(n);
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      for (auto n : ids)
        ++impl_->point_offsets[n + 1];
    }
    for (std::size_t i = 1; i < impl_->point_offsets.size(); ++i)
      impl_->point_offsets[i] += impl_->point_offsets[i - 1];
    impl_->point_cells.resize(impl_->point_offsets.back());
    auto cursor = impl_->point_offsets;
    for (Index c = 0; c < m.n_cells; ++c) {
      ids.clear();
      for (Index j = cell_faces().offsets[c]; j < cell_faces().offsets[c + 1];
           ++j)
        for (auto n : nodes(m, cell_faces().faces[j]))
          ids.push_back(n);
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      for (auto n : ids)
        impl_->point_cells[cursor[n]++] = c;
    }
  }
  Sample out;
  out.components = info(name).components;
  double sum = 0;
  bool missing = false;
  for (Index j = impl_->point_offsets[node]; j < impl_->point_offsets[node + 1];
       ++j) {
    const auto c = impl_->point_cells[j];
    const auto value = sample(name, c, m.points[node]);
    const auto delta = m.points[node] - geometry().cell_centroid[c];
    const double weight =
        1 / std::max(dot(delta, delta), tolerance() * tolerance());
    for (std::uint32_t k = 0; k < out.components; ++k) {
      out.value[k] += weight * value.value[k];
      missing |= !std::isfinite(value.value[k]);
    }
    sum += weight;
  }
  for (std::uint32_t k = 0; k < out.components; ++k)
    out.value[k] = missing || sum == 0 ? nan : out.value[k] / sum;
  return out;
}
bool Volume::contains(Index cell, Vec3d p) const {
  if (cell < 0)
    return false;
  const auto &s = *impl_->spatial;
  const auto &m = s.mesh;
  const auto &cf = s.faces;
  const auto &g = s.geometry;
  if (s.convex[cell] < 0) {
    bool convex = true;
    for (Index i = cf.offsets[cell]; i < cf.offsets[cell + 1] && convex; ++i) {
      const auto face = cf.faces[i];
      const auto normal =
          g.face_area[face] * (cf.signs[i] / norm(g.face_area[face]));
      for (Index j = cf.offsets[cell]; j < cf.offsets[cell + 1] && convex; ++j)
        for (auto n : nodes(m, cf.faces[j]))
          if (dot(m.points[n] - g.face_centroid[face], normal) > s.eps) {
            convex = false;
            break;
          }
    }
    s.convex[cell] = convex ? 1 : 0;
  }
  if (s.convex[cell]) {
    for (Index i = cf.offsets[cell]; i < cf.offsets[cell + 1]; ++i) {
      const auto f = cf.faces[i];
      if (dot(p - g.face_centroid[f], g.face_area[f]) * cf.signs[i] >
          s.eps * norm(g.face_area[f]))
        return false;
    }
    return true;
  }
  double angle = 0;
  for (Index i = cf.offsets[cell]; i < cf.offsets[cell + 1]; ++i) {
    const auto f = cf.faces[i];
    const auto ns = nodes(m, f);
    const auto triangles = face_triangles(f);
    for (std::size_t t = 0; t < triangles.size(); t += 3) {
      const auto a = m.points[ns[triangles[t]]] - p,
                 b = m.points[ns[triangles[t + 1]]] - p,
                 c = m.points[ns[triangles[t + 2]]] - p;
      if (std::min({norm(a), norm(b), norm(c)}) <= s.eps)
        return true;
      const auto normal = cross(b - a, c - a);
      const double normal2 = dot(normal, normal);
      if (normal2 > 0 &&
          std::abs(dot(a, normal)) <= s.eps * std::sqrt(normal2)) {
        const double edge = std::max({norm(b - a), norm(c - b), norm(a - c)}),
                     e = s.eps / edge;
        if (dot(cross(b, c), normal) / normal2 >= -e &&
            dot(cross(c, a), normal) / normal2 >= -e &&
            dot(cross(a, b), normal) / normal2 >= -e)
          return true;
      }
      const double numerator = dot(a, cross(b, c));
      const double denominator = norm(a) * norm(b) * norm(c) +
                                 dot(a, b) * norm(c) + dot(b, c) * norm(a) +
                                 dot(c, a) * norm(b);
      angle += double(cf.signs[i]) * 2 * std::atan2(numerator, denominator);
    }
  }
  return std::abs(angle) > 2 * std::acos(-1.0);
}
Index Volume::locate(Vec3d p) const {
  const auto &s = *impl_->spatial;
  s.index_cells();
  std::vector<Index> stack{0};
  Index found = -1;
  while (!stack.empty()) {
    const auto id = stack.back();
    stack.pop_back();
    const auto &n = s.tree[id];
    if (!n.box.contains(p, s.eps))
      continue;
    if (n.left >= 0) {
      stack.push_back(n.left);
      stack.push_back(n.right);
    } else
      for (Index i = n.begin; i < n.end; ++i) {
        const auto c = s.order[i];
        if ((found < 0 || c < found) && s.boxes[c].contains(p, s.eps) &&
            contains(c, p))
          found = c;
      }
  }
  return found;
}
Walk Volume::walk(Vec3d from, Index cell, Vec3d to) const {
  const auto &s = *impl_->spatial;
  const auto &m = s.mesh;
  const auto delta = to - from;
  const double length = norm(delta);
  if (length == 0)
    return {cell, -1, to, false};
  double travelled = 0;
  Index entered = -1;
  for (int crossing = 0; crossing < 1024; ++crossing) {
    const bool target_inside = contains(cell, to);
    if (target_inside && s.convex[cell])
      return {cell, -1, to, false};
    double best = INFINITY;
    Index face = -1;
    for (Index j = s.faces.offsets[cell]; j < s.faces.offsets[cell + 1]; ++j) {
      const auto f = s.faces.faces[j];
      if (f == entered)
        continue;
      std::vector<Vec3d> p;
      for (auto n : nodes(m, f))
        p.push_back(m.points[n]);
      const auto triangles = face_triangles(f);
      for (std::size_t t = 0; t < triangles.size(); t += 3) {
        const auto a = p[triangles[t]], b = p[triangles[t + 1]],
                   c = p[triangles[t + 2]];
        if (dot(delta, cross(b - a, c - a)) * s.faces.signs[j] <= 0)
          continue;
        const auto hit = hit_triangle(from, delta, a, b, c, 1e-10);
        if (hit && *hit >= travelled - 1e-10 &&
            (*hit < best - 1e-10 ||
             (std::abs(*hit - best) <= 1e-10 && f >= m.n_internal_faces()))) {
          best = *hit;
          face = f;
        }
      }
    }
    if (face < 0)
      return {target_inside ? cell : Index(-1), -1, to, false};
    const auto hit = from + delta * best;
    if (face >= m.n_internal_faces())
      return {cell, face, hit, true};
    travelled = best;
    entered = face;
    cell = m.owner[face] == cell ? m.neighbour[face] : m.owner[face];
  }
  throw Error(
      "post-processing cell walk exceeded 1024 crossings; reduce step_length");
}
} // namespace kutta::post
