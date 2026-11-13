// Snap: each wall point moves to the least-squares intersection of the STL
// tangent planes seen by its wall faces (implicit feature capture: one plane
// on a smooth surface, a line at a feature edge, a point at a corner),
// regularized toward its own nearest surface point. Cells that the snap would
// crush (typically a fine cell protruding past a coarser staircase at a
// refinement transition) are removed first and the snap is recomputed. The
// displacement then decays into the volume, retreating where a cell would
// fall below snap_min_quality (see displace).
#include "internal.hpp"
#include <numeric>

namespace kutta::mesher::internal {
using namespace detail;
namespace {
constexpr int snap_iterations = 4, removal_passes = 12, max_retries = 4;
constexpr double regularization = 0.05, crushed = 0.3, coincident = 0.05;

bool solve3(std::array<std::array<double, 3>, 3> a, std::array<double, 3> b, std::array<double, 3>& x) {
  const auto det3 = [](const std::array<std::array<double, 3>, 3>& m) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  };
  const auto det = det3(a);
  if (!(std::abs(det) > 1e-300)) return false;
  for (int c = 0; c < 3; ++c) {
    auto m = a;
    for (int r = 0; r < 3; ++r) m[r][c] = b[r];
    x[c] = det3(m) / det;
  }
  return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
}
double& at(Vec3d& v, int a) { return a == 0 ? v.x : a == 1 ? v.y : v.z; }
double at(const Vec3d& v, int a) { return a == 0 ? v.x : a == 1 ? v.y : v.z; }
Vec3d unit_normal(const Surface& sf, Index t) {
  const auto& n = sf.triangles[t].nodes;
  const auto v = cross(sf.points[n[1]] - sf.points[n[0]], sf.points[n[2]] - sf.points[n[0]]);
  return v / norm(v);
}
Vec3d closest(const Surface& sf, const SurfaceQuery& q, Vec3d p, Index* triangle = nullptr) {
  const auto s = q.nearest(p);
  const auto& n = sf.triangles[s.triangle].nodes;
  if (triangle) *triangle = s.triangle;
  return closest_on_triangle(p, sf.points[n[0]], sf.points[n[1]], sf.points[n[2]]);
}
bool is_wall(const Face& f) { return f.neighbour < 0 && f.patch >= first_wall; }

struct Targets {
  std::vector<Index> wall_points;
  std::vector<double> h;            // per point: finest adjacent cell size
  std::vector<Vec3d> displacement;  // per point (wall points only)
  std::vector<std::uint8_t> feature;
};

// Iterated: wall faces take the tangent plane at the nearest STL point of
// their current centre; each wall point solves least squares over its faces'
// planes, regularized toward its own nearest STL point. A point whose faces
// see one plane lands on the surface; two distinct planes put it on their
// intersection (a feature edge), three on a corner. Inconsistent planes (for
// example opposite walls of a slot) fall back to the nearest point.
Targets targets(const Work& w, const PointFaces& pf, const Surface& sf, const SurfaceQuery& query) {
  const auto np = index(w.points.size());
  Targets t;
  for (Index p = 0; p < np; ++p)
    for (auto i = pf.offsets[p]; i < pf.offsets[p + 1]; ++i)
      if (is_wall(w.faces[pf.faces[i]])) { t.wall_points.push_back(p); break; }
  t.h.assign(np, std::numeric_limits<double>::infinity());
  for (const auto& f : w.faces)
    for (auto n : f.nodes) { t.h[n] = std::min(t.h[n], w.size[f.owner]); if (f.neighbour >= 0) t.h[n] = std::min(t.h[n], w.size[f.neighbour]); }
  std::vector<Index> wall_faces;
  for (Index f = 0; f < index(w.faces.size()); ++f) if (is_wall(w.faces[f])) wall_faces.push_back(f);
  t.displacement.assign(np, {});
  t.feature.assign(np, 0);
  std::vector<Vec3d> plane_point(w.faces.size()), plane_normal(w.faces.size());
  for (int iteration = 0; iteration < snap_iterations; ++iteration) {
    auto current = [&](Index p) { return w.points[p] + t.displacement[p]; };
    parallel_for(index(wall_faces.size()), [&](Index i) {
      const auto f = wall_faces[i];
      Vec3d c{};
      for (auto n : w.faces[f].nodes) c = c + current(n);
      c = c / double(w.faces[f].nodes.size());
      Index tri = -1;
      plane_point[f] = closest(sf, query, c, &tri);
      plane_normal[f] = unit_normal(sf, tri);
    });
    std::vector<Vec3d> next(np);
    parallel_for(index(t.wall_points.size()), [&](Index i) {
      const auto p = t.wall_points[i];
      const auto origin = w.points[p];
      const auto q = closest(sf, query, current(p));
      std::array<std::array<double, 3>, 3> a{};
      std::array<double, 3> b{};
      std::vector<Index> faces;
      for (auto k = pf.offsets[p]; k < pf.offsets[p + 1]; ++k) if (is_wall(w.faces[pf.faces[k]])) faces.push_back(pf.faces[k]);
      for (int r = 0; r < 3; ++r) { a[r][r] = regularization; b[r] = regularization * at(q, r); }
      double spread = 1;
      for (auto f : faces) {
        const auto n = plane_normal[f];
        const auto d = dot(n, plane_point[f]);
        spread = std::min(spread, dot(n, plane_normal[faces[0]]));
        for (int r = 0; r < 3; ++r) { b[r] += at(n, r) * d; for (int c = 0; c < 3; ++c) a[r][c] += at(n, r) * at(n, c); }
      }
      Vec3d result = q;
      for (int axis = 0; axis < 3; ++axis)
        if (w.planes[p] & (3 << (2 * axis))) { a[axis] = {0, 0, 0}; a[axis][axis] = 1; b[axis] = at(origin, axis); at(result, axis) = at(origin, axis); }
      std::array<double, 3> x{};
      bool feature = false;
      if (solve3(a, b, x)) {
        const Vec3d candidate{x[0], x[1], x[2]};
        bool consistent = norm(candidate - q) < .5 * t.h[p];
        for (auto f : faces) consistent = consistent && std::abs(dot(plane_normal[f], candidate - plane_point[f])) < .1 * t.h[p];
        if (consistent) { result = candidate; feature = spread < .9; }
      }
      t.feature[p] = feature;
      next[p] = result - origin;
    });
    for (auto p : t.wall_points) t.displacement[p] = next[p];
  }
  return t;
}

double volume(const Work& w, const Adjacency& a, const std::vector<Vec3d>& points, Index c) {
  std::vector<FaceGeom> faces;
  Vec3d estimate{};
  for (auto i = a.offsets[c]; i < a.offsets[c + 1]; ++i) {
    faces.push_back(polygon(points, w.faces[a.faces[i]].nodes));
    estimate = estimate + faces.back().centroid;
  }
  estimate = estimate / double(faces.size());
  double v = 0;
  for (std::size_t k = 0; k < faces.size(); ++k) v += a.signs[a.offsets[c] + Index(k)] * dot(faces[k].area, faces[k].centroid - estimate) / 3;
  return v;
}

// Removes cells; faces they shared with kept cells become walls of the zone
// of the nearest facet, oriented out of the kept cell. A removal can cut off
// a pocket of fluid; every cell not connected to the inlet plane (or, without
// one, to the largest region) is removed too, and their count returned.
Index remove_cells(Work& w, const std::vector<char>& removed, const Surface& sf, const SurfaceQuery& query) {
  std::vector<Index> id(w.cells, -1);
  Index cells = 0;
  std::vector<CellKind> kind;
  std::vector<double> size;
  for (Index c = 0; c < w.cells; ++c)
    if (!removed[c]) { id[c] = cells++; kind.push_back(w.kind[c]); size.push_back(w.size[c]); }
  std::vector<Face> faces;
  for (auto& f : w.faces) {
    const bool owner_gone = removed[f.owner], neighbour_gone = f.neighbour >= 0 && removed[f.neighbour];
    if (f.neighbour < 0) { if (!owner_gone) faces.push_back({std::move(f.nodes), id[f.owner], -1, f.patch}); continue; }
    if (owner_gone && neighbour_gone) continue;
    if (!owner_gone && !neighbour_gone) { faces.push_back({std::move(f.nodes), id[f.owner], id[f.neighbour], -1}); continue; }
    Face wall{std::move(f.nodes), id[owner_gone ? f.neighbour : f.owner], -1, -1};
    if (owner_gone) std::reverse(wall.nodes.begin(), wall.nodes.end());
    Vec3d c{};
    for (auto n : wall.nodes) c = c + w.points[n];
    wall.patch = first_wall + int(sf.triangles[query.nearest(c / double(wall.nodes.size())).triangle].zone);
    faces.push_back(std::move(wall));
  }
  w.faces = std::move(faces);
  w.cells = cells;
  w.kind = std::move(kind);
  w.size = std::move(size);

  std::vector<Index> parent(w.cells);
  std::iota(parent.begin(), parent.end(), 0);
  const auto find = [&](Index c) { while (parent[c] != c) c = parent[c] = parent[parent[c]]; return c; };
  for (const auto& f : w.faces) if (f.neighbour >= 0) parent[find(f.owner)] = find(f.neighbour);
  std::vector<Index> count(w.cells);
  for (Index c = 0; c < w.cells; ++c) ++count[find(c)];
  Index keep = -1;
  for (const auto& f : w.faces) if (f.neighbour < 0 && f.patch == 0) { keep = find(f.owner); break; }
  if (keep < 0) keep = Index(std::max_element(count.begin(), count.end()) - count.begin());
  if (count[keep] == w.cells) return 0;
  std::vector<char> cut(w.cells);
  for (Index c = 0; c < w.cells; ++c) cut[c] = find(c) != keep;
  return w.cells - count[keep] + remove_cells(w, cut, sf, query);
}
// Joins wall points that share an edge and snap to (nearly) the same place,
// such as the two edges of a step face between a fine cell and a coarser
// staircase. Degenerate faces are deleted; a merge that would pinch a face
// (repeat a point non-consecutively) is refused.
Index merge_coincident(Work& w, const Targets& t, const PointFaces& pf, const std::vector<char>& keep_apart) {
  const auto np = index(w.points.size());
  auto blocked = keep_apart;
  for (int attempt = 0; attempt < 8; ++attempt) {
    std::vector<Index> parent(np);
    std::iota(parent.begin(), parent.end(), 0);
    const auto find = [&](Index p) { while (parent[p] != p) p = parent[p] = parent[parent[p]]; return p; };
    const auto incident = [&](Index p) { return pf.offsets[p + 1] - pf.offsets[p]; };
    Index merged = 0;
    for (const auto& f : w.faces) {
      if (!is_wall(f)) continue;
      for (std::size_t j = 0; j < f.nodes.size(); ++j) {
        const auto a = f.nodes[j], b = f.nodes[(j + 1) % f.nodes.size()];
        if (blocked[a] || blocked[b]) continue;
        if (norm(w.points[a] + t.displacement[a] - w.points[b] - t.displacement[b]) >= coincident * std::min(t.h[a], t.h[b])) continue;
        auto ra = find(a), rb = find(b);
        if (ra == rb) continue;
        if (incident(rb) > incident(ra)) std::swap(ra, rb);  // keep the better connected point
        parent[rb] = ra;
        ++merged;
      }
    }
    if (!merged) return 0;
    std::vector<std::vector<Index>> mapped(w.faces.size());
    bool pinched = false;
    for (std::size_t i = 0; i < w.faces.size(); ++i) {
      auto& m = mapped[i];
      for (auto n : w.faces[i].nodes) { const auto r = find(n); if (m.empty() || m.back() != r) m.push_back(r); }
      while (m.size() > 1 && m.front() == m.back()) m.pop_back();
      auto sorted = m;
      std::sort(sorted.begin(), sorted.end());
      if (m.size() >= 3 && std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
        pinched = true;
        for (auto n : w.faces[i].nodes) blocked[n] = 1;
      }
    }
    if (pinched) continue;
    std::vector<Face> faces;
    for (std::size_t i = 0; i < w.faces.size(); ++i)
      if (mapped[i].size() >= 3) { w.faces[i].nodes = std::move(mapped[i]); faces.push_back(std::move(w.faces[i])); }
    w.faces = std::move(faces);
    return merged;
  }
  return 0;
}
}  // namespace

Json::Object snap(Work& w, const Surface& sf, const SurfaceQuery& query, const Settings&, const Progress& progress) {
  // A merge can leave the unsnapped start invalid at a corner (a cell still
  // below quality with every scale at zero). The snap is then redone from the
  // castellated mesh without merges at that cell's points; only if that makes
  // no progress are such cells removed.
  const Work castellated = w;
  std::vector<char> keep_apart(w.points.size());
  Targets t;
  Index removed_cells = 0, merged_points = 0, disconnected = 0, stuck_removed = 0;
  int passes = 0, retries = 0;
  std::vector<double> scale;
  DisplaceResult recovery;
  for (;;) {
    w = castellated;
    removed_cells = merged_points = disconnected = 0;
    for (passes = 0;; ++passes) {
      const auto pf = point_faces(w);
      t = targets(w, pf, sf, query);
      if (progress) progress("trim", .55, "Snapping wall points to surface and feature edges");
      if (passes == removal_passes) {
        // Out of passes: still drop any cell a merge reduced below four faces.
        const auto adjacency = cell_faces(w);
        std::vector<char> removed(w.cells);
        Index count = 0;
        for (Index c = 0; c < w.cells; ++c) if (adjacency.offsets[c + 1] - adjacency.offsets[c] < 4) { removed[c] = 1; ++count; }
        if (!count) break;
        removed_cells += count;
        disconnected += remove_cells(w, removed, sf, query);
        t = targets(w, point_faces(w), sf, query);
        break;
      }
      const auto merged = merge_coincident(w, t, pf, keep_apart);
      merged_points += merged;
      if (merged) continue;  // faces changed: recompute targets first
      const auto adjacency = cell_faces(w);
      std::vector<Vec3d> trial = w.points;
      std::vector<char> touched(w.cells), removed(w.cells);
      for (auto p : t.wall_points) {
        trial[p] = w.points[p] + t.displacement[p];
        for (auto i = pf.offsets[p]; i < pf.offsets[p + 1]; ++i) {
          const auto& f = w.faces[pf.faces[i]];
          touched[f.owner] = 1;
          if (f.neighbour >= 0) touched[f.neighbour] = 1;
        }
      }
      parallel_for(w.cells, [&](Index c) {
        if (adjacency.offsets[c + 1] - adjacency.offsets[c] < 4) removed[c] = 1;
        else if (touched[c]) removed[c] = volume(w, adjacency, trial, c) < crushed * volume(w, adjacency, w.points, c);
      });
      const auto count = std::count(removed.begin(), removed.end(), 1);
      if (!count) break;
      removed_cells += index(std::size_t(count));
      disconnected += remove_cells(w, removed, sf, query);
    }
    const auto adjacency = cell_faces(w);
    const auto pf = point_faces(w);
    scale.assign(w.points.size(), 1.);
    const auto start = w.points;
    recovery = displace(w, adjacency, pf, t.wall_points, t.displacement, scale, snap_min_quality, progress, "trim");
    if (recovery.stuck.empty()) break;
    Index fresh = 0;
    for (auto c : recovery.stuck)
      for (auto i = adjacency.offsets[c]; i < adjacency.offsets[c + 1]; ++i)
        for (auto p : w.faces[adjacency.faces[i]].nodes) if (!keep_apart[p]) { keep_apart[p] = 1; ++fresh; }
    if (fresh && retries < max_retries) { ++retries; continue; }
    // Last resort: remove the cells and snap what remains.
    for (int attempt = 0; !recovery.stuck.empty(); ++attempt) {
      KUTTA_REQUIRE(attempt < 10, "mesh: snap could not recover cell quality; refine the surface or report the geometry");
      w.points = start;
      std::vector<char> removed(w.cells);
      for (auto c : recovery.stuck) removed[c] = 1;
      stuck_removed += index(recovery.stuck.size());
      disconnected += remove_cells(w, removed, sf, query);
      t = targets(w, point_faces(w), sf, query);
      const auto a = cell_faces(w);
      scale.assign(w.points.size(), 1.);
      recovery = displace(w, a, point_faces(w), t.wall_points, t.displacement, scale, snap_min_quality, progress, "trim");
    }
    break;
  }

  // Geometric deviation of the snapped wall.
  std::vector<double> deviation, relative;
  Index full = 0, partial = 0, none = 0, features = 0;
  for (auto p : t.wall_points) {
    const auto d = query.nearest(w.points[p]).distance;
    deviation.push_back(d);
    relative.push_back(d / t.h[p]);
    full += scale[p] == 1; none += scale[p] == 0; partial += scale[p] > 0 && scale[p] < 1;
    features += t.feature[p] && scale[p] == 1;
  }
  auto stats = [](std::vector<double> v) {
    if (v.empty()) return Json::Object{};
    std::sort(v.begin(), v.end());
    const auto pick = [&](double q) { return v[std::min(v.size() - 1, std::size_t(q * double(v.size() - 1) + .5))]; };
    return Json::Object{{"mean", std::accumulate(v.begin(), v.end(), 0.) / double(v.size())}, {"p95", pick(.95)}, {"p99", pick(.99)}, {"max", v.back()}};
  };
  return Json::Object{
    {"wall_points", double(t.wall_points.size())}, {"fully_snapped", double(full)}, {"partially_snapped", double(partial)}, {"unsnapped", double(none)},
    {"feature_points", double(features)}, {"crushed_cells_removed", double(removed_cells)}, {"stuck_cells_removed", double(stuck_removed)}, {"merge_retries", retries}, {"disconnected_cells_removed", double(disconnected)}, {"merged_points", double(merged_points)}, {"removal_passes", passes},
    {"cells_below_quality_before_recovery", double(recovery.first_bad)}, {"recovery_rounds", recovery.rounds},
    {"min_quality", snap_min_quality},
    {"point_distance_m", stats(deviation)}, {"point_distance_over_cell_size", stats(relative)}};
}
}  // namespace kutta::mesher::internal
