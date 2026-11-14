// External-aero volume mesher: balanced octree -> castellated hex mesh ->
// snapped wall -> optional prism layers (snappyHexMesh-style). Cells cut by
// the STL are removed rather than clipped, so no cut-cell fragments exist;
// the snap stage moves the staircase wall onto the surface.
#include "internal.hpp"
#include "kutta/mesh/topology.hpp"
#include "kutta/core/memory.hpp"
#include <chrono>
#include <map>
#include <numeric>
#include <unordered_map>

namespace kutta::mesher {
using namespace detail;
using namespace internal;
namespace {
using Grid = std::array<int, 3>;
struct GridHash {
  std::size_t operator()(const Grid& k) const {
    std::size_t h = 0;
    for (auto v : k) h = (h * 1000003) ^ std::hash<int>{}(v);
    return h;
  }
};
struct Octant { Grid lo; int size, level, child = -1; Sample sample{}; };

// Per-facet target size: curvature (radius * angle), and gaps to any facet or
// the ground along the outward normal (slots, ground clearance).
std::vector<double> facet_sizes(const Surface& sf, const SurfaceQuery& query, const Settings& st) {
  const double angle = st.curvature_degrees * std::acos(-1.) / 180, floor = st.proximity_floor();
  std::vector<double> out(sf.triangles.size());
  parallel_for(index(sf.triangles.size()), [&](Index t) {
    const auto& tri = sf.triangles[t];
    const auto a = sf.points[tri.nodes[0]], b = sf.points[tri.nodes[1]], c = sf.points[tri.nodes[2]];
    auto n = cross(b - a, c - a); n = n / norm(n);
    const auto centroid = (a + b + c) / 3.;
    double size = st.max_size;
    // Facets buried inside another shell (a wing tip inside its endplate)
    // are not walls and do not drive sizing.
    if (query.sample(centroid + n * (1e-6 * norm(sf.bounds.max - sf.bounds.min))).inside) { out[t] = size; return; }
    if (std::isfinite(tri.radius)) size = std::min(size, std::max(st.min_size, tri.radius * angle));
    auto gap = query.ray(centroid, n, t);
    if (st.ground && n.z < -1e-3) gap = std::min(gap, std::max(0., centroid.z) / -n.z);
    // Castellation removes about one cell at each wall of a gap.
    if (std::isfinite(gap)) size = std::min(size, std::max(floor, gap / (st.proximity_cells + 1.5)));
    out[t] = size;
  });
  return out;
}

struct Tree {
  const Surface& surface; const Settings& settings; const SurfaceQuery& query;
  Box bounds; double unit; int root_size; Grid counts;
  std::vector<Octant> nodes;
  Index leaves = 0;
  std::vector<char> layered;  // per STL zone
  Tree(const Surface& sf, const Settings& st, const SurfaceQuery& q) : surface(sf), settings(st), query(q), layered(sf.zones.size()) {
    if (st.layer_count > 0)
      for (const auto& name : st.layer_zones)
        for (std::size_t z = 0; z < sf.zones.size(); ++z) if (name == "*" || name == sf.zones[z]) layered[z] = 1;
    const auto length = sf.bounds.max.x - sf.bounds.min.x;
    KUTTA_REQUIRE(length > 0, "mesh: body must have positive x length");
    bounds = {{sf.bounds.min.x - st.upstream * length, st.symmetry_y ? 0 : sf.bounds.min.y - st.side * length, st.ground ? 0 : sf.bounds.min.z - st.side * length},
              {sf.bounds.max.x + st.downstream * length, sf.bounds.max.y + st.side * length, sf.bounds.max.z + st.top * length}};
    KUTTA_REQUIRE(!st.ground || sf.bounds.min.z >= -1e-7 * length, "mesh: ground is z=0 but geometry extends below it; translate the STL or disable ground");
    KUTTA_REQUIRE(!st.symmetry_y || sf.bounds.max.y > 0, "mesh: y symmetry retains y>=0, but geometry is wholly outside it");
    auto minimum = std::min(st.min_size, st.proximity_floor());
    if (st.layer_count > 0) minimum = std::min(minimum, layer_transition * st.last_layer());
    for (const auto& b : st.boxes) minimum = std::min(minimum, b.size);
    int depth = 0;
    while (st.base_size / std::pow(2., depth) > minimum * (1 + 1e-12)) ++depth;
    KUTTA_REQUIRE(depth <= 16, "mesh: requested size range exceeds 16 octree levels");
    root_size = 1 << (depth + 1);
    unit = st.base_size / root_size;  // extra bit for leaf centres
    const auto extent = (bounds.max - bounds.min) / st.base_size;
    KUTTA_REQUIRE(extent.x > 0 && extent.y > 0 && extent.z > 0 && extent.x < 100000 && extent.y < 100000 && extent.z < 100000, "mesh: invalid tunnel dimensions");
    counts = {int(std::ceil(extent.x)), int(std::ceil(extent.y)), int(std::ceil(extent.z))};
    for (int a = 0; a < 3; ++a)
      KUTTA_REQUIRE(std::int64_t(counts[a]) * root_size < std::numeric_limits<int>::max() - 1, "mesh: requested tunnel/size range exceeds integer lattice capacity");
    const auto roots = std::int64_t(counts[0]) * counts[1] * counts[2];
    KUTTA_REQUIRE(roots <= st.max_cells, "mesh/max_cells: tunnel base grid alone exceeds the cell budget; increase base_size");
    bounds.max = bounds.min + Vec3d{double(counts[0]), double(counts[1]), double(counts[2])} * st.base_size;
    nodes.reserve(std::size_t(roots));
    for (int k = 0; k < counts[2]; ++k)
      for (int j = 0; j < counts[1]; ++j)
        for (int i = 0; i < counts[0]; ++i) nodes.push_back({{i * root_size, j * root_size, k * root_size}, root_size, 0, -1, {}});
    leaves = Index(roots);
  }
  Vec3d position(Grid g) const { return bounds.min + Vec3d{double(g[0]), double(g[1]), double(g[2])} * unit; }
  Vec3d center(const Octant& n) const { return position({n.lo[0] + n.size / 2, n.lo[1] + n.size / 2, n.lo[2] + n.size / 2}); }
  int locate(Grid p) const {
    for (int a = 0; a < 3; ++a) if (p[a] < 0 || p[a] >= counts[a] * root_size) return -1;
    int id = (p[2] / root_size * counts[1] + p[1] / root_size) * counts[0] + p[0] / root_size;
    while (nodes[id].child >= 0) {
      const auto& n = nodes[id];
      int c = 0;
      for (int a = 0; a < 3; ++a) if (p[a] >= n.lo[a] + n.size / 2) c |= 1 << a;
      id = n.child + c;
    }
    return id;
  }
  void split(int id) {
    const auto n = nodes[id];
    KUTTA_REQUIRE(n.size > 2, "mesh: balance exceeded minimum octree size");
    KUTTA_REQUIRE(leaves <= settings.max_cells - 7, "mesh/max_cells: octree budget exceeded; increase sizes or max_cells");
    const auto first = int(nodes.size());
    nodes[id].child = first;
    for (int i = 0; i < 8; ++i) {
      auto p = n.lo;
      for (int a = 0; a < 3; ++a) p[a] += (i >> a & 1) * n.size / 2;
      nodes.push_back({p, n.size / 2, n.level + 1, -1, {}});
    }
    leaves += 7;
  }
  double target(const Octant& n, const std::vector<double>& facet) const {
    const auto p = center(n);
    const auto h = n.size * unit;
    double t = settings.base_size;
    if (n.sample.distance < 1.75 * h) {
      const auto& nearest = surface.triangles[n.sample.triangle];
      t = facet[n.sample.triangle];
      if (layered[std::size_t(nearest.zone)]) t = std::min(t, layer_transition * settings.last_layer());
      for (int e = 0; e < 3; ++e)
        if (nearest.feature_edges[e]) {
          const auto a = surface.points[nearest.nodes[e]], b = surface.points[nearest.nodes[(e + 1) % 3]], edge = b - a;
          const auto s = std::clamp(dot(p - a, edge) / dot(edge, edge), 0., 1.);
          if (norm(p - (a + edge * s)) < 1.75 * h) t = std::min(t, settings.min_size);
        }
    }
    const auto a = position(n.lo), b = a + Vec3d{h, h, h};
    for (const auto& box : settings.boxes)
      if (a.x < box.bounds.max.x && b.x > box.bounds.min.x && a.y < box.bounds.max.y && b.y > box.bounds.min.y && a.z < box.bounds.max.z && b.z > box.bounds.min.z)
        t = std::min(t, box.size);
    // Grow at most about one cell size per cell away from the surface.
    return std::min(t, std::max(settings.max_size, n.sample.distance * .6));
  }
  void refine(const std::vector<double>& facet, const Progress& progress) {
    // Breadth-first by level so each level's samples run in parallel.
    std::vector<int> level;
    for (int i = 0; i < int(nodes.size()); ++i) level.push_back(i);
    while (!level.empty()) {
      parallel_for(index(level.size()), [&](Index k) { auto& n = nodes[level[k]]; n.sample = query.nearest(center(n)); });
      std::vector<int> next;
      for (auto id : level) {
        if (nodes[id].size <= 2 || nodes[id].size * unit <= target(nodes[id], facet) * (1 + 1e-12)) continue;
        split(id);
        for (int c = 0; c < 8; ++c) next.push_back(nodes[id].child + c);
      }
      level = std::move(next);
      if (progress) progress("octree", .2, "Refining surface, gaps and volume boxes");
    }
    // Each leaf asks every face/edge/corner neighbour for at most a 2:1 jump.
    for (bool changed = true; changed;) {
      changed = false;
      for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i].child >= 0) continue;
        const auto n = nodes[i];
        for (int x = -1; x <= 1; ++x) for (int y = -1; y <= 1; ++y) for (int z = -1; z <= 1; ++z) {
          if (x == 0 && y == 0 && z == 0) continue;
          const Grid direction{x, y, z};
          auto p = n.lo;
          for (int a = 0; a < 3; ++a) p[a] += direction[a] < 0 ? -1 : direction[a] > 0 ? n.size : n.size / 2;
          const auto j = locate(p);
          if (j >= 0 && nodes[j].size > 2 * n.size) { split(j); changed = true; }
        }
      }
    }
    std::vector<int> all;
    for (int i = 0; i < int(nodes.size()); ++i) if (nodes[i].child < 0) all.push_back(i);
    parallel_for(index(all.size()), [&](Index k) { auto& n = nodes[all[k]]; n.sample = query.sample(center(n)); });
  }
  int max_level_jump() const {
    int jump = 0;
    for (const auto& n : nodes) {
      if (n.child >= 0) continue;
      for (int x = -1; x <= 1; ++x) for (int y = -1; y <= 1; ++y) for (int z = -1; z <= 1; ++z) {
        if (x == 0 && y == 0 && z == 0) continue;
        const Grid direction{x, y, z};
        auto p = n.lo;
        for (int a = 0; a < 3; ++a) p[a] += direction[a] < 0 ? -1 : direction[a] > 0 ? n.size : n.size / 2;
        const auto j = locate(p);
        if (j >= 0) jump = std::max(jump, std::abs(n.level - nodes[j].level));
      }
    }
    return jump;
  }
  // Leaves across one side at its four quadrant centres (-1 outside).
  std::array<int, 4> across(const Octant& n, int axis, int side) const {
    std::array<int, 4> out{};
    const int u = (axis + 1) % 3, v = (axis + 2) % 3;
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 2; ++i) {
        Grid p = n.lo;
        p[axis] = side ? n.lo[axis] + n.size : n.lo[axis] - 1;
        p[u] = n.lo[u] + (2 * i + 1) * n.size / 4;
        p[v] = n.lo[v] + (2 * j + 1) * n.size / 4;
        out[j * 2 + i] = locate(p);
      }
    return out;
  }
};

enum Leaf : std::int8_t { Interior = -1, Removed = 0, Kept = 1, Fluid = 2 };

// Fluid cells crossed by each facet's outward normal ray across a narrow gap
// (narrower than four max surface sizes) to another shell or to the ground,
// summarized per zone pair. Zero cells means the gap closed.
Json::Array gap_report(const Tree& tree, const SurfaceQuery& query, const Surface& sf, const std::vector<Leaf>& state) {
  const auto& st = tree.settings;
  const auto limit = 4 * st.max_size;
  struct Gap { int cells; double width; Index other; };  // other: zone, or -1 ground
  std::vector<Gap> gaps(sf.triangles.size(), {-1, 0, -2});
  const auto scale = 1e-6 * norm(sf.bounds.max - sf.bounds.min);
  parallel_for(index(sf.triangles.size()), [&](Index t) {
    const auto& tri = sf.triangles[t];
    const auto a = sf.points[tri.nodes[0]], b = sf.points[tri.nodes[1]], c = sf.points[tri.nodes[2]];
    auto n = cross(b - a, c - a); n = n / norm(n);
    const auto centroid = (a + b + c) / 3.;
    if (query.sample(centroid + n * scale).inside) return;
    Index hit = -1;
    auto width = query.ray(centroid, n, t, &hit);
    Index other = -2;
    if (hit >= 0 && sf.triangles[hit].component != tri.component) other = sf.triangles[hit].zone;
    else width = std::numeric_limits<double>::infinity();
    if (st.ground && n.z < -1e-3 && centroid.z / -n.z < width) { width = centroid.z / -n.z; other = -1; }
    if (!(width < limit) || other == -2) return;
    // Distinct fluid leaves along the ray at quarter finest-cell steps.
    int cells = 0, last = -1;
    const auto step = tree.unit * .5;
    for (double s = step; s < width; s += step) {
      const auto p = centroid + n * s - tree.bounds.min;
      const auto id = tree.locate({int(std::floor(p.x / tree.unit)), int(std::floor(p.y / tree.unit)), int(std::floor(p.z / tree.unit))});
      if (id >= 0 && id != last && state[id] == Fluid) ++cells;
      if (id >= 0) last = id;
    }
    gaps[t] = {cells, width, other};
  });
  std::map<std::pair<std::string, std::string>, Gap> pairs;
  for (std::size_t t = 0; t < gaps.size(); ++t) {
    if (gaps[t].other == -2) continue;
    auto key = std::make_pair(sf.zones[sf.triangles[t].zone], gaps[t].other < 0 ? std::string("ground") : sf.zones[gaps[t].other]);
    if (key.second != "ground" && key.second < key.first) std::swap(key.first, key.second);
    const auto it = pairs.find(key);
    if (it == pairs.end() || gaps[t].cells < it->second.cells || (gaps[t].cells == it->second.cells && gaps[t].width < it->second.width)) pairs[key] = gaps[t];
  }
  Json::Array out;
  for (const auto& [key, g] : pairs) out.push_back(Json::Object{{"between", Json::Array{key.first, key.second}}, {"min_cells_across", g.cells}, {"width_m", g.width}});
  return out;
}

// Keeps leaves that are outside the body and not crossed by any facet, and
// connected to the far field. Cells are octree cubes; a side facing a finer
// neighbour is split into four quads (2:1 hanging faces), and edge midpoints
// used by any finer cell are inserted into coarser polygons.
Work castellate(const Tree& tree, const SurfaceQuery& query, const Surface& surface, Json::Object& stats) {
  const auto& nodes = tree.nodes;
  std::vector<Leaf> state(nodes.size(), Interior);
  std::vector<int> leaves;
  for (int i = 0; i < int(nodes.size()); ++i) if (nodes[i].child < 0) leaves.push_back(i);
  parallel_for(index(leaves.size()), [&](Index k) {
    const auto& n = nodes[leaves[k]];
    const double h = n.size * tree.unit;
    bool removed = n.sample.inside;
    if (!removed && n.sample.distance < 0.8660254037844387 * h * (1 + 1e-9)) {
      const auto lo = tree.position(n.lo);
      const Vec3d e{h * 1e-6, h * 1e-6, h * 1e-6};
      removed = query.crosses({lo + e, lo + Vec3d{h, h, h} - e});
    }
    state[leaves[k]] = removed ? Removed : Kept;
  });
  const auto start = tree.locate({0, 0, 0});
  KUTTA_REQUIRE(start >= 0 && state[start] == Kept, "mesh: tunnel corner is inside the body; enlarge the tunnel");
  std::vector<int> queue{start};
  state[start] = Fluid;
  for (std::size_t q = 0; q < queue.size(); ++q)
    for (int a = 0; a < 3; ++a)
      for (int side = 0; side < 2; ++side)
        for (auto id : tree.across(nodes[queue[q]], a, side))
          if (id >= 0 && state[id] == Kept) { state[id] = Fluid; queue.push_back(id); }
  Index isolated = 0, removed = 0;
  for (auto i : leaves) { if (state[i] == Kept) { state[i] = Removed; ++isolated; } removed += state[i] == Removed; }
  stats["removed_cells"] = double(removed);
  stats["isolated_cells_removed"] = double(isolated);

  stats["gaps"] = gap_report(tree, query, surface, state);
  Work w;
  w.bounds = tree.bounds;
  std::vector<Index> cell(nodes.size(), -1);
  for (auto i : leaves)
    if (state[i] == Fluid) { cell[i] = w.cells++; w.size.push_back(nodes[i].size * tree.unit); w.kind.push_back(Hex); }
  KUTTA_REQUIRE(w.cells > 0, "mesh: no fluid cells remain");

  // patch: -1 internal, [0,6) tunnel plane, -2 wall (zone resolved below).
  auto visit = [&](auto&& emit) {
    for (auto i : leaves) {
      if (state[i] != Fluid) continue;
      const auto& n = nodes[i];
      for (int axis = 0; axis < 3; ++axis)
        for (int side = 0; side < 2; ++side) {
          const auto ids = tree.across(n, axis, side);
          const bool split = ids[0] >= 0 && nodes[ids[0]].size < n.size;
          const int count = split ? 2 : 1, s = n.size / count, u = (axis + 1) % 3, v = (axis + 2) % 3;
          for (int j = 0; j < count; ++j)
            for (int k = 0; k < count; ++k) {
              const auto id = ids[j * 2 + k];
              Index neighbour = -1;
              int patch = -1;
              if (id < 0) patch = 2 * axis + side;
              else if (state[id] == Removed) patch = -2;
              else {
                if (nodes[id].size < n.size || (nodes[id].size == n.size && id < i)) continue;
                neighbour = cell[id];
              }
              Grid a = n.lo;
              a[axis] += side * n.size; a[u] += k * s; a[v] += j * s;
              Grid b = a, c = a, d = a;
              b[u] += s; c[u] += s; c[v] += s; d[v] += s;
              std::array<Grid, 4> quad = side ? std::array<Grid, 4>{a, b, c, d} : std::array<Grid, 4>{a, d, c, b};
              emit(quad, cell[i], neighbour, patch);
            }
        }
    }
  };
  std::unordered_map<Grid, Index, GridHash> vertex;
  visit([&](const std::array<Grid, 4>& q, Index, Index, int) { for (const auto& g : q) vertex.emplace(g, -1); });
  vertex.rehash(vertex.size());
  auto id = [&](const Grid& g) {
    auto& v = vertex.at(g);
    if (v < 0) {
      v = index(w.points.size());
      w.points.push_back(tree.position(g));
      std::uint8_t planes = 0;
      for (int a = 0; a < 3; ++a) {
        if (g[a] == 0) planes |= 1 << (2 * a);
        if (g[a] == tree.counts[a] * tree.root_size) planes |= 1 << (2 * a + 1);
      }
      w.planes.push_back(planes);
    }
    return v;
  };
  visit([&](const std::array<Grid, 4>& q, Index owner, Index neighbour, int patch) {
    Face f{{}, owner, neighbour, patch};
    for (int j = 0; j < 4; ++j) {
      f.nodes.push_back(id(q[j]));
      Grid mid;
      for (int a = 0; a < 3; ++a) mid[a] = (q[j][a] + q[(j + 1) % 4][a]) / 2;
      if (vertex.contains(mid)) f.nodes.push_back(id(mid));
    }
    w.faces.push_back(std::move(f));
  });
  parallel_for(index(w.faces.size()), [&](Index i) {
    auto& f = w.faces[i];
    if (f.patch != -2) return;
    Vec3d c{};
    for (auto n : f.nodes) c = c + w.points[n];
    f.patch = first_wall + int(surface.triangles[query.nearest(c / double(f.nodes.size())).triangle].zone);
  });
  return w;
}

std::vector<Patch> patch_list(const Surface& surface, const Settings& st) {
  std::vector<Patch> patches{{"inlet", PatchKind::VelocityInlet, "velocity-inlet", 3}, {"outlet", PatchKind::PressureOutlet, "pressure-outlet", 4},
    {st.symmetry_y ? "symmetry" : "side-low", PatchKind::Symmetry, "symmetry", 5}, {"side-high", PatchKind::Symmetry, "symmetry", 6},
    {st.ground ? "ground" : "bottom", st.ground ? PatchKind::Wall : PatchKind::Symmetry, st.ground ? "wall" : "symmetry", 7}, {"top", PatchKind::Symmetry, "symmetry", 8}};
  for (const auto& name : surface.zones) {
    KUTTA_REQUIRE(std::none_of(patches.begin(), patches.end(), [&](const Patch& p) { return p.name == name; }), "STL zone name collides with tunnel boundary: " + name);
    patches.push_back({name, PatchKind::Wall, "wall", int(patches.size() + 3)});
  }
  return patches;
}

// Drops empty zones (for example a wheel wholly outside a symmetry model) and
// unused points, then restores PolyMesh face ordering.
PolyMesh finish(Work w, std::vector<Patch> patches) {
  std::vector<Index> counts(patches.size());
  for (const auto& f : w.faces) if (f.neighbour < 0) ++counts[f.patch];
  std::vector<Index> map(patches.size(), -1);
  std::vector<Patch> used;
  for (std::size_t i = 0; i < patches.size(); ++i) if (counts[i]) { map[i] = index(used.size()); used.push_back(patches[i]); }
  for (auto& f : w.faces) if (f.neighbour < 0) f.patch = map[f.patch];
  std::vector<Index> remap(w.points.size(), -1);
  std::vector<Vec3d> points;
  for (auto& f : w.faces)
    for (auto& n : f.nodes) {
      if (remap[n] < 0) { remap[n] = index(points.size()); points.push_back(w.points[n]); }
      n = remap[n];
    }
  return pack(std::move(points), w.cells, std::move(w.faces), std::move(used));
}

Json vec(Vec3d p) { return Json::Array{p.x, p.y, p.z}; }
double percentile(std::vector<double> v, double p) {
  if (v.empty()) return 0;
  const auto k = std::min(v.size() - 1, std::size_t(p * double(v.size() - 1) + .5));
  std::nth_element(v.begin(), v.begin() + std::ptrdiff_t(k), v.end());
  return v[k];
}
Json::Array histogram(const std::vector<double>& values, const std::vector<double>& edges) {
  Json::Array out;
  for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
    const auto n = std::count_if(values.begin(), values.end(), [&](double v) { return v >= edges[i] && (v < edges[i + 1] || (i + 2 == edges.size() && v <= edges[i + 1])); });
    out.push_back(Json::Object{{"from", edges[i]}, {"to", edges[i + 1]}, {"cells", double(n)}});
  }
  return out;
}

Json::Object quality_report(const PolyMesh& mesh, const MeshGeometry& geometry, const std::vector<double>& quality, const std::vector<CellKind>& kind) {
  const auto worst = std::min_element(quality.begin(), quality.end());
  const auto worst_cell = Index(worst - quality.begin());
  double max_skew = 0;
  std::vector<double> ratio_core, ratio_layer;
  for (Index f = 0; f < mesh.n_internal_faces(); ++f) {
    const auto o = mesh.owner[f], n = mesh.neighbour[f];
    const auto link = geometry.cell_centroid[n] - geometry.cell_centroid[o];
    const auto t = dot(geometry.face_centroid[f] - geometry.cell_centroid[o], geometry.face_area[f]) / dot(link, geometry.face_area[f]);
    max_skew = std::max(max_skew, norm(geometry.face_centroid[f] - (geometry.cell_centroid[o] + link * t)) / std::max(norm(link), 1e-30));
    const auto r = std::max(geometry.cell_volume[o] / geometry.cell_volume[n], geometry.cell_volume[n] / geometry.cell_volume[o]);
    (kind[o] == Prism || kind[n] == Prism ? ratio_layer : ratio_core).push_back(r);
  }
  double wall_ratio = std::numeric_limits<double>::infinity(), area_ratio = std::numeric_limits<double>::infinity();
  for (const auto& p : mesh.patches)
    if (p.kind == PatchKind::Wall)
      for (Index f = p.start; f < p.start + p.size; ++f) {
        const auto o = mesh.owner[f];
        const auto r = std::abs(dot(geometry.face_centroid[f] - geometry.cell_centroid[o], geometry.face_area[f])) / norm(geometry.face_area[f]) / std::cbrt(geometry.cell_volume[o]);
        wall_ratio = std::min(wall_ratio, r);
      }
  for (Index f = 0; f < mesh.n_faces(); ++f) {
    const auto r = norm(geometry.face_area[f]) / std::pow(geometry.cell_volume[mesh.owner[f]], 2. / 3);
    area_ratio = std::min(area_ratio, r);
  }
  const auto below = [&](double q) { return double(std::count_if(quality.begin(), quality.end(), [&](double v) { return v < q; })); };
  Index nonconvex = 0;
  const auto cell_faces = build_cell_faces(mesh);
  for (Index c = 0; c < mesh.n_cells; ++c) {
    bool convex = true;
    for (auto i = cell_faces.offsets[c]; i < cell_faces.offsets[c + 1] && convex; ++i) {
      const auto f = cell_faces.faces[i];
      const auto normal = geometry.face_area[f] * cell_faces.signs[i];
      const auto tolerance = 1e-3 * std::cbrt(geometry.cell_volume[c]) * norm(normal);
      for (auto j = cell_faces.offsets[c]; j < cell_faces.offsets[c + 1] && convex; ++j)
        for (auto v : nodes(mesh, cell_faces.faces[j]))
          if (dot(normal, mesh.points[v] - geometry.face_centroid[f]) > tolerance) { convex = false; break; }
    }
    nonconvex += !convex;
  }
  const auto ratios = [&](const std::vector<double>& r) {
    return Json::Object{{"faces", double(r.size())}, {"max", r.empty() ? 1. : *std::max_element(r.begin(), r.end())},
      {"p50", percentile(r, .5)}, {"p99", percentile(r, .99)}, {"p999", percentile(r, .999)},
      {"histogram", histogram(r, {1, 2, 4, 8, 16, 20, 50, 1e300})}};
  };
  auto all = ratio_core;
  all.insert(all.end(), ratio_layer.begin(), ratio_layer.end());
  return Json::Object{
    {"orthogonal_min", *worst},
    {"orthogonal_mean", std::accumulate(quality.begin(), quality.end(), 0.) / double(quality.size())},
    {"cells_below_0_05", below(.05)}, {"cells_below_0_15", below(.15)},
    {"fraction_below_0_15", below(.15) / double(quality.size())},
    {"orthogonal_histogram", histogram(quality, {0, .05, .1, .15, .2, .3, .5, .7, .9, 1})},
    {"worst_cell", Json::Object{{"index", double(worst_cell)}, {"centroid", vec(geometry.cell_centroid[worst_cell])}, {"volume", geometry.cell_volume[worst_cell]}}},
    {"skewness_max", max_skew},
    {"adjacent_volume_ratio_max", all.empty() ? 1. : *std::max_element(all.begin(), all.end())},
    {"adjacent_volume_ratio_core", ratios(ratio_core)},
    {"adjacent_volume_ratio_layers", ratios(ratio_layer)},
    {"nonconvex_cells", double(nonconvex)}, {"nonconvex_tolerance", "1e-3 of cell size"},
    {"negative_cells", 0}, {"wall_distance_over_cell_size_min", wall_ratio}, {"face_area_over_cell_area_min", area_ratio},
    {"volume_min", *std::min_element(geometry.cell_volume.begin(), geometry.cell_volume.end())}};
}
}  // namespace

Result generate(const Surface& surface, const Settings& settings, double speed, const Progress& progress) {
  const auto start = std::chrono::steady_clock::now();
  Json::Object timing;
  auto lap = [&, last = start](const char* name) mutable {
    const auto now = std::chrono::steady_clock::now();
    timing[name] = std::chrono::duration<double>(now - last).count();
    last = now;
  };
  const auto emit = [&](std::string_view stage, double f, std::string_view m) { if (progress) progress(stage, f, m); };
  KUTTA_REQUIRE(std::isfinite(speed) && speed > 0, "mesh: speed must be positive and finite");
  SurfaceQuery query(surface);
  emit("octree", .15, "Sizing facets by curvature, gaps and features");
  const auto facet = facet_sizes(surface, query, settings);
  Tree tree(surface, settings, query);
  tree.refine(facet, progress);
  const auto level_jump = tree.max_level_jump();
  KUTTA_REQUIRE(level_jump <= 1, "mesh: octree 2:1 balance invariant failed");
  lap("octree");
  emit("trim", .4, "Castellating: removing cells inside or crossed by the surface");
  Json::Object castellation;
  auto work = castellate(tree, query, surface, castellation);
  const auto octree_cells = work.cells;
  lap("castellate");
  emit("trim", .5, "Snapping the wall onto the surface and its feature edges");
  auto snap_report = snap(work, surface, query, settings, progress);
  lap("snap");
  emit("layers", .75, settings.layer_count ? "Inserting prism layers" : "No prism layers requested");
  Json::Object layer_report{{"requested", settings.layer_count}, {"generated", 0}, {"status", "no layers requested; wall-function sizing"}};
  if (settings.layer_count) layer_report = add_layers(work, surface, settings, progress);
  KUTTA_REQUIRE(work.cells <= settings.max_cells, "mesh/max_cells: cell budget exceeded after layers; increase sizes or max_cells");
  lap("layers");
  emit("quality", .85, "Checking closure, volumes and orthogonal quality");
  const auto kind = work.kind;
  auto mesh = finish(std::move(work), patch_list(surface, settings));
  validate(mesh);
  const auto geometry = compute_geometry(mesh);
  const auto quality = orthogonal_quality(mesh, geometry);
  auto quality_json = quality_report(mesh, geometry, quality, kind);
  const auto min_quality = *std::min_element(quality.begin(), quality.end());
  KUTTA_REQUIRE(min_quality > 1e-8, "mesh: orthogonal quality below hard floor (1e-8); nothing written");
  Index hexes = 0, pure = 0, prisms = 0;
  {
    const auto cf = build_cell_faces(mesh);
    for (Index c = 0; c < mesh.n_cells; ++c) {
      prisms += kind[c] == Prism;
      // An octree cell is a hex (possibly with split hanging faces) unless a
      // snap merge collapsed an edge, leaving fewer faces or a triangle.
      if (kind[c] != Hex || cf.offsets[c + 1] - cf.offsets[c] < 6) continue;
      bool quads = true, simple = cf.offsets[c + 1] - cf.offsets[c] == 6;
      for (auto i = cf.offsets[c]; i < cf.offsets[c + 1]; ++i) {
        const auto n = nodes(mesh, cf.faces[i]).size();
        quads = quads && n >= 4;
        simple = simple && n == 4;
      }
      hexes += quads;
      pure += simple;
    }
  }
  Json::Array zones, missing;
  for (const auto& p : mesh.patches) zones.push_back(Json::Object{{"name", p.name}, {"faces", double(p.size)}, {"type", p.source_type}});
  for (const auto& name : surface.zones)
    if (std::none_of(mesh.patches.begin(), mesh.patches.end(), [&](const Patch& p) { return p.name == name; })) missing.push_back(name);
  // Area-weighted first-cell y+ estimate on STL walls (centroid distance).
  const double nu = 1.7894e-5 / 1.225, length = surface.bounds.max.x - surface.bounds.min.x;
  const double cf = .026 / std::pow(speed * length / nu, 1. / 7), ustar = speed * std::sqrt(cf / 2);
  std::vector<std::pair<double, double>> yplus;  // (y+, area)
  double wall_area = 0;
  for (const auto& p : mesh.patches) {
    if (p.kind != PatchKind::Wall || p.name == "ground") continue;
    for (Index f = p.start; f < p.start + p.size; ++f) {
      const auto a = norm(geometry.face_area[f]);
      const auto d = std::abs(dot(geometry.face_centroid[f] - geometry.cell_centroid[mesh.owner[f]], geometry.face_area[f])) / a;
      yplus.push_back({d * ustar / nu, a});
      wall_area += a;
    }
  }
  std::sort(yplus.begin(), yplus.end());
  auto yplus_at = [&](double q) { double acc = 0; for (const auto& [y, a] : yplus) { acc += a; if (acc >= q * wall_area) return y; } return yplus.empty() ? 0. : yplus.back().first; };
  auto yplus_fraction = [&](double lo, double hi) { double acc = 0; for (const auto& [y, a] : yplus) if (y >= lo && y < hi) acc += a; return wall_area > 0 ? acc / wall_area : 0.; };
  lap("quality");
  const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  const auto usage = memory::usage();
  castellation["octree_cells_kept"] = double(octree_cells);
  Json::Object report{{"v", 1}, {"cells", double(mesh.n_cells)}, {"faces", double(mesh.n_faces())}, {"nodes", double(mesh.n_points())},
    {"cell_types", Json::Object{{"hex", double(hexes)}, {"hex_without_hanging_nodes", double(pure)}, {"poly", double(mesh.n_cells - hexes - prisms)}, {"prism", double(prisms)},
                                {"hex_fraction_of_core", double(hexes) / double(std::max<Index>(1, mesh.n_cells - prisms))}}},
    {"octree_leaves", double(tree.leaves)}, {"max_face_level_difference", level_jump},
    {"bounds", Json::Array{vec(tree.bounds.min), vec(tree.bounds.max)}}, {"body_bounds", Json::Array{vec(surface.bounds.min), vec(surface.bounds.max)}},
    {"zones", std::move(zones)}, {"missing_surface_zones", std::move(missing)}, {"repairs", surface.repairs},
    {"castellation", std::move(castellation)}, {"snap", std::move(snap_report)},
    {"quality", std::move(quality_json)}, {"layers", std::move(layer_report)},
    {"estimated_first_cell_y_plus", yplus.empty() ? Json(nullptr) : Json(Json::Object{{"speed_m_s", speed}, {"min", yplus.front().first}, {"p05", yplus_at(.05)}, {"p50", yplus_at(.5)}, {"p95", yplus_at(.95)}, {"max", yplus.back().first},
      {"area_fraction_below_5", yplus_fraction(0, 5)}, {"area_fraction_5_to_30", yplus_fraction(5, 30)}, {"area_fraction_30_to_300", yplus_fraction(30, 300)}, {"area_fraction_above_300", yplus_fraction(300, 1e300)},
      {"method", "flat-plate Cf=0.026/Re_L^(1/7) at the reference speed; owner-centroid distance on STL walls, area-weighted; estimate, not solved y+"}})},
    {"stage_seconds", std::move(timing)},
    {"wall_seconds", elapsed}, {"peak_memory_bytes", usage.peak ? Json(double(usage.peak)) : Json(nullptr)},
    {"geometry_method", "castellated octree (cells inside or crossed by the STL removed), wall snapped to the STL with implicit feature capture"},
    {"geometry_resolution_m", settings.min_size},
    {"warnings", Json::Array{"Snapping retreats where a cell would fall below orthogonal quality 0.15; see snap/point_distance_m for the remaining geometric error.",
                             "Slots narrower than about two cells close; see castellation/gaps.",
                             "Octree cubes with 2:1 hanging faces count as hex; Fluent stores every cell as a polyhedron."}}};
  return {std::move(mesh), std::move(report)};
}
}  // namespace kutta::mesher
