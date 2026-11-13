// Prism layers (snappyHexMesh addLayers style): the snapped mesh is pushed
// away from the selected walls along point normals, the push decaying over a
// few cell layers, and the gap is filled with one prism stack per wall face.
// Where layers end, the stack's side faces close the neighbouring cell
// (addPatchCellLayer style: points on the termination split into a wall copy
// and a lifted copy), so the unlayered wall stays on the surface. Where that
// split would fail quality, the termination is kept lifted instead: the
// unlayered face beside it tilts up to the stack and the stack ends in a
// step wall. (A wedge taper to zero thickness has orthogonal quality ~2t/L,
// far below the gate for thin first layers on coarse faces.) A face keeps its
// layers only while every cell stays above layer_min_quality; elsewhere
// layers collapse locally.
#include "internal.hpp"
#include <unordered_map>

namespace kutta::mesher::internal {
using namespace detail;
namespace {
constexpr double visibility = 0.5;  // each wall face sees the point normal within 60 degrees
constexpr double collapse = 0.25;   // a stack thinner than this fraction of its target collapses
constexpr double push_limit = 1.5;  // push at most 1.5 adjacent wall-cell heights (spread over the band)
constexpr int max_steps = 3;        // a stack may end on at most two sides (one for triangles)
double& at(Vec3d& v, int a) { return a == 0 ? v.x : a == 1 ? v.y : v.z; }

// The edges of the STL wall faces, numbered once.
struct WallEdges {
  std::vector<Index> start;            // per face: first entry in `ids`, -1 off the STL walls
  std::vector<Index> ids;              // edge j of face f is ids[start[f] + j], from nodes[j] to nodes[j+1]
  std::vector<std::array<Index, 2>> ends;
  std::vector<int> uses;               // wall faces sharing the edge
  std::unordered_map<std::uint64_t, Index> lookup;
  Index find(Index a, Index b) const {
    const auto it = lookup.find(key(a, b));
    return it == lookup.end() ? -1 : it->second;
  }
  static std::uint64_t key(Index a, Index b) { return std::uint64_t(std::uint32_t(std::min(a, b))) << 32 | std::uint32_t(std::max(a, b)); }
};
WallEdges wall_edges(const Work& w) {
  WallEdges e;
  e.start.assign(w.faces.size(), -1);
  std::vector<std::pair<std::uint64_t, Index>> keyed;  // (key, slot)
  for (std::size_t f = 0; f < w.faces.size(); ++f) {
    const auto& face = w.faces[f];
    if (face.neighbour >= 0 || face.patch < first_wall) continue;
    e.start[f] = index(e.ids.size());
    for (std::size_t j = 0; j < face.nodes.size(); ++j) {
      keyed.emplace_back(WallEdges::key(face.nodes[j], face.nodes[(j + 1) % face.nodes.size()]), index(e.ids.size()));
      e.ids.push_back(-1);
    }
  }
  std::sort(keyed.begin(), keyed.end());
  for (std::size_t i = 0; i < keyed.size(); ++i) {
    if (i == 0 || keyed[i].first != keyed[i - 1].first) {
      e.ends.push_back({Index(keyed[i].first >> 32), Index(keyed[i].first & 0xffffffffu)});
      e.uses.push_back(0);
    }
    e.ids[keyed[i].second] = index(e.ends.size()) - 1;
    ++e.uses.back();
  }
  e.lookup.reserve(e.ends.size());
  for (std::size_t i = 0; i < e.ends.size(); ++i) e.lookup.emplace(WallEdges::key(e.ends[i][0], e.ends[i][1]), index(i));
  return e;
}
}  // namespace

Json::Object add_layers(Work& w, const Surface& sf, const Settings& st, const Progress& progress) {
  const int layers = st.layer_count;
  const double g = st.growth, h1 = st.first_height;
  const double total = h1 * (std::pow(g, layers) - 1) / (g - 1);
  std::vector<double> fraction(layers + 1);
  for (int k = 0; k <= layers; ++k) fraction[k] = h1 * (std::pow(g, k) - 1) / (g - 1) / total;
  std::vector<char> zone_on(sf.zones.size());
  for (const auto& name : st.layer_zones) {
    if (name == "*") { std::fill(zone_on.begin(), zone_on.end(), 1); continue; }
    const auto it = std::find(sf.zones.begin(), sf.zones.end(), name);
    KUTTA_REQUIRE(it != sf.zones.end(), "mesh/layers/zones: unknown STL zone " + name);
    zone_on[std::size_t(it - sf.zones.begin())] = 1;
  }
  const Work base = w;
  const auto adjacency = cell_faces(base);
  const auto pf = point_faces(base);
  const auto base_geometry = geometry(base, adjacency);
  const auto edges = wall_edges(base);
  const auto np = index(base.points.size());
  const auto nf = index(base.faces.size());
  const auto ne = index(edges.ends.size());
  const auto selected = [&](Index f) {
    const auto& face = base.faces[f];
    return face.neighbour < 0 && face.patch >= first_wall && zone_on[std::size_t(face.patch - first_wall)];
  };
  const auto is_wall = [&](Index f) { return edges.start[f] >= 0; };
  std::vector<Vec3d> normal(np);
  std::vector<double> limit(np, std::numeric_limits<double>::infinity());
  std::vector<Index> sources;
  for (Index f = 0; f < nf; ++f) {
    const auto& face = base.faces[f];
    if (!is_wall(f)) continue;
    for (auto p : face.nodes) sources.push_back(p);
    if (!selected(f)) continue;
    const auto area = base_geometry.face[f].area;
    const auto height = base_geometry.volume[face.owner] / norm(area);
    for (auto p : face.nodes) { normal[p] = normal[p] - area; limit[p] = std::min(limit[p], height); }
  }
  std::sort(sources.begin(), sources.end());
  sources.erase(std::unique(sources.begin(), sources.end()), sources.end());
  for (Index p = 0; p < np; ++p) {
    for (int axis = 0; axis < 3; ++axis) if (base.planes[p] & (3 << (2 * axis))) at(normal[p], axis) = 0;
    const auto n = norm(normal[p]);
    normal[p] = n > 0 ? normal[p] / n : Vec3d{};
  }
  std::vector<char> extrude(nf);
  for (Index f = 0; f < nf; ++f) {
    if (!selected(f)) continue;
    const auto area = base_geometry.face[f].area;
    const auto unit = area / -norm(area);
    extrude[f] = std::all_of(base.faces[f].nodes.begin(), base.faces[f].nodes.end(), [&](Index p) { return dot(normal[p], unit) >= visibility; });
  }

  // Settling drops faces whose stacks cannot be built or would not be
  // stable, to a fixed point: an edge joining more than two extruded faces
  // cannot pair its side faces; a wall edge of more than two wall faces (a
  // thin feature after merges) leaves the termination ambiguous; and a stack
  // ending on all but one side (an island or a peninsula tip) made a
  // near-closed pocket the coupled solver diverged in.
  std::vector<int> uses(ne);
  const auto settle = [&] {
    for (bool dropped = true; dropped;) {
      std::fill(uses.begin(), uses.end(), 0);
      for (Index f = 0; f < nf; ++f)
        if (extrude[f]) for (std::size_t j = 0; j < base.faces[f].nodes.size(); ++j) ++uses[edges.ids[edges.start[f] + index(j)]];
      std::vector<Index> drop;
      for (Index f = 0; f < nf; ++f) {
        if (!extrude[f]) continue;
        const auto& ns = base.faces[f].nodes;
        int ends = 0;
        for (std::size_t j = 0; j < ns.size(); ++j) {
          const auto e = edges.ids[edges.start[f] + index(j)];
          if (uses[e] > 2 || edges.uses[e] > 2) ends = max_steps;
          else if (uses[e] == 1 && !(base.planes[ns[j]] & base.planes[ns[(j + 1) % ns.size()]])) ++ends;
        }
        if (ends >= std::min(max_steps, int(ns.size()) - 1)) drop.push_back(f);
      }
      for (auto f : drop) extrude[f] = 0;
      dropped = !drop.empty();
    }
  };

  std::vector<double> thickness(np);
  std::vector<Index> source_face;  // per new layer cell
  std::vector<char> lift(np);      // termination points kept lifted after a split failed quality
  Index split_points = 0, lifted_points = 0;
  int rounds = 0;
  bool built = true;  // w holds layer topology (or is not yet the base mesh)
  for (;; ++rounds) {
    KUTTA_REQUIRE(rounds < 50, "mesh: layer collapse did not terminate");
    settle();
    if (built) w = base;
    else { w.points = base.points; w.planes = base.planes; }
    built = false;
    std::vector<char> pushed(np);
    for (Index f = 0; f < nf; ++f) if (extrude[f]) for (auto p : base.faces[f].nodes) pushed[p] = 1;
    std::vector<Vec3d> target(np);
    for (auto p : sources) if (pushed[p]) target[p] = normal[p] * std::min(total, push_limit * limit[p]);
    std::vector<double> scale(np, 1.);
    displace(w, adjacency, pf, sources, target, scale, layer_min_quality, progress, "layers");
    for (Index p = 0; p < np; ++p) thickness[p] = norm(target[p]) * scale[p];
    bool changed = false;
    for (Index f = 0; f < nf; ++f)
      if (extrude[f] && std::any_of(base.faces[f].nodes.begin(), base.faces[f].nodes.end(), [&](Index p) { return thickness[p] < collapse * total; })) { extrude[f] = 0; changed = true; }
    if (changed) continue;

    // Topology: point p keeps its index as the stack top; levels 0..N-1 are new.
    built = true;
    std::vector<Index> level0(np, -1);
    for (Index f = 0; f < nf; ++f)
      if (extrude[f])
        for (auto p : base.faces[f].nodes)
          if (level0[p] < 0) {
            level0[p] = index(w.points.size());
            for (int k = 0; k < layers; ++k) {
              w.points.push_back(base.points[p] + (w.points[p] - base.points[p]) * fraction[k]);
              w.planes.push_back(base.planes[p]);
            }
          }
    const auto level = [&](Index p, int k) { return k == layers ? p : level0[p] + k; };
    // Per wall edge: on a stack, and the owner of an unlayered wall face on it.
    std::vector<char> stacked(ne);
    std::vector<Index> bare(ne, -1);
    for (Index f = 0; f < nf; ++f) {
      if (!is_wall(f)) continue;
      for (std::size_t j = 0; j < base.faces[f].nodes.size(); ++j) {
        const auto e = edges.ids[edges.start[f] + index(j)];
        if (extrude[f]) stacked[e] = 1; else bare[e] = base.faces[f].owner;
      }
    }
    const auto termination = [&](Index e) { return e >= 0 && stacked[e] && bare[e] >= 0; };
    // A lifted termination stays whole as stacks retreat: both ends of a
    // termination edge split, or both stay lifted.
    for (bool grew = true; grew;) {
      grew = false;
      for (Index e = 0; e < ne; ++e)
        if (termination(e) && lift[edges.ends[e][0]] != lift[edges.ends[e][1]]) { lift[edges.ends[e][0]] = lift[edges.ends[e][1]] = 1; grew = true; }
    }
    std::vector<char> split(np);
    for (Index f = 0; f < nf; ++f)
      if (is_wall(f) && !extrude[f])
        for (auto p : base.faces[f].nodes) if (level0[p] >= 0 && !lift[p]) split[p] = 1;
    split_points = index(std::size_t(std::count(split.begin(), split.end(), 1)));
    lifted_points = 0;
    for (Index p = 0; p < np; ++p) lifted_points += lift[p] && level0[p] >= 0;
    // A split point's unlayered wall faces take its wall copy. In any other
    // face, an edge running along the unlayered wall (on an unlayered wall
    // face, on no stack) ends at the wall copy, and the face climbs to the
    // lifted point through every level, as the stack side faces do.
    const auto along_wall = [&](Index a, Index b) { const auto e = edges.find(a, b); return e >= 0 && bare[e] >= 0 && !stacked[e]; };
    for (Index f = 0; f < nf; ++f) {
      auto& face = w.faces[f];
      if (extrude[f] || std::none_of(face.nodes.begin(), face.nodes.end(), [&](Index p) { return split[p]; })) continue;
      if (is_wall(f)) {
        for (auto& p : face.nodes) if (split[p]) p = level0[p];
        continue;
      }
      const auto& ns = base.faces[f].nodes;
      std::vector<Index> out;
      for (std::size_t j = 0; j < ns.size(); ++j) {
        const auto p = ns[j], before = ns[(j + ns.size() - 1) % ns.size()], after = ns[(j + 1) % ns.size()];
        if (!split[p]) { out.push_back(p); continue; }
        const bool back = along_wall(p, before), ahead = along_wall(p, after);
        if (back && ahead) out.push_back(level0[p]);
        else if (back) for (int k = 0; k <= layers; ++k) out.push_back(level(p, k));
        else if (ahead) for (int k = layers; k >= 0; --k) out.push_back(level(p, k));
        else out.push_back(p);
      }
      face.nodes = std::move(out);
    }
    source_face.clear();
    std::vector<Index> side(ne, -1);  // first of an edge's unpaired side faces (one per level)
    for (Index f = 0; f < nf; ++f) {
      if (!extrude[f]) continue;
      const auto first = w.cells;
      for (int k = 0; k < layers; ++k) { w.kind.push_back(Prism); w.size.push_back(base.size[base.faces[f].owner]); source_face.push_back(f); }
      w.cells += layers;
      const auto ns = base.faces[f].nodes;
      const auto patch = base.faces[f].patch;
      w.faces[f].neighbour = first + layers - 1;
      w.faces[f].patch = -1;
      for (int k = 0; k < layers; ++k) {
        Face face{{}, first + k, k ? first + k - 1 : -1, k ? -1 : patch};
        for (auto p : ns) face.nodes.push_back(level(p, k));
        w.faces.push_back(std::move(face));
      }
      for (std::size_t j = 0; j < ns.size(); ++j) {
        const auto a = ns[j], b = ns[(j + 1) % ns.size()];
        const auto e = edges.ids[edges.start[f] + index(j)];
        if (side[e] >= 0) {
          for (int k = 0; k < layers; ++k) w.faces[side[e] + k].neighbour = first + k;
          side[e] = -2;
          continue;
        }
        side[e] = index(w.faces.size());
        for (int k = 0; k < layers; ++k) w.faces.push_back({{level(a, k), level(a, k + 1), level(b, k + 1), level(b, k)}, first + k, -1, patch});
      }
    }
    // Unpaired side faces close the cell on the unlayered wall beside a split
    // termination, lie on a tunnel plane (the stack slides along it), or end
    // a lifted termination as a step wall of the stack's zone.
    for (Index e = 0; e < ne; ++e) {
      if (side[e] < 0) continue;
      const auto a = edges.ends[e][0], b = edges.ends[e][1];
      const auto shared = base.planes[a] & base.planes[b];
      for (int k = 0; k < layers; ++k) {
        auto& face = w.faces[side[e] + k];
        if (bare[e] >= 0 && split[a] && split[b]) { face.neighbour = bare[e]; face.patch = -1; continue; }
        for (int bit = 0; bit < 6; ++bit) if (shared & (1 << bit)) { face.patch = bit; break; }
      }
    }
    const auto layered = cell_faces(w);
    const auto check = geometry(w, layered);
    for (Index c = 0; c < w.cells; ++c) {
      if (check.quality[c] >= layer_min_quality) continue;
      changed = true;
      if (c >= base.cells) { extrude[source_face[c - base.cells]] = 0; continue; }
      // A bad octree cell drops its own stacks.
      bool own = false;
      for (auto i = adjacency.offsets[c]; i < adjacency.offsets[c + 1]; ++i)
        if (extrude[adjacency.faces[i]]) { extrude[adjacency.faces[i]] = 0; own = true; }
      if (own) continue;
      // Beside a split termination it keeps that whole termination (connected
      // along its edges) lifted: the form the displacement already kept above
      // quality. Dropping the stacks instead would only move the failing end
      // of the termination along by one cell.
      std::vector<Index> front;
      for (auto i = adjacency.offsets[c]; i < adjacency.offsets[c + 1]; ++i)
        for (auto p : base.faces[adjacency.faces[i]].nodes) if (split[p] && !lift[p]) { lift[p] = 1; front.push_back(p); }
      if (!front.empty()) {
        while (!front.empty()) {
          const auto p = front.back();
          front.pop_back();
          point_neighbours(base, pf, p, [&](Index q) {
            if (split[q] && !lift[q] && termination(edges.find(p, q))) { lift[q] = 1; front.push_back(q); }
          });
        }
        continue;
      }
      for (auto i = adjacency.offsets[c]; i < adjacency.offsets[c + 1]; ++i)
        for (auto p : base.faces[adjacency.faces[i]].nodes)
          for (auto j = pf.offsets[p]; j < pf.offsets[p + 1]; ++j) extrude[pf.faces[j]] = 0;
    }
    if (!changed) break;
    if (progress) progress("layers", .8, "Collapsing layers where cell quality would fail");
  }

  double area_selected = 0, area_full = 0, area_thin = 0, first_min = std::numeric_limits<double>::infinity(), first_sum = 0;
  Index stacks = 0;
  for (Index f = 0; f < nf; ++f) {
    if (!selected(f)) continue;
    const auto a = norm(base_geometry.face[f].area);
    area_selected += a;
    if (!extrude[f]) continue;
    ++stacks;
    const bool full = std::all_of(base.faces[f].nodes.begin(), base.faces[f].nodes.end(), [&](Index p) { return thickness[p] >= .99 * total; });
    (full ? area_full : area_thin) += a;
    double mean = 0;
    for (auto p : base.faces[f].nodes) { first_min = std::min(first_min, thickness[p] * fraction[1]); mean += thickness[p] * fraction[1]; }
    first_sum += mean / double(base.faces[f].nodes.size());
  }
  return Json::Object{
    {"requested", layers}, {"generated", stacks ? layers : 0}, {"first_height", h1}, {"growth", g}, {"total_thickness", total},
    {"stacks", double(stacks)}, {"cells", double(stacks) * layers},
    {"coverage", Json::Object{{"full", area_selected > 0 ? area_full / area_selected : 0.}, {"thinned", area_selected > 0 ? area_thin / area_selected : 0.},
                              {"none", area_selected > 0 ? 1 - (area_full + area_thin) / area_selected : 0.}, {"wall_area_m2", area_selected}}},
    {"first_height_achieved", Json::Object{{"min", stacks ? first_min : 0.}, {"mean", stacks ? first_sum / double(stacks) : 0.}}},
    {"collapse_rounds", rounds}, {"min_quality", layer_min_quality},
    {"terminations", Json::Object{{"split_points", double(split_points)}, {"lifted_points", double(lifted_points)}}},
    {"status", "prism layers by pushing the snapped mesh; terminations split so unlayered walls stay on the surface, kept lifted (a step) where a split fails quality; collapsed locally where quality fails"}};
}
}  // namespace kutta::mesher::internal
