#include "kutta/mesh/mesh.hpp"

#include <map>
#include <set>

#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta {
using namespace detail;
namespace {
void structure(const PolyMesh& m) {
  index(m.points.size());
  index(m.owner.size());
  index(m.face_nodes.size());
  KUTTA_REQUIRE(m.n_cells > 0, "mesh must contain at least one cell");
  KUTTA_REQUIRE(m.neighbour.size() <= m.owner.size(),
                "internal face count exceeds total face count");
  KUTTA_REQUIRE(m.face_offsets.size() == m.owner.size() + 1,
                "face_offsets must have n_faces+1 entries");
  KUTTA_REQUIRE(m.face_offsets.front() == 0, "face_offsets must begin at zero");
  KUTTA_REQUIRE(m.face_offsets.back() == index(m.face_nodes.size()),
                "last face offset must equal face_nodes size");
  for (Index p = 0; p < m.n_points(); ++p)
    KUTTA_REQUIRE(finite(m.points[p]), "non-finite point " + std::to_string(p));
  for (Index f = 0; f < m.n_faces(); ++f) {
    const auto label = "face " + std::to_string(f) + ": ";
    KUTTA_REQUIRE(m.face_offsets[f] >= 0 &&
                      m.face_offsets[f + 1] >= m.face_offsets[f] &&
                      m.face_offsets[f + 1] <= index(m.face_nodes.size()),
                  label + "invalid CSR offsets");
    KUTTA_REQUIRE(m.face_offsets[f + 1] - m.face_offsets[f] >= 3,
                  label + "requires at least three nodes");
    std::set<Index> unique;
    for (auto p : nodes(m, f)) {
      KUTTA_REQUIRE(p >= 0 && p < m.n_points(),
                    label + "node index outside points");
      KUTTA_REQUIRE(unique.insert(p).second,
                    label + "repeated node " + std::to_string(p));
    }
    KUTTA_REQUIRE(m.owner[f] >= 0 && m.owner[f] < m.n_cells,
                  label + "owner outside cell range");
    if (f < m.n_internal_faces())
      KUTTA_REQUIRE(
          m.neighbour[f] > m.owner[f] && m.neighbour[f] < m.n_cells,
          label + "internal face requires 0 <= owner < neighbour < n_cells");
  }
  Index next = m.n_internal_faces();
  std::set<std::string> patch_names;
  for (const auto& p : m.patches) {
    KUTTA_REQUIRE(patch_names.insert(p.name).second,
                  "duplicate patch name: " + p.name);
    KUTTA_REQUIRE(
        p.start == next && p.size >= 0 && p.size <= m.n_faces() - next,
        "patch '" + p.name +
            "' must contiguously cover boundary faces in patch order");
    next += p.size;
  }
  KUTTA_REQUIRE(next == m.n_faces(), "patches do not cover all boundary faces");
  for (const auto& z : m.cell_zones) {
    std::set<Index> seen;
    for (auto c : z.cells) {
      KUTTA_REQUIRE(c >= 0 && c < m.n_cells,
                    "cell zone '" + z.name + "' has out-of-range cell");
      KUTTA_REQUIRE(seen.insert(c).second,
                    "cell zone '" + z.name + "' repeats a cell");
    }
  }
}
MeshGeometry geometry(const PolyMesh& m) {
  MeshGeometry g;
  g.face_centroid.resize(m.n_faces());
  g.face_area.resize(m.n_faces());
  g.cell_centroid.resize(m.n_cells);
  g.cell_volume.resize(m.n_cells);
  std::vector<Vec3d> estimate(m.n_cells), anchor(m.n_cells);
  std::vector<Index> count(m.n_cells);
  // Work with differences from a local anchor, including for translated meshes.
  for (Index f = 0; f < m.n_faces(); ++f) {
    const Vec3d center = face_center(m, f);
    Vec3d moment{};
    double area = 0;
    const auto ns = nodes(m, f);
    for (std::size_t j = 0; j < ns.size(); ++j)
      g.face_area[f] = g.face_area[f] +
          cross(m.points[ns[j]] - center,
                m.points[ns[(j + 1) % ns.size()]] - center) * 0.5;
    const auto magnitude = norm(g.face_area[f]);
    KUTTA_REQUIRE(magnitude > 0 && std::isfinite(magnitude),
                  "face " + std::to_string(f) + " has zero/invalid area");
    const auto normal = g.face_area[f] / magnitude;
    for (std::size_t j = 0; j < ns.size(); ++j) {
      auto a = m.points[ns[j]] - center,
           b = m.points[ns[(j + 1) % ns.size()]] - center;
      const auto w = dot(cross(a, b) * 0.5, normal);
      area += w;
      moment = moment + (a + b) * (w / 3);
    }
    KUTTA_REQUIRE(area > 0 && std::isfinite(area),
                  "face " + std::to_string(f) + " has invalid projected area");
    g.face_centroid[f] = center + moment / area;
    auto add = [&](Index c) {
      if (!count[c]) anchor[c] = center;
      estimate[c] = estimate[c] + (center - anchor[c]);
      ++count[c];
    };
    add(m.owner[f]);
    if (f < m.n_internal_faces()) add(m.neighbour[f]);
  }
  for (Index c = 0; c < m.n_cells; ++c) {
    KUTTA_REQUIRE(count[c] >= 4,
                  "cell " + std::to_string(c) + " has fewer than four faces");
    estimate[c] = anchor[c] + estimate[c] / count[c];
  }
  // Signed tetrahedral pyramids use the same center fan as the face geometry.
  // This integrates first moments of warped faces without flattening them.
  for (Index f = 0; f < m.n_faces(); ++f) {
    auto ns = nodes(m, f);
    auto center = face_center(m, f);
    auto add = [&](Index c, double sign) {
      auto r = center - estimate[c];
      for (std::size_t j = 0; j < ns.size(); ++j) {
        auto a = m.points[ns[j]] - estimate[c],
             b = m.points[ns[(j + 1) % ns.size()]] - estimate[c];
        double v = sign * dot(r, cross(a, b)) / 6;
        g.cell_volume[c] += v;
        g.cell_centroid[c] = g.cell_centroid[c] + (r + a + b) * (v / 4);
      }
    };
    add(m.owner[f], 1);
    if (f < m.n_internal_faces()) add(m.neighbour[f], -1);
  }
  for (Index c = 0; c < m.n_cells; ++c) {
    KUTTA_REQUIRE(g.cell_volume[c] > 0 && std::isfinite(g.cell_volume[c]),
                  "cell " + std::to_string(c) +
                      " has non-positive/invalid volume; check face winding "
                      "and connectivity");
    g.cell_centroid[c] = estimate[c] + g.cell_centroid[c] / g.cell_volume[c];
    KUTTA_REQUIRE(finite(g.cell_centroid[c]),
                  "invalid centroid for cell " + std::to_string(c));
  }
  return g;
}
}  // namespace
MeshGeometry compute_geometry(const PolyMesh& m) {
  structure(m);
  return geometry(m);
}
void validate(const PolyMesh& m) {
  structure(m);
  auto g = geometry(m);
  std::vector<Vec3d> closure(m.n_cells);
  std::vector<double> areas(m.n_cells);
  for (Index f = 0; f < m.n_faces(); ++f) {
    auto add = [&](Index c, int sign) {
      closure[c] = closure[c] + g.face_area[f] * sign;
      areas[c] += norm(g.face_area[f]);
    };
    add(m.owner[f], 1);
    if (f < m.n_internal_faces()) add(m.neighbour[f], -1);
    auto d = g.face_centroid[f] - g.cell_centroid[m.owner[f]];
    KUTTA_REQUIRE(
        dot(g.face_area[f], d) > 0,
        "face " + std::to_string(f) +
            " normal does not point out of owner; reverse node winding");
    if (f < m.n_internal_faces()) {
      auto dn = g.cell_centroid[m.neighbour[f]] - g.face_centroid[f];
      KUTTA_REQUIRE(dot(g.face_area[f], dn) > 0,
                    "face " + std::to_string(f) +
                        " normal does not point into neighbour");
    }
  }
  auto adjacency = build_cell_faces(m);
  for (Index c = 0; c < m.n_cells; ++c) {
    KUTTA_REQUIRE(norm(closure[c]) <= 1e-10 * areas[c],
                  "cell " + std::to_string(c) +
                      " is not closed: outward face areas do not sum to zero");
    // Keep the edge map local to one cell rather than retaining millions of
    // tree nodes for the entire mesh. Area closure alone can hide holes.
    std::map<std::pair<Index, Index>, std::pair<int, int>> edges;
    for (auto i = adjacency.offsets[c]; i < adjacency.offsets[c + 1]; ++i) {
      auto ns = nodes(m, adjacency.faces[i]);
      for (std::size_t j = 0; j < ns.size(); ++j) {
        auto a = ns[j], b = ns[(j + 1) % ns.size()];
        if (adjacency.signs[i] < 0) std::swap(a, b);
        auto& edge = edges[std::minmax(a, b)];
        ++edge.first;
        edge.second += (a < b ? 1 : -1);
      }
    }
    for (const auto& [edge, n] : edges) {
      (void)edge;
      KUTTA_REQUIRE(
          n.first == 2 && n.second == 0,
          "cell " + std::to_string(c) +
              " is not a closed oriented manifold: unmatched/repeated edge");
    }
  }
}
}  // namespace kutta
