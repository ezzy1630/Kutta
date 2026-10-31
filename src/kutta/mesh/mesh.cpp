#include "kutta/mesh/mesh.hpp"

#include <algorithm>
#include <set>
#include <tuple>
#include <vector>

#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta {
using namespace detail;
FaceGeometry face_geometry(const PolyMesh& m, Index f) {
  const Vec3d center = face_center(m, f);
  Vec3d moment{}, area_vector{};
  double area = 0;
  const auto ns = nodes(m, f);
  for (std::size_t j = 0; j < ns.size(); ++j)
    area_vector = area_vector +
        cross(m.points[ns[j]] - center,
              m.points[ns[(j + 1) % ns.size()]] - center) * 0.5;
  const auto magnitude = norm(area_vector);
  KUTTA_REQUIRE(magnitude > 0 && std::isfinite(magnitude),
                "face " + std::to_string(f) + " has zero/invalid area");
  const auto normal = area_vector / magnitude;
  for (std::size_t j = 0; j < ns.size(); ++j) {
    auto a = m.points[ns[j]] - center,
         b = m.points[ns[(j + 1) % ns.size()]] - center;
    const auto w = dot(cross(a, b) * 0.5, normal);
    area += w;
    moment = moment + (a + b) * (w / 3);
  }
  KUTTA_REQUIRE(area > 0 && std::isfinite(area),
                "face " + std::to_string(f) + " has invalid projected area");
  return {center + moment / area, area_vector};
}
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
    // Built only on failure: millions of faces make eager labels a startup cost.
    const auto label = [f] { return "face " + std::to_string(f) + ": "; };
    KUTTA_REQUIRE(m.face_offsets[f] >= 0 &&
                      m.face_offsets[f + 1] >= m.face_offsets[f] &&
                      m.face_offsets[f + 1] <= index(m.face_nodes.size()),
                  label() + "invalid CSR offsets");
    KUTTA_REQUIRE(m.face_offsets[f + 1] - m.face_offsets[f] >= 3,
                  label() + "requires at least three nodes");
    const auto ns = nodes(m, f);
    for (std::size_t j = 0; j < ns.size(); ++j) {
      const auto p = ns[j];
      KUTTA_REQUIRE(p >= 0 && p < m.n_points(),
                    label() + "node index outside points");
      KUTTA_REQUIRE(std::find(ns.begin(), ns.begin() + j, p) == ns.begin() + j,
                    label() + "repeated node " + std::to_string(p));
    }
    KUTTA_REQUIRE(m.owner[f] >= 0 && m.owner[f] < m.n_cells,
                  label() + "owner outside cell range");
    if (f < m.n_internal_faces())
      KUTTA_REQUIRE(
          m.neighbour[f] > m.owner[f] && m.neighbour[f] < m.n_cells,
          label() + "internal face requires 0 <= owner < neighbour < n_cells");
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
    std::vector<bool> seen(m.n_cells);
    for (auto c : z.cells) {
      KUTTA_REQUIRE(c >= 0 && c < m.n_cells,
                    "cell zone '" + z.name + "' has out-of-range cell");
      KUTTA_REQUIRE(!seen[c], "cell zone '" + z.name + "' repeats a cell");
      seen[c] = true;
    }
  }
}
// Cell centroids and volumes from `face(f)`, which may recompute each face or
// read stored arrays: face_geometry gives the same doubles either way, so
// validate() needs no per-face storage.
template <class Faces>
void cell_geometry(const PolyMesh& m, const Faces& face,
                   std::vector<Vec3d>& centroid, std::vector<double>& volume) {
  centroid.assign(m.n_cells, {});
  volume.assign(m.n_cells, 0);
  std::vector<Vec3d> estimate(m.n_cells);
  {
    std::vector<Vec3d> anchor(m.n_cells);
    std::vector<Index> count(m.n_cells);
    // Work with differences from a local anchor, including for translated meshes.
    for (Index f = 0; f < m.n_faces(); ++f) {
      const auto fc = face(f).centroid;
      auto add = [&](Index c) {
        if (!count[c]) anchor[c] = fc;
        estimate[c] = estimate[c] + (fc - anchor[c]);
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
  }
  // Fluent/OpenFOAM convention: one pyramid per face from the mean of the
  // cell's face centroids to the face's (flattened) area vector and centroid,
  // with the pyramid centroid 3/4 of the way to the face. Saved Fluent wall
  // distances, y+ and wall omega reproduce only with this convention (exact
  // warped-face tetrahedra shift wall-normal distances by up to 0.6%; see
  // docs/research/fluent-sst-forensics.md and fv/turbulence/README.md).
  for (Index f = 0; f < m.n_faces(); ++f) {
    const auto g = face(f);
    auto add = [&](Index c, double sign) {
      const auto r = g.centroid - estimate[c];
      const double v = sign * dot(g.area, r) / 3;
      volume[c] += v;
      centroid[c] = centroid[c] + r * (0.75 * v);
    };
    add(m.owner[f], 1);
    if (f < m.n_internal_faces()) add(m.neighbour[f], -1);
  }
  for (Index c = 0; c < m.n_cells; ++c) {
    KUTTA_REQUIRE(volume[c] > 0 && std::isfinite(volume[c]),
                  "cell " + std::to_string(c) +
                      " has non-positive/invalid volume; check face winding "
                      "and connectivity");
    centroid[c] = estimate[c] + centroid[c] / volume[c];
    KUTTA_REQUIRE(finite(centroid[c]),
                  "invalid centroid for cell " + std::to_string(c));
  }
}
MeshGeometry geometry(const PolyMesh& m) {
  MeshGeometry g;
  g.face_centroid.resize(m.n_faces());
  g.face_area.resize(m.n_faces());
  for (Index f = 0; f < m.n_faces(); ++f) {
    const auto face = face_geometry(m, f);
    g.face_centroid[f] = face.centroid;
    g.face_area[f] = face.area;
  }
  cell_geometry(m, [&](Index f) { return FaceGeometry{g.face_centroid[f], g.face_area[f]}; },
                g.cell_centroid, g.cell_volume);
  return g;
}
}  // namespace
MeshGeometry compute_geometry(const PolyMesh& m) {
  structure(m);
  return geometry(m);
}
std::vector<double> orthogonal_quality(const PolyMesh& m, const MeshGeometry& g) {
  std::vector<double> quality(m.n_cells, 1);
  const auto cosine = [](Vec3d a, Vec3d b) { return dot(a, b) / (norm(a) * norm(b)); };
  for (Index f = 0; f < m.n_faces(); ++f) {
    const auto o = m.owner[f];
    quality[o] = std::min(quality[o], cosine(g.face_area[f], g.face_centroid[f] - g.cell_centroid[o]));
    if (f >= m.n_internal_faces()) continue;
    const auto n = m.neighbour[f];
    const auto link = g.cell_centroid[n] - g.cell_centroid[o];
    quality[o] = std::min({quality[o], cosine(g.face_area[f], link)});
    quality[n] = std::min({quality[n], cosine(g.face_area[f], g.cell_centroid[n] - g.face_centroid[f]),
                           cosine(g.face_area[f], link)});
  }
  for (auto& q : quality) q = std::max(q, 0.);
  return quality;
}
void validate(const PolyMesh& m) {
  structure(m);
  // Same checks as on compute_geometry's arrays, with face geometry
  // recomputed where needed instead of stored for every face.
  std::vector<Vec3d> cell_centroid;
  {
    std::vector<double> volume;
    cell_geometry(m, [&](Index f) { return face_geometry(m, f); }, cell_centroid, volume);
  }
  std::vector<Vec3d> closure(m.n_cells);
  std::vector<double> areas(m.n_cells);
  for (Index f = 0; f < m.n_faces(); ++f) {
    const auto g = face_geometry(m, f);
    auto add = [&](Index c, int sign) {
      closure[c] = closure[c] + g.area * sign;
      areas[c] += norm(g.area);
    };
    add(m.owner[f], 1);
    if (f < m.n_internal_faces()) add(m.neighbour[f], -1);
    auto d = g.centroid - cell_centroid[m.owner[f]];
    KUTTA_REQUIRE(
        dot(g.area, d) > 0,
        "face " + std::to_string(f) +
            " normal does not point out of owner; reverse node winding");
    if (f < m.n_internal_faces()) {
      auto dn = cell_centroid[m.neighbour[f]] - g.centroid;
      KUTTA_REQUIRE(dot(g.area, dn) > 0,
                    "face " + std::to_string(f) +
                        " normal does not point into neighbour");
    }
  }
  std::vector<Vec3d>().swap(cell_centroid);
  auto adjacency = build_cell_faces(m);
  // One cell's directed edges, reused across cells. Area closure alone can
  // hide holes: every undirected edge must appear exactly twice, once in
  // each direction.
  std::vector<std::tuple<Index, Index, int>> edges;
  for (Index c = 0; c < m.n_cells; ++c) {
    KUTTA_REQUIRE(norm(closure[c]) <= 1e-10 * areas[c],
                  "cell " + std::to_string(c) +
                      " is not closed: outward face areas do not sum to zero");
    edges.clear();
    for (auto i = adjacency.offsets[c]; i < adjacency.offsets[c + 1]; ++i) {
      auto ns = nodes(m, adjacency.faces[i]);
      for (std::size_t j = 0; j < ns.size(); ++j) {
        auto a = ns[j], b = ns[(j + 1) % ns.size()];
        if (adjacency.signs[i] < 0) std::swap(a, b);
        const auto [lo, hi] = std::minmax(a, b);
        edges.emplace_back(lo, hi, a < b ? 1 : -1);
      }
    }
    std::sort(edges.begin(), edges.end());
    for (std::size_t j = 0; j < edges.size();) {
      std::size_t k = j;
      int count = 0, orientation = 0;
      for (; k < edges.size() && std::get<0>(edges[k]) == std::get<0>(edges[j]) &&
             std::get<1>(edges[k]) == std::get<1>(edges[j]); ++k) {
        ++count;
        orientation += std::get<2>(edges[k]);
      }
      KUTTA_REQUIRE(
          count == 2 && orientation == 0,
          "cell " + std::to_string(c) +
              " is not a closed oriented manifold: unmatched/repeated edge");
      j = k;
    }
  }
}
}  // namespace kutta
