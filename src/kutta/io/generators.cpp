#include <map>
#include <numeric>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

namespace kutta {
using namespace detail;
namespace {
PolyMesh box(int nx, int ny, int nz, Vec3d lo, Vec3d hi, bool tet) {
  KUTTA_REQUIRE(nx > 0 && ny > 0 && nz > 0, "box dimensions must be positive");
  KUTTA_REQUIRE(
      finite(lo) && finite(hi) && hi.x > lo.x && hi.y > lo.y && hi.z > lo.z,
      "box requires finite lo < hi in every coordinate");
  const auto capacity =
      static_cast<std::uint64_t>(std::numeric_limits<Index>::max()) /
      (tet ? 36 : 6);
  auto count = static_cast<std::uint64_t>(nx);
  KUTTA_REQUIRE(
      count <= capacity && static_cast<std::uint64_t>(ny) <= capacity / count,
      "box exceeds 32-bit mesh capacity");
  count *= ny;
  KUTTA_REQUIRE(static_cast<std::uint64_t>(nz) <= capacity / count,
                "box exceeds 32-bit mesh capacity");
  count *= nz;
  auto np = (static_cast<std::uint64_t>(nx) + 1) *
            (static_cast<std::uint64_t>(ny) + 1) *
            (static_cast<std::uint64_t>(nz) + 1);
  KUTTA_REQUIRE(
      np <= static_cast<std::uint64_t>(std::numeric_limits<Index>::max()),
      "box exceeds 32-bit point capacity");
  auto node = [&](int i, int j, int k) {
    return static_cast<Index>(
        (static_cast<std::uint64_t>(k) * (ny + 1) + j) * (nx + 1) + i);
  };
  std::vector<Vec3d> points;
  points.reserve(static_cast<std::size_t>(np));
  for (int k = 0; k <= nz; ++k)
    for (int j = 0; j <= ny; ++j)
      for (int i = 0; i <= nx; ++i)
        points.push_back({lo.x + (hi.x - lo.x) * i / nx,
                          lo.y + (hi.y - lo.y) * j / ny,
                          lo.z + (hi.z - lo.z) * k / nz});
  std::vector<detail::CellFaces> cells;
  cells.reserve(static_cast<std::size_t>(count) * (tet ? 6 : 1));
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        std::array<Index, 8> v{node(i, j, k),
                               node(i + 1, j, k),
                               node(i + 1, j + 1, k),
                               node(i, j + 1, k),
                               node(i, j, k + 1),
                               node(i + 1, j, k + 1),
                               node(i + 1, j + 1, k + 1),
                               node(i, j + 1, k + 1)};
        if (!tet)
          cells.push_back({{v[0], v[3], v[2], v[1]},
                           {v[4], v[5], v[6], v[7]},
                           {v[0], v[1], v[5], v[4]},
                           {v[1], v[2], v[6], v[5]},
                           {v[2], v[3], v[7], v[6]},
                           {v[3], v[0], v[4], v[7]}});
        else {
          constexpr int ring[6] = {1, 2, 3, 7, 4, 5};
          for (int t = 0; t < 6; ++t) {
            auto a = v[0], b = v[ring[t]], c = v[ring[(t + 1) % 6]], d = v[6];
            cells.push_back({{a, b, c}, {a, d, b}, {b, d, c}, {c, d, a}});
          }
        }
      }
  std::vector<Patch> patches;
  for (auto name : {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"})
    patches.push_back({name, PatchKind::Wall, ""});
  auto classify_ids = [&](const std::vector<Index>& ns) -> Index {
    for (Index side = 0; side < 6; ++side) {
      bool all = true;
      for (auto n : ns) {
        auto ii = n % (nx + 1), jj = (n / (nx + 1)) % (ny + 1),
             kk = n / ((nx + 1) * (ny + 1));
        all = all && (side == 0   ? ii == 0
                      : side == 1 ? ii == nx
                      : side == 2 ? jj == 0
                      : side == 3 ? jj == ny
                      : side == 4 ? kk == 0
                                  : kk == nz);
      }
      if (all) return side;
    }
    throw Error("box boundary classification failed");
  };
  auto m =
      from_cells(std::move(points), cells, std::move(patches), classify_ids);
  validate(m);
  return m;
}
}  // namespace
PolyMesh make_box_hex(int nx, int ny, int nz, Vec3d lo, Vec3d hi) {
  return box(nx, ny, nz, lo, hi, false);
}
PolyMesh make_box_tet(int nx, int ny, int nz, Vec3d lo, Vec3d hi) {
  return box(nx, ny, nz, lo, hi, true);
}

PolyMesh extrude_2d(const PlanarMesh& p, double thickness) {
  KUTTA_REQUIRE(std::isfinite(thickness) && thickness > 0,
                "2D extrusion thickness must be positive and finite");
  KUTTA_REQUIRE(p.n_cells > 0, "2D mesh has no cells");
  KUTTA_REQUIRE(p.edges.size() == p.owner.size() &&
                    p.edges.size() == p.neighbour.size() &&
                    p.edges.size() == p.edge_patch.size(),
                "2D edge adjacency lengths disagree");
  auto n = index(p.points.size());
  KUTTA_REQUIRE(n <= std::numeric_limits<Index>::max() / 2,
                "extrusion exceeds point capacity");
  std::vector<Vec3d> points = p.points;
  for (auto a : p.points) {
    KUTTA_REQUIRE(finite(a) && a.z == 0,
                  "2D points must be finite and lie in xy plane");
    a.z = thickness;
    points.push_back(a);
  }
  std::vector<Face> faces;
  std::vector<std::map<Index, Index>> loops(p.n_cells);
  for (std::size_t i = 0; i < p.edges.size(); ++i) {
    auto [a, b] = p.edges[i];
    auto o = p.owner[i], nb = p.neighbour[i];
    KUTTA_REQUIRE(a >= 0 && a < n && b >= 0 && b < n && a != b,
                  "2D edge has invalid nodes");
    KUTTA_REQUIRE(
        o >= 0 && o < p.n_cells && nb >= -1 && nb < p.n_cells && nb != o,
        "2D edge has invalid adjacent cells");
    KUTTA_REQUIRE(loops[o].emplace(a, b).second,
                  "2D owner has branched/repeated oriented edges");
    if (nb >= 0)
      KUTTA_REQUIRE(loops[nb].emplace(b, a).second,
                    "2D neighbour has branched/repeated oriented edges");
    faces.push_back({{a, b, b + n, a + n}, o, nb, p.edge_patch[i]});
  }
  auto patches = p.patches;
  auto front = index(patches.size()), back = front + 1;
  patches.push_back({"front", PatchKind::Symmetry, "symmetry"});
  patches.push_back({"back", PatchKind::Symmetry, "symmetry"});
  for (Index c = 0; c < p.n_cells; ++c) {
    auto& loop = loops[c];
    KUTTA_REQUIRE(loop.size() >= 3, "2D cell has fewer than three edges");
    std::vector<Index> polygon;
    auto start = loop.begin()->first, cur = start;
    do {
      KUTTA_REQUIRE(polygon.size() < loop.size(),
                    "2D cell edges do not form one closed loop");
      polygon.push_back(cur);
      auto it = loop.find(cur);
      KUTTA_REQUIRE(it != loop.end(), "2D cell has an open edge loop");
      cur = it->second;
    } while (cur != start);
    KUTTA_REQUIRE(polygon.size() == loop.size(),
                  "2D cell contains disconnected edge loops");
    auto top = polygon;
    for (auto& v : top) v += n;
    std::reverse(polygon.begin(), polygon.end());
    faces.push_back({std::move(polygon), c, -1, front});
    faces.push_back({std::move(top), c, -1, back});
  }
  auto m = pack(std::move(points), p.n_cells, std::move(faces),
                std::move(patches), p.cell_zones);
  validate(m);
  return m;
}
Index split_patch(PolyMesh& m, Index patch, const std::string& name,
                  const std::function<bool(Vec3d)>& predicate) {
  validate(m);
  KUTTA_REQUIRE(patch >= 0 && patch < index(m.patches.size()),
                "split_patch: patch id out of range");
  auto g = compute_geometry(m);
  auto patches = m.patches;
  auto added = patches[patch];
  added.name = name;
  patches.insert(patches.begin() + patch + 1, added);
  std::vector<Face> faces;
  for (Index f = 0; f < m.n_internal_faces(); ++f) {
    auto ns = nodes(m, f);
    faces.push_back({{ns.begin(), ns.end()}, m.owner[f], m.neighbour[f], -1});
  }
  for (Index i = 0; i < index(m.patches.size()); ++i)
    for (Index f = m.patches[i].start;
         f < m.patches[i].start + m.patches[i].size; ++f) {
      auto pi = i > patch ? i + 1 : i;
      if (i == patch && predicate(g.face_centroid[f])) ++pi;
      auto ns = nodes(m, f);
      faces.push_back({{ns.begin(), ns.end()}, m.owner[f], -1, pi});
    }
  auto result = pack(m.points, m.n_cells, std::move(faces), std::move(patches),
                     m.cell_zones);
  validate(result);
  m = std::move(result);
  return patch + 1;
}
}  // namespace kutta
