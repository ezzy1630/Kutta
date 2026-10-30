#include "kutta/mesh/topology.hpp"

#include <algorithm>
#include <map>
#include <numeric>
#include <queue>
#include <tuple>

#include "kutta/mesh/detail.hpp"

namespace kutta::detail {
PolyMesh pack(std::vector<Vec3d> points, Index cells, std::vector<Face> faces,
              std::vector<Patch> patches, std::vector<CellZone> zones) {
  index(faces.size());
  for (auto& f : faces) {
    KUTTA_REQUIRE(f.owner >= 0 && f.owner < cells,
                  "construction: invalid face owner");
    KUTTA_REQUIRE(f.neighbour < cells && f.neighbour >= -1,
                  "construction: invalid face neighbour");
    if (f.neighbour >= 0) {
      KUTTA_REQUIRE(f.owner != f.neighbour,
                    "construction: face has identical adjacent cells");
      if (f.owner > f.neighbour) {
        std::swap(f.owner, f.neighbour);
        std::reverse(f.nodes.begin(), f.nodes.end());
      }
    } else
      KUTTA_REQUIRE(f.patch >= 0 && f.patch < index(patches.size()),
                    "construction: boundary face has invalid patch");
  }
  std::stable_sort(
      faces.begin(), faces.end(), [](const Face& a, const Face& b) {
        return std::tuple(a.neighbour < 0, a.neighbour < 0 ? a.patch : a.owner,
                          a.neighbour) <
               std::tuple(b.neighbour < 0, b.neighbour < 0 ? b.patch : b.owner,
                          b.neighbour);
      });
  PolyMesh m;
  m.points = std::move(points);
  m.n_cells = cells;
  m.cell_zones = std::move(zones);
  for (auto& p : patches) p.size = 0;
  for (const auto& f : faces) {
    m.owner.push_back(f.owner);
    if (f.neighbour >= 0)
      m.neighbour.push_back(f.neighbour);
    else
      ++patches[f.patch].size;
    m.face_nodes.insert(m.face_nodes.end(), f.nodes.begin(), f.nodes.end());
    m.face_offsets.push_back(index(m.face_nodes.size()));
  }
  auto start = m.n_internal_faces();
  for (auto& p : patches) {
    p.start = start;
    start += p.size;
  }
  m.patches = std::move(patches);
  return m;
}
PolyMesh from_cells(
    std::vector<Vec3d> points, const std::vector<CellFaces>& cells,
    std::vector<Patch> patches,
    const std::function<Index(const std::vector<Index>&)>& boundary_patch) {
  std::map<std::vector<Index>, std::size_t> lookup;
  std::vector<Face> faces;
  for (Index c = 0; c < index(cells.size()); ++c) {
    std::vector<Index> vertices;
    for (const auto& ns : cells[c])
      vertices.insert(vertices.end(), ns.begin(), ns.end());
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()),
                   vertices.end());
    KUTTA_REQUIRE(!vertices.empty(), "construction: cell has no vertices");
    auto origin = points[vertices[0]];
    Vec3d center{};
    for (auto n : vertices) center = center + (points[n] - origin);
    center = origin + center / static_cast<double>(vertices.size());
    for (auto ns : cells[c]) {
      auto key = ns;
      std::sort(key.begin(), key.end());
      auto [it, added] = lookup.emplace(key, faces.size());
      if (added) {
        Vec3d fc{}, av{};
        auto r = points[ns[0]];
        for (auto n : ns) fc = fc + (points[n] - r);
        fc = r + fc / static_cast<double>(ns.size());
        for (std::size_t i = 0; i < ns.size(); ++i)
          av = av +
               cross(points[ns[i]] - fc, points[ns[(i + 1) % ns.size()]] - fc);
        if (dot(av, fc - center) < 0) std::reverse(ns.begin(), ns.end());
        faces.push_back({std::move(ns), c, -1, -1});
      } else {
        auto& f = faces[it->second];
        KUTTA_REQUIRE(f.neighbour < 0,
                      "construction: more than two cells share a face");
        f.neighbour = c;
      }
    }
  }
  for (auto& f : faces)
    if (f.neighbour < 0) f.patch = boundary_patch(f.nodes);
  return pack(std::move(points), index(cells.size()), std::move(faces),
              std::move(patches));
}
}  // namespace kutta::detail

namespace kutta {
CellFaces build_cell_faces(const PolyMesh& m) {
  detail::index(m.owner.size() + m.neighbour.size());
  CellFaces a;
  a.offsets.resize(static_cast<std::size_t>(m.n_cells) + 1);
  for (auto c : m.owner) ++a.offsets[c + 1];
  for (auto c : m.neighbour) ++a.offsets[c + 1];
  for (Index c = 0; c < m.n_cells; ++c) a.offsets[c + 1] += a.offsets[c];
  a.faces.resize(a.offsets.back());
  a.signs.resize(a.faces.size());
  auto cursor = a.offsets;
  for (Index f = 0; f < m.n_faces(); ++f) {
    auto i = cursor[m.owner[f]]++;
    a.faces[i] = f;
    a.signs[i] = 1;
    if (f < m.n_internal_faces()) {
      i = cursor[m.neighbour[f]]++;
      a.faces[i] = f;
      a.signs[i] = -1;
    }
  }
  return a;
}
CellNeighbours build_cell_neighbours(const PolyMesh& m) {
  auto faces = build_cell_faces(m);
  CellNeighbours a;
  a.offsets.push_back(0);
  for (Index c = 0; c < m.n_cells; ++c) {
    std::vector<Index> row;
    for (auto i = faces.offsets[c]; i < faces.offsets[c + 1]; ++i) {
      auto f = faces.faces[i];
      if (f < m.n_internal_faces())
        row.push_back(m.owner[f] == c ? m.neighbour[f] : m.owner[f]);
    }
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
    a.cells.insert(a.cells.end(), row.begin(), row.end());
    a.offsets.push_back(detail::index(a.cells.size()));
  }
  return a;
}
MeshPermutation renumber_cells_rcm(PolyMesh& m) {
  validate(m);
  auto adj = build_cell_neighbours(m);
  MeshPermutation p;
  std::vector<bool> seen(m.n_cells);
  auto degree = [&](Index c) { return adj.offsets[c + 1] - adj.offsets[c]; };
  std::vector<Index> seeds(m.n_cells);
  std::iota(seeds.begin(), seeds.end(), 0);
  std::stable_sort(seeds.begin(), seeds.end(),
                   [&](Index a, Index b) { return degree(a) < degree(b); });
  for (auto seed : seeds)
    if (!seen[seed]) {
      std::vector<Index> component;
      std::queue<Index> q;
      q.push(seed);
      seen[seed] = true;
      while (!q.empty()) {
        auto c = q.front();
        q.pop();
        component.push_back(c);
        std::vector<Index> row;
        for (auto i = adj.offsets[c]; i < adj.offsets[c + 1]; ++i)
          if (!seen[adj.cells[i]]) row.push_back(adj.cells[i]);
        std::stable_sort(row.begin(), row.end(), [&](Index a, Index b) {
          return degree(a) < degree(b);
        });
        for (auto n : row)
          if (!seen[n]) {
            seen[n] = true;
            q.push(n);
          }
      }
      p.cells_new_to_old.insert(p.cells_new_to_old.end(), component.rbegin(),
                                component.rend());
    }
  std::vector<Index> old_to_new(m.n_cells);
  for (Index c = 0; c < m.n_cells; ++c) old_to_new[p.cells_new_to_old[c]] = c;
  // Faces are permuted in place of the CSR arrays: one detail::Face (and node
  // vector) per face would cost millions of small allocations that macOS
  // keeps resident after they are freed. The result equals detail::pack of
  // the same faces: it would keep this order and these patch ranges.
  const auto faces = m.n_faces(), internal = m.n_internal_faces();
  std::vector<Index> patch(faces, -1), owner(faces), neighbour(faces, -1);
  std::vector<bool> flipped(faces);
  for (Index i = 0; i < detail::index(m.patches.size()); ++i)
    for (Index f = m.patches[i].start;
         f < m.patches[i].start + m.patches[i].size; ++f)
      patch[f] = i;
  for (Index f = 0; f < faces; ++f) {
    owner[f] = old_to_new[m.owner[f]];
    if (f < internal) neighbour[f] = old_to_new[m.neighbour[f]];
    if (neighbour[f] >= 0 && owner[f] > neighbour[f]) {
      std::swap(owner[f], neighbour[f]);
      flipped[f] = true;
    }
  }
  p.faces_new_to_old.resize(faces);
  std::iota(p.faces_new_to_old.begin(), p.faces_new_to_old.end(), 0);
  std::stable_sort(
      p.faces_new_to_old.begin(), p.faces_new_to_old.end(),
      [&](Index a, Index b) {
        return std::tuple(neighbour[a] < 0, patch[a], owner[a], neighbour[a]) <
               std::tuple(neighbour[b] < 0, patch[b], owner[b], neighbour[b]);
      });
  auto zones = m.cell_zones;
  for (auto& z : zones)
    for (auto& c : z.cells) c = old_to_new[c];
  auto face_zones=std::move(m.internal_face_zones);
  std::vector<Index> inverse_faces(m.n_faces());
  for(Index f=0;f<m.n_faces();++f) inverse_faces[p.faces_new_to_old[f]]=f;
  for(auto& zone:face_zones) for(auto& face:zone.faces) {
    const auto f=face.face;
    KUTTA_REQUIRE(f>=0 && f<m.n_internal_faces(),"internal face-zone index outside internal faces");
    if(old_to_new[m.owner[f]]>old_to_new[m.neighbour[f]]) face.direction=-face.direction;
    face.face=inverse_faces[f];
  }
  PolyMesh r;
  r.points = std::move(m.points);
  r.n_cells = m.n_cells;
  r.cell_zones = std::move(zones);
  r.owner.resize(faces);
  r.neighbour.resize(internal);
  r.face_offsets.resize(std::size_t(faces) + 1);
  r.face_nodes.resize(m.face_nodes.size());
  Index cursor = 0;
  for (Index i = 0; i < faces; ++i) {
    const auto f = p.faces_new_to_old[i];
    r.owner[i] = owner[f];
    if (i < internal) r.neighbour[i] = neighbour[f];
    const auto ns = detail::nodes(m, f);
    if (flipped[f]) std::reverse_copy(ns.begin(), ns.end(), r.face_nodes.begin() + cursor);
    else std::copy(ns.begin(), ns.end(), r.face_nodes.begin() + cursor);
    cursor += detail::index(ns.size());
    r.face_offsets[i + 1] = cursor;
  }
  r.patches = std::move(m.patches);
  m = std::move(r);
  m.internal_face_zones=std::move(face_zones);
  validate(m);
  return p;
}
}  // namespace kutta
