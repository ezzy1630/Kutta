#include "kutta/mesh/polyhedral.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>

#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta {
namespace {
using namespace detail;

std::uint64_t edge_key(Index a, Index b) {
  if (a > b) std::swap(a, b);
  return (std::uint64_t(std::uint32_t(a)) << 32) | std::uint32_t(b);
}

// Area vector of a polygon as a fan about its node mean.
Vec3d polygon_area(const std::vector<Vec3d>& points, const std::vector<Index>& polygon) {
  Vec3d centre{};
  for (auto n : polygon) centre = centre + points[n];
  centre = centre / double(polygon.size());
  Vec3d area{};
  for (std::size_t k = 0; k < polygon.size(); ++k)
    area = area + cross(points[polygon[k]] - centre, points[polygon[(k + 1) % polygon.size()]] - centre) * 0.5;
  return area;
}

class Dual {
 public:
  explicit Dual(const PolyMesh& mesh) : mesh_(mesh), cells_(build_cell_faces(mesh)), tet_(mesh.n_cells, 0) {
    const Index n_cells = mesh.n_cells;
    for (Index c = 0; c < n_cells; ++c) {
      bool all_triangles = cells_.offsets[c + 1] - cells_.offsets[c] == 4;
      for (Index j = cells_.offsets[c]; all_triangles && j < cells_.offsets[c + 1]; ++j)
        all_triangles = mesh.face_offsets[cells_.faces[j] + 1] - mesh.face_offsets[cells_.faces[j]] == 3;
      tet_[c] = all_triangles;
      tets += all_triangles;
    }
    if (!tets) return;
    zone_of_.assign(n_cells, -1);
    for (Index z = 0; z < index(mesh.cell_zones.size()); ++z)
      for (auto c : mesh.cell_zones[z].cells) zone_of_[c] = z;
    for (const auto& zone : mesh.internal_face_zones)
      for (const auto& f : zone.faces)
        KUTTA_REQUIRE(!tet_[mesh.owner[f.face]] && !tet_[mesh.neighbour[f.face]],
                      "polyhedral conversion: internal face zone '" + zone.name + "' touches tetrahedra");
    face_patch_.assign(mesh.n_faces(), -1);
    for (Index p = 0; p < index(mesh.patches.size()); ++p)
      for (Index f = mesh.patches[p].start; f < mesh.patches[p].start + mesh.patches[p].size; ++f) face_patch_[f] = p;

    // All new points up front, so rebuilding with more split vertices adds none.
    points_ = mesh.points;
    vertex_zone_.assign(mesh.n_points(), -2);
    tet_centre_.assign(n_cells, -1);
    face_centre_.assign(mesh.n_faces(), -1);
    std::vector<Index> count(mesh.n_points() + 1, 0);
    for (Index c = 0; c < n_cells; ++c) {
      if (!tet_[c]) continue;
      std::array<Index, 4> vertices{-1, -1, -1, -1};
      int n_vertices = 0;
      for (Index j = cells_.offsets[c]; j < cells_.offsets[c + 1]; ++j) {
        const auto f = cells_.faces[j];
        const auto n = nodes(mesh, f);
        if (face_centre_[f] < 0) face_centre_[f] = add((mesh.points[n[0]] + mesh.points[n[1]] + mesh.points[n[2]]) / 3.0);
        for (int k = 0; k < 3; ++k) {
          const auto key = edge_key(n[k], n[(k + 1) % 3]);
          auto [it, inserted] = ring_.try_emplace(key);
          if (inserted) midpoint_[key] = add((mesh.points[n[k]] + mesh.points[n[(k + 1) % 3]]) * 0.5);
          it->second.push_back({c, f});
          if (std::find(vertices.begin(), vertices.begin() + n_vertices, n[k]) == vertices.begin() + n_vertices) vertices[n_vertices++] = n[k];
        }
      }
      KUTTA_REQUIRE(n_vertices == 4, "polyhedral conversion: tetrahedron without four distinct vertices");
      Vec3d centre{};
      for (auto v : vertices) {
        centre = centre + mesh.points[v];
        KUTTA_REQUIRE(vertex_zone_[v] == -2 || vertex_zone_[v] == zone_of_[c],
                      "polyhedral conversion: tetrahedra of different cell zones share a vertex");
        vertex_zone_[v] = zone_of_[c];
        ++count[v + 1];
      }
      tet_centre_[c] = add(centre / 4.0);
    }
    // vertex -> tets around it (CSR), in cell order.
    for (Index v = 0; v < mesh.n_points(); ++v) count[v + 1] += count[v];
    vertex_tets_offsets_ = count;
    vertex_tets_.resize(count.back());
    for (Index c = 0; c < n_cells; ++c) {
      if (!tet_[c]) continue;
      for (auto v : tet_vertices(c)) vertex_tets_[count[v]++] = c;
    }
    keys_.reserve(ring_.size());
    for (const auto& [key, _] : ring_) keys_.push_back(key);
    std::sort(keys_.begin(), keys_.end());
  }

  Index tets = 0;

  // Dual vertices whose cell fails the orientation checks become one hex
  // quarter per tetrahedron, as in a tet-to-hex split.
  PolyhedralConversion build(const std::vector<char>& split, std::vector<Index>& bad_vertices) const {
    const Index n_cells = mesh_.n_cells, n_internal = mesh_.n_internal_faces();
    std::vector<Index> cell_id(n_cells, -1), base(mesh_.n_points(), -1), vertex_of;
    Index next = 0;
    for (Index c = 0; c < n_cells; ++c)
      if (!tet_[c]) cell_id[c] = next++;
    const Index kept = next;
    vertex_of.assign(kept, -1);
    Index polyhedra = 0;
    for (Index v = 0; v < mesh_.n_points(); ++v) {
      if (vertex_zone_[v] == -2) continue;
      base[v] = next;
      const auto n = split[v] ? vertex_tets_offsets_[v + 1] - vertex_tets_offsets_[v] : 1;
      polyhedra += !split[v];
      next += n;
      vertex_of.insert(vertex_of.end(), std::size_t(n), v);
    }
    const auto cell_of = [&](Index v, Index t) {
      if (!split[v]) return base[v];
      const auto first = vertex_tets_.begin() + vertex_tets_offsets_[v], last = vertex_tets_.begin() + vertex_tets_offsets_[v + 1];
      return base[v] + Index(std::lower_bound(first, last, t) - first);
    };
    const auto mid = [&](Index a, Index b) { return midpoint_.at(edge_key(a, b)); };

    std::vector<Face> out;
    out.reserve(std::size_t(mesh_.n_faces()) + ring_.size());
    for (Index f = 0; f < mesh_.n_faces(); ++f) {
      const auto owner = mesh_.owner[f];
      const auto neighbour = f < n_internal ? mesh_.neighbour[f] : -1;
      const bool owner_tet = tet_[owner], neighbour_tet = neighbour >= 0 && tet_[neighbour];
      const auto n = nodes(mesh_, f);
      if (!owner_tet && !neighbour_tet) {
        // A kept face along a tet edge gains its midpoint, as the dual faces
        // there do, so the mesh stays conforming.
        std::vector<Index> polygon;
        for (std::size_t k = 0; k < n.size(); ++k) {
          polygon.push_back(n[k]);
          const auto b = n[(k + 1) % n.size()];
          if (ring_.contains(edge_key(n[k], b))) polygon.push_back(mid(n[k], b));
        }
        out.push_back({std::move(polygon), cell_id[owner], neighbour >= 0 ? cell_id[neighbour] : -1, face_patch_[f]});
        continue;
      }
      // Each vertex's quarter of a tet triangle, wound as the triangle, so it
      // still points out of the input owner. Between two tets it separates
      // two quarters of one split vertex; elsewhere it is a boundary or
      // interface face of the vertex's cell.
      for (int k = 0; k < 3; ++k) {
        const auto a = n[k], b = n[(k + 1) % 3], c = n[(k + 2) % 3];
        if (owner_tet && neighbour_tet && !split[a]) continue;
        Face q{{a, mid(a, b), face_centre_[f], mid(c, a)}, -1, -1, face_patch_[f]};
        q.owner = owner_tet ? cell_of(a, owner) : cell_id[owner];
        if (neighbour >= 0) q.neighbour = neighbour_tet ? cell_of(a, neighbour) : cell_id[neighbour];
        out.push_back(std::move(q));
      }
    }

    // One face per tet edge: around the edge, alternate triangle centroids and
    // tet centroids. An edge on the region's surface opens the ring there and
    // the polygon passes through the edge midpoint. A split vertex takes one
    // quad per tetrahedron instead.
    for (auto key : keys_) {
      const auto& entries = ring_.at(key);
      const Index v = Index(key >> 32), w = Index(key & 0xffffffffu);
      const auto oriented = [&](std::vector<Index> polygon, Index owner, Index neighbour) {
        if (dot(polygon_area(points_, polygon), mesh_.points[w] - mesh_.points[v]) < 0) std::reverse(polygon.begin(), polygon.end());
        out.push_back({std::move(polygon), owner, neighbour, -1});
      };
      std::vector<Index> ring_tets, ring_faces;
      for (const auto& [c, f] : entries) {
        if (std::find(ring_tets.begin(), ring_tets.end(), c) == ring_tets.end()) ring_tets.push_back(c);
        if (std::find(ring_faces.begin(), ring_faces.end(), f) == ring_faces.end()) ring_faces.push_back(f);
      }
      const auto other_face = [&](Index c, Index f) {
        for (const auto& [d, g] : entries)
          if (d == c && g != f) return g;
        throw Error("polyhedral conversion: tetrahedron edge without two faces");
      };
      if (split[v] || split[w]) {
        for (auto c : ring_tets) {
          Index f1 = -1;
          for (const auto& [d, g] : entries)
            if (d == c) f1 = g;
          oriented({mid(v, w), face_centre_[f1], tet_centre_[c], face_centre_[other_face(c, f1)]}, cell_of(v, c), cell_of(w, c));
        }
        continue;
      }
      std::vector<char> used(ring_tets.size(), 0);
      const auto free_tets = [&](Index f) {
        std::vector<Index> t;
        for (const auto& [c, g] : entries)
          if (g == f && !used[std::find(ring_tets.begin(), ring_tets.end(), c) - ring_tets.begin()]) t.push_back(c);
        return t;
      };
      const auto face_tets = [&](Index f) {
        return std::count_if(entries.begin(), entries.end(), [&](const auto& e) { return e[1] == f; });
      };
      for (;;) {
        Index start = -1;
        bool open = false;
        for (auto f : ring_faces) {
          if (free_tets(f).empty()) continue;
          if (face_tets(f) == 1) {
            start = f;
            open = true;
            break;
          }
          if (start < 0) start = f;
        }
        if (start < 0) break;
        std::vector<Index> polygon;
        if (open) polygon.push_back(mid(v, w));
        polygon.push_back(face_centre_[start]);
        for (Index f = start;;) {
          const auto t = free_tets(f);
          if (t.empty()) break;
          used[std::find(ring_tets.begin(), ring_tets.end(), t.front()) - ring_tets.begin()] = 1;
          polygon.push_back(tet_centre_[t.front()]);
          f = other_face(t.front(), f);
          if (f == start) break;
          polygon.push_back(face_centre_[f]);
        }
        oriented(std::move(polygon), base[v], base[w]);
      }
    }

    std::vector<CellZone> zones = mesh_.cell_zones;
    for (Index z = 0; z < index(zones.size()); ++z) {
      auto& zc = zones[z].cells;
      zc.clear();
      for (Index c = 0; c < next; ++c)
        if (c < kept ? false : vertex_zone_[vertex_of[c]] == z) zc.push_back(c);
      for (Index c = 0; c < n_cells; ++c)
        if (!tet_[c] && zone_of_[c] == z) zc.push_back(cell_id[c]);
      std::sort(zc.begin(), zc.end());
    }
    auto converted = pack(points_, next, std::move(out), mesh_.patches, std::move(zones));

    // Every face must lie in front of its owner's centroid and behind its
    // neighbour's (validate()). A dual cell wrapped around a sharp convex edge
    // of the region can put its centroid behind one of its faces.
    bad_vertices.clear();
    const auto g = compute_geometry(converted);
    const auto mark = [&](Index c) {
      if (c >= kept && !split[vertex_of[c]]) bad_vertices.push_back(vertex_of[c]);
    };
    for (Index f = 0; f < converted.n_faces(); ++f) {
      const auto o = converted.owner[f];
      bool ok = dot(g.face_area[f], g.face_centroid[f] - g.cell_centroid[o]) > 0;
      if (f < converted.n_internal_faces()) {
        const auto n = converted.neighbour[f];
        ok = ok && dot(g.face_area[f], g.cell_centroid[n] - g.face_centroid[f]) > 0;
        if (!ok) mark(n);
      }
      if (!ok) mark(o);
    }
    std::sort(bad_vertices.begin(), bad_vertices.end());
    bad_vertices.erase(std::unique(bad_vertices.begin(), bad_vertices.end()), bad_vertices.end());
    return {std::move(converted), tets, polyhedra, kept};
  }

  Index n_points() const { return mesh_.n_points(); }

 private:
  Index add(Vec3d p) {
    points_.push_back(p);
    return index(points_.size() - 1);
  }
  std::array<Index, 4> tet_vertices(Index c) const {
    std::array<Index, 4> out{-1, -1, -1, -1};
    int n = 0;
    for (Index j = cells_.offsets[c]; j < cells_.offsets[c + 1]; ++j)
      for (auto v : nodes(mesh_, cells_.faces[j]))
        if (std::find(out.begin(), out.begin() + n, v) == out.begin() + n) out[n++] = v;
    return out;
  }

  const PolyMesh& mesh_;
  kutta::CellFaces cells_;
  std::vector<char> tet_;
  std::vector<Index> zone_of_, face_patch_, vertex_zone_, tet_centre_, face_centre_, vertex_tets_offsets_, vertex_tets_;
  std::vector<Vec3d> points_;
  std::unordered_map<std::uint64_t, std::vector<std::array<Index, 2>>> ring_;
  std::unordered_map<std::uint64_t, Index> midpoint_;
  std::vector<std::uint64_t> keys_;
};

PolyMesh compact_points(PolyMesh mesh) {
  std::vector<Index> id(mesh.points.size(), -1);
  std::vector<Vec3d> points;
  for (auto& n : mesh.face_nodes) {
    if (id[n] < 0) {
      id[n] = index(points.size());
      points.push_back(mesh.points[n]);
    }
    n = id[n];
  }
  mesh.points = std::move(points);
  return mesh;
}
}  // namespace

PolyhedralConversion convert_tets_to_polyhedra(const PolyMesh& mesh) {
  const Dual dual(mesh);
  if (!dual.tets) return {mesh, 0, 0, mesh.n_cells};
  std::vector<char> split(mesh.n_points(), 0);
  std::vector<Index> bad;
  auto result = dual.build(split, bad);
  // Splitting changes only the split cells and their neighbours' shared faces.
  for (int round = 0; !bad.empty() && round < 8; ++round) {
    for (auto v : bad) split[v] = 1;
    result = dual.build(split, bad);
  }
  for (auto s : split) result.split_vertices += s;
  result.mesh = compact_points(std::move(result.mesh));
  validate(result.mesh);
  return result;
}

}  // namespace kutta
