#include <numeric>
#include <tuple>

#include "kutta/io/fluent_detail.hpp"

namespace kutta::detail {
std::string fluent_type(int id) {
  switch (id) {
    case 2:
      return "interior";
    case 3:
      return "wall";
    case 4:
      return "pressure-inlet";
    case 5:
      return "pressure-outlet";
    case 7:
      return "symmetry";
    case 8:
      return "periodic-shadow";
    case 9:
      return "pressure-far-field";
    case 10:
      return "velocity-inlet";
    case 12:
      return "periodic";
    case 14:
      return "fan";
    case 20:
      return "mass-flow-inlet";
    case 24:
      return "interface";
    case 31:
      return "parent";
    case 36:
      return "outflow";
    case 37:
      return "axis";
    default:
      return "fluent-type-" + std::to_string(id);
  }
}
PatchKind patch_kind(const std::string& t) {
  if (t == "wall") return PatchKind::Wall;
  if (t == "velocity-inlet") return PatchKind::VelocityInlet;
  if (t == "pressure-inlet") return PatchKind::PressureInlet;
  if (t == "pressure-outlet" || t == "outflow")
    return PatchKind::PressureOutlet;
  if (t == "symmetry" || t == "symmetry-plane") return PatchKind::Symmetry;
  if (t == "pressure-far-field") return PatchKind::Farfield;
  return PatchKind::Other;
}
FluentMeshData assemble_fluent(const FluentRaw& r, double thickness) {
  KUTTA_REQUIRE(r.dimension == 2 || r.dimension == 3,
                "Fluent dimension must be 2 or 3");
  KUTTA_REQUIRE(!r.nodes.empty() && !r.cells.empty() && !r.faces.empty(),
                "Fluent mesh requires nodes, active cells and faces "
                "(surface-only files are unsupported)");
  // Points and cells are numbered in ascending source id, so a source id's
  // index is its rank in the table.
  std::vector<Vec3d> points;
  points.reserve(r.nodes.size());
  r.nodes.each([&](Index id, Vec3d p) {
    KUTTA_REQUIRE(id > 0 && finite(p), "Fluent: invalid source node");
    points.push_back(p);
  });
  r.cells.each([&](Index id, int) {
    KUTTA_REQUIRE(id > 0, "Fluent: invalid source cell id");
  });
  auto lookup = [](const auto& table, Index id, const char* what) -> Index {
    KUTTA_REQUIRE(table.find(id) != nullptr, std::string("Fluent: unknown ") + what +
                                       " id " + std::to_string(id));
    return index(table.rank(id));
  };
  std::vector<Patch> patches;
  std::map<int, Index> patch_map;
  r.faces.each([&](Index, const FluentFace& f) {
    if (f.c0 == 0 || f.c1 == 0) {
      if (!patch_map.count(f.zone)) {
        auto z = r.zones.find(f.zone);
        KUTTA_REQUIRE(z != r.zones.end(),
                      "Fluent: boundary has no zone metadata");
        patch_map.emplace(f.zone, index(patches.size()));
        patches.push_back({z->second.name.empty()
                               ? "zone" + std::to_string(f.zone)
                               : z->second.name,
                           patch_kind(z->second.type), z->second.type, f.zone});
      }
    }
  });
  std::vector<CellZone> zones;
  std::map<int, std::size_t> zone_map;
  Index cell = 0;
  r.cells.each([&](Index, int z) {
    if (!zone_map.count(z)) {
      auto it = r.zones.find(z);
      FluentZone meta = it == r.zones.end()
                            ? FluentZone{"zone" + std::to_string(z), "fluid"}
                            : it->second;
      zone_map.emplace(z, zones.size());
      zones.push_back({meta.name, meta.type, z, {}});
    }
    zones[zone_map[z]].cells.push_back(cell++);
  });
  // Faces in source-id order, as flat arrays; nodes already mapped to point
  // indices and wound so the normal leaves the owner.
  const auto n_faces = r.faces.size();
  std::vector<Index> owner(n_faces), neighbour(n_faces), patch(n_faces),
      node_offsets{0}, face_nodes;
  node_offsets.reserve(n_faces + 1);
  face_nodes.reserve(r.face_nodes.size());
  std::size_t k = 0;
  r.faces.each([&](Index id, const FluentFace& f) {
    KUTTA_REQUIRE(f.c0 != 0 || f.c1 != 0, "Fluent face " + std::to_string(id) +
                                              " has no adjacent cell");
    auto o =
        f.c0 ? lookup(r.cells, f.c0, "cell") : lookup(r.cells, f.c1, "cell");
    auto nb = f.c0 && f.c1 ? lookup(r.cells, f.c1, "cell") : -1;
    bool owner_is_c0 = f.c0 != 0;
    if (nb >= 0 && o > nb) {
      std::swap(o, nb);
      owner_is_c0 = false;
    }
    const auto begin = face_nodes.size();
    for (auto n : r.nodes_of(f)) face_nodes.push_back(lookup(r.nodes, n, "node"));
    // Legacy/CFF file winding points into c0 in 3D. In 2D the
    // directed edge leaves c0 on its left, as required by extrude_2d.
    if ((r.dimension == 3 && owner_is_c0) || (r.dimension == 2 && !owner_is_c0))
      std::reverse(face_nodes.begin() + std::ptrdiff_t(begin), face_nodes.end());
    node_offsets.push_back(index(face_nodes.size()));
    owner[k] = o;
    neighbour[k] = nb;
    patch[k] = nb < 0 ? patch_map.at(f.zone) : -1;
    ++k;
  });
  // detail::pack's construction checks and stable face order (internal faces
  // by owner and neighbour, then boundary faces by patch), applied to the
  // flat arrays; source ids are tracked through it for periodic pairs and
  // internal face zones.
  for (std::size_t f = 0; f < n_faces; ++f)
    KUTTA_REQUIRE(neighbour[f] < 0 || owner[f] != neighbour[f],
                  "construction: face has identical adjacent cells");
  std::vector<Index> order(n_faces);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](Index a, Index b) {
    return std::tuple(neighbour[a] < 0, neighbour[a] < 0 ? patch[a] : owner[a],
                      neighbour[a]) <
           std::tuple(neighbour[b] < 0, neighbour[b] < 0 ? patch[b] : owner[b],
                      neighbour[b]);
  });
  std::vector<Index> sorted(n_faces);  // source-order position -> face index
  for (std::size_t i = 0; i < n_faces; ++i) sorted[std::size_t(order[i])] = index(i);
  FluentMeshData result;
  if (r.dimension == 2) {
    PlanarMesh p;
    p.points = std::move(points);
    p.n_cells = index(r.cells.size());
    p.patches = patches;
    p.cell_zones = zones;
    // Extrusion appends caps after all original side faces. Its stable pack
    // places caps after the original boundary patches, preserving this map.
    for (std::size_t f = 0; f < n_faces; ++f) {
      KUTTA_REQUIRE(node_offsets[f + 1] - node_offsets[f] == 2,
                    "Fluent 2D face must have two nodes");
      p.edges.push_back({face_nodes[node_offsets[f]], face_nodes[node_offsets[f] + 1]});
      p.owner.push_back(owner[f]);
      p.neighbour.push_back(neighbour[f]);
      p.edge_patch.push_back(patch[f]);
    }
    result.mesh = extrude_2d(p, thickness);
  } else {
    auto& m = result.mesh;
    m.points = std::move(points);
    m.n_cells = index(r.cells.size());
    m.cell_zones = std::move(zones);
    for (auto& q : patches) q.size = 0;
    m.owner.reserve(n_faces);
    m.face_offsets.reserve(n_faces + 1);
    m.face_nodes.reserve(face_nodes.size());
    for (auto f : order) {
      m.owner.push_back(owner[f]);
      if (neighbour[f] >= 0)
        m.neighbour.push_back(neighbour[f]);
      else
        ++patches[patch[f]].size;
      m.face_nodes.insert(m.face_nodes.end(), face_nodes.begin() + node_offsets[f],
                          face_nodes.begin() + node_offsets[f + 1]);
      m.face_offsets.push_back(index(m.face_nodes.size()));
    }
    auto start = m.n_internal_faces();
    for (auto& q : patches) {
      q.start = start;
      start += q.size;
    }
    m.patches = std::move(patches);
    std::vector<Index>().swap(face_nodes);
    std::vector<Index>().swap(node_offsets);
    validate(m);
  }
  std::vector<bool> paired(n_faces);
  // Preserve internal face-zone identity; face numbering alone loses fan and
  // porous-jump BCs after owner/neighbour sorting. File winding points into c0.
  std::map<int,std::size_t> internal_zones;
  k = 0;
  r.faces.each([&](Index, const FluentFace& f) {
    const auto position = k++;
    if (!f.c0 || !f.c1) return;
    const auto z=r.zones.find(f.zone);
    KUTTA_REQUIRE(z!=r.zones.end(),"Fluent: internal face has no zone metadata");
    if(z->second.type=="interior") return; // ordinary faces need no duplicated identity arrays
    if(!internal_zones.contains(f.zone)) {
      internal_zones.emplace(f.zone,result.mesh.internal_face_zones.size());
      result.mesh.internal_face_zones.push_back({z->second.name,z->second.type,f.zone,{}});
    }
    const bool owner_is_c0=r.cells.rank(f.c0)<r.cells.rank(f.c1);
    result.mesh.internal_face_zones[internal_zones.at(f.zone)].faces.push_back({sorted[position],owner_is_c0?-1:1});
  });
  for (auto pair : r.periodic) {
    const auto *f = r.faces.find(pair.face), *s = r.faces.find(pair.shadow);
    KUTTA_REQUIRE(f && s, "Fluent periodic pair references missing face");
    KUTTA_REQUIRE(
        f->zone == pair.zone && s->zone == pair.shadow_zone,
        "Fluent periodic pair zone mismatch");
    KUTTA_REQUIRE(
        (!f->c0 || !f->c1) && (!s->c0 || !s->c1),
        "Fluent periodic pair must reference boundary faces");
    const auto a = r.faces.rank(pair.face), b = r.faces.rank(pair.shadow);
    KUTTA_REQUIRE(!paired[a] && !paired[b] && a != b,
                  "Fluent periodic face paired more than once");
    paired[a] = paired[b] = true;
    pair.face = sorted[a];
    pair.shadow = sorted[b];
    result.periodic_faces.push_back(pair);
  }
  return result;
}
}  // namespace kutta::detail
