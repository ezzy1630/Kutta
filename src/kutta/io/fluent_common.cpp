#include <numeric>
#include <set>
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
  std::map<Index, Index> point_map, cell_map;
  std::vector<Vec3d> points;
  for (auto [id, p] : r.nodes) {
    KUTTA_REQUIRE(id > 0 && finite(p), "Fluent: invalid source node");
    point_map.emplace(id, index(points.size()));
    points.push_back(p);
  }
  for (auto [id, z] : r.cells) {
    (void)z;
    KUTTA_REQUIRE(id > 0, "Fluent: invalid source cell id");
    cell_map.emplace(id, index(cell_map.size()));
  }
  auto lookup = [](const auto& map, Index id, const char* what) -> Index {
    auto it = map.find(id);
    KUTTA_REQUIRE(it != map.end(), std::string("Fluent: unknown ") + what +
                                       " id " + std::to_string(id));
    return it->second;
  };
  std::vector<Patch> patches;
  std::map<int, Index> patch_map;
  for (const auto& [id, f] : r.faces)
    if (f.c0 == 0 || f.c1 == 0) {
      (void)id;
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
  std::vector<CellZone> zones;
  std::map<int, std::size_t> zone_map;
  for (auto [id, z] : r.cells) {
    if (!zone_map.count(z)) {
      auto it = r.zones.find(z);
      FluentZone meta = it == r.zones.end()
                            ? FluentZone{"zone" + std::to_string(z), "fluid"}
                            : it->second;
      zone_map.emplace(z, zones.size());
      zones.push_back({meta.name, meta.type, z, {}});
    }
    zones[zone_map[z]].cells.push_back(cell_map.at(id));
  }
  std::vector<Face> faces;
  std::vector<Index> source_ids;
  for (const auto& [id, f] : r.faces) {
    KUTTA_REQUIRE(f.c0 != 0 || f.c1 != 0, "Fluent face " + std::to_string(id) +
                                              " has no adjacent cell");
    auto o =
        f.c0 ? lookup(cell_map, f.c0, "cell") : lookup(cell_map, f.c1, "cell");
    auto nb = f.c0 && f.c1 ? lookup(cell_map, f.c1, "cell") : -1;
    bool owner_is_c0 = f.c0 != 0;
    if (nb >= 0 && o > nb) {
      std::swap(o, nb);
      owner_is_c0 = false;
    }
    std::vector<Index> ns;
    for (auto n : f.nodes) ns.push_back(lookup(point_map, n, "node"));
    // Legacy/CFF file winding points into c0 in 3D. In 2D the
    // directed edge leaves c0 on its left, as required by extrude_2d.
    if ((r.dimension == 3 && owner_is_c0) || (r.dimension == 2 && !owner_is_c0))
      std::reverse(ns.begin(), ns.end());
    faces.push_back({std::move(ns), o, nb, nb < 0 ? patch_map.at(f.zone) : -1});
    source_ids.push_back(id);
  }
  // Track source ids through the exact stable ordering used by pack, so
  // periodic pairs survive both the internal-face sort and boundary patch
  // grouping.
  std::vector<std::size_t> order(faces.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) {
    const auto& x = faces[a];
    const auto& y = faces[b];
    return std::tuple(x.neighbour < 0, x.neighbour < 0 ? x.patch : x.owner,
                      x.neighbour) <
           std::tuple(y.neighbour < 0, y.neighbour < 0 ? y.patch : y.owner,
                      y.neighbour);
  });
  std::map<Index, Index> face_map;
  for (std::size_t i = 0; i < order.size(); ++i)
    face_map.emplace(source_ids[order[i]], index(i));
  FluentMeshData result;
  if (r.dimension == 2) {
    PlanarMesh p;
    p.points = std::move(points);
    p.n_cells = index(r.cells.size());
    p.patches = patches;
    p.cell_zones = zones;
    // Extrusion appends caps after all original side faces. Its stable pack
    // places caps after the original boundary patches, preserving this map.
    for (const auto& f : faces) {
      KUTTA_REQUIRE(f.nodes.size() == 2, "Fluent 2D face must have two nodes");
      p.edges.push_back({f.nodes[0], f.nodes[1]});
      p.owner.push_back(f.owner);
      p.neighbour.push_back(f.neighbour);
      p.edge_patch.push_back(f.patch);
    }
    result.mesh = extrude_2d(p, thickness);
  } else {
    result.mesh = pack(std::move(points), index(r.cells.size()),
                       std::move(faces), std::move(patches), std::move(zones));
    validate(result.mesh);
  }
  std::set<Index> paired;
  for (auto pair : r.periodic) {
    auto f = r.faces.find(pair.face), s = r.faces.find(pair.shadow);
    KUTTA_REQUIRE(f != r.faces.end() && s != r.faces.end(),
                  "Fluent periodic pair references missing face");
    KUTTA_REQUIRE(
        f->second.zone == pair.zone && s->second.zone == pair.shadow_zone,
        "Fluent periodic pair zone mismatch");
    KUTTA_REQUIRE(
        (!f->second.c0 || !f->second.c1) && (!s->second.c0 || !s->second.c1),
        "Fluent periodic pair must reference boundary faces");
    KUTTA_REQUIRE(
        paired.insert(pair.face).second && paired.insert(pair.shadow).second,
        "Fluent periodic face paired more than once");
    pair.face = face_map.at(pair.face);
    pair.shadow = face_map.at(pair.shadow);
    result.periodic_faces.push_back(pair);
  }
  return result;
}
}  // namespace kutta::detail
