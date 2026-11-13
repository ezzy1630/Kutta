// Fluent legacy ASCII mesh writer (sections 0, 2, 10, 12, 13, 45).
#include <fstream>
#include <iomanip>
#include <numeric>

#include "kutta/io/fluent_detail.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta {
namespace {
int boundary_code(const Patch& patch) {
  for (int code = 2; code <= 37; ++code)
    if (detail::fluent_type(code) == patch.source_type) return code;
  switch (patch.kind) {
    case PatchKind::Wall: return 3;
    case PatchKind::VelocityInlet: return 10;
    case PatchKind::PressureOutlet: return 5;
    case PatchKind::PressureInlet: return 4;
    case PatchKind::Symmetry: return 7;
    case PatchKind::Farfield: return 9;
    case PatchKind::Other: break;
  }
  throw Error("write_fluent_msh: patch '" + patch.name + "' has no Fluent boundary type");
}

// Fluent's element types: 2 tet, 4 hex, 5 pyramid, 6 wedge, 7 polyhedron.
int element_type(const PolyMesh& mesh, const CellFaces& cells, Index cell) {
  Index triangles = 0, quads = 0;
  const auto count = cells.offsets[cell + 1] - cells.offsets[cell];
  for (Index j = cells.offsets[cell]; j < cells.offsets[cell + 1]; ++j) {
    const auto size = mesh.face_offsets[cells.faces[j] + 1] - mesh.face_offsets[cells.faces[j]];
    triangles += size == 3;
    quads += size == 4;
  }
  if (count == 4 && triangles == 4) return 2;
  if (count == 5 && triangles == 2 && quads == 3) return 6;
  if (count == 5 && triangles == 4 && quads == 1) return 5;
  if (count == 6 && quads == 6) return 4;
  return 7;
}
}  // namespace

void write_fluent_msh(const std::filesystem::path& path, const PolyMesh& mesh, std::string_view comment) {
  validate(mesh);
  KUTTA_REQUIRE(mesh.internal_face_zones.empty(), "write_fluent_msh: internal face zones are not supported");
  // Fluent cell zones are contiguous id ranges: order cells by zone.
  std::vector<Index> zone_of(mesh.n_cells, Index(mesh.cell_zones.size()));
  for (Index z = 0; z < detail::index(mesh.cell_zones.size()); ++z)
    for (auto c : mesh.cell_zones[z].cells) zone_of[c] = z;
  std::vector<Index> order(mesh.n_cells), id(mesh.n_cells);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](Index a, Index b) { return zone_of[a] < zone_of[b]; });
  for (Index i = 0; i < mesh.n_cells; ++i) id[order[i]] = i;
  KUTTA_REQUIRE(mesh.cell_zones.empty() || std::all_of(zone_of.begin(), zone_of.end(), [&](Index z) { return z < Index(mesh.cell_zones.size()); }),
                "write_fluent_msh: every cell needs a cell zone when zones are present");
  const auto cells = build_cell_faces(mesh);
  std::ofstream file(path);
  KUTTA_REQUIRE(bool(file), "write_fluent_msh: cannot open " + path.string());
  std::string text(comment);
  for (auto& c : text)
    if (c == '"' || c == '\n') c = ' ';
  file << std::setprecision(17) << "(0 \"" << text << "\")\n(2 3)\n";
  file << std::hex << "(10 (0 1 " << mesh.points.size() << " 0 3))\n(10 (1 1 " << mesh.points.size() << " 1 3)(\n" << std::dec;
  for (auto p : mesh.points) file << p.x << ' ' << p.y << ' ' << p.z << '\n';
  file << "))\n" << std::hex << "(12 (0 1 " << mesh.n_cells << " 0))\n";
  const auto zone_count = std::max<Index>(1, detail::index(mesh.cell_zones.size()));
  int next_id = 2;
  std::vector<int> cell_zone_ids;
  for (Index z = 0, first = 0; z < zone_count; ++z) {
    Index last = first;
    while (last < mesh.n_cells && (mesh.cell_zones.empty() || zone_of[order[last]] == z)) ++last;
    cell_zone_ids.push_back(next_id);
    if (last > first) {
      file << "(12 (" << next_id << ' ' << first + 1 << ' ' << last << " 1 0)(\n" << std::dec;
      for (Index i = first; i < last; ++i) file << element_type(mesh, cells, order[i]) << '\n';
      file << "))\n" << std::hex;
    }
    ++next_id;
    first = last;
  }
  file << "(13 (0 1 " << mesh.n_faces() << " 0))\n";
  const auto faces = [&](int zone, Index start, Index count, int type) {
    if (!count) return;
    file << std::hex << "(13 (" << zone << ' ' << start + 1 << ' ' << start + count << ' ' << type << " 0)(\n";
    for (Index f = start; f < start + count; ++f) {
      // Legacy Fluent winding points into c0; PolyMesh winding points out of the owner.
      const auto n = detail::nodes(mesh, f);
      Index c0 = id[mesh.owner[f]], c1 = f < mesh.n_internal_faces() ? id[mesh.neighbour[f]] : -1;
      file << n.size();
      for (auto j = n.size(); j-- > 0;) file << ' ' << n[j] + 1;
      file << ' ' << c0 + 1 << ' ' << c1 + 1 << '\n';
    }
    file << "))\n";
  };
  const int interior_id = next_id++;
  faces(interior_id, 0, mesh.n_internal_faces(), 2);
  std::vector<int> patch_ids;
  for (const auto& p : mesh.patches) {
    patch_ids.push_back(next_id);
    faces(next_id++, p.start, p.size, boundary_code(p));
  }
  file << std::dec;
  for (Index z = 0; z < zone_count; ++z) {
    if (!mesh.cell_zones.empty() && mesh.cell_zones[z].cells.empty()) continue;
    const auto name = mesh.cell_zones.empty() ? std::string("fluid") : mesh.cell_zones[z].name;
    const auto type = mesh.cell_zones.empty() || mesh.cell_zones[z].source_type.empty() ? std::string("fluid") : mesh.cell_zones[z].source_type;
    file << "(45 (" << cell_zone_ids[z] << ' ' << type << ' ' << name << ")())\n";
  }
  file << "(45 (" << interior_id << " interior interior)())\n";
  for (std::size_t i = 0; i < mesh.patches.size(); ++i)
    if (mesh.patches[i].size) file << "(45 (" << patch_ids[i] << ' ' << detail::fluent_type(boundary_code(mesh.patches[i])) << ' ' << mesh.patches[i].name << ")())\n";
  file.close();
  KUTTA_REQUIRE(bool(file), "write_fluent_msh: write failed for " + path.string());
}

}  // namespace kutta
