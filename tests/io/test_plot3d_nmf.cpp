#include <doctest/doctest.h>

#include <array>
#include <map>
#include <set>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "test_support.hpp"
using namespace kutta;
using namespace kutta::detail;
namespace {
void boundary_checks(const PolyMesh& m) {
  validate(m);
  const auto g = compute_geometry(m);
  Vec3d closure{};
  double area = 0;
  std::set<std::vector<std::array<double, 3>>> boundary;
  for (Index f = m.n_internal_faces(); f < m.n_faces(); ++f) {
    closure = closure + g.face_area[f];
    area += norm(g.face_area[f]);
    std::vector<std::array<double, 3>> key;
    for (auto n : nodes(m, f)) {
      const auto p = m.points[n];
      // TMR wake coordinates may differ by roundoff across the cut.
      key.push_back({std::round(p.x * 1e10) / 1e10,
                     std::round(p.y * 1e10) / 1e10,
                     std::round(p.z * 1e10) / 1e10});
    }
    std::sort(key.begin(), key.end());
    CHECK(boundary.insert(std::move(key)).second);
  }
  CHECK(norm(closure) < 1e-10 * area);
}
}  // namespace
TEST_CASE("NASA NMF ranges supply exact BC patches and stitch the NACA wake") {
  const std::filesystem::path root = KUTTA_SOURCE_DIR "/validation/data/nasa_tmr";
  if (!std::filesystem::exists(root)) {
    MESSAGE("SKIP: external NASA TMR fixtures unavailable");
    return;
  }
  struct Case {
    const char* file;
    Index cells, internal, points;
    std::vector<std::pair<std::string, Index>> patches;
  };
  const std::vector<Case> cases{
      {"2DZP_flatplate/grids/flatplate_clust2_4levelsdown_35x25.p2dfmt", 816, 1574, 1750,
       {{"front", 816},
        {"back", 816},
        {"subsonic_inflow_pt_block1_f5_1-25_1-2", 24},
        {"back_pressure_block1_f6_1-25_1-2", 24},
        {"symmetry_z_strong_block1_f1_1-2_1-7", 6},
        {"viscous_solid_block1_f1_1-2_7-35", 28},
        {"farfield_riem_block1_f2_1-2_1-35", 34}}},
      {"2DZP_flatplate/grids/flatplate_clust2_3levelsdown_69x49.p2dfmt", 3264, 6412, 6762,
       {{"front", 3264},
        {"back", 3264},
        {"subsonic_inflow_pt_block1_f5_1-49_1-2", 48},
        {"back_pressure_block1_f6_1-49_1-2", 48},
        {"symmetry_z_strong_block1_f1_1-2_1-13", 12},
        {"viscous_solid_block1_f1_1-2_13-69", 56},
        {"farfield_riem_block1_f2_1-2_1-69", 68}}},
      {"2DZP_flatplate/grids/flatplate_clust2_2levelsdown_137x97.p2dfmt", 13056, 25880, 26578,
       {{"front", 13056},
        {"back", 13056},
        {"subsonic_inflow_pt_block1_f5_1-97_1-2", 96},
        {"back_pressure_block1_f6_1-97_1-2", 96},
        {"symmetry_z_strong_block1_f1_1-2_1-25", 24},
        {"viscous_solid_block1_f1_1-2_25-137", 112},
        {"farfield_riem_block1_f2_1-2_1-137", 136}}},
      {"2DBC_bump/grids/bump_4levelsdown_89x41.p2dfmt", 3520, 6912, 7298,
       {{"front", 3520},
        {"back", 3520},
        {"subsonic_inflow_pt_block1_f5_1-41_1-2", 40},
        {"back_pressure_block1_f6_1-41_1-2", 40},
        {"symmetry_z_strong_block1_f1_1-2_1-25", 24},
        {"viscous_solid_block1_f1_1-2_25-65", 40},
        {"symmetry_z_strong_block1_f1_1-2_65-89", 24},
        {"symmetry_z_strong_block1_f2_1-2_1-89", 88}}},
      {"2DBC_bump/grids/bump_3levelsdown_177x81.p2dfmt", 14080, 27904, 28674,
       {{"front", 14080},
        {"back", 14080},
        {"subsonic_inflow_pt_block1_f5_1-81_1-2", 80},
        {"back_pressure_block1_f6_1-81_1-2", 80},
        {"symmetry_z_strong_block1_f1_1-2_1-49", 48},
        {"viscous_solid_block1_f1_1-2_49-129", 80},
        {"symmetry_z_strong_block1_f1_1-2_129-177", 48},
        {"symmetry_z_strong_block1_f2_1-2_1-177", 176}}},
      {"2DN00_naca0012/grids/n0012_113-33.p2dfmt", 3584, 7048, 7408,
       {{"front", 3584},
        {"back", 3584},
        {"farfield_riem_block1_f5_1-33_1-2", 32},
        {"farfield_riem_block1_f6_1-33_1-2", 32},
        {"viscous_solid_block1_f1_1-2_25-89", 64},
        {"farfield_riem_block1_f2_1-2_1-113", 112}}},
      {"2DN00_naca0012/grids/n0012_225-65.p2dfmt", 14336, 28432, 29152,
       {{"front", 14336},
        {"back", 14336},
        {"farfield_riem_block1_f5_1-65_1-2", 64},
        {"farfield_riem_block1_f6_1-65_1-2", 64},
        {"viscous_solid_block1_f1_1-2_49-177", 128},
        {"farfield_riem_block1_f2_1-2_1-225", 224}}},
  };
  for (const auto& c : cases) {
    INFO(std::string(c.file));
    auto m = read_plot3d(root / c.file, .01);
    boundary_checks(m);
    CHECK(m.n_cells == c.cells);
    CHECK(m.n_internal_faces() == c.internal);
    CHECK(m.n_points() == c.points);
    REQUIRE(m.patches.size() == c.patches.size());
    for (std::size_t i = 0; i < c.patches.size(); ++i) {
      const auto& p = m.patches[i];
      CHECK(p.name == c.patches[i].first);
      CHECK(p.size == c.patches[i].second);
      if (p.name == "front" || p.name == "back") {
        CHECK(p.kind == PatchKind::Symmetry);
        CHECK(p.source_type == "symmetry_y_strong");
      } else if (p.source_type == "viscous_solid") CHECK(p.kind == PatchKind::Wall);
      else if (p.source_type.starts_with("farfield")) CHECK(p.kind == PatchKind::Farfield);
      else if (p.source_type.starts_with("symmetry")) CHECK(p.kind == PatchKind::Symmetry);
      else if (p.source_type == "subsonic_inflow_pt") CHECK(p.kind == PatchKind::PressureInlet);
      else if (p.source_type == "back_pressure") CHECK(p.kind == PatchKind::PressureOutlet);
      else FAIL("unexpected BC type");
      MESSAGE(std::string(c.file) << ": " << p.name << " = " << p.size);
    }
  }
}

TEST_CASE("NMF explicit maps, discovery, strict ranges and complete coverage") {
  MeshTemp temp;
  const auto grid = temp.path("grid.xyz"), map = temp.path("boundary.map");
  mesh_text(grid, "1\n3 2\n0 1 2 0 1 2\n0 0 0 1 1 1\n");
  const std::string header = "# comment\n1 \\\n1 2 3 2\n";
  const std::string caps = "'symmetry_y_strong' 1 3 1 3 1 2\n'symmetry_y_strong' 1 4 1 3 1 2\n";
  const std::string sides = "WALL 1 1 1 2 1 3\noutflow 1 2 1 2 1 3\nInflow 1 5 1 2 1 2\ncustom_bc 1 6 1 2 1 2 # trailing comment\n";
  std::string text = header + caps + sides;
  bool bad = false;
  SUBCASE("valid explicit map with custom BC") {}
  SUBCASE("invalid block count") { text.replace(text.find("1 \\"), 1, "2"); bad = true; }
  SUBCASE("mismatched dimensions") { text.replace(text.find("1 2 3 2"), 7, "1 2 4 2"); bad = true; }
  SUBCASE("out of range block") { text += "WALL 2 1 1 2 1 3\n"; bad = true; }
  SUBCASE("out of range face") { text += "WALL 1 7 1 2 1 3\n"; bad = true; }
  SUBCASE("out of range index") { text = header + caps + "WALL 1 1 1 2 0 3\n"; bad = true; }
  SUBCASE("integer overflow") { text += "WALL 1 1 1 9999999999999999999 1 3\n"; bad = true; }
  SUBCASE("zero extent") { text = header + caps + "WALL 1 1 1 1 1 3\n"; bad = true; }
  SUBCASE("uncovered face") { text = header + caps; bad = true; }
  SUBCASE("overlapping ranges") { text += "WALL 1 1 1 2 1 3\n"; bad = true; }
  SUBCASE("trailing garbage") { text += "WALL 1 1 1 2 1 3 garbage\n"; bad = true; }
  SUBCASE("unsupported topology") { text = header + "PATCHED 1 1 1 2 1 3\n"; bad = true; }
  SUBCASE("truncated interface") { text = header + "ONE_TO_ONE 1 1 1 2 1 3\n"; bad = true; }
  SUBCASE("invalid swap") { text = header + "ONE_TO_ONE 1 1 1 2 1 3 1 2 1 2 1 3 yes\n"; bad = true; }
  SUBCASE("noncoincident interface") { text = header + "ONE_TO_ONE 1 1 1 2 1 3 1 2 1 2 1 3 FALSE\n"; bad = true; }
  SUBCASE("interface size mismatch") { text = header + "ONE_TO_ONE 1 1 1 2 1 3 1 2 1 2 1 2 FALSE\n"; bad = true; }
  SUBCASE("interface to itself") { text = header + "ONE_TO_ONE 1 1 1 2 1 3 1 1 1 2 1 3 FALSE\n"; bad = true; }
  SUBCASE("non-symmetry cap") { text = header + "WALL 1 3 1 3 1 2\n"; bad = true; }
  mesh_text(map, text);
  if (bad) CHECK_THROWS_AS(read_plot3d(grid, map, .2), Error);
  else {
    auto m = read_plot3d(grid, map, .2);
    boundary_checks(m);
    REQUIRE(m.patches.size() == 6);
    CHECK(m.patches[2].name == "wall_block1_f1_1-2_1-3");
    CHECK(m.patches[2].source_type == "WALL");
    CHECK(m.patches[2].kind == PatchKind::Wall);
    CHECK(m.patches[5].kind == PatchKind::Other);
  }
  // Discovery must neither ignore a malformed sibling nor load an unrelated map.
  CHECK(read_plot3d(grid).patches[0].name == "block1_imin");
  mesh_text(temp.path("grid.nmf"), text);
  if (bad) CHECK_THROWS_AS(read_plot3d(grid), Error);
  else CHECK(read_plot3d(grid).patches[2].kind == PatchKind::Wall);
  CHECK(read_plot3d(grid, std::filesystem::path{}).patches[0].name == "block1_imin");
  CHECK_THROWS_AS(read_plot3d(grid, temp.path("missing.nmf")), Error);
}

TEST_CASE("NMF multi-block extrusion shares caps and stitches reversed ranges") {
  MeshTemp temp;
  auto grid = temp.path("grid.xyz");
  // Second block runs downward in j, so the interface's primary range reverses.
  mesh_text(grid, "2\n3 2 3 2\n0 1 2 0 1 2\n0 0 0 1 1 1\n2 3 4 2 3 4\n1 1 1 0 0 0\n");
  std::string map = "2\n2 2 3 2\n1 2 3 2\n";
  for (int b : {1, 2}) {
    const auto id = std::to_string(b);
    map += "Symmetry-Y " + id + " 3 1 3 1 2\nSymmetry-Y " + id + " 4 1 3 1 2\n";
    map += "WALL " + id + " 1 1 2 1 3\nWALL " + id + " 2 1 2 1 3\n";
  }
  map += "inflow 1 5 1 2 1 2\noutflow 2 6 1 2 1 2\nONE_TO_ONE 1 6 1 2 1 2 2 5 2 1 1 2 FALSE\n";
  mesh_text(temp.path("grid.nmf"), map);
  auto m = read_plot3d(grid, .3);
  boundary_checks(m);
  CHECK(m.n_cells == 4);
  CHECK(m.n_internal_faces() == 3);
  CHECK(m.n_points() == 20);
  REQUIRE(m.patches.size() == 8);
  CHECK(m.patches[0].name == "front");
  CHECK(m.patches[0].size == 4);
  CHECK(m.patches[1].name == "back");
  CHECK(m.patches[1].size == 4);
}

TEST_CASE("NMF native 3D swapped primary directions stitch conforming blocks") {
  MeshTemp temp;
  auto grid = temp.path("grid.xyz");
  std::ostringstream data;
  data << "2\n2 3 4 2 4 3\n";
  for (int block : {0, 1}) {
    const int nj = block ? 4 : 3, nk = block ? 3 : 4;
    for (int dim = 0; dim < 3; ++dim)
      for (int k = 0; k < nk; ++k)
        for (int j = 0; j < nj; ++j)
          for (int i = 0; i < 2; ++i)
            data << (dim == 0 ? i + block : dim == 1 ? (block ? k : j) : (block ? j : k)) << ' ';
    data << '\n';
  }
  mesh_text(grid, data.str());
  std::string map = "2\n1 2 3 4\n2 2 4 3\n";
  for (int block : {1, 2}) {
    const auto id = std::to_string(block);
    const int nj = block == 1 ? 3 : 4, nk = block == 1 ? 4 : 3;
    map += "wall " + id + " 1 1 2 1 " + std::to_string(nj) + "\n";
    map += "wall " + id + " 2 1 2 1 " + std::to_string(nj) + "\n";
    map += "wall " + id + (block == 1 ? " 3" : " 4") + " 1 " + std::to_string(nj) + " 1 " + std::to_string(nk) + "\n";
    for (int face : {5, 6})
      map += "wall " + id + " " + std::to_string(face) + " 1 " + std::to_string(nk) + " 1 2\n";
  }
  map += "'one-to-one' 1 4 1 3 1 4 2 3 1 4 1 3 TRUE\n";
  mesh_text(temp.path("grid.nmf"), map);
  auto m = read_plot3d(grid);
  boundary_checks(m);
  CHECK(m.n_cells == 12);
  CHECK(m.n_internal_faces() == 20);
  CHECK(m.n_points() == 36);
  CHECK(m.patches.size() == 10);
  auto g = compute_geometry(m);
  for (auto v : g.cell_volume) CHECK(v == doctest::Approx(1));
}
