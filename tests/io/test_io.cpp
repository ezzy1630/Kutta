#include <doctest/doctest.h>

#include <bit>
#include <numeric>
#include <regex>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
#include "test_support.hpp"
using namespace kutta;
using namespace kutta::detail;
namespace {
std::filesystem::path fixture(const char* name) {
  return std::filesystem::path(KUTTA_SOURCE_DIR) / "tests/data/meshes" / name;
}
// Load the generated XML tree sufficiently to check nesting, then inspect the
// typed VTK arrays below. No rendering or optional XML/VTK dependency is
// needed.
void check_xml_tree(const std::string& text) {
  std::regex tags(R"(<([^>]+)>)");
  std::vector<std::string> stack;
  for (std::sregex_iterator i(text.begin(), text.end(), tags), end; i != end;
       ++i) {
    auto tag = (*i)[1].str();
    if (tag.starts_with('?')) continue;
    if (tag.starts_with('/')) {
      REQUIRE(!stack.empty());
      CHECK(stack.back() == tag.substr(1));
      stack.pop_back();
    } else if (!tag.ends_with('/'))
      stack.push_back(tag.substr(0, tag.find_first_of(" \t\n")));
  }
  CHECK(stack.empty());
}
std::vector<long long> xml_array(const std::string& xml,
                                 const std::string& name) {
  auto pos = xml.find("Name=\"" + name + "\"");
  CHECK(pos != std::string::npos);
  if (pos == std::string::npos) return {};
  pos = xml.find('>', pos) + 1;
  auto end = xml.find("</DataArray>", pos);
  std::istringstream data(xml.substr(pos, end - pos));
  std::vector<long long> out;
  long long v;
  while (data >> v) out.push_back(v);
  return out;
}
void binary_word(std::string& data, std::uint64_t word, int width,
                 bool little) {
  for (int i = 0; i < width; ++i)
    data +=
        static_cast<char>((word >> (8 * (little ? i : width - 1 - i))) & 255);
}
std::string binary_cube(bool wide, bool little) {
  auto m = make_box_hex(1, 1, 1, {0, 0, 0}, {1, 1, 1});
  auto section = wide ? 3010 : 2010;
  std::string s = "(2 3)\n(4 (" + std::string(little ? "60" : "61") + "))\n(" +
                  std::to_string(section) + " (1 1 8 1 3)(";
  for (auto p : m.points)
    for (auto x : {p.x, p.y, p.z}) {
      if (wide)
        binary_word(s, std::bit_cast<std::uint64_t>(x), 8, little);
      else
        binary_word(s, std::bit_cast<std::uint32_t>(static_cast<float>(x)), 4,
                    little);
    }
  s += ") End of Binary Section " + std::to_string(section) +
       ")\n(2012 (14 1 1 1 0)(";
  binary_word(s, 4, 4, little);
  s += ") End of Binary Section 2012)\n(3013 (a 1 6 3 0)(";
  for (Index f = 0; f < m.n_faces(); ++f) {
    auto ns = nodes(m, f);
    binary_word(s, ns.size(), 4, little);
    for (auto it = ns.rbegin(); it != ns.rend(); ++it)
      binary_word(s, *it + 1, 4, little);
    binary_word(s, 1, 4, little);
    binary_word(s, 0, 4, little);
  }
  s += ") End of Binary Section 3013)\n(45 (10 wall cube 1)())\n(45 (20 fluid "
       "water 1)())";
  return s;
}
}  // namespace
TEST_CASE("Fluent hand-made ASCII 2D quad, 3D hex and mixed cells") {
  auto q = read_mesh(fixture("quad2d.msh"));
  auto g = compute_geometry(q);
  CHECK(q.n_cells == 1);
  CHECK(q.n_faces() == 6);
  CHECK(g.cell_volume[0] == doctest::Approx(1));
  CHECK(q.patches[0].name == "three walls");
  CHECK(q.patches[0].source_zone_id == 10);
  CHECK(q.patches[1].source_zone_id == 11);
  CHECK(q.patches[1].kind == PatchKind::VelocityInlet);
  CHECK(q.patches[1].source_type == "velocity-inlet");
  CHECK(q.patches[2].name == "front");
  CHECK(q.patches[2].kind == PatchKind::Symmetry);
  CHECK(q.patches[3].name == "back");
  CHECK(q.cell_zones[0].name == "water");
  CHECK(q.cell_zones[0].source_zone_id == 20);
  auto h = read_fluent_msh(fixture("hex3d.msh"));
  validate(h);
  CHECK(compute_geometry(h).cell_volume[0] == doctest::Approx(1));
  auto mixed = read_fluent_msh(fixture("mixed2d.msh"), .25);
  g = compute_geometry(mixed);
  CHECK(mixed.n_cells == 3);
  CHECK(mixed.n_internal_faces() == 2);
  CHECK(g.cell_volume[0] == doctest::Approx(.25));
  CHECK(g.cell_volume[1] == doctest::Approx(.125));
  CHECK(g.cell_volume[2] == doctest::Approx(.125));
}
TEST_CASE("legacy binary floats/doubles and both byte orders") {
  MeshTemp temp;
  for (bool wide : {false, true})
    for (bool little : {false, true}) {
      auto path = temp.path("binary.msh");
      mesh_text(path, binary_cube(wide, little));
      auto m = read_mesh(path);
      validate(m);
      CHECK(m.n_cells == 1);
      CHECK(m.n_faces() == 6);
      CHECK(compute_geometry(m).cell_volume[0] == doctest::Approx(1));
    }
  auto truncated = binary_cube(true, true);
  truncated.resize(80);
  auto path = temp.path("short.msh");
  mesh_text(path, truncated);
  CHECK_THROWS_AS(read_mesh(path), Error);
}
TEST_CASE("legacy periodic pairs survive boundary grouping and 2D extrusion") {
  MeshTemp temp;
  auto text = mesh_text(fixture("quad2d.msh"));
  text += "\n(18 (1 1 a b)(1 4))\n";
  auto path = temp.path("periodic.msh");
  mesh_text(path, text);
  auto result = read_fluent_msh_with_metadata(path);
  REQUIRE(result.periodic_faces.size() == 1);
  auto pair = result.periodic_faces[0];
  CHECK(pair.zone == 10);
  CHECK(pair.shadow_zone == 11);
  CHECK(pair.face == result.mesh.patches[0].start);
  CHECK(pair.shadow == result.mesh.patches[1].start);
  text += "(18 (2 2 a b)(1 4))";
  mesh_text(path, text);
  CHECK_THROWS_WITH_AS(read_fluent_msh(path),
                       doctest::Contains("more than once"), Error);
}
TEST_CASE(
    "Fluent malformed references and unsupported adaptations fail clearly") {
  MeshTemp temp;
  auto text = mesh_text(fixture("quad2d.msh"));
  auto path = temp.path("bad.msh");
  SUBCASE("missing node") {
    auto pos = text.find("1 2 1 0");
    text.replace(pos, 7, "1 f 1 0");
    mesh_text(path, text);
    CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains("unknown node"),
                         Error);
  }
  SUBCASE("unclosed cell") {
    auto pos = text.find("1 4 0 1");
    text.replace(pos, 7, "2 4 0 1");
    mesh_text(path, text);
    CHECK_THROWS_AS(read_mesh(path), Error);
  }
  SUBCASE("AMR") {
    mesh_text(path, text + "\n(58 (1 1 1 1)(0))");
    CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains("adapted"), Error);
  }
  SUBCASE("unsupported extension") {
    CHECK_THROWS_AS(read_mesh(temp.path("x.foo")), Error);
  }
}
TEST_CASE("Plot3D formatted 2D multiblock, extrusion and boundary splitting") {
  MeshTemp temp;
  auto path = temp.path("two.p2dfmt");
  mesh_text(
      path,
      "2\n3 2\n2 2\n0 1 2 0 1 2\n0 0 0 1 1 1\n3D0 4D0 3D0 4D0\n0 0 1 1\n");
  auto m = read_mesh(path);
  validate(m);
  CHECK(m.n_cells == 3);
  CHECK(m.n_internal_faces() == 1);
  auto g = compute_geometry(m);
  CHECK(std::accumulate(g.cell_volume.begin(), g.cell_volume.end(), 0.) ==
        doctest::Approx(3));
  REQUIRE(m.patches.size() == 10);
  const std::vector<std::string> names{
      "block1_imin", "block1_imax", "block1_jmin", "block1_jmax", "front", "back",
      "block2_imin", "block2_imax", "block2_jmin", "block2_jmax"};
  for (std::size_t i = 0; i < names.size(); ++i) CHECK(m.patches[i].name == names[i]);
  CHECK(m.patches[4].size == 3);
  CHECK(m.patches[5].size == 3);
  CHECK(m.patches[0].name == "block1_imin");
  CHECK(m.patches[2].name == "block1_jmin");
  auto id = split_patch(m, 2, "plate", [](Vec3d p) { return p.x >= 1; });
  CHECK(id == 3);
  CHECK(m.patches[2].size == 1);
  CHECK(m.patches[3].size == 1);
  validate(m);
  auto before = m.face_nodes;
  CHECK_THROWS_AS(
      split_patch(m, 3, "bad", [](Vec3d) -> bool { throw Error("test"); }),
      Error);
  CHECK(m.face_nodes == before);
  CHECK_THROWS_WITH_AS(read_plot3d(temp.path("grid.p2dfmt.gz")),
                       doctest::Contains("gunzip"), Error);
}
TEST_CASE("Plot3D formatted 3D and single k-plane") {
  MeshTemp temp;
  auto path = temp.path("cube.p3dfmt");
  mesh_text(path,
            "1\n2 2 2\n0 1 0 1 0 1 0 1\n0 0 1 1 0 0 1 1\n0 0 0 0 1 1 1 1\n");
  CHECK(compute_geometry(read_mesh(path)).cell_volume[0] == doctest::Approx(1));
  mesh_text(path, "1\n2 2 1\n0 1 0 1\n0 0 1 1\n2 2 2 2\n");
  auto m = read_plot3d(path, .2);
  CHECK(compute_geometry(m).cell_volume[0] == doctest::Approx(.2));
  CHECK(m.patches[4].name == "front");
}
TEST_CASE(
    "NASA TMR real formatted grids (optional external "
    "fixtures)") {
  const std::filesystem::path root =
      KUTTA_SOURCE_DIR "/validation/data/nasa_tmr";
  if (!std::filesystem::exists(root)) {
    MESSAGE("SKIP: external NASA TMR fixtures unavailable");
    return;
  }
  int checked = 0;
  for (const auto& entry :
       std::filesystem::recursive_directory_iterator(root)) {
    if (entry.path().parent_path().filename() != "grids" ||
        entry.path().extension() != ".p2dfmt")
      continue;
    // Exercise real coarsest grids without copying or decompressing large
    // inputs.
    if (entry.path().filename().string().find("4levelsdown") ==
        std::string::npos)
      continue;
    INFO(entry.path().string());
    auto m = read_plot3d(entry.path(), .01);
    validate(m);
    CHECK(m.n_cells > 0);
    ++checked;
  }
  CHECK(checked >= 2);
  MESSAGE("NASA TMR real grids checked: " << checked);
}
TEST_CASE(
    "VTU ASCII structure, hex/tet fast paths, polyhedron stream and fields") {
  MeshTemp temp;
  for (bool tet : {false, true}) {
    auto m = tet ? make_box_tet(1, 1, 1, {0, 0, 0}, {1, 1, 1})
                 : make_box_hex(1, 1, 1, {0, 0, 0}, {1, 1, 1});
    std::vector<double> p(m.n_cells, 2.5);
    CellField field{"p & \"gauge\"", 1, p};
    auto path = temp.path("box.vtu");
    std::vector<double> velocity(static_cast<std::size_t>(m.n_cells) * 3, 1.5);
    std::array<CellField, 2> fields{field, CellField{"velocity", 3, velocity}};
    write_vtu(path, m, fields);
    auto text = mesh_text(path);
    check_xml_tree(text);
    CHECK(text.find("Name=\"velocity\" NumberOfComponents=\"3\"") !=
          std::string::npos);
    CHECK(text.find("<UnstructuredGrid>") != std::string::npos);
    CHECK(text.find("p &amp; &quot;gauge&quot;") != std::string::npos);
    auto types = xml_array(text, "types"),
         conn = xml_array(text, "connectivity"),
         offsets = xml_array(text, "offsets");
    REQUIRE(types.size() == static_cast<std::size_t>(m.n_cells));
    for (auto t : types) CHECK(t == (tet ? 10 : 12));
    CHECK(offsets.back() == static_cast<long long>(conn.size()));
    // The emitted VTK node order must have positive orientation.
    for (Index c = 0; c < m.n_cells; ++c) {
      auto off = c ? offsets[c - 1] : 0;
      auto a = m.points[conn[off]], b = m.points[conn[off + 1]],
           d = m.points[conn[off + (tet ? 2 : 3)]],
           e = m.points[conn[off + (tet ? 3 : 4)]];
      CHECK(dot(b - a, cross(d - a, e - a)) > 0);
    }
    CHECK(text.find("Name=\"faces\"") == std::string::npos);
  }
  PlanarMesh p;
  p.n_cells = 1;
  p.patches = {{"wall", PatchKind::Wall}};
  for (Index i = 0; i < 5; ++i) {
    double a = 2 * std::acos(-1.) * i / 5;
    p.points.push_back({std::cos(a), std::sin(a), 0});
    p.edges.push_back({i, (i + 1) % 5});
    p.owner.push_back(0);
    p.neighbour.push_back(-1);
    p.edge_patch.push_back(0);
  }
  auto m = extrude_2d(p);
  auto path = temp.path("poly.vtu");
  write_vtu(path, m);
  auto text = mesh_text(path);
  check_xml_tree(text);
  auto faces = xml_array(text, "faces"),
       offsets = xml_array(text, "faceoffsets");
  REQUIRE(offsets.size() == 1);
  CHECK(offsets[0] == static_cast<long long>(faces.size()));
  REQUIRE(faces.size() > 1);
  CHECK(faces[0] == 7);
  CHECK(xml_array(text, "types")[0] == 42);
  std::size_t i = 1;
  for (int f = 0; f < 7; ++f) {
    auto count = faces[i++];
    CHECK((count == 4 || count == 5));
    i += static_cast<std::size_t>(count);
  }
  CHECK(i == faces.size());
  std::vector<double> bad;
  CellField field{"bad", 1, bad};
  CHECK_THROWS_AS(write_vtu(path, m, std::span(&field, 1)), Error);
}

TEST_CASE("legacy 3D mixed and polygonal cells honor file c0/c1 winding") {
  auto m = read_fluent_msh(fixture("mixed3d.msh"));
  validate(m);
  auto g = compute_geometry(m);
  REQUIRE(m.n_cells == 2);
  REQUIRE(m.n_internal_faces() == 1);
  CHECK(g.cell_volume[0] == doctest::Approx(1));
  CHECK(g.cell_volume[1] == doctest::Approx(1. / 3));
  CHECK(g.face_area[0].z > 0);
  CHECK(g.cell_centroid[1].z == doctest::Approx(1.25));
  CHECK(m.patches[1].kind == PatchKind::Symmetry);
  auto poly = read_fluent_msh(fixture("poly3d.msh"));
  validate(poly);
  g = compute_geometry(poly);
  CHECK(poly.n_faces() == 7);
  CHECK(g.cell_volume[0] == doctest::Approx(2.5));
  CHECK(g.cell_centroid[0].z == doctest::Approx(.5));
}
TEST_CASE("binary periodic records and 2D nodes") {
  MeshTemp temp;
  std::string text = "(2 2)(4 (60))(2010 (1 1 4 1 2)(";
  for (auto point :
       std::array<std::array<float, 2>, 4>{{{0, 0}, {1, 0}, {1, 1}, {0, 1}}})
    for (auto x : point)
      binary_word(text, std::bit_cast<std::uint32_t>(x), 4, true);
  text += ") End of Binary Section 2010)(12 (14 1 1 1 3))(2013 (a 1 3 c 2)(";
  for (auto x : {1, 2, 1, 0, 2, 3, 1, 0, 3, 4, 1, 0})
    binary_word(text, x, 4, true);
  text += ") End of Binary Section 2013)(3013 (b 4 4 8 2)(";
  for (auto x : {4, 1, 1, 0}) binary_word(text, x, 4, true);
  text += ") End of Binary Section 3013)(2018 (1 1 a b)(";
  binary_word(text, 1, 4, true);
  binary_word(text, 4, 4, true);
  text += ") End of Binary Section 2018)";
  auto path = temp.path("periodic2d.msh");
  mesh_text(path, text);
  auto result = read_fluent_msh_with_metadata(path, .25);
  REQUIRE(result.periodic_faces.size() == 1);
  CHECK(result.mesh.patches[0].source_type == "periodic");
  CHECK(result.mesh.patches[1].source_type == "periodic-shadow");
  CHECK(compute_geometry(result.mesh).cell_volume[0] == doctest::Approx(.25));
}

TEST_CASE("Fluent optional dimension section is inferred from node ND") {
  MeshTemp temp;
  auto text = mesh_text(fixture("quad2d.msh"));
  auto pos = text.find("(2 2)");
  text.erase(pos, 5);
  auto path = temp.path("no-dimension.msh");
  mesh_text(path, text);
  auto m = read_fluent_msh(path);
  CHECK(compute_geometry(m).cell_volume[0] == doctest::Approx(1));
}

TEST_CASE("Fluent face sections may precede later coordinate chunks") {
  MeshTemp temp;
  auto text = mesh_text(fixture("hex3d.msh"));
  auto start = text.find("(10 (1 1 8");
  auto end = text.find("))", start) + 2;
  text.replace(start, end - start, "(10 (1 1 3 1 3)(0 0 0 1 0 0 1 1 0))");
  text += "\n(10 (1 4 8 1 3)(0 1 0 0 0 1 1 0 1 1 1 1 0 1 1))";
  auto path = temp.path("chunks.msh");
  mesh_text(path, text);
  auto m = read_fluent_msh(path);
  CHECK(compute_geometry(m).cell_volume[0] == doctest::Approx(1));
}

TEST_CASE("VTU mixed fast cells and polyhedra use correct face offsets") {
  std::vector<Vec3d> points;
  std::vector<kutta::detail::CellFaces> cells;
  for (const auto& part : {make_box_hex(1, 1, 1, {0, 0, 0}, {1, 1, 1}),
                           read_fluent_msh(fixture("poly3d.msh")),
                           make_box_tet(1, 1, 1, {3, 0, 0}, {4, 1, 1})}) {
    auto offset = detail::index(points.size());
    points.insert(points.end(), part.points.begin(), part.points.end());
    auto adj = build_cell_faces(part);
    for (Index c = 0; c < part.n_cells; ++c) {
      kutta::detail::CellFaces faces;
      for (auto i = adj.offsets[c]; i < adj.offsets[c + 1]; ++i) {
        auto ns = nodes(part, adj.faces[i]);
        std::vector<Index> face(ns.begin(), ns.end());
        for (auto& p : face) p += offset;
        faces.push_back(std::move(face));
      }
      cells.push_back(std::move(faces));
    }
  }
  auto mixed =
      from_cells(std::move(points), cells, {{"walls", PatchKind::Wall}},
                 [](const auto&) { return 0; });
  MeshTemp temp;
  auto path = temp.path("mixed.vtu");
  write_vtu(path, mixed);
  auto text = mesh_text(path);
  check_xml_tree(text);
  auto types = xml_array(text, "types"),
       faceoffsets = xml_array(text, "faceoffsets"),
       faces = xml_array(text, "faces");
  REQUIRE(types.size() == 8);
  REQUIRE(faceoffsets.size() == 8);
  CHECK(types[0] == 12);
  CHECK(types[1] == 42);
  CHECK(faceoffsets[0] == -1);
  CHECK(faceoffsets[1] == static_cast<long long>(faces.size()));
  for (std::size_t i = 2; i < 8; ++i) {
    CHECK(types[i] == 10);
    CHECK(faceoffsets[i] == -1);
  }
}

TEST_CASE("legacy uniform binary cells have a trailer without a payload") {
  MeshTemp temp;
  for (int section : {2012, 3012}) {
    auto text = mesh_text(fixture("hex3d.msh"));
    const std::string cells = "(12 (14 1 1 1 4))";
    auto start = text.find(cells);
    REQUIRE(start != std::string::npos);
    text.replace(start, cells.size(),
                 "(" + std::to_string(section) +
                     " (14 1 1 1 4)\nEnd of Binary Section " +
                     std::to_string(section) + ")");
    auto path = temp.path("uniform.cas");
    mesh_text(path, text);
    auto m = read_mesh(path);
    CHECK(compute_geometry(m).cell_volume[0] == doctest::Approx(1));
  }
}

TEST_CASE("legacy settings skip balanced lists and escaped quoted parentheses") {
  MeshTemp temp;
  auto text = mesh_text(fixture("hex3d.msh"));
  text += R"settings(
(37 ((nested (list () ())) (label "(unbalanced \" quote") (more 1)))
(38 ((label "closing ) and opening ( parentheses") (nested (() ()))))
(45 (99 wall "unused (zone) name")())
)settings";
  auto path = temp.path("settings.cas");
  mesh_text(path, text);
  auto m = read_mesh(path);
  CHECK(compute_geometry(m).cell_volume[0] == doctest::Approx(1));
  mesh_text(path, text + "(37 ((truncated \"string)");
  CHECK_THROWS_AS(read_mesh(path), Error);
  // Never scan unknown binary payloads looking for parentheses/terminators.
  mesh_text(path, text + "(2999 (1 2)(arbitrary bytes))");
  CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains("unsupported binary"),
                       Error);
}

TEST_CASE("binary node flags skip exact counts even with parentheses in raw bytes") {
  MeshTemp temp;
  for (int section : {2041, 3041}) {
    auto text = mesh_text(fixture("hex3d.msh"));
    text += "(" + std::to_string(section) + " (1 2)(";
    binary_word(text, 0x29222828, 4, true);
    binary_word(text, 0xffffffff, 4, true);
    auto truncated = text.substr(0, text.size() - 1);
    text += ") End of Binary Section " + std::to_string(section) + ")";
    auto path = temp.path("flags.cas");
    mesh_text(path, text);
    CHECK(compute_geometry(read_mesh(path)).cell_volume[0] == doctest::Approx(1));
    mesh_text(path, truncated);
    CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains("truncated binary"),
                         Error);
  }
}

TEST_CASE("Plot3D multi-block 2D dimensions on one line without an extension hint") {
  MeshTemp temp;
  auto path = temp.path("two.xyz");
  mesh_text(path, "2\n2 2 2 2\n0 1 0 1\n0 0 1 1\n2 3 2 3\n0 0 1 1\n");
  auto m = read_plot3d(path);
  CHECK(m.n_cells == 2);
  CHECK(m.patches.size() == 10);
}
