// Legacy (non-HDF5) Fluent .dat files: synthetic format rules, then the
// corpus pairs (cylinder, Borealis, VKI) when KUTTA_FLUENT_EXAMPLES has them.
#include "../../bench/fluent_restore_cr.hpp"
#include "../fluent_examples.hpp"
#include "kutta/case/fluent_data.hpp"
#include "kutta/core/sha256.hpp"
#include "kutta/io/mesh_io.hpp"
#include <bit>
#include <chrono>
#include <cmath>
#include <doctest/doctest.h>
#include <fstream>
#include <iterator>
#include <random>

using namespace kutta::fluent;
namespace {
std::filesystem::path temporary(const std::string& stem) {
  return std::filesystem::temp_directory_path() /
         (stem + "-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + "-" +
          std::to_string(std::random_device{}()) + ".dat");
}
struct File {
  std::filesystem::path path;
  explicit File(const std::string& bytes, const std::string& stem = "kutta-legacy-data") : path(temporary(stem)) {
    std::ofstream(path, std::ios::binary).write(bytes.data(), std::streamsize(bytes.size()));
  }
  ~File() {std::error_code ec;std::filesystem::remove(path, ec);}
};
std::string doubles(std::initializer_list<double> values) {
  std::string out;
  for (double v : values) {
    const auto bits = std::bit_cast<std::uint64_t>(v);
    for (int i = 0; i < 8; ++i) out.push_back(char((bits >> (8 * i)) & 0xff));
  }
  return out;
}
// Two fluid cells (zone 2, IDs 1-2), a wall (zone 3, faces 1-2) and an
// interior face (zone 4, face 3); residual names as Fluent's rpvar.
FluentCase topology() {
  FluentCase c;
  c.mesh.dimension = 2;c.mesh.cells = 2;c.mesh.faces = 3;c.mesh.nodes = 4;
  c.mesh.cell_zones = {{2, "fluid", 1, 2, 1}};
  c.mesh.face_zones = {{3, "wall", 1, 2, 3}, {4, "interior", 3, 3, 2}};
  CellZone fluid;fluid.id = 2;fluid.name = "fluid";fluid.type = ZoneType::fluid;c.cell_zones.push_back(fluid);
  c.variables.alist = parse_sexpr("((residuals/settings ((continuity #t 0.) (x-velocity #t 0.))))");
  return c;
}
std::string header() {
  return "(0 \"fluent18.0.0  build-id: 311\")\n(4 (60 0 0 1 2 4 4 4 8 8 8 4))\n(33 (2 3 4))\n"
         "(37 (\n(flow-time 0)\n(number-of-iterations 7)))\n";
}
std::string sample() {
  return header() +
         "(0 \"SV_P, domain 1, cell zone 2, 2 cells:\")\n(300 (1 2 1 0 1 1 2)\n(\n1.5\n-2.25e1\n))\n" +
         "(0 \"SV_U, domain 1, cell zone 2, 2 cells:\")\n(3300 (111 2 1 0 1 1 2)\n(" + doubles({3, 4}) +
         "\nEnd of Binary Section   3300))\n" +
         "(0 \"SV_WALL_SHEAR, domain 1, wall zone 3, 2 faces:\")\n(3300 (19 3 3 0 1 1 2)\n(" +
         doubles({1, 2, 0, 5, 6, 0}) + "\nEnd of Binary Section   3300))\n" +
         "(0 \"Integer data:\")\n(3314 (1 3 1 0 1 1 2)\n(" + doubles({7, 8}) + "\nEnd of Binary Section   3314))\n" +
         "(0 \"Residuals:\")\n(3302 (2 1 1 1)\n(" + doubles({1, 10, 10, 2, 5, 10}) +
         ")\n\nEnd of Binary Section   3302)\n(3302 (2 14 1 1)\n(" + doubles({1, 3, 6, 2, 1, 6}) +
         ")\n\nEnd of Binary Section   3302)\n";
}
}  // namespace

TEST_CASE("Legacy data reads ASCII and binary fields, residual rows and the grid size") {
  const File file(sample());
  const auto d = read_data(file.path, topology());
  CHECK(d.cell_fields.at("SV_P").at(1)[0] == 1.5);
  CHECK(d.cell_fields.at("SV_P").at(2)[0] == -22.5);
  CHECK(d.cell_fields.at("SV_U").at(2)[0] == 4);
  const auto shear = d.face_fields.at("SV_WALL_SHEAR").at(2);
  REQUIRE(shear.size() == 3);
  CHECK(shear[0] == 5);CHECK(shear[1] == 6);CHECK(shear[2] == 0);
  CHECK(d.iteration_count == 7);
  CHECK(d.last_residual_iteration == 2);
  // Equations follow the case's residuals/settings order.
  const auto& continuity = d.residuals.at("continuity");
  REQUIRE(continuity.rows.size() == 2);
  CHECK(continuity.iterations == std::vector<std::uint64_t>{1, 2});
  CHECK(continuity.rows[1][0] == 5);CHECK(continuity.rows[1][1] == 10);
  CHECK(std::isnan(continuity.rows[1][2]));
  CHECK(d.residuals.at("x-velocity").rows[0][0] == 3);
  CHECK(std::any_of(d.issues.begin(), d.issues.end(), [](auto& i) {return i.path == "data/section/3314";}));
}

TEST_CASE("Legacy data rejects mismatched grids, unknown zones, short bodies and renamed IDs") {
  const auto c = topology();
  auto grid = sample();grid.replace(grid.find("(33 (2 3 4))"), 12, "(33 (3 3 4))");
  CHECK_THROWS_WITH(read_data(File(grid).path, c), doctest::Contains("grid size disagrees"));
  auto zone = sample();zone.replace(zone.find("(3300 (111 2 "), 13, "(3300 (111 9 ");
  CHECK_THROWS_WITH(read_data(File(zone).path, c), doctest::Contains("unknown zone 9"));
  auto range = sample();range.replace(range.find("(300 (1 2 1 0 1 1 2)"), 20, "(300 (1 2 1 0 1 2 3)");
  CHECK_THROWS_WITH(read_data(File(range).path, c), doctest::Contains("outside matching case topology"));
  // A text-mode transfer deleting a CR before LF shortens a binary body.
  auto stripped = sample();
  const auto body = stripped.find("(3300 (111 2 1 0 1 1 2)\n(") + 25;
  stripped.replace(body, 16, doubles({3, 4}).substr(0, 15));
  stripped += std::string("\r");  // a lone CR, but no CR LF pair
  CHECK_THROWS_WITH(read_data(File(stripped).path, c), doctest::Contains("text-mode transfer"));
  // The same loss in a residual body, whose last read byte is then its ')'.
  stripped = sample();
  const auto rows = stripped.find("(3302 (2 1 1 1)\n(") + 17;
  stripped.erase(rows, 1);stripped += std::string("\r");
  CHECK_THROWS_WITH(read_data(File(stripped).path, c), doctest::Contains("text-mode transfer"));
  // One numeric ID under two names means the comments are not the fields'.
  auto renamed = sample();renamed.replace(renamed.find("(3300 (19 3 "), 12, "(3300 (111 3 ");
  CHECK_THROWS_WITH(read_data(File(renamed).path, c), doctest::Contains("names both"));
  // A cell field must cover every fluid cell.
  const std::string full = "(300 (1 2 1 0 1 1 2)\n(\n1.5\n-2.25e1\n))";
  auto partial = sample();partial.replace(partial.find(full), full.size(), "(300 (1 2 1 0 1 1 1)\n(\n1.5\n))");
  CHECK_THROWS_WITH(read_data(File(partial).path, c), doctest::Contains("does not cover"));
}

TEST_CASE("Legacy cylinder data: upstream CR loss is detected, and the checked restoration reads") {
  const auto cas = fluent_examples() / "adjoint_cylinder_tutorial/cylinder_tutorial.cas";
  const auto dat = fluent_examples() / "adjoint_cylinder_tutorial/cylinder_tutorial.dat";
  if (!std::filesystem::exists(cas) || !std::filesystem::exists(dat)) {
    MESSAGE("SKIP: external Fluent cylinder pair absent (set KUTTA_FLUENT_EXAMPLES)");
    return;
  }
  const auto c = read_case(cas);
  CHECK(c.dimension == 2);
  CHECK_THROWS_WITH(read_data(dat, c), doctest::Contains("text-mode transfer"));
  std::ifstream in(dat, std::ios::binary);
  const std::string bytes{std::istreambuf_iterator<char>(in), {}};
  const auto r = kutta::forensics::restore_stripped_cr(bytes);
  // validation/results/cylinder_2d/REPORT.md lists the 15 restored bytes.
  CHECK(r.inserted.size() == 15);
  CHECK(r.bytes.size() == bytes.size() + 15);
  const File restored(r.bytes, "kutta-cylinder-restored");
  CHECK(kutta::sha256_file(restored.path) == "cb7dce25827ca4aea45f6fc4bd9936d4ff55fa84ee14d361cec65c7d90ec57cb");
  const auto d = read_data(restored.path, c);
  // Independently decoded from the restored bytes with Python struct.
  CHECK(d.cell_fields.at("SV_P").at(1)[0] == -393.29536858061044);
  CHECK(d.cell_fields.at("SV_P").at(25740)[0] == -2.456773394238133);
  CHECK(d.cell_fields.at("SV_V").at(4044)[0] == 13.970779809010404);
  CHECK(d.face_fields.at("SV_WALL_SHEAR").at(1)[1] == 0.17156632703069016);
  CHECK_FALSE(d.cell_fields.contains("SV_W"));
  // Fluent stopped at its 1e-5 continuity criterion.
  const auto& continuity = d.residuals.at("continuity");
  CHECK(continuity.iterations.back() == 66);
  CHECK(continuity.rows.back()[0] / continuity.rows.back()[1] < 1e-5);
}

TEST_CASE("Legacy Borealis and VKI data read against their cases") {
  const auto root = fluent_examples();
  for (const auto* stem : {"borealis_fracture/26 Mpa-solution-final", "borealis_fracture/final solution-before-26 MPa",
                           "borealis_fracture/final solution-before-94 MPa", "vki_turbine/vki_turbine"}) {
    INFO(stem);
    const auto cas = root / (std::string(stem) + ".cas"), dat = root / (std::string(stem) + ".dat");
    if (!std::filesystem::exists(cas) || !std::filesystem::exists(dat)) {MESSAGE("SKIP: absent ", stem);continue;}
    const auto c = read_case(cas);
    const auto d = read_data(dat, c);
    CHECK(d.cell_fields.at("SV_P").size() == *c.mesh.cells);
    CHECK(d.residuals.contains("continuity"));
  }
  const auto borealis = root / "borealis_fracture/26 Mpa-solution-final";
  if (std::filesystem::exists(borealis.string() + ".dat")) {
    const auto c = read_case(borealis.string() + ".cas");
    const auto d = read_data(borealis.string() + ".dat", c);
    CHECK(d.cell_fields.at("SV_P").at(1)[0] == 8330.398965126962);
    CHECK(d.cell_fields.at("SV_P").at(511466)[0] == 14688.538710315972);
    // Section 2040 (cell partitions) no longer blocks the mesh.
    CHECK(kutta::read_mesh(borealis.string() + ".cas").n_cells == 511466);
  }
  const auto vki = root / "vki_turbine/vki_turbine.cas";
  if (std::filesystem::exists(vki))
    CHECK_THROWS_WITH(kutta::read_mesh(vki), doctest::Contains("nonconformal"));
}
