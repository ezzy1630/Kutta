#include <doctest/doctest.h>
#include <algorithm>
#include <iterator>

#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <set>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
using namespace kutta;
using namespace kutta::detail;

namespace {
struct ExpectedZone {
  int id;
  const char* name;
  Index size;
};
struct Example {
  const char* file;
  // Source counts from h5dump /meshes/1 attributes and zoneTopology ranges,
  // or legacy declaration headers 10/12/13. 2D counts precede extrusion.
  Index cells, faces, points, internal;
  bool planar;
  ExpectedZone cell_zone;
  std::vector<ExpectedZone> patches;
};
const std::vector<ExpectedZone> mixer_patches{
    {20, "wall", 3873}, {19, "outlet", 104}, {18, "inlet2", 103},
    {17, "inlet1", 101}};
const std::vector<Example> examples{
    {"flow_over_cylinder/flow_over_cylinder.msh", 46949, 94494, 47545, 93281,
     true, {2, "surface_body", 46949},
     {{5, "inlet", 150}, {6, "outlet", 150}, {7, "walls", 600},
      {8, "cylinder", 313}}},
    {"lid_driven_cavity/lid_driven_cavity.msh", 3550, 7221, 3672, 6979,
     true, {2, "surface_body", 3550},
     {{5, "top_moving_wall", 71}, {6, "bottom_wall", 71},
      {7, "stationary_side_wall", 100}}},
    {"simple_duct/simple-duct.msh.h5", 10, 51, 44, 9, false, {2, "fluid", 10},
     {{5, "symm", 40}, {4, "inlet", 1}, {3, "outflow", 1}}},
    {"box/poly.msh.h5", 86803, 591932, 505132, 575870, false,
     {2, "fluid", 86803},
     {{8, "symmetry", 2907}, {9, "top", 2587}, {7, "right", 2545},
      {6, "left", 2532}, {4, "front", 2907}, {3, "bottom", 2584}}},
    {"aero/wing.msh.h5", 114219, 654343, 478848, 641410, false,
     {240, "solid", 114219},
     {{86, "symmetry", 1895}, {90, "farfield", 877}, {89, "wing_te", 415},
      {87, "wing_top", 4886}, {88, "wing_bot", 4860}}},
    {"file_session/elbow1.cas.h5", 13852, 44500, 16968, 38612, false,
     {2, "fluid", 13852},
     {{3, "wall", 3630}, {4, "symmetry", 2018}, {5, "pressure-outlet-7", 100},
      {6, "velocity-inlet-6", 40}, {7, "velocity-inlet-5", 100}}},
    {"mixing_elbow/mixing_elbow.cas.h5", 17822, 91581, 66417, 86782, false,
     {87, "elbow-fluid", 17822},
     {{34, "wall-elbow", 2168}, {33, "wall-inlet", 268}, {32, "outlet", 155},
      {31, "cold-inlet", 152}, {30, "hot-inlet", 55},
      {29, "symmetry-xyplane", 2001}}},
    {"static_mixer/Static_Mixer_Parameters.cas.h5", 22771, 119140, 82247, 114959,
     false, {97, "fluid", 22771}, mixer_patches},
    {"static_mixer/Static_Mixer_Parameters_legacy_text.cas", 22771, 119140,
     82247, 114959, false, {97, "fluid", 22771}, mixer_patches},
    {"static_mixer/Static_Mixer_Parameters_legacy_binary.cas", 22771, 119140,
     82247, 114959, false, {97, "fluid", 22771}, mixer_patches},
};
const std::set<std::string> data_files{
    "file_session/elbow1.dat.h5", "mixing_elbow/mixing_elbow.dat.h5"};

void same_mesh(const PolyMesh& a, const PolyMesh& b) {
  CHECK(a.n_cells == b.n_cells);
  CHECK(a.face_offsets == b.face_offsets);
  CHECK(a.face_nodes == b.face_nodes);
  CHECK(a.owner == b.owner);
  CHECK(a.neighbour == b.neighbour);
  REQUIRE(a.points.size() == b.points.size());
  double scale = 0, displacement = 0;
  for (std::size_t i = 0; i < a.points.size(); ++i) {
    scale = std::max(scale, norm(a.points[i]));
    displacement = std::max(displacement, norm(a.points[i] - b.points[i]));
  }
  CHECK(displacement <= 1e-12 * scale);
  REQUIRE(a.patches.size() == b.patches.size());
  for (std::size_t i = 0; i < a.patches.size(); ++i) {
    CHECK(a.patches[i].name == b.patches[i].name);
    CHECK(a.patches[i].source_zone_id == b.patches[i].source_zone_id);
    CHECK(a.patches[i].kind == b.patches[i].kind);
    CHECK(a.patches[i].start == b.patches[i].start);
    CHECK(a.patches[i].size == b.patches[i].size);
  }
}
}  // namespace

TEST_CASE("real Fluent external mesh acceptance") {
  const std::filesystem::path root =
      KUTTA_SOURCE_DIR "/validation/fluent_examples";
  if (!std::filesystem::exists(root)) {
    MESSAGE("SKIP: external Fluent examples directory is absent");
    return;
  }
  // Inventory prevents an added fixture from silently escaping the concrete
  // acceptance expectations. Only two fixtures are tracked; the rest are
  // gitignored downloads (SOURCE.md), so absent ones are reported and skipped.
  // Data-only files cannot supply a mesh.
  std::set<std::string> expected = data_files, found;
  for (const auto& e : examples) expected.insert(e.file);
  for (const auto& entry : std::filesystem::recursive_directory_iterator(root))
    if (entry.is_regular_file() && entry.path().extension() != ".md" &&
        !entry.path().filename().string().starts_with("LICENSE"))
      found.insert(entry.path().lexically_relative(root).generic_string());
  std::vector<std::string> unexpected, absent;
  std::set_difference(found.begin(), found.end(), expected.begin(), expected.end(),
                      std::back_inserter(unexpected));
  std::set_difference(expected.begin(), expected.end(), found.begin(), found.end(),
                      std::back_inserter(absent));
  for (const auto& file : unexpected) FAIL_CHECK("fixture without expectations: " << file);
  for (const auto& file : absent) MESSAGE("SKIP: fixture absent (see SOURCE.md): " << file);

  std::optional<PolyMesh> mixer;
  double mixer_volume = 0;
  int checked = 0, readable = 0;
  for (const auto& e : examples) {
    const auto path = root / e.file;
    if (!found.contains(e.file)) continue;
#if !KUTTA_HAS_HDF5
    if (path.extension() == ".h5") continue;
#endif
    ++readable;
    INFO(e.file);
    try {
      auto m = read_mesh(path);
      validate(m);
      REQUIRE(m.n_cells == e.cells);
      CHECK(m.n_faces() == e.faces + (e.planar ? 2 * e.cells : 0));
      CHECK(m.n_points() == e.points * (e.planar ? 2 : 1));
      CHECK(m.n_internal_faces() == e.internal);
      auto patches = e.patches;
      if (e.planar) {
        patches.push_back({-1, "front", e.cells});
        patches.push_back({-1, "back", e.cells});
      }
      REQUIRE(m.patches.size() == patches.size());
      for (std::size_t i = 0; i < patches.size(); ++i) {
        CHECK(m.patches[i].source_zone_id == patches[i].id);
        CHECK(m.patches[i].name == patches[i].name);
        CHECK(m.patches[i].size == patches[i].size);
      }
      REQUIRE(m.cell_zones.size() == 1);
      CHECK(m.cell_zones[0].source_zone_id == e.cell_zone.id);
      CHECK(m.cell_zones[0].name == e.cell_zone.name);
      CHECK(m.cell_zones[0].cells.size() == static_cast<std::size_t>(e.cell_zone.size));
      const auto g = compute_geometry(m);
      auto volume = std::accumulate(g.cell_volume.begin(), g.cell_volume.end(), 0.);
      CHECK(std::isfinite(volume));
      CHECK(volume > 0);
      CHECK(*std::min_element(g.cell_volume.begin(), g.cell_volume.end()) > 0);
      Vec3d closure;
      double area = 0;
      Index wrong_orientation = 0;
      double minimum_orientation = std::numeric_limits<double>::infinity();
      for (Index f = 0; f < m.n_faces(); ++f) {
        auto delta = g.face_centroid[f] - g.cell_centroid[m.owner[f]];
        auto projection = dot(g.face_area[f], delta);
        if (!(projection > 0)) ++wrong_orientation;
        minimum_orientation = std::min(minimum_orientation,
                                      projection / (norm(g.face_area[f]) * norm(delta)));
        if (f >= m.n_internal_faces()) {
          closure = closure + g.face_area[f];
          area += norm(g.face_area[f]);
        }
      }
      // Every face points away from its owner's center, including boundaries.
      CHECK(wrong_orientation == 0);
      CHECK(minimum_orientation > 0);
      CHECK(norm(closure) <= 1e-10 * area);
      std::cout << std::setprecision(15) << "REAL " << e.file
                << " cells=" << m.n_cells << " faces=" << m.n_faces()
                << " points=" << m.n_points() << " patches=" << m.patches.size()
                << " volume=" << volume << " min-normal-cos=" << minimum_orientation
                << " boundary-closure=" << norm(closure) / area << '\n';
      if (path.parent_path().filename() == "static_mixer") {
        if (mixer) {
          same_mesh(*mixer, m);
          CHECK(std::abs(volume - mixer_volume) <= 1e-12 * mixer_volume);
        } else {
          mixer = std::move(m);
          mixer_volume = volume;
        }
      }
      ++checked;
    } catch (const Error& error) {
      FAIL_CHECK(std::string(error.what()));
    }
  }
  CHECK(checked == readable);
#if KUTTA_HAS_HDF5
  for (const auto& file : data_files) {
    if (!found.contains(file)) continue;
    INFO(file);
    CHECK_THROWS_WITH_AS(read_mesh(root / file),
                         doctest::Contains("Unsupported mesh extension"), Error);
    CHECK_THROWS_WITH_AS(read_fluent_cff(root / file),
                         doctest::Contains("no /meshes/1 mesh"), Error);
  }
#else
  MESSAGE("SKIP: CFF mesh cases and data-only files require HDF5");
#endif
}
