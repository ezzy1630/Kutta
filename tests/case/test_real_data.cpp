#include "kutta/case/fluent_data.hpp"
#include <cstdlib>
#include <algorithm>
#include <doctest/doctest.h>
using namespace kutta::fluent;
#if KUTTA_HAS_HDF5
namespace {
std::filesystem::path examples() {
  auto env = std::getenv("KUTTA_FLUENT_EXAMPLES");
  return env ? env : KUTTA_SOURCE_DIR "/validation/fluent_examples";
}
bool available(const std::filesystem::path &path) {
  if (std::filesystem::exists(path))
    return true;
  MESSAGE("Fluent fixture absent, optional real-file test skipped: ", path.string());
  return false;
}
} // namespace
TEST_CASE("Real mixing elbow data fields match mesh zones") {
  auto path = examples() / "mixing_elbow/mixing_elbow.cas.h5";
  if (!available(path))
    return;
  auto c = read_case(path);
  auto data = examples() / "mixing_elbow/mixing_elbow.dat.h5";
  if (!available(data))
    return;
  auto d = read_data(data, c);
  CHECK(d.case_file == "mixing_elbow.cas.h5");
  CHECK(c.source_path->filename() == "mixing_elbow.cas.h5");
  CHECK(std::none_of(d.issues.begin(), d.issues.end(), [](auto &i) { return i.path == "data/Case File"; }));
  CHECK(d.cell_fields.at("SV_P").size() == 17822);
  CHECK(d.cell_fields.at("SV_O").size() == 17822);
  CHECK(d.face_fields.at("SV_P").size() == 91581);
  CHECK(d.cell_fields.contains("SV_U"));
  CHECK(d.cell_fields.contains("SV_V"));
  CHECK(d.cell_fields.contains("SV_W"));
  CHECK(d.cell_fields.contains("SV_K"));
  CHECK(d.cell_fields.contains("SV_MU_T"));
  CHECK(d.cell_fields.contains("SV_DENSITY"));
  CHECK(d.face_fields.at("SV_WALL_SHEAR").size() == 2436);
  CHECK(d.face_fields.at("SV_WALL_SHEAR").sections[0].components == 3);
  // Independently read with h5dump -m '%.17g' -s 0,0 -c 2,3.
  auto shear = d.face_fields.at("SV_WALL_SHEAR").at(1);
  REQUIRE(shear.size() == 3);
  CHECK(shear[0] == -2.1333901821377996e-08);
  CHECK(shear[1] == -2.4499093077299929e-08);
  CHECK(shear[2] == -3.5408833061117729e-09);
  CHECK(d.face_fields.at("SV_WALL_SHEAR").at(2)[1] == -3.0258976277082841e-07);

  CHECK(d.zone_field(FieldLocation::face, 34, "SV_WALL_SHEAR")[0].values.size() == 2168 * 3);
  CHECK(d.zone_field(FieldLocation::face, 33, "SV_WALL_SHEAR")[0].values.size() == 268 * 3);
  CHECK(d.zone_field(FieldLocation::face, 32, "SV_WALL_SHEAR").empty());
  CHECK(d.last_residual_iteration == 294);
  CHECK(d.residuals.at("omega").rows.size() == 294);
  CHECK_FALSE(d.iteration_count);
  CHECK(d.cell_fields.at("SV_P").at(0).empty());
  CHECK(d.cell_fields.at("SV_P").at(17823).empty());
  for (auto &z : c.mesh.face_zones) {
    auto slices = d.zone_field(FieldLocation::face, z.zone_id, "SV_FLUX");
    std::uint64_t n = 0;
    for (auto &s : slices)
      n += s.last - s.first + 1;
    CHECK(n == z.size());
  }
}

TEST_CASE("Real mismatched case data pairs fail primary field coverage") {
  auto case_path = examples() / "mixing_elbow/mixing_elbow.cas.h5";
  auto data_path = examples() / "file_session/elbow1.dat.h5";
  if (!available(case_path) || !available(data_path))
    return;
  auto c = read_case(case_path);
  CHECK_THROWS_WITH(read_data(data_path, c),
                    doctest::Contains("does not cover expected fluid cell IDs"));
}
TEST_CASE("Real elbow1 data fields match mesh zones") {
  auto path = examples() / "file_session/elbow1.cas.h5";
  if (!available(path))
    return;
  auto c = read_case(path);
  auto data = examples() / "file_session/elbow1.dat.h5";
  if (!available(data))
    return;
  auto d = read_data(data, c);
  CHECK(d.cell_fields.at("SV_P").size() == 13852);
  CHECK(d.cell_fields.at("SV_D").size() == 13852);
  CHECK(d.face_fields.at("SV_WALL_SHEAR").size() == 3630);
  CHECK(d.face_fields.at("SV_FLUX").size() == 44500);
  CHECK(d.residuals.at("continuity").rows.size() == 198);
  CHECK(d.last_residual_iteration == 198);
  CHECK(d.zone_field(FieldLocation::face, 3, "SV_WALL_SHEAR")[0].values.size() == 3630 * 3);
}
#endif
