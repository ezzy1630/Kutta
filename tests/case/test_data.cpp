#include "kutta/case/fluent_data.hpp"
#include <doctest/doctest.h>
#include <algorithm>
#if KUTTA_HAS_HDF5
#include "kutta/case/hdf5.hpp"
#include <chrono>
#include <random>
using namespace kutta::fluent;
namespace {
struct Fixture {
  std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      ("kutta-e-data-" +
       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + "-" +
       std::to_string(std::random_device{}()) + ".h5");
  h5::Handle file{H5Fcreate(path.string().c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT),
                  H5Fclose};
  ~Fixture() {
    std::error_code ec;
    std::filesystem::remove(path, ec);
  }
  void group(const std::string &path) {
    h5::Handle prop(H5Pcreate(H5P_LINK_CREATE), H5Pclose);
    h5::check(H5Pset_create_intermediate_group(prop, 1), path);
    h5::Handle group(H5Gcreate2(file, path.c_str(), prop, H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
  }
  void str(const std::string &path, const std::string &value) {
    h5::Handle space(H5Screate(H5S_SCALAR), H5Sclose), type(H5Tcopy(H5T_C_S1), H5Tclose);
    h5::check(H5Tset_size(type, value.size() + 1), path);
    h5::Handle ds(
        H5Dcreate2(file, path.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
        H5Dclose);
    h5::check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value.c_str()), path);
  }
  void attr(hid_t ds, const char *name, std::uint64_t n) {
    h5::Handle space(H5Screate(H5S_SCALAR), H5Sclose);
    h5::Handle attr(H5Acreate2(ds, name, H5T_NATIVE_UINT64, space, H5P_DEFAULT, H5P_DEFAULT),
                    H5Aclose);
    h5::check(H5Awrite(attr, H5T_NATIVE_UINT64, &n), name);
  }
  void numbers(const std::string &path, std::vector<double> values, std::vector<hsize_t> dims,
               std::optional<std::uint64_t> first = {}, std::optional<std::uint64_t> last = {}) {
    h5::Handle space(H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr),
                     H5Sclose);
    h5::Handle ds(H5Dcreate2(file, path.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT,
                             H5P_DEFAULT),
                  H5Dclose);
    REQUIRE(values.size() == h5::count(dims));
    h5::check(H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()), path);
    if (first)
      attr(ds, "minId", *first);
    if (last)
      attr(ds, "maxId", *last);
  }
  void integers(const std::string &path, std::vector<std::int64_t> values) {
    const hsize_t size = values.size();
    h5::Handle space(H5Screate_simple(1, &size, nullptr), H5Sclose);
    h5::Handle ds(H5Dcreate2(file, path.c_str(), H5T_STD_I64LE, space, H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT), H5Dclose);
    h5::check(H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()), path);
  }
  void base() {
    group("/settings");
    str("/settings/Data Variables", "(37 ((number-of-iterations 20)))");
    group("/results/1/phase-1/cells/SV_P");
    group("/results/1/phase-1/faces/SV_WALL_SHEAR");
  }
  void flush() { h5::check(H5Fflush(file, H5F_SCOPE_GLOBAL), "flush fixture"); }
};
FluentCase matching() {
  FluentCase c;
  c.mesh.cells = 3;
  c.mesh.faces = 5;
  c.mesh.cell_zones = {{8, "fluid", 101, 103, 1}};
  c.mesh.face_zones = {{4, "wall", 11, 12, 3}, {5, "outlet", 21, 23, 5}};
  return c;
}
} // namespace
TEST_CASE("Data chunks use Fluent IDs not ordinal positions and vectors are row major") {
  Fixture f;
  f.base();
  f.numbers("/results/1/phase-1/cells/SV_P/1", {30}, {1}, 103, 103);
  f.numbers("/results/1/phase-1/cells/SV_P/10", {10, 20}, {2}, 101, 102);
  f.numbers("/results/1/phase-1/faces/SV_WALL_SHEAR/2", {1, 2, 3, 4, 5, 6}, {2, 3}, 11, 12);
  f.group("/results/residuals/phase-1/continuity");
  f.numbers("/results/residuals/phase-1/continuity/iterations", {1, 11}, {2});
  f.numbers("/results/residuals/phase-1/continuity/data", {3, 5, 7, 9, 2, 4, 6, 8}, {2, 4});
  f.flush();
  auto d = read_data(f.path, matching());
  CHECK(d.cell_fields.at("SV_P").at(101)[0] == 10);
  CHECK(d.cell_fields.at("SV_P").at(103)[0] == 30);
  CHECK(d.cell_fields.at("SV_P").at(1).empty());
  CHECK(d.cell_fields.at("SV_P").size() == 3);
  auto v = d.face_fields.at("SV_WALL_SHEAR").at(12);
  REQUIRE(v.size() == 3);
  CHECK(v[0] == 4);
  CHECK(v[1] == 5);
  CHECK(v[2] == 6);
  CHECK(d.zone_field(FieldLocation::face, 4, "SV_WALL_SHEAR")[0].values.size() == 6);
  CHECK(d.zone_field(FieldLocation::face, 5, "SV_WALL_SHEAR").empty());
  CHECK_THROWS_AS(d.zone_field(FieldLocation::face, 999, "SV_WALL_SHEAR"), std::out_of_range);
  CHECK_THROWS_AS(d.zone_field(FieldLocation::face, 4, "missing"), std::out_of_range);
  CHECK(d.iteration_count == 20);
  CHECK(d.last_residual_iteration == 11);
  CHECK(d.residuals.at("continuity").rows[1][3] == 8);
}
TEST_CASE("Data rejects invalid ranges shapes overlaps and external links") {
  for (int mode = 0; mode < 8; ++mode) {
    INFO("malformed fixture mode ", mode);
    Fixture f;
    f.base();
    const std::string cell = "/results/1/phase-1/cells/SV_P/";
    if (mode == 0)
      f.numbers(cell + "1", {1}, {1}, 100, 100);
    if (mode == 1)
      f.numbers(cell + "1", {1}, {1}, 101, 102);
    if (mode == 2) {
      f.numbers(cell + "1", {1}, {1}, 101, 101);
      f.numbers(cell + "2", {2}, {1}, 101, 101);
    }
    if (mode == 3) {
      f.numbers(cell + "1", {1, 2, 3}, {3}, 101, 103);
      f.numbers("/results/1/phase-1/faces/SV_WALL_SHEAR/1", std::vector<double>(10), {10}, 12, 21);
    }
    if (mode == 4)
      f.numbers(cell + "1", {1}, {1}, 101, {});
    if (mode == 5)
      h5::check(H5Lcreate_external("/untrusted/other.h5", "/data", f.file, (cell + "1").c_str(),
                                   H5P_DEFAULT, H5P_DEFAULT),
                "external link fixture");
    if (mode == 6) {
      f.numbers(cell + "1", {1}, {1}, 101, 101);
      f.numbers(cell + "2", {2, 3}, {1, 2}, 102, 102);
    }
    if (mode == 7) {
      f.numbers(cell + "1", {1, 2, 3}, {3}, 101, 103);
      f.numbers("/results/1/phase-1/faces/SV_WALL_SHEAR/1", {1, 2, 3, 4}, {2, 2}, 11, 12);
    }
    f.flush();
    const char *errors[] = {"outside matching", "shape disagrees", "overlapping field",
                            "outside matching", "missing/zero ID", "non-local HDF5 link",
                            "scalar field",     "wall shear"};
    CHECK_THROWS_WITH(read_data(f.path, matching()), doctest::Contains(errors[mode]));
  }
}
#endif
#if KUTTA_HAS_HDF5
TEST_CASE("Data Case File provenance is compared by basename without claiming mesh identity") {
  for (const auto &saved : {"expected.cas.h5", "C:\\run\\expected.cas.h5", "different.cas.h5"}) {
    Fixture f;
    f.base();
    f.str("/settings/Case File", saved);
    f.numbers("/results/1/phase-1/cells/SV_P/1", {1, 2, 3}, {3}, 101, 103);
    f.numbers("/results/1/phase-1/faces/SV_WALL_SHEAR/1", {1, 2, 3, 4, 5, 6}, {2, 3}, 11, 12);
    f.flush();
    auto c = matching();
    c.source_path = "/fixtures/expected.cas.h5";
    auto d = read_data(f.path, c);
    const bool mismatch = std::any_of(d.issues.begin(), d.issues.end(), [](auto &i) {
      return i.path == "data/Case File" && i.reason.find("differs") != std::string::npos;
    });
    CHECK(mismatch == (std::string_view(saved) == "different.cas.h5"));
    CHECK(d.case_file == saved);
    for (const auto &i : d.issues)
      if (i.path == "data/Case File")
        CHECK(i.classification == IssueClass::informational);
  }
}
TEST_CASE("Case topology requires matching thread records names and types") {
  for (bool cell : {false, true})
  for (int mode = 0; mode < 4; ++mode) {
    Fixture f;
    f.group("/settings");
    f.str("/settings/Rampant Variables", "(37 ((case-config ((rp-3d? . #t)))))");
    const auto type = cell ? "fluid" : "wall";
    const auto name = mode == 2 ? "other" : "body";
    const auto actual_type = mode == 3 ? (cell ? "solid" : "symmetry") : type;
    auto thread = mode == 1 ? std::string{} : "(39 (8 " + std::string(actual_type) + " " + name + " 1) ())";
    f.str("/settings/Thread Variables", thread);
    const std::string root = cell ? "/meshes/1/cells/zoneTopology" : "/meshes/1/faces/zoneTopology";
    f.group(root);
    h5::Handle mesh(H5Gopen2(f.file, "/meshes/1", H5P_DEFAULT), H5Gclose);
    f.attr(mesh, "dimension", 3);
    f.attr(mesh, cell ? "cellCount" : "faceCount", 1);
    const std::string p = root + "/";
    f.integers(p + "id", {8});
    f.integers(p + "minId", {1});
    f.integers(p + "maxId", {1});
    f.integers(p + "zoneType", {cell ? 1 : 3});
    f.str(p + "name", "body");
    f.flush();
    auto c = read_case(f.path);
    const bool mismatch = std::any_of(c.issues.begin(), c.issues.end(), [](auto &i) {
      return i.path.starts_with("mesh/zone/8");
    });
    CHECK(mismatch == (mode != 0));
  }
}
TEST_CASE("CFF reads large complete Thread Variables strings with a bounded parser") {
  Fixture f;
  f.group("/settings");
  f.str("/settings/Rampant Variables", "(37 ((case-config ((rp-3d? . #t)))))");
  std::string threads;
  for (int zone = 1; zone <= 1200; ++zone) {
    threads += "(39 (" + std::to_string(zone) + " wall body 1) (";
    for (int key = 0; key < 560; ++key)
      threads += "(switch . #f) ";
    threads += "))\n";
  }
  f.str("/settings/Thread Variables", threads);
  f.flush();
  CHECK(read_case(f.path).boundaries.size() == 1200);
}
#endif
