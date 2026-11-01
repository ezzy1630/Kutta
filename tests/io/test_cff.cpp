#include "../fluent_examples.hpp"
#include "kutta/core/path.hpp"
#include <doctest/doctest.h>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
#include "test_support.hpp"
#if KUTTA_HAS_HDF5
#include "cff_fixture.hpp"
#include <hdf5.h>
#include <zlib.h>

#include <bit>

#include <numeric>
using namespace kutta;
using namespace kutta::detail;
namespace {
using namespace cff_fixture;
}  // namespace
TEST_CASE(
    "synthetic Fluent CFF HDF5: 3D tet/hex, polygon faces, section ranges and "
    "zone names") {
  MeshTemp temp;
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
  for (const auto& m :
       {make_box_hex(2, 1, 1, {0, 0, 0}, {2, 1, 1}),
        make_box_tet(2, 1, 1, {0, 0, 0}, {2, 1, 1}), extrude_2d(p)})
    for (bool variable : {false, true}) {
      auto path = temp.path("synthetic.cas.h5");
      synthetic(path, m, variable);
      auto got = read_mesh(path);
      validate(got);
      REQUIRE(got.n_cells == m.n_cells);
      CHECK(got.n_faces() == m.n_faces());
      CHECK(got.n_internal_faces() == m.n_internal_faces());
      auto before = compute_geometry(m), after = compute_geometry(got);
      for (Index c = 0; c < m.n_cells; ++c) {
        CHECK(after.cell_volume[c] == doctest::Approx(before.cell_volume[c]));
        CHECK(norm(after.cell_centroid[c] - before.cell_centroid[c]) < 1e-12);
      }
      CHECK(got.patches[0].source_zone_id == 10);
      CHECK(got.patches[0].name == "shell");
      CHECK(got.patches[0].kind == PatchKind::Wall);
      CHECK(got.patches[0].source_type == "wall");
      CHECK(got.cell_zones[0].name == "water");
    }
}
TEST_CASE(
    "synthetic Fluent CFF HDF5: 2D extrusion and malformed dataset rejection") {
  MeshTemp temp;
  auto path = temp.path("synthetic.msh.h5");
  auto m = make_box_hex(1, 1, 1, {0, 0, 0}, {1, 1, 1});
  synthetic(path, m, false, true);
  auto q = read_fluent_cff(path, .4);
  CHECK(compute_geometry(q).cell_volume[0] == doctest::Approx(.4));
  CHECK(q.patches[1].kind == PatchKind::Symmetry);
  CHECK(q.patches[2].name == "back");
  auto f = H5Fopen(kutta::utf8_path(path).c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  REQUIRE(f >= 0);
  REQUIRE(H5Ldelete(f, "/meshes/1/faces/nodes/1/nnodes", H5P_DEFAULT) >= 0);
  dataset(f, "/meshes/1/faces/nodes/1/nnodes", {99});
  H5Fclose(f);
  auto handles = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
  CHECK_THROWS_WITH_AS(read_fluent_cff(path),
                       doctest::Contains("nnodes length"), Error);
  CHECK(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) == handles);
}

TEST_CASE("CFF cell type sections use ranges, with complete coverage and count checks") {
  MeshTemp temp;
  auto path = temp.path("sections.msh.h5");
  synthetic(path, make_box_hex(2, 1, 1, {0, 0, 0}, {2, 1, 1}), false);
  auto f = H5Fopen(kutta::utf8_path(path).c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  REQUIRE(f >= 0);
  const std::string base = "/meshes/1/cells/ctype";
  REQUIRE(H5Ldelete(f, (base + "/1").c_str(), H5P_DEFAULT) >= 0);
  auto replace_attr = [&](const std::string& object, const char* name,
                          std::int64_t value) {
    auto o = H5Oopen(f, object.c_str(), H5P_DEFAULT);
    auto a = H5Aopen(o, name, H5P_DEFAULT);
    REQUIRE(a >= 0);
    REQUIRE(H5Awrite(a, H5T_NATIVE_INT64, &value) >= 0);
    H5Aclose(a);
    H5Oclose(o);
  };
  replace_attr(base, "nSections", 2);
  for (int part = 1; part <= 2; ++part) {
    auto section = base + "/" + std::to_string(part);
    group(f, section);
    attr(f, section, "minId", 3 - part);
    attr(f, section, "maxId", 3 - part);
    attr(f, section, "elementType", 4);
  }
  REQUIRE(H5Fflush(f, H5F_SCOPE_GLOBAL) >= 0);
  SUBCASE("valid sections in reversed range order") {
    auto m = read_mesh(path);
    CHECK(m.n_cells == 2);
    CHECK(compute_geometry(m).cell_volume[1] == doctest::Approx(1));
  }
  SUBCASE("overlapping sections") {
    replace_attr(base + "/2", "minId", 2);
    replace_attr(base + "/2", "maxId", 2);
    H5Fflush(f, H5F_SCOPE_GLOBAL);
    CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains("overlapping"), Error);
  }
  SUBCASE("missing section coverage") {
    replace_attr(base, "nSections", 1);
    H5Fflush(f, H5F_SCOPE_GLOBAL);
    CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains("missing ctype"), Error);
  }
  SUBCASE("invalid cell element type") {
    replace_attr(base + "/2", "elementType", 32);
    H5Fflush(f, H5F_SCOPE_GLOBAL);
    CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains("invalid cell element"), Error);
  }
  SUBCASE("uniform zone and section types disagree") {
    REQUIRE(H5Ldelete(f, "/meshes/1/cells/zoneTopology/cellType", H5P_DEFAULT) >= 0);
    dataset(f, "/meshes/1/cells/zoneTopology/cellType", {2});
    H5Fflush(f, H5F_SCOPE_GLOBAL);
    CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains("ctype disagrees"), Error);
  }
  SUBCASE("invalid mixed cell element type") {
    replace_attr(base + "/1", "elementType", 0);
    dataset(f, base + "/1/cell-types", {8});
    H5Fflush(f, H5F_SCOPE_GLOBAL);
    CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains("invalid cell element"),
                         Error);
  }
  SUBCASE("mixed cell type payload length") {
    replace_attr(base + "/1", "elementType", 0);
    dataset(f, base + "/1/cell-types", {});
    H5Fflush(f, H5F_SCOPE_GLOBAL);
    CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains("cell-types length"), Error);
  }
  SUBCASE("declared mesh counts must agree") {
    for (const auto& [name, actual] :
         {std::pair{"cellCount", 2}, std::pair{"faceCount", 11},
          std::pair{"nodeCount", 12}}) {
      replace_attr("/meshes/1", name, 999);
      H5Fflush(f, H5F_SCOPE_GLOBAL);
      auto message = std::string(name) + " disagrees with topology";
      CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains(message.c_str()), Error);
      replace_attr("/meshes/1", name, actual);
    }
  }
  SUBCASE("node coordinate coverage must match topology") {
    REQUIRE(H5Ldelete(f, "/meshes/1/nodes/coords/1", H5P_DEFAULT) >= 0);
    H5Fflush(f, H5F_SCOPE_GLOBAL);
    CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains("missing section data"), Error);
  }
  REQUIRE(H5Fclose(f) >= 0);
}

TEST_CASE("CFF legacy deflate partial chunks preserve numeric conversion and reject corruption") {
  MeshTemp temp;
  auto path = temp.path("partial.msh.h5");
  auto original = make_box_hex(2, 1, 1, {0, 0, 0}, {2, 1, 1});
  synthetic(path, original, false);
  auto f = H5Fopen(kutta::utf8_path(path).c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  REQUIRE(f >= 0);
  const std::string coords = "/meshes/1/nodes/coords/1";
  REQUIRE(H5Ldelete(f, coords.c_str(), H5P_DEFAULT) >= 0);
  enum class Stream {
    Partial, Full, Corrupt, Short, Overlong, Trailing, ShortInterior, Oversized,
    Truncated
  };
  auto stream = Stream::Partial;
  bool flat = false, narrow = false;
  SUBCASE("valid short final chunk") {}
  SUBCASE("full padded final chunk") { stream = Stream::Full; }
  SUBCASE("rank one chunks") { flat = true; }
  SUBCASE("float32 widens to native double") { narrow = true; }
  SUBCASE("corrupt deflate stream bypasses compatibility fallback") {
    stream = Stream::Corrupt;
  }
#if H5_VERSION_GE(2, 2, 0)
  // Earlier HDF5 versions accept undersized streams without invoking fallback.
  SUBCASE("complete stream with missing numeric values") {
    stream = Stream::Short;
  }
  SUBCASE("complete stream with excess numeric values") {
    stream = Stream::Overlong;
  }
  SUBCASE("trailing bytes after complete stream") { stream = Stream::Trailing; }
  SUBCASE("short interior chunk") { stream = Stream::ShortInterior; }
  SUBCASE("decoded stream exceeds full buffer") { stream = Stream::Oversized; }
  SUBCASE("truncated compressed stream") { stream = Stream::Truncated; }
#endif
  hsize_t dims[2]{flat ? 18u : 6u, 3}, chunk[2]{flat ? 12u : 4u, 3};
  const int rank = flat ? 1 : 2;
  auto space = H5Screate_simple(rank, dims, nullptr);
  auto props = H5Pcreate(H5P_DATASET_CREATE);
  REQUIRE(H5Pset_chunk(props, rank, chunk) >= 0);
  REQUIRE(H5Pset_deflate(props, 1) >= 0);
  auto d = H5Dcreate2(f, coords.c_str(), narrow ? H5T_IEEE_F32BE : H5T_IEEE_F64BE,
                     space, H5P_DEFAULT, props, H5P_DEFAULT);
  REQUIRE(d >= 0);
  attr(f, coords, "minId", 7);
  attr(f, coords, "maxId", 12);
  const std::size_t width = narrow ? 4 : 8;
  for (hsize_t row : {0, 4}) {
    std::vector<unsigned char> bytes;
    for (hsize_t i = row; i < std::min<hsize_t>(row + 4, 6); ++i) {
      auto point = original.points[6 + i];
      for (auto x : {point.x, point.y, point.z}) {
        const std::uint64_t word = narrow
            ? std::bit_cast<std::uint32_t>(static_cast<float>(x))
            : std::bit_cast<std::uint64_t>(x);
        for (int j = static_cast<int>(width) - 1; j >= 0; --j)
          bytes.push_back((word >> (8 * j)) & 255);
      }
    }
    if ((stream == Stream::Short && row == 4) ||
        (stream == Stream::ShortInterior && row == 0))
      bytes.resize(bytes.size() - width);
    if (stream == Stream::Overlong && row == 4) bytes.resize(bytes.size() + width);
    if (stream == Stream::Full && row == 4) bytes.resize(12 * width);
    if (stream == Stream::Oversized && row == 4) bytes.resize(13 * width);
    auto length = compressBound(bytes.size());
    std::vector<unsigned char> compressed(length);
    REQUIRE(compress2(compressed.data(), &length, bytes.data(), bytes.size(), 1) == Z_OK);
    if (stream == Stream::Corrupt && row == 4) compressed[length - 1] ^= 255;
    if (stream == Stream::Truncated && row == 4) --length;
    if (stream == Stream::Trailing && row == 4) {
      compressed.resize(length + 1);
      compressed[length++] = 42;
    }
    hsize_t offset[2]{flat ? row * 3 : row, 0};
    REQUIRE(H5Dwrite_chunk(d, H5P_DEFAULT, 0, offset, length, compressed.data()) >= 0);
  }
  H5Dclose(d);
  H5Pclose(props);
  H5Sclose(space);
  REQUIRE(H5Fclose(f) >= 0);
  if (stream == Stream::Corrupt) {
    const auto expected = "CFF: HDF5 failure reading " + coords;
    CHECK_THROWS_WITH_AS(read_mesh(path), doctest::Contains(expected.c_str()), Error);
  } else if (stream != Stream::Partial && stream != Stream::Full) {
    CHECK_THROWS_AS(read_mesh(path), Error);
  } else {
    auto mesh = read_mesh(path);
    CHECK(mesh.n_points() == original.n_points());
    CHECK(compute_geometry(mesh).cell_volume[1] == doctest::Approx(1));
    for (std::size_t i = 0; i < mesh.points.size(); ++i)
      CHECK(norm(mesh.points[i] - original.points[i]) == 0);
  }
}
TEST_CASE("CFF reads a file whose name is not ASCII") {
  // Windows passes HDF5 the UTF-8 bytes, not the ANSI code page conversion.
  const auto source = fluent_examples() / "simple_duct/simple-duct.msh.h5";
  const auto dir = std::filesystem::temp_directory_path() / std::filesystem::path(u8"kutta-ü-网格");
  std::filesystem::create_directories(dir);
  const auto copy = dir / std::filesystem::path(u8"düct-网格.msh.h5");
  std::filesystem::copy_file(source, copy, std::filesystem::copy_options::overwrite_existing);
  const auto mesh = read_mesh(copy);
  CHECK(mesh.n_cells == read_mesh(source).n_cells);
  std::filesystem::remove_all(dir);
}
#else
TEST_CASE("CFF disabled build explains missing HDF5") {
  CHECK_THROWS_WITH_AS(kutta::read_fluent_cff("absent.msh.h5"),
                       doctest::Contains("KUTTA_HAS_HDF5=0"), kutta::Error);
}
#endif
