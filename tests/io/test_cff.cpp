#include <doctest/doctest.h>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
#include "test_support.hpp"
#if KUTTA_HAS_HDF5
#include <hdf5.h>
#include <zlib.h>

#include <bit>

#include <numeric>
using namespace kutta;
using namespace kutta::detail;
namespace {
// Independent test writer for the dataset layout described by VTK's CFF
// reader. This is synthetic format coverage, not a real Fluent acceptance file.
void group(hid_t f, const std::string& path) {
  auto lc = H5Pcreate(H5P_LINK_CREATE);
  H5Pset_create_intermediate_group(lc, 1);
  auto g = H5Gcreate2(f, path.c_str(), lc, H5P_DEFAULT, H5P_DEFAULT);
  REQUIRE(g >= 0);
  H5Gclose(g);
  H5Pclose(lc);
}
void attr(hid_t f, const std::string& object, const char* name,
          std::int64_t value) {
  auto o = H5Oopen(f, object.c_str(), H5P_DEFAULT);
  REQUIRE(o >= 0);
  auto s = H5Screate(H5S_SCALAR);
  auto a = H5Acreate2(o, name, H5T_STD_I64LE, s, H5P_DEFAULT, H5P_DEFAULT);
  REQUIRE(a >= 0);
  REQUIRE(H5Awrite(a, H5T_NATIVE_INT64, &value) >= 0);
  H5Aclose(a);
  H5Sclose(s);
  H5Oclose(o);
}
void dataset(hid_t f, const std::string& path,
             const std::vector<std::int64_t>& values) {
  hsize_t dim = values.size();
  auto s = H5Screate_simple(1, &dim, nullptr);
  auto d = H5Dcreate2(f, path.c_str(), H5T_STD_I64LE, s, H5P_DEFAULT,
                      H5P_DEFAULT, H5P_DEFAULT);
  REQUIRE(d >= 0);
  if (!values.empty())
    REQUIRE(H5Dwrite(d, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     values.data()) >= 0);
  H5Dclose(d);
  H5Sclose(s);
}
void labels(hid_t f, const std::string& path,
            const std::vector<std::string>& names, bool variable) {
  auto t = H5Tcopy(H5T_C_S1);
  std::string joined;
  for (const auto& n : names) joined += n + ";";
  H5Tset_size(t, variable ? H5T_VARIABLE : joined.size() + 1);
  hsize_t dim = names.size();
  auto s =
      variable ? H5Screate_simple(1, &dim, nullptr) : H5Screate(H5S_SCALAR);
  auto d =
      H5Dcreate2(f, path.c_str(), t, s, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  REQUIRE(d >= 0);
  if (variable) {
    std::vector<const char*> ptr;
    for (const auto& n : names) ptr.push_back(n.c_str());
    REQUIRE(H5Dwrite(d, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr.data()) >= 0);
  } else
    REQUIRE(H5Dwrite(d, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, joined.c_str()) >= 0);
  H5Dclose(d);
  H5Sclose(s);
  H5Tclose(t);
}
void topology(hid_t f, const std::string& base,
              const std::vector<std::int64_t>& ids,
              const std::vector<std::int64_t>& first,
              const std::vector<std::int64_t>& last,
              const std::vector<std::int64_t>& types,
              const std::vector<std::string>& names, bool variable,
              int dimension = 3) {
  group(f, base);
  attr(f, base, "nZones", ids.size());
  dataset(f, base + "/id", ids);
  dataset(f, base + "/minId", first);
  dataset(f, base + "/maxId", last);
  dataset(f, base + "/zoneType", types);
  dataset(f, base + "/dimension",
          std::vector<std::int64_t>(ids.size(), dimension));
  if (base.find("/faces/") != std::string::npos) {
    dataset(f, base + "/childZoneId", std::vector<std::int64_t>(ids.size(), 0));
    dataset(f, base + "/shadowZoneId",
            std::vector<std::int64_t>(ids.size(), 0));
    dataset(f, base + "/flags", std::vector<std::int64_t>(ids.size(), 0));
  } else if (base.find("/cells/") != std::string::npos) {
    dataset(f, base + "/childZoneId", std::vector<std::int64_t>(ids.size(), 0));
  }

  labels(f, base + "/name", names, variable);
}
void synthetic(const std::filesystem::path& path, const PolyMesh& m,
               bool variable, bool planar = false) {
  auto f =
      H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  REQUIRE(f >= 0);
  const std::string root = "/meshes/1";
  group(f, root);
  attr(f, root, "dimension", planar ? 2 : 3);
  auto points =
      planar ? std::vector<Vec3d>{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}
             : m.points;
  auto nf = planar ? 4 : m.n_faces(), ni = planar ? 0 : m.n_internal_faces();
  attr(f, root, "nodeOffset", 0);
  attr(f, root, "nodeCount", points.size());
  attr(f, root, "faceOffset", 0);
  attr(f, root, "faceCount", nf);
  attr(f, root, "cellOffset", 0);
  attr(f, root, "cellCount", m.n_cells);
  auto split = static_cast<std::int64_t>(points.size() / 2);
  topology(f, root + "/nodes/zoneTopology", {8, 41}, {1, split + 1},
           {split, static_cast<std::int64_t>(points.size())}, {0, 0},
           {"nodes-a", "nodes-b"}, variable, planar ? 2 : 3);
  group(f, root + "/nodes/coords");
  std::vector<double> coords;
  for (auto p : points) {
    coords.push_back(p.x);
    coords.push_back(p.y);
    if (!planar) coords.push_back(p.z);
  }
  // Both layouts occur in Fluent: zone-id names and numbered sections.
  // Reverse the numbered section ranges to ensure placement uses minId/maxId.
  for (int part = 0; part < 2; ++part) {
    auto first = part == 0 ? 0 : split;
    auto last = part == 0 ? split : static_cast<std::int64_t>(points.size());
    auto name = variable ? (part == 0 ? "8" : "41") : (part == 0 ? "2" : "1");
    auto path_coords = root + "/nodes/coords/" + name;
    hsize_t dims[2] = {static_cast<hsize_t>(last - first),
                       static_cast<hsize_t>(planar ? 2 : 3)};
    auto s = H5Screate_simple(2, dims, nullptr);
    auto d = H5Dcreate2(f, path_coords.c_str(), H5T_IEEE_F64LE, s,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    REQUIRE(d >= 0);
    REQUIRE(H5Dwrite(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     coords.data() + first * (planar ? 2 : 3)) >= 0);
    H5Dclose(d);
    H5Sclose(s);
    attr(f, path_coords, "minId", first + 1);
    attr(f, path_coords, "maxId", last);
  }
  topology(f, root + "/cells/zoneTopology", {20}, {1}, {m.n_cells}, {1},
           {"water"}, variable, planar ? 2 : 3);
  dataset(f, root + "/cells/zoneTopology/cellType",
          {0});  // mixed, connectivity is face-based

  group(f, root + "/cells/ctype");
  attr(f, root + "/cells/ctype", "nSections", 1);
  group(f, root + "/cells/ctype/1");
  attr(f, root + "/cells/ctype/1", "minId", 1);
  attr(f, root + "/cells/ctype/1", "maxId", m.n_cells);
  attr(f, root + "/cells/ctype/1", "elementType", 0);
  auto cell_faces = kutta::build_cell_faces(m);
  std::vector<std::int64_t> cell_types;
  for (Index c = 0; c < m.n_cells; ++c) {
    auto n = cell_faces.offsets[c + 1] - cell_faces.offsets[c];
    cell_types.push_back(planar ? 3 : n == 4 ? 2 : n == 6 ? 4 : 7);
  }
  dataset(f, root + "/cells/ctype/1/cell-types", cell_types);
  if (ni)
    topology(f, root + "/faces/zoneTopology", {2, 10}, {1, ni + 1}, {ni, nf},
             {2, 3}, {"interior", "shell"}, variable, planar ? 2 : 3);
  else
    topology(f, root + "/faces/zoneTopology", {10}, {1}, {nf}, {3}, {"shell"},
             variable, planar ? 2 : 3);
  // Element topology differs from the boundary code; this guards their mapping.
  dataset(f, root + "/faces/zoneTopology/faceType",
          ni ? std::vector<std::int64_t>{3, 0} : std::vector<std::int64_t>{0});
  group(f, root + "/faces/nodes");
  attr(f, root + "/faces/nodes", "nSections", 2);
  for (int section = 1; section <= 2; ++section) {
    auto base = root + "/faces/nodes/" + std::to_string(section);
    group(f, base);
    Index begin = section == 1 ? 0 : nf / 2, end = section == 1 ? nf / 2 : nf;
    attr(f, base, "minId", begin + 1);
    attr(f, base, "maxId", end);
    std::vector<std::int64_t> counts, nodes_out;
    for (Index i = begin; i < end; ++i) {
      if (planar) {
        counts.push_back(2);
        nodes_out.push_back(i + 1);
        nodes_out.push_back((i + 1) % 4 + 1);
      } else {
        auto ns = nodes(m, i);
        counts.push_back(ns.size());
        for (auto n = ns.rbegin(); n != ns.rend(); ++n)
          nodes_out.push_back(*n + 1);
      }
    }
    dataset(f, base + "/nnodes", counts);
    dataset(f, base + "/nodes", nodes_out);
  }
  group(f, root + "/faces/c0");
  attr(f, root + "/faces/c0", "nSections", 1);
  std::vector<std::int64_t> owners;
  for (Index i = 0; i < nf; ++i) owners.push_back(planar ? 1 : m.owner[i] + 1);
  dataset(f, root + "/faces/c0/1", owners);
  attr(f, root + "/faces/c0/1", "minId", 1);
  attr(f, root + "/faces/c0/1", "maxId", nf);
  group(f, root + "/faces/c1");
  auto c1_count = variable ? nf : ni;
  attr(f, root + "/faces/c1", "nSections", c1_count ? 1 : 0);
  if (c1_count) {
    std::vector<std::int64_t> nb;
    for (Index face = 0; face < c1_count; ++face)
      nb.push_back(face < ni ? m.neighbour[face] + 1 : 0);
    dataset(f, root + "/faces/c1/1", nb);
    attr(f, root + "/faces/c1/1", "minId", 1);
    attr(f, root + "/faces/c1/1", "maxId", c1_count);
  }
  REQUIRE(H5Fclose(f) >= 0);
}
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
  auto f = H5Fopen(path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
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
  auto f = H5Fopen(path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
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
  auto f = H5Fopen(path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  REQUIRE(f >= 0);
  const std::string coords = "/meshes/1/nodes/coords/1";
  REQUIRE(H5Ldelete(f, coords.c_str(), H5P_DEFAULT) >= 0);
  hsize_t dims[2]{6, 3}, chunk[2]{4, 3};
  auto space = H5Screate_simple(2, dims, nullptr);
  auto props = H5Pcreate(H5P_DATASET_CREATE);
  REQUIRE(H5Pset_chunk(props, 2, chunk) >= 0);
  REQUIRE(H5Pset_deflate(props, 1) >= 0);
  auto d = H5Dcreate2(f, coords.c_str(), H5T_IEEE_F64BE, space, H5P_DEFAULT,
                     props, H5P_DEFAULT);
  REQUIRE(d >= 0);
  attr(f, coords, "minId", 7);
  attr(f, coords, "maxId", 12);
  bool corrupt = false, short_data = false;
  SUBCASE("valid short final chunk") {}
  SUBCASE("corrupt deflate stream") { corrupt = true; }
#if H5_VERSION_GE(2, 2, 0)
  // Older HDF5 accepts undersized streams itself; this exercises our fallback
  // on versions that reject them and must never manufacture missing values.
  SUBCASE("complete stream with missing numeric values") { short_data = true; }
#endif
  for (hsize_t row : {0, 4}) {
    std::vector<unsigned char> bytes;
    for (hsize_t i = row; i < std::min(row + 4, dims[0]); ++i) {
      auto point = original.points[6 + i];
      for (auto x : {point.x, point.y, point.z}) {
        auto word = std::bit_cast<std::uint64_t>(x);
        for (int j = 7; j >= 0; --j) bytes.push_back((word >> (8 * j)) & 255);
      }
    }
    if (short_data && row == 4) bytes.resize(bytes.size() - sizeof(double));
    auto length = compressBound(bytes.size());
    std::vector<unsigned char> compressed(length);
    REQUIRE(compress2(compressed.data(), &length, bytes.data(), bytes.size(), 1) == Z_OK);
    if (corrupt && row == 4) compressed[length - 1] ^= 255;
    hsize_t offset[2]{row, 0};
    REQUIRE(H5Dwrite_chunk(d, H5P_DEFAULT, 0, offset, length, compressed.data()) >= 0);
  }
  H5Dclose(d);
  H5Pclose(props);
  H5Sclose(space);
  REQUIRE(H5Fclose(f) >= 0);
  if (corrupt || short_data) {
    H5E_BEGIN_TRY { CHECK_THROWS_AS(read_mesh(path), Error); }
    H5E_END_TRY;
  } else {
    auto mesh = read_mesh(path);
    CHECK(mesh.n_points() == original.n_points());
    CHECK(compute_geometry(mesh).cell_volume[1] == doctest::Approx(1));
    for (std::size_t i = 0; i < mesh.points.size(); ++i)
      CHECK(norm(mesh.points[i] - original.points[i]) == 0);
  }
}
#else
TEST_CASE("CFF disabled build explains missing HDF5") {
  CHECK_THROWS_WITH_AS(kutta::read_fluent_cff("absent.msh.h5"),
                       doctest::Contains("KUTTA_HAS_HDF5=0"), kutta::Error);
}
#endif
