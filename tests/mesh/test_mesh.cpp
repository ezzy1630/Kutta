#include <doctest/doctest.h>

#include <numeric>
#include <set>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
using namespace kutta;
using namespace kutta::detail;
namespace {
void volume_and_closure(const PolyMesh& m, double expected) {
  validate(m);
  auto g = compute_geometry(m);
  auto a = build_cell_faces(m);
  CHECK(std::accumulate(g.cell_volume.begin(), g.cell_volume.end(), 0.0) ==
        doctest::Approx(expected).epsilon(1e-11));
  for (Index c = 0; c < m.n_cells; ++c) {
    Vec3d sum{};
    double scale = 0;
    for (auto i = a.offsets[c]; i < a.offsets[c + 1]; ++i) {
      sum = sum + g.face_area[a.faces[i]] * a.signs[i];
      scale += norm(g.face_area[a.faces[i]]);
    }
    CHECK(norm(sum) < 1e-11 * scale);
  }
}
}  // namespace
TEST_CASE(
    "hex and conforming six-tet boxes: volumes, centroids and linear moments") {
  for (bool tet : {false, true}) {
    auto m = tet ? make_box_tet(3, 2, 2, {0, 0, 0}, {3, 2, 2})
                 : make_box_hex(3, 2, 2, {0, 0, 0}, {3, 2, 2});
    volume_and_closure(m, 12);
    auto g = compute_geometry(m);
    Vec3d moment{};
    for (Index c = 0; c < m.n_cells; ++c)
      moment = moment + g.cell_centroid[c] * g.cell_volume[c];
    CHECK(moment.x == doctest::Approx(18));
    CHECK(moment.y == doctest::Approx(12));
    CHECK(moment.z == doctest::Approx(12));
    CHECK(m.n_cells == (tet ? 72 : 12));
  }
  auto m = make_box_hex(1, 1, 1, {0, 0, 0}, {1, 1, 1});
  auto g = compute_geometry(m);
  CHECK(g.cell_centroid[0].x == doctest::Approx(.5));
  CHECK(g.cell_centroid[0].y == doctest::Approx(.5));
  CHECK(g.cell_centroid[0].z == doctest::Approx(.5));
}
TEST_CASE("shear, skew, scale and translation preserve affine volume moments") {
  for (bool tet : {false, true}) {
    auto m = tet ? make_box_tet(2, 2, 2, {0, 0, 0}, {1, 1, 1})
                 : make_box_hex(2, 2, 2, {0, 0, 0}, {1, 1, 1});
    for (auto& p : m.points)
      p = {1e6 + 2 * p.x + .7 * p.y + .3 * p.z, -2e6 + 3 * p.y + .2 * p.z,
           3e6 + 4 * p.z};
    volume_and_closure(m, 24);
    auto g = compute_geometry(m);
    Vec3d moment{};
    for (Index c = 0; c < m.n_cells; ++c)
      moment = moment +
               (g.cell_centroid[c] - Vec3d{1e6, -2e6, 3e6}) * g.cell_volume[c];
    CHECK(moment.x / 24 == doctest::Approx(1.5));
    CHECK(moment.y / 24 == doctest::Approx(1.6));
    CHECK(moment.z / 24 == doctest::Approx(2));
  }
}
TEST_CASE("prism, polygonal prism and warped face geometry") {
  for (int sides : {3, 5}) {
    PlanarMesh p;
    p.n_cells = 1;
    p.patches = {{"wall", PatchKind::Wall}};
    for (int i = 0; i < sides; ++i) {
      double t = 2 * std::acos(-1.) * i / sides;
      p.points.push_back({std::cos(t), std::sin(t), 0});
      p.edges.push_back({i, (i + 1) % sides});
      p.owner.push_back(0);
      p.neighbour.push_back(-1);
      p.edge_patch.push_back(0);
    }
    auto m = extrude_2d(p, 2);
    volume_and_closure(m, sides * std::sin(2 * std::acos(-1.) / sides));
    auto g = compute_geometry(m);
    CHECK(g.cell_centroid[0].x == doctest::Approx(0).epsilon(1e-12));
    CHECK(g.cell_centroid[0].y == doctest::Approx(0).epsilon(1e-12));
    CHECK(g.cell_centroid[0].z == doctest::Approx(1));
  }
  auto m = make_box_hex(1, 1, 1, {0, 0, 0}, {1, 1, 1});
  m.points[7].z = 1.4;
  // Bilinear top z=1+.4*x*y: exact integral of its center fan is 1.1.
  volume_and_closure(m, 1.1);
  auto g = compute_geometry(m);
  CHECK(g.cell_centroid[0].x == doctest::Approx((.5 + .4 / 6) / 1.1));
  CHECK(g.cell_centroid[0].y == doctest::Approx((.5 + .4 / 6) / 1.1));
  // Area centroid must differ from the unweighted vertex average on a warped
  // quad.
  auto f = m.patches[5].start;
  CHECK(norm(g.face_centroid[f] - face_center(m, f)) > 1e-5);
  CHECK(g.face_centroid[f].x == doctest::Approx(41. / 81));
  CHECK(g.face_centroid[f].y == doctest::Approx(41. / 81));
  CHECK(g.face_centroid[f].z == doctest::Approx(893. / 810));
  CHECK(g.cell_centroid[0].z == doctest::Approx((14.62 / 24) / 1.1));
}
TEST_CASE("validation rejects malformed topology and winding") {
  auto good = make_box_hex(2, 1, 1, {0, 0, 0}, {2, 1, 1});
  SUBCASE("CSR size") {
    auto m = good;
    m.face_offsets.pop_back();
    CHECK_THROWS_WITH_AS(validate(m), doctest::Contains("face_offsets"), Error);
  }
  SUBCASE("node index") {
    auto m = good;
    m.face_nodes[0] = -1;
    CHECK_THROWS_WITH_AS(validate(m), doctest::Contains("node index"), Error);
  }
  SUBCASE("owner order") {
    auto m = good;
    m.neighbour[0] = m.owner[0];
    CHECK_THROWS_WITH_AS(validate(m), doctest::Contains("owner < neighbour"),
                         Error);
  }
  SUBCASE("duplicate patch name") {
    auto m = good;
    m.patches[1].name = m.patches[0].name;
    CHECK_THROWS_WITH_AS(validate(m), doctest::Contains("duplicate patch name"), Error);
  }
  SUBCASE("patch gap") {
    auto m = good;
    ++m.patches[0].start;
    CHECK_THROWS_WITH_AS(validate(m), doctest::Contains("contiguously"), Error);
  }
  SUBCASE("winding") {
    auto m = good;
    std::reverse(m.face_nodes.begin() + m.face_offsets[1],
                 m.face_nodes.begin() + m.face_offsets[2]);
    CHECK_THROWS_AS(validate(m), Error);
  }
  SUBCASE("non-finite") {
    auto m = good;
    m.points[0].x = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(validate(m), Error);
  }
}
TEST_CASE("gather adjacency and RCM preserve topology, zones and geometry") {
  auto m = make_box_tet(4, 2, 1, {0, 0, 0}, {4, 2, 1});
  m.cell_zones = {{"fluid", "fluid", 7, {0, 2, 5}}};
  auto a = build_cell_faces(m);
  auto n = build_cell_neighbours(m);
  auto g = compute_geometry(m);
  auto old = m;
  CHECK(a.faces.size() == m.owner.size() + m.neighbour.size());
  for (Index c = 0; c < m.n_cells; ++c)
    for (auto i = a.offsets[c]; i < a.offsets[c + 1]; ++i) {
      auto f = a.faces[i];
      CHECK((a.signs[i] == 1 ? m.owner[f] : m.neighbour[f]) == c);
    }
  for (Index c = 0; c < m.n_cells; ++c)
    for (auto i = n.offsets[c]; i < n.offsets[c + 1]; ++i) {
      auto b = n.cells[i];
      CHECK(std::binary_search(n.cells.begin() + n.offsets[b],
                               n.cells.begin() + n.offsets[b + 1], c));
    }
  auto p = renumber_cells_rcm(m);
  validate(m);
  auto after = compute_geometry(m);
  CHECK(p.cells_new_to_old.size() == static_cast<std::size_t>(m.n_cells));
  for (Index c = 0; c < m.n_cells; ++c) {
    CHECK(after.cell_volume[c] ==
          doctest::Approx(g.cell_volume[p.cells_new_to_old[c]]));
    CHECK(norm(after.cell_centroid[c] -
               g.cell_centroid[p.cells_new_to_old[c]]) < 1e-12);
  }
  for (Index f = 0; f < m.n_faces(); ++f) {
    auto src = p.faces_new_to_old[f];
    CHECK(norm(after.face_centroid[f] - g.face_centroid[src]) < 1e-12);
    const auto sign = p.cells_new_to_old[m.owner[f]] == old.owner[src] ? 1 : -1;
    CHECK(norm(after.face_area[f] - g.face_area[src] * sign) < 1e-12);
  }
  REQUIRE(m.patches.size() == old.patches.size());
  for (std::size_t i = 0; i < m.patches.size(); ++i) {
    const auto& patch = m.patches[i];
    const auto& original = old.patches[i];
    CHECK(patch.name == original.name);
    CHECK(patch.kind == original.kind);
    CHECK(patch.size == original.size);
    for (Index f = patch.start; f < patch.start + patch.size; ++f) {
      CHECK(p.faces_new_to_old[f] >= original.start);
      CHECK(p.faces_new_to_old[f] < original.start + original.size);
    }
  }
  for (std::size_t i = 0; i < m.cell_zones[0].cells.size(); ++i)
    CHECK(p.cells_new_to_old[m.cell_zones[0].cells[i]] ==
          old.cell_zones[0].cells[i]);
  auto again = old;
  CHECK(renumber_cells_rcm(again).cells_new_to_old == p.cells_new_to_old);
}

TEST_CASE(
    "individual tet centroid is the vertex average; trapezoid face moments") {
  auto m = make_box_tet(2, 1, 1, {0, 0, 0}, {2, 1, 1});
  auto g = compute_geometry(m);
  auto a = build_cell_faces(m);
  for (Index c = 0; c < m.n_cells; ++c) {
    std::set<Index> ns;
    for (auto i = a.offsets[c]; i < a.offsets[c + 1]; ++i)
      for (auto v : nodes(m, a.faces[i])) ns.insert(v);
    REQUIRE(ns.size() == 4);
    Vec3d mean{};
    for (auto v : ns) mean = mean + m.points[v] / 4;
    CHECK(norm(g.cell_centroid[c] - mean) < 1e-12);
  }
  PlanarMesh p;
  p.n_cells = 1;
  p.points = {{0, 0, 0}, {2, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  p.edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  p.owner = {0, 0, 0, 0};
  p.neighbour = {-1, -1, -1, -1};
  p.edge_patch = {0, 0, 0, 0};
  p.patches = {{"wall", PatchKind::Wall}};
  m = extrude_2d(p, 2);
  g = compute_geometry(m);
  CHECK(g.cell_volume[0] == doctest::Approx(3));
  CHECK(g.cell_centroid[0].x == doctest::Approx(7. / 9));
  CHECK(g.cell_centroid[0].y == doctest::Approx(4. / 9));
  for (auto id : {1, 2}) {
    auto fc = g.face_centroid[m.patches[id].start];
    CHECK(fc.x == doctest::Approx(7. / 9));
    CHECK(fc.y == doctest::Approx(4. / 9));
  }
}
TEST_CASE(
    "RCM reduces a scrambled chain bandwidth and handles disconnected cells") {
  auto m = make_box_hex(32, 1, 1, {0, 0, 0}, {32, 1, 1});
  for (auto& c : m.owner) c = (c * 17) % 32;
  for (auto& c : m.neighbour) c = (c * 17) % 32;
  for (Index f = 0; f < m.n_internal_faces(); ++f)
    if (m.owner[f] > m.neighbour[f]) {
      std::swap(m.owner[f], m.neighbour[f]);
      std::reverse(m.face_nodes.begin() + m.face_offsets[f],
                   m.face_nodes.begin() + m.face_offsets[f + 1]);
    }
  auto bandwidth = [](const PolyMesh& mesh) {
    Index maximum = 0;
    for (Index f = 0; f < mesh.n_internal_faces(); ++f)
      maximum = std::max(maximum, mesh.neighbour[f] - mesh.owner[f]);
    return maximum;
  };
  CHECK(bandwidth(m) > 1);
  renumber_cells_rcm(m);
  CHECK(bandwidth(m) == 1);
  validate(m);
  std::vector<Vec3d> points = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
                               {3, 0, 0}, {4, 0, 0}, {3, 1, 0}, {3, 0, 1}};
  std::vector<kutta::detail::CellFaces> cells = {
      {{0, 1, 2}, {0, 3, 1}, {1, 3, 2}, {2, 3, 0}},
      {{4, 5, 6}, {4, 7, 5}, {5, 7, 6}, {6, 7, 4}}};
  auto disconnected = from_cells(points, cells, {{"wall", PatchKind::Wall}},
                                 [](const auto&) { return 0; });
  auto perm = renumber_cells_rcm(disconnected);
  CHECK(perm.cells_new_to_old == std::vector<Index>{0, 1});
  validate(disconnected);
}
TEST_CASE("oversized box dimensions fail before arithmetic overflows") {
  auto n = std::numeric_limits<int>::max();
  CHECK_THROWS_AS(make_box_hex(n, n, n, {0, 0, 0}, {1, 1, 1}), Error);
  CHECK_THROWS_AS(make_box_tet(n, n, n, {0, 0, 0}, {1, 1, 1}), Error);
}

TEST_CASE("signed face moments handle a dart quad outside the vertex-average kernel") {
  std::vector<Vec3d> points{{-1, 0, 0}, {0, 1.5, 0}, {1, 0, 0}, {0, 3, 0},
                            {-1, 0, 1}, {0, 1.5, 1}, {1, 0, 1}, {0, 3, 1}};
  auto m = pack(points, 1,
                {{{3, 2, 1, 0}, 0, -1, 0}, {{4, 5, 6, 7}, 0, -1, 0},
                 {{0, 1, 5, 4}, 0, -1, 0}, {{1, 2, 6, 5}, 0, -1, 0},
                 {{2, 3, 7, 6}, 0, -1, 0}, {{3, 0, 4, 7}, 0, -1, 0}},
                {{"shell", PatchKind::Wall}});
  auto g = compute_geometry(m);
  for (Index f : {0, 1}) {
    CHECK(g.face_centroid[f].x == doctest::Approx(0));
    CHECK(g.face_centroid[f].y == doctest::Approx(1.5));
  }
  CHECK(g.cell_volume[0] == doctest::Approx(1.5));
}
