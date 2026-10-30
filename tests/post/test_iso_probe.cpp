#include "fixture.hpp"
#include <doctest/doctest.h>
using namespace post_fixture;

TEST_CASE("Marching tetrahedra exports an affine iso plane with correct area "
          "fields and winding") {
  auto s = uniform(4);
  post::Volume volume(s);
  post::Iso iso{"plane", "phi", .43, -1, {"pressure", "velocity"}};
  const auto g = post::isosurface(volume, iso);
  REQUIRE(!g.indices.empty());
  double area = 0;
  for (std::size_t i = 0; i < g.positions.size() / 3; ++i) {
    const auto p = point(g, i);
    CHECK(p.x == doctest::Approx(.43).epsilon(2e-6));
    CHECK(g.fields[0].values[i] ==
          doctest::Approx(2 + 3 * p.x - 2 * p.y + p.z).epsilon(3e-6));
  }
  for (std::size_t i = 0; i < g.indices.size(); i += 3) {
    const auto a = point(g, g.indices[i]), b = point(g, g.indices[i + 1]),
               c = point(g, g.indices[i + 2]);
    const auto n = cross(b - a, c - a);
    CHECK(n.x > 0);
    area += norm(n) * .5;
  }
  CHECK(area == doctest::Approx(1).epsilon(2e-6));
  iso.value = 2;
  CHECK(post::isosurface(volume, iso).indices.empty());
}
TEST_CASE(
    "Line probes preserve affine values and mask samples outside the volume") {
  auto s = uniform(3);
  post::Volume volume(s);
  post::LineProbe spec{
      "line", {{-.1, .5, .5}, {1.1, .5, .5}, 13}, {"phi", "pressure"}};
  const auto result = post::probe(volume, spec);
  const auto &g = result.geometry;
  REQUIRE(g.positions.size() == 39);
  CHECK(g.fields.back().values.front() == 0);
  CHECK(g.fields.back().values.back() == 0);
  for (std::size_t i = 1; i < 12; ++i) {
    const auto p = point(g, i);
    CHECK(g.fields.back().values[i] == 1);
    CHECK(g.fields[0].values[i] == doctest::Approx(p.x).epsilon(2e-6));
    CHECK(g.fields[1].values[i] ==
          doctest::Approx(1.5 + 3 * p.x).epsilon(2e-6));
  }
}
TEST_CASE("Surface-plane section probes produce ordered closed Cp curves") {
  auto s = snapshot(
      2, 2, 2, {0, 0, 0}, {1, 1, 1}, [](Vec3d) { return Vec3d{1, 0, 0}; },
      [](Vec3d) { return 10; });
  // One surface containing the whole boundary, like a Fluent wall zone.
  const auto begin = s.mesh.n_internal_faces();
  s.mesh.patches = {
      {"body", PatchKind::Wall, "wall", 10, begin, s.mesh.n_faces() - begin}};
  post::Volume volume(s);
  post::SectionProbe spec{
      "section", {"body"}, {{0, 0, .5}, {0, 0, 1}}, {"pressure_coefficient"}};
  const auto result = post::probe(volume, spec);
  REQUIRE(result.lines.size() == 1);
  CHECK(std::get<bool>(object(result.lines[0]).at("closed").value));
  CHECK(number(object(result.lines[0]).at("length")) == doctest::Approx(4));
  for (float value : result.geometry.fields[0].values)
    CHECK(value == doctest::Approx(2.25));
  for (std::size_t i = 0; i < result.geometry.positions.size() / 3; ++i)
    CHECK(point(result.geometry, i).z == doctest::Approx(.5));
}

TEST_CASE("Nonconvex cell iso decomposition finds a positive kernel away from "
          "the Fluent centroid") {
  view::Snapshot s;
  const std::vector<Vec3d> base{{0, 0, 0},     {2, 0, 0},   {2, .25, 0},
                                {.25, .25, 0}, {.25, 2, 0}, {0, 2, 0}};
  auto points = base;
  for (auto p : base)
    points.push_back({p.x, p.y, 1});
  std::vector<kutta::detail::CellFaces> cells(1);
  cells[0] = {{5, 4, 3, 2, 1, 0}, {6, 7, 8, 9, 10, 11}};
  for (Index i = 0; i < 6; ++i) {
    const auto next = (i + 1) % 6;
    cells[0].push_back({i, next, next + 6, i + 6});
  }
  // Preserve the known outward winding. The generic convex-cell constructor
  // orients faces against a vertex average, which lies outside this L prism.
  std::vector<kutta::detail::Face> faces;
  for (const auto &ns : cells[0])
    faces.push_back({ns, 0, -1, 0});
  s.mesh = pack(points, 1, faces, {{"shell", PatchKind::Wall, "wall", 1}});
  const auto geometry = compute_geometry(s.mesh);
  s.cells = {{"phi", "m", 1, {.5f}}};
  view::Field boundary{"phi", "m", 1, {}};
  for (auto p : geometry.face_centroid)
    boundary.values.push_back(float(p.z));
  s.boundary.push_back(boundary);
  post::Volume volume(s);
  CHECK(geometry.cell_centroid[0].x > .25);
  CHECK(volume.decomposition_center(0).x < .25);
  CHECK(volume.decomposition_center(0).y < .25);
  CHECK(volume.locate({.1, 1, .5}) == 0);
  CHECK(volume.locate({1, 1, .5}) == -1);
  CHECK(volume.locate({.25, 1, .5}) == 0);
  const auto crossing = volume.walk({.1, 1, .5}, 0, {1, .1, .5});
  CHECK(crossing.boundary);
  CHECK(crossing.point.x == doctest::Approx(.25));
  const auto result =
      post::isosurface(volume, {"section", "phi", .5, -1, {"phi"}});
  double area = 0;
  for (std::size_t i = 0; i < result.indices.size(); i += 3) {
    const auto a = point(result, result.indices[i]),
               b = point(result, result.indices[i + 1]),
               c = point(result, result.indices[i + 2]);
    const auto n = cross(b - a, c - a);
    CHECK(n.z > 0);
    area += norm(n) * .5;
  }
  CHECK(area == doctest::Approx(.9375).epsilon(1e-6));
}
