#include "fixture.hpp"
#include <doctest/doctest.h>
using namespace post_fixture;

TEST_CASE("Post reconstruction is affine exact and walks real face neighbours "
          "to a boundary") {
  auto s = uniform(4);
  post::Volume v(s);
  const Vec3d p{.23, .37, .41};
  const auto cell = v.locate(p);
  REQUIRE(cell >= 0);
  CHECK(v.sample("pressure", cell, p).value[0] ==
        doctest::Approx(2 + 3 * p.x - 2 * p.y + p.z).epsilon(2e-6));
  CHECK(v.sample("velocity_x", cell, p).value[0] == doctest::Approx(1));
  const auto walked = v.walk(p, cell, {.93, .37, .41});
  CHECK_FALSE(walked.boundary);
  CHECK(walked.cell != cell);
  CHECK(v.contains(walked.cell, walked.point));
  const auto exit = v.walk(p, cell, {1.5, .37, .41});
  CHECK(exit.boundary);
  CHECK(exit.point.x == doctest::Approx(1));
  CHECK(v.locate({-1, .5, .5}) == -1);
  CHECK_THROWS(v.info("invented-field"));
  CHECK_THROWS(v.sample("wall_shear", cell, p));
}
TEST_CASE("Solid-body rotation gives analytic vorticity Q helicity and "
          "pressure coefficients") {
  auto s = snapshot(
      6, 6, 2, {-1, -1, -.1}, {1, 1, .1},
      [](Vec3d p) { return Vec3d{-2 * p.y, 2 * p.x, 1}; },
      [](Vec3d) { return 5; });
  post::Volume v(s);
  const Vec3d p{.3, .2, .01};
  const auto cell = v.locate(p);
  REQUIRE(cell >= 0);
  const auto omega = v.sample("vorticity", cell, p);
  CHECK(omega.value[0] == doctest::Approx(0).epsilon(1e-5));
  CHECK(omega.value[1] == doctest::Approx(0).epsilon(1e-5));
  CHECK(omega.value[2] == doctest::Approx(4).epsilon(1e-5));
  CHECK(v.sample("q_criterion", cell, p).value[0] ==
        doctest::Approx(4).epsilon(1e-5));
  CHECK(v.sample("helicity", cell, p).value[0] ==
        doctest::Approx(4).epsilon(1e-5));
  CHECK(v.sample("pressure_coefficient", cell, p).value[0] ==
        doctest::Approx(1));
  CHECK(v.sample("total_pressure_coefficient", cell, p).value[0] ==
        doctest::Approx((5 + (.16 + .36 + 1) - 1) / 4).epsilon(1e-5));
  CHECK(v.boundary("skin_friction", s.mesh.patches[0].start, {}).value[0] ==
        doctest::Approx(.75));
  CHECK(v.boundary("y_plus", s.mesh.patches[0].start, {}).value[0] ==
        doctest::Approx(7));
}
TEST_CASE("Post reference overrides recompute existing coefficients and "
          "preserve solution fields") {
  auto s = uniform(2);
  const auto nc = s.mesh.n_cells,
             nb = s.mesh.n_faces() - s.mesh.n_internal_faces();
  s.cells.push_back(
      {"pressure_coefficient", "1", 1, std::vector<float>(nc, 123)});
  s.boundary.push_back(
      {"pressure_coefficient", "1", 1, std::vector<float>(nb, 123)});
  const auto pressure = s.cells[1].values;
  post::Request r;
  r.reference = {{"pressure", 3}, {"reference_density", 4}, {"velocity", 4}};
  r.fields = {"pressure_coefficient"};
  post::enrich(s, r);
  CHECK(s.cells[1].values == pressure);
  post::Volume v(s);
  for (Index f = s.mesh.n_internal_faces(); f < s.mesh.n_faces(); ++f) {
    const auto p = v.boundary("pressure", f, {}).value[0];
    CHECK(v.boundary("pressure_coefficient", f, {}).value[0] ==
          doctest::Approx((p - 3) / 32).epsilon(1e-6));
  }
}

TEST_CASE("Nonlinear Fluent deltas subtract independently reconstructed "
          "diagnostics") {
  auto s = snapshot(
      4, 4, 2, {-1, -1, -.1}, {1, 1, .1},
      [](Vec3d p) { return Vec3d{-2 * p.y, 2 * p.x, 1}; },
      [](Vec3d) { return 5; });
  for (auto *fields : {&s.cells, &s.boundary}) {
    auto u = (*fields)[0], p = (*fields)[1];
    u.name = "velocity@fluent";
    p.name = "pressure@fluent";
    for (std::size_t i = 0; i < u.values.size(); i += 3) {
      u.values[i] *= .5f;
      u.values[i + 1] *= .5f;
    }
    for (auto &value : p.values)
      value = 3;
    fields->push_back(std::move(u));
    fields->push_back(std::move(p));
  }
  post::Volume volume(s);
  const Vec3d point{.3, .2, 0};
  const auto cell = volume.locate(point);
  REQUIRE(cell >= 0);
  CHECK(volume.sample("q_criterion@fluent", cell, point).value[0] ==
        doctest::Approx(1).epsilon(1e-5));
  CHECK(volume.sample("q_criterion@delta", cell, point).value[0] ==
        doctest::Approx(3).epsilon(1e-5));
  CHECK(
      volume.sample("total_pressure_coefficient@delta", cell, point).value[0] ==
      doctest::Approx(.5975).epsilon(1e-5));
  post::Request request;
  request.fields = {"q_criterion", "q_criterion@fluent", "q_criterion@delta"};
  post::enrich(s, request);
  view::append_deltas(s.cells);
  CHECK(std::count_if(s.cells.begin(), s.cells.end(), [](const auto &f) {
          return f.name == "q_criterion@delta";
        }) == 1);
}
