#include "fixture.hpp"
#include <doctest/doctest.h>
using namespace post_fixture;
namespace {
Vec3d vector(const Json &j) {
  const auto &a = array(j);
  return {number(a[0]), number(a[1]), number(a[2])};
}
} // namespace
TEST_CASE("Aerodynamic surface reports reproduce analytic area flow force "
          "moment and constrained pressure centre") {
  auto s = snapshot(
      1, 1, 1, {0, 0, 0}, {1, 1, 1}, [](Vec3d) { return Vec3d{2, 0, 0}; },
      [](Vec3d) { return 10; });
  post::Volume v(s);
  post::Integral spec;
  spec.zones = {"xmax"};
  spec.fields = {"pressure", "velocity"};
  spec.pressure_plane = view::Plane{{1, 0, 0}, {1, 0, 0}};
  const auto report = post::integrate(v, spec);
  const auto &total = object(object(report).at("total"));
  CHECK(number(total.at("area")) == doctest::Approx(1));
  CHECK(number(total.at("mass_flow")) == doctest::Approx(4));
  const auto &force = object(total.at("force"));
  CHECK(norm(vector(force.at("pressure")) - Vec3d{9, 0, 0}) < 1e-12);
  CHECK(norm(vector(force.at("viscous")) - Vec3d{0, 3, 0}) < 1e-12);
  const auto &moment = object(total.at("moment"));
  CHECK(norm(vector(moment.at("pressure")) - Vec3d{0, 4.5, -4.5}) < 1e-12);
  CHECK(norm(vector(moment.at("viscous")) - Vec3d{-1.5, 0, 3}) < 1e-12);
  const auto &coefficient =
      object(object(total.at("coefficients")).at("total"));
  CHECK(number(coefficient.at("drag")) == doctest::Approx(.75));
  CHECK(number(coefficient.at("lift")) == doctest::Approx(.25));
  const auto &center = object(total.at("centre_of_pressure"));
  CHECK(norm(vector(center.at("point")) - Vec3d{1, .5, .5}) < 1e-12);
  CHECK(std::get<bool>(center.at("zero_moment").value));
  CHECK(number(object(object(total.at("averages")).at("pressure"))
                   .at("area_weighted_average")) == doctest::Approx(10));
}
TEST_CASE("Mass-weighted averages use absolute face flux while mass flow "
          "retains sign") {
  auto s = snapshot(
      1, 2, 1, {0, 0, 0}, {1, 1, 1}, [](Vec3d) { return Vec3d{1, 0, 0}; },
      [](Vec3d) { return 0; });
  const auto &p = s.mesh.patches[1];
  const auto b = p.start - s.mesh.n_internal_faces();
  for (auto &f : s.boundary)
    if (f.name == "pressure") {
      f.values[b] = 10;
      f.values[b + 1] = 20;
    } else if (f.name == "mass_flux") {
      f.values[b] = -1;
      f.values[b + 1] = 3;
    }
  s.properties.pressure = 0;
  post::Volume v(s);
  post::Integral spec;
  spec.zones = {"xmax"};
  spec.pressure_plane = view::Plane{{1, 0, 0}, {1, 0, 0}};
  const auto result = post::integrate(v, spec);
  const auto &total = object(object(result).at("total"));
  const auto &average = object(object(total.at("averages")).at("pressure"));
  CHECK(number(total.at("mass_flow")) == doctest::Approx(2));
  CHECK(number(total.at("absolute_mass_flow")) == doctest::Approx(4));
  CHECK(number(average.at("area_weighted_average")) == doctest::Approx(15));
  CHECK(number(average.at("mass_flow_weighted_average")) ==
        doctest::Approx(17.5));
  CHECK(vector(object(total.at("centre_of_pressure")).at("point")).y ==
        doctest::Approx(7.0 / 12));
}
TEST_CASE("Surface reports return null for undefined zero-flow averages and "
          "zero-force pressure centres") {
  auto s = snapshot(
      1, 1, 1, {0, 0, 0}, {1, 1, 1}, [](Vec3d) { return Vec3d{}; },
      [](Vec3d) { return 1; });
  post::Volume v(s);
  post::Integral spec;
  spec.zones = {"xmax"};
  const auto result = post::integrate(v, spec);
  const auto &total = object(object(result).at("total"));
  CHECK(object(object(total.at("averages")).at("pressure"))
            .at("mass_flow_weighted_average") == Json(nullptr));
  CHECK(object(total.at("centre_of_pressure")).at("point") == Json(nullptr));
}
