#include "fixture.hpp"
#include <doctest/doctest.h>
#include <numbers>
using namespace post_fixture;

TEST_CASE("Uniform-flow adaptive streamlines stop on a boundary and honour "
          "maximum length") {
  auto s = uniform(4);
  post::Volume volume(s);
  post::Trace spec;
  spec.seeds = post::Rake{{.1, .2, .5}, {.1, .8, .5}, 4};
  spec.step_length = .2;
  spec.max_length = 2;
  spec.max_steps = 200;
  auto result = post::trace(volume, spec);
  REQUIRE(result.lines.size() == 4);
  REQUIRE(result.geometry.line_offsets.size() == 5);
  for (std::size_t i = 0; i < 4; ++i) {
    const auto &report = object(result.lines[i]);
    CHECK(std::get<std::string>(report.at("stop_reason").value) == "boundary");
    const auto last =
        result.geometry.indices[result.geometry.line_offsets[i + 1] - 1];
    const auto p = point(result.geometry, last);
    CHECK(p.x == doctest::Approx(1).epsilon(1e-6));
    CHECK(p.y == doctest::Approx(.2 + .2 * i).epsilon(1e-6));
    CHECK(p.z == doctest::Approx(.5));
    CHECK(number(report.at("length")) == doctest::Approx(.9).epsilon(1e-6));
  }
  spec.seeds = post::Rake{{.1, .5, .5}, {.1, .5, .5}, 1};
  spec.max_length = .3;
  result = post::trace(volume, spec);
  CHECK(std::get<std::string>(
            object(result.lines[0]).at("stop_reason").value) == "max_length");
  CHECK(point(result.geometry, result.geometry.positions.size() / 3 - 1).x ==
        doctest::Approx(.4).epsilon(1e-6));
}
TEST_CASE("RK4 reconstructed solid-body rotation follows a circular "
          "closed-form path") {
  auto s = snapshot(
      8, 8, 2, {-1, -1, -.1}, {1, 1, .1},
      [](Vec3d p) { return Vec3d{-2 * p.y, 2 * p.x, 0}; },
      [](Vec3d) { return 0; });
  post::Volume volume(s);
  post::Trace spec;
  spec.seeds = post::Rake{{.5, 0, 0}, {.5, 0, 0}, 1};
  spec.step_length = .02;
  spec.max_time = std::numbers::pi;
  spec.max_length = 10;
  spec.relative_tolerance = 1e-8;
  spec.absolute_tolerance = 1e-10;
  const auto result = post::trace(volume, spec);
  REQUIRE(result.lines.size() == 1);
  CHECK(std::get<std::string>(
            object(result.lines[0]).at("stop_reason").value) == "max_time");
  const auto end =
      point(result.geometry, result.geometry.positions.size() / 3 - 1);
  CHECK(norm(end - Vec3d{.5, 0, 0}) < 2e-5);
  for (std::size_t i = 0; i < result.geometry.positions.size() / 3; ++i) {
    const auto p = point(result.geometry, i);
    CHECK(std::hypot(p.x, p.y) == doctest::Approx(.5).epsilon(2e-5));
  }
}
TEST_CASE("Time-dependent pathlines interpolate saved frame velocity with "
          "physical time") {
  auto first = snapshot(
      4, 2, 2, {0, 0, 0}, {10, 1, 1}, [](Vec3d) { return Vec3d{1, 0, 0}; },
      [](Vec3d) { return 0; });
  auto second = first;
  for (auto *fields : {&second.cells, &second.boundary})
    for (auto &f : *fields)
      if (f.name == "velocity")
        for (std::size_t i = 0; i < f.values.size(); i += 3)
          f.values[i] = 3;
  post::Volume a(first), b(second, &a);
  std::array<const post::Volume *, 2> frames{&a, &b};
  post::Trace spec;
  spec.pathline = true;
  spec.frames = {{0, {}}, {2, {}}};
  spec.seeds = post::Rake{{.5, .5, .5}, {.5, .5, .5}, 1};
  spec.max_time = 1;
  spec.max_length = 5;
  spec.step_length = .2;
  const auto result = post::trace(a, spec, frames);
  CHECK(point(result.geometry, result.geometry.positions.size() / 3 - 1).x ==
        doctest::Approx(2).epsilon(2e-6));
  CHECK(number(object(result.lines[0]).at("time")) == doctest::Approx(1));
}
TEST_CASE("Plane and deterministic area-distributed surface seeds respect "
          "their specifications") {
  auto s = uniform(2);
  post::Volume volume(s);
  const auto grid = post::seeds(
      volume, post::PlaneGrid{{.1, .2, .2}, {0, .6, 0}, {0, 0, .6}, 3, 2});
  CHECK(grid.size() == 6);
  CHECK(grid[0].x == doctest::Approx(.1));
  CHECK(grid.back().z == doctest::Approx(.8));
  for (bool random : {false, true}) {
    post::SurfaceSeeds spec{{"xmax"}, 100, random, 42, .01};
    const auto a = post::seeds(volume, spec), b = post::seeds(volume, spec);
    REQUIRE(a.size() == 100);
    for (std::size_t i = 0; i < a.size(); ++i) {
      CHECK(a[i].x == doctest::Approx(.99));
      CHECK(a[i].y >= 0);
      CHECK(a[i].y <= 1);
      CHECK(a[i].z >= 0);
      CHECK(a[i].z <= 1);
      CHECK(a[i].x == b[i].x);
      CHECK(a[i].y == b[i].y);
      CHECK(a[i].z == b[i].z);
    }
  }
}

TEST_CASE("Pathline steps honour saved time knots in forward and backward "
          "integration") {
  auto first = snapshot(
      4, 2, 2, {0, 0, 0}, {10, 1, 1}, [](Vec3d) { return Vec3d{1, 0, 0}; },
      [](Vec3d) { return 0; });
  auto peak = first;
  for (auto *fields : {&peak.cells, &peak.boundary})
    for (auto &f : *fields)
      if (f.name == "velocity")
        for (std::size_t i = 0; i < f.values.size(); i += 3)
          f.values[i] = 3;
  post::Volume a(first), b(peak, &a);
  std::array<const post::Volume *, 3> frames{&a, &b, &a};
  post::Trace spec;
  spec.pathline = true;
  spec.frames = {{0, {}}, {1, {}}, {2, {}}};
  spec.max_time = 2;
  spec.max_length = 8;
  spec.step_length = 8;
  for (int direction : {-1, 1}) {
    spec.direction = direction;
    spec.start_time = direction > 0 ? 0 : 2;
    spec.seeds = post::Rake{{5, .5, .5}, {5, .5, .5}, 1};
    const auto result = post::trace(a, spec, frames);
    const auto last = result.geometry.positions.size() / 3 - 1;
    CHECK(point(result.geometry, last).x ==
          doctest::Approx(5 + direction * 4).epsilon(2e-6));
    CHECK(number(object(result.lines[0]).at("time")) ==
          doctest::Approx(direction > 0 ? 2 : 0).epsilon(1e-8));
  }
}

TEST_CASE("Unsteady pathlines can leave rest and retain physical time at a "
          "length limit") {
  auto first = snapshot(
      2, 2, 2, {0, 0, 0}, {2, 1, 1}, [](Vec3d) { return Vec3d{}; },
      [](Vec3d) { return 0; });
  auto second = first;
  for (auto *fields : {&second.cells, &second.boundary})
    for (auto &f : *fields)
      if (f.name == "velocity")
        for (std::size_t i = 0; i < f.values.size(); i += 3)
          f.values[i] = 1;
  post::Volume a(first), b(second, &a);
  std::array<const post::Volume *, 2> frames{&a, &b};
  post::Trace spec;
  spec.pathline = true;
  spec.frames = {{0, {}}, {1, {}}};
  spec.max_time = 1;
  spec.max_length = .125;
  spec.step_length = 1;
  spec.absolute_tolerance = 1e-9;
  spec.seeds = post::Rake{{.5, .5, .5}, {.5, .5, .5}, 1};
  const auto result = post::trace(a, spec, frames);
  const auto last = result.geometry.positions.size() / 3 - 1;
  CHECK(point(result.geometry, last).x == doctest::Approx(.625).epsilon(1e-6));
  CHECK(number(object(result.lines[0]).at("time")) ==
        doctest::Approx(.5).epsilon(1e-6));
  CHECK(std::get<std::string>(
            object(result.lines[0]).at("stop_reason").value) == "max_length");
}
