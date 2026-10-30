#include "../io/test_support.hpp"
#include "fixture.hpp"
#include <cstdlib>
#include <doctest/doctest.h>

namespace {
std::string quote(const std::filesystem::path &path) {
#ifdef _WIN32
  return '"' + std::filesystem::path(path).make_preferred().string() + '"';
#else
  std::string out = "'";
  for (char c : path.string())
    out += c == '\'' ? "'\\''" : std::string(1, c);
  return out + "'";
#endif
}
} // namespace
TEST_CASE("Post CLI publishes independently decoded compact traces iso "
          "surfaces probes and reports") {
  MeshTemp t;
  auto s = post_fixture::uniform(3);
  kutta::view::write_snapshot(t.path("view-input.bin"), s);
  mesh_text(t.path("summary.json"), "{\"v\":1}");
  std::filesystem::create_directory(t.path("frame1"));
  auto later = s;
  for (auto *fields : {&later.cells, &later.boundary})
    for (auto &f : *fields)
      if (f.name == "velocity")
        for (std::size_t i = 0; i < f.values.size(); i += 3)
          f.values[i] = 3;
  kutta::view::write_snapshot(t.path("frame1/view-input.bin"), later);
  mesh_text(t.path("frame1/summary.json"), "{\"v\":1}");
  auto command = quote(KUTTA_PYTHON_PATH) + " " +
                 quote(std::filesystem::path(KUTTA_SOURCE_DIR) /
                       "tests/post/check_post.py") +
                 " " + quote(KUTTA_CLI_PATH) + " " + quote(t.directory);
#ifdef _WIN32
  command = '"' + command + '"';
#endif
  CHECK(std::system(command.c_str()) == 0);
}
