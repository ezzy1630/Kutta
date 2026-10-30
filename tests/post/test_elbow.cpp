#include "../fluent_examples.hpp"
#include "../io/test_support.hpp"
#include <cstdlib>
#include <doctest/doctest.h>
namespace {
std::string quote_elbow(const std::filesystem::path &path) {
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
TEST_CASE("Mixing elbow post-processing smoke exercises every requested "
          "operation through the CLI") {
#if KUTTA_HAS_HDF5
  const auto fixture = fluent_examples() / "mixing_elbow";
  if (!std::filesystem::exists(fixture / "mixing_elbow.cas.h5") ||
      !std::filesystem::exists(fixture / "mixing_elbow.dat.h5")) {
    MESSAGE("SKIP: optional Fluent mixing elbow case/data absent");
    return;
  }
  MeshTemp t;
  auto command = quote_elbow(KUTTA_PYTHON_PATH) + " " +
                 quote_elbow(std::filesystem::path(KUTTA_SOURCE_DIR) /
                             "tests/post/check_elbow.py") +
                 " " + quote_elbow(KUTTA_CLI_PATH) + " " +
                 quote_elbow(fixture) + " " + quote_elbow(t.directory);
#ifdef _WIN32
  command = '"' + command + '"';
#endif
  CHECK(std::system(command.c_str()) == 0);
#else
  MESSAGE("SKIP: HDF5 disabled");
#endif
}
