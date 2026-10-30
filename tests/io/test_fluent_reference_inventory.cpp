#include <doctest/doctest.h>
#include "../fluent_reference_files.hpp"
#include "kutta/core/sha256.hpp"

TEST_CASE("Fluent research downloads retain published byte identities" * doctest::test_suite("slow")) {
  const std::filesystem::path root=KUTTA_SOURCE_DIR "/validation/fluent_examples";
  for(const auto& file:fluent_reference_files()) {
    const auto path=root/file.path;INFO(file.path.string());
    if(!std::filesystem::exists(path)) {MESSAGE("SKIP: external research file absent (see SOURCE.md)");continue;}
    REQUIRE(std::filesystem::file_size(path)==file.bytes);
    CHECK(kutta::sha256_file(path)==file.sha256);
  }
}
