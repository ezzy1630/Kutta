#include <doctest/doctest.h>

#include <clocale>
#include <string>

#include "kutta/core/parse.hpp"
#include "kutta/core/types.hpp"

using kutta::parse_double;

TEST_CASE("Error locations expose a source filename without the build directory") {
  CHECK(kutta::source_filename("/build/checkout/src/file.cpp") == "file.cpp");
  CHECK(kutta::source_filename(R"(C:\build\checkout\src\file.cpp)") == "file.cpp");
  CHECK(kutta::source_filename("file.cpp") == "file.cpp");
  try {
    KUTTA_REQUIRE(false, "diagnostic detail");
    FAIL("the failed invariant must throw");
  } catch (const kutta::Error& e) {
    const std::string message(e.what());
    CHECK(message.starts_with("test_parse.cpp:"));
    CHECK(message.ends_with(": diagnostic detail"));
    CHECK(message.find(KUTTA_SOURCE_DIR) == std::string::npos);
  }
}

TEST_CASE("parse_double accepts whole decimal and exponent forms only") {
  CHECK(parse_double("1.5") == 1.5);
  CHECK(parse_double("-2.5e-3") == -2.5e-3);
  CHECK(parse_double("6.02E23") == 6.02e23);
  CHECK(parse_double("0") == 0.0);
  CHECK_FALSE(parse_double(""));
  CHECK_FALSE(parse_double("+1"));
  CHECK_FALSE(parse_double(" 1"));
  CHECK_FALSE(parse_double("1.5x"));
  CHECK_FALSE(parse_double("1,5"));
}

TEST_CASE("parse_double ignores the process locale") {
  // A comma-decimal locale must not change parsing; skip if not installed.
  const char* old = std::setlocale(LC_NUMERIC, nullptr);
  const std::string saved = old ? old : "C";
  if (std::setlocale(LC_NUMERIC, "de_DE.UTF-8") || std::setlocale(LC_NUMERIC, "de_DE")) {
    CHECK(parse_double("1.25") == 1.25);
    CHECK_FALSE(parse_double("1,25"));
  }
  std::setlocale(LC_NUMERIC, saved.c_str());
}
