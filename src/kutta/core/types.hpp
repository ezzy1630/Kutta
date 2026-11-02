#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

namespace kutta {

// Cell, face and node indices. 32-bit keeps device arrays small and matches the
// kernel-side `int`; 2^31 cells is far beyond any single-GPU mesh we target.
using Index = std::int32_t;

// Device field storage and arithmetic precision. Geometry is computed on the
// host in double and stored on the device as float *differences* (face-to-cell
// vectors, volumes), never as absolute coordinates, so float stays accurate.
using Real = float;

struct Vec3d {
  double x = 0, y = 0, z = 0;
};

class Error : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// Diagnostics keep their source location without publishing the build host's
// directory layout. Both separators matter for cross-compiled Windows builds.
constexpr std::string_view source_filename(std::string_view path) {
  const auto slash = path.find_last_of("/\\");
  return slash == std::string_view::npos ? path : path.substr(slash + 1);
}

#define KUTTA_REQUIRE(cond, msg)                                              \
  do {                                                                        \
    if (!(cond)) throw ::kutta::Error(std::string(::kutta::source_filename(__FILE__)) + ":" + \
                                      std::to_string(__LINE__) + ": " + (msg)); \
  } while (0)

}  // namespace kutta
