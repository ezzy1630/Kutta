#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

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

#define KUTTA_REQUIRE(cond, msg)                                              \
  do {                                                                        \
    if (!(cond)) throw ::kutta::Error(std::string(__FILE__) + ":" +           \
                                      std::to_string(__LINE__) + ": " + (msg)); \
  } while (0)

}  // namespace kutta
