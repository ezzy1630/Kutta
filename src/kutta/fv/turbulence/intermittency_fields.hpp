#pragma once
#include "kutta/fv/device_mesh.hpp"

namespace kutta {
// Per-cell gamma-equation coefficients and diagnostics of the intermittency
// transition model, S = source - sink*gamma per volume (see sst.hpp).
struct IntermittencyFields {
  Buffer<float> source, sink, face_diffusivity;
  Buffer<float> normal_derivative;  // dV/dy = n . grad(n . U)
  Buffer<float> re_theta_c, onset;  // Re_theta_c(Tu_L, lambda_thetaL) and F_onset
  DeviceVec3 wall_direction;        // unit n, from the nearest wall point
};
}  // namespace kutta
