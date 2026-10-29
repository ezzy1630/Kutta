#pragma once

#include "kutta/fv/device_mesh.hpp"
#include "kutta/linalg/csr.hpp"

namespace kutta::fv {

DeviceVec3 vector_field(Device& dev, std::size_t n, float value = 0);

// Geometry-only work, built once. Incidence arrays follow cell_face_offsets.
// Negative adjacent ids encode boundary faces as -(boundary_index + 1).
struct Stencil {
  CsrPattern pattern;
  Buffer<Index> adjacent, matrix_entry;
  DeviceVec3 lsq_delta, face_delta;
  Buffer<float> lsq_inverse;  // six SoA planes: xx, xy, xz, yy, yz, zz
  Buffer<float> diffusion;
  DeviceVec3 nonorthogonal;
};
Stencil make_stencil(Device& dev, const PolyMesh& mesh, const MeshGeometry& geometry);

// Boundary entries are actual face values, including extrapolated Neumann data.
// LSQ uses inverse-distance-squared weights and includes boundary faces.
void gradient(Device& dev, const DeviceMesh& mesh, const Stencil& stencil,
              const Buffer<float>& cell, const Buffer<float>& boundary, DeviceVec3& result);
void limit_gradient(Device& dev, const DeviceMesh& mesh, const Stencil& stencil,
                    const Buffer<float>& cell, const Buffer<float>& boundary,
                    const DeviceVec3& gradient, Buffer<float>& limiter);

}  // namespace kutta::fv
