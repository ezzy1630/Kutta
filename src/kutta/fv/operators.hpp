#pragma once

#include <span>

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

struct BoundaryLeastSquares {
  DeviceVec3 delta;
  Buffer<float> inverse;
};
enum class BoundaryObservation { Omitted, FaceValue, NormalFootValue };
// A projected cell value at symmetry is located at the normal foot, rather
// than at the tangentially displaced face centroid. Extrapolated data is omitted.
// Missing directions receive zero-difference rows in the null space only.
BoundaryLeastSquares boundary_least_squares(Device& dev,const PolyMesh& mesh,
    const MeshGeometry& geometry,std::span<const BoundaryObservation> observations);
void gradient(Device& dev,const DeviceMesh& mesh,const Stencil& stencil,
    const BoundaryLeastSquares& lsq,const Buffer<float>& cell,
    const Buffer<float>& boundary,DeviceVec3& result);
void vector_gradient(Device& dev,const DeviceMesh& mesh,const Stencil& stencil,
    const BoundaryLeastSquares& lsq,const DeviceVec3& cell,const DeviceVec3& boundary,
    DeviceVec3& u,DeviceVec3& v,DeviceVec3& w);

// Face-centroid rows (1/|r_f|^2) fit the difference between the two side
// reconstructions. This derivative correction is separate from convection.
BoundaryLeastSquares face_jump_least_squares(Device& dev,const PolyMesh& mesh,
    const MeshGeometry& geometry);
void jump_corrected_gradient(Device& dev,const DeviceMesh& mesh,const Stencil& stencil,
    const BoundaryLeastSquares& fit,const Buffer<float>& cell,const Buffer<float>& boundary,
    const DeviceVec3& base,DeviceVec3& result);

// Boundary entries are actual face values, including extrapolated Neumann data.
// LSQ uses inverse-distance-squared weights and includes boundary faces.
void gradient(Device& dev, const DeviceMesh& mesh, const Stencil& stencil,
              const Buffer<float>& cell, const Buffer<float>& boundary, DeviceVec3& result);
void vector_gradient(Device& dev,const DeviceMesh& mesh,const Stencil& stencil,
                     const DeviceVec3& cell,const DeviceVec3& boundary,
                     DeviceVec3& u,DeviceVec3& v,DeviceVec3& w);
void limit_gradient(Device& dev, const DeviceMesh& mesh, const Stencil& stencil,
                    const Buffer<float>& cell, const Buffer<float>& boundary,
                    const DeviceVec3& gradient, Buffer<float>& limiter);

}  // namespace kutta::fv
