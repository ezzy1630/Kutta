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
  DeviceVec3 boundary_lsq_delta;  // lsq_delta's boundary-face rows, by boundary face
  Buffer<float> lsq_inverse;  // six SoA planes: xx, xy, xz, yy, yz, zz
  Buffer<float> diffusion;
  DeviceVec3 nonorthogonal;
};
Stencil make_stencil(Device& dev, const PolyMesh& mesh, const MeshGeometry& geometry);

// Fit rows per cell-face entry. Internal-face rows depend on geometry alone
// and equal Stencil::lsq_delta, so a fit built with a stencil uses the
// stencil's rows and stores only its boundary-face rows; fits differing only
// in boundary observations no longer each hold every row.
struct BoundaryLeastSquares {
  const DeviceVec3& rows() const { return stencil_rows ? *stencil_rows : delta; }
  DeviceVec3 delta;           // owned entry rows; empty when stencil_rows is set
  DeviceVec3 boundary_delta;  // rows of boundary faces, by boundary face
  Buffer<float> inverse;
  const DeviceVec3* stencil_rows = nullptr;  // that stencil must outlive the fit
};
enum class BoundaryObservation { Omitted, FaceValue, NormalFootValue };
// A projected cell value at symmetry is located at the normal foot, rather
// than at the tangentially displaced face centroid. Extrapolated data is omitted.
// Missing directions receive zero-difference rows in the null space only.
BoundaryLeastSquares boundary_least_squares(Device& dev,const PolyMesh& mesh,
    const MeshGeometry& geometry,std::span<const BoundaryObservation> observations,
    const Stencil* stencil=nullptr);
// Per boundary face, sum_j |r_b . M^-1 w_j d_j| over its cell's rows: how far
// that fit's extrapolation to the face centroid can amplify the neighbour
// differences. At most 1, the extrapolated value stays within the cell value
// plus or minus the largest difference.
std::vector<double> boundary_extrapolation_gain(const PolyMesh& mesh,const MeshGeometry& geometry,
    std::span<const BoundaryObservation> observations);
// Fluent's Green-Gauss cell-based gradient (Theory Guide 2024 R2 23.3.3.1,
// Eq. 23-30, with the arithmetic face mean of Eq. 23-31) in the same linear
// form: since a closed cell's area vectors sum to zero,
//   (1/V) sum_f phi_f A_f = (1/V) sum_f (phi_f - phi_c) A_f,
// so the rows are A_f/(2V) at internal faces, A_f/V at boundary faces (the
// actual boundary value is the face value) and the inverse is the identity.
// Every boundary face is observed. The operator suits any consumer of a fit,
// including the Coupled block system.
BoundaryLeastSquares green_gauss_cell(Device& dev,const PolyMesh& mesh,const MeshGeometry& geometry);
void gradient(Device& dev,const DeviceMesh& mesh,const Stencil& stencil,
    const BoundaryLeastSquares& lsq,const Buffer<float>& cell,
    const Buffer<float>& boundary,DeviceVec3& result);
void vector_gradient(Device& dev,const DeviceMesh& mesh,const Stencil& stencil,
    const BoundaryLeastSquares& lsq,const DeviceVec3& cell,const DeviceVec3& boundary,
    DeviceVec3& u,DeviceVec3& v,DeviceVec3& w);

// Face-centroid rows (1/|r_f|^2) fit the difference between the two side
// reconstructions. This derivative correction is separate from convection.
struct FaceJumpFit {
  DeviceVec3 delta;  // per cell-face entry
  Buffer<float> inverse;
};
FaceJumpFit face_jump_least_squares(Device& dev,const PolyMesh& mesh,
    const MeshGeometry& geometry);
void jump_corrected_gradient(Device& dev,const DeviceMesh& mesh,const Stencil& stencil,
    const FaceJumpFit& fit,const Buffer<float>& cell,const Buffer<float>& boundary,
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
