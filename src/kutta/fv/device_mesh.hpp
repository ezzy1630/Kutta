#pragma once

// Device-resident mesh for finite-volume kernels: the contract between the
// flow solver (fv/) and the turbulence model (fv/turbulence/).
//
// Precision rule (docs/ARCHITECTURE.md): everything is computed from
// MeshGeometry in double on the host, then stored as float *differences*
// (face-to-cell and cell-to-cell vectors), never absolute coordinates.
//
// Indexing:
// - Faces keep PolyMesh order: internal faces [0, n_internal_faces), then
//   boundary faces grouped by patch. Boundary face b = f - n_internal_faces
//   indexes the per-boundary-face arrays (boundary values, wall data).
// - Cell -> face lists come from build_cell_faces(); kernels gather over them
//   (no scatter, no atomics).

#include <string>
#include <vector>

#include "kutta/mesh/mesh.hpp"
#include "kutta/runtime/device.hpp"

namespace kutta {

struct DeviceVec3 {
  Buffer<float> x, y, z;
};

struct DevicePatch {
  std::string name;
  PatchKind kind = PatchKind::Other;
  Index start = 0;  // first face index (global face numbering)
  Index size = 0;
};

struct DeviceMesh {
  Index n_cells = 0;
  Index n_faces = 0;
  Index n_internal_faces = 0;
  Index n_boundary_faces() const { return n_faces - n_internal_faces; }

  // Faces (size n_faces unless noted).
  Buffer<std::int32_t> owner;
  Buffer<std::int32_t> neighbour;  // size n_internal_faces
  DeviceVec3 area;                 // area vector S_f, out of owner
  DeviceVec3 owner_to_face;        // x_f - C_owner
  DeviceVec3 neighbour_to_face;    // x_f - C_neighbour, size n_internal_faces
  DeviceVec3 owner_to_neighbour;   // C_neighbour - C_owner, size n_internal_faces;
                                   // for boundary faces use owner_to_face
  Buffer<float> weight;            // owner interpolation weight w, phi_f = w*phi_P + (1-w)*phi_N,
                                   // size n_internal_faces

  // Cells (size n_cells unless noted).
  Buffer<float> volume;
  Buffer<std::int32_t> cell_face_offsets;  // size n_cells+1
  Buffer<std::int32_t> cell_faces;         // CSR into faces
  Buffer<float> cell_face_sign;            // +1 owner, -1 neighbour (float: kernels multiply)

  // Wall data. wall_distance: exact centroid-to-nearest-wall-surface distance
  // (matches Fluent SV_WALL_DIST). wall_normal_distance: per boundary face, the
  // owner centroid's distance to that face's plane (Fluent's wall-function dy);
  // zero on non-wall faces.
  Buffer<float> wall_distance;              // size n_cells
  Buffer<float> wall_normal_distance;       // size n_boundary_faces

  std::vector<DevicePatch> patches;  // host copy of patch ranges
};

// Builds and uploads all of the above. wall_distance is computed exactly on the
// host in double (BVH over faces of PatchKind::Wall patches).
DeviceMesh upload_mesh(Device& dev, const PolyMesh& mesh, const MeshGeometry& geometry);

}  // namespace kutta
