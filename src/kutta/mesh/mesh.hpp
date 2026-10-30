#pragma once

// Host-side polyhedral mesh in face-based (owner/neighbour) form, the layout
// Fluent, OpenFOAM and CGNS NGON all reduce to. Everything here is double
// precision; device data is derived from it (see mesh/device_mesh.hpp, later).
//
// Invariants (checked by validate()):
// - Faces [0, n_internal_faces) are internal; boundary faces follow, grouped
//   contiguously by patch in patch order.
// - Patch names are unique.
// - Internal faces: owner[f] < neighbour[f].
// - Face nodes are ordered so the right-hand normal points out of owner[f]
//   (into neighbour[f] for internal faces, out of the domain on boundaries).
// - Every cell is closed: the sum of its outward face area vectors is ~0.

#include <string>
#include <vector>

#include "kutta/core/types.hpp"
#include "kutta/mesh/face_zones.hpp"

namespace kutta {

// Generic role of a boundary patch. The case setup assigns the actual boundary
// condition; this is only the hint the mesh file carried (Fluent zone type).
enum class PatchKind : std::uint8_t {
  Wall,
  VelocityInlet,
  PressureOutlet,
  PressureInlet,
  Symmetry,
  Farfield,  // pressure far-field
  Other,
};

struct Patch {
  std::string name;
  PatchKind kind = PatchKind::Other;
  std::string source_type;  // type string from the file, e.g. "velocity-inlet"
  int source_zone_id = -1;  // Fluent zone id or equivalent, -1 if none
  Index start = 0;          // first face index
  Index size = 0;           // number of faces
};

// Optional cell grouping (fluid/solid regions, MRF zones, porous zones).
struct CellZone {
  std::string name;
  std::string source_type;  // e.g. "fluid", "solid"
  int source_zone_id = -1;
  std::vector<Index> cells;
};

struct PolyMesh {
  std::vector<Vec3d> points;
  std::vector<Index> face_offsets{0};  // size n_faces()+1, CSR into face_nodes
  std::vector<Index> face_nodes;
  std::vector<Index> owner;      // size n_faces()
  std::vector<Index> neighbour;  // size n_internal_faces()
  Index n_cells = 0;
  std::vector<Patch> patches;
  std::vector<CellZone> cell_zones;
  // Interior BC identity and source orientation; optional on generated meshes.
  std::vector<InternalFaceZone> internal_face_zones;

  Index n_faces() const { return static_cast<Index>(owner.size()); }
  Index n_internal_faces() const { return static_cast<Index>(neighbour.size()); }
  Index n_points() const { return static_cast<Index>(points.size()); }
};

// Throws kutta::Error describing the first violated invariant.
void validate(const PolyMesh& mesh);

// Double-precision geometry in the Fluent/OpenFOAM convention. Faces are
// triangle fans about the node average, with centroids weighted by triangle
// area projected on the resultant normal (exact for planar simple polygons).
// Cells are pyramids from the mean of their face centroids to each flattened
// face, with the centroid 3/4 of the way to the face. Saved Fluent fields
// (wall distance, y+, wall omega) reproduce only with this convention.
struct MeshGeometry {
  std::vector<Vec3d> face_centroid;  // per face
  std::vector<Vec3d> face_area;      // area-weighted normal, points out of owner
  std::vector<Vec3d> cell_centroid;  // per cell
  std::vector<double> cell_volume;   // per cell, > 0
};

MeshGeometry compute_geometry(const PolyMesh& mesh);

}  // namespace kutta
