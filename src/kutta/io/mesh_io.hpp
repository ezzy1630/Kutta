#pragma once

// Mesh readers, writers and generators. The solver is 3D only: 2D meshes are
// extruded one cell thick, and the two extrusion faces become Symmetry patches
// named "front" and "back".

#include <array>
#include <filesystem>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "kutta/mesh/mesh.hpp"

namespace kutta {

// Fluent legacy mesh/case (.msh, .cas), ASCII or binary sections; 2D or 3D.
PolyMesh read_fluent_msh(const std::filesystem::path& path,
                         double extrude_2d = 1.0);

// Optional source metadata for periodic coupling; face ids are Kutta ids after
// reordering/extrusion. read_fluent_msh() returns only this result's mesh.
struct PeriodicFacePair {
  Index face = -1, shadow = -1;
  int zone = -1, shadow_zone = -1;
};
struct FluentMeshData {
  PolyMesh mesh;
  std::vector<PeriodicFacePair> periodic_faces;
};
FluentMeshData read_fluent_msh_with_metadata(const std::filesystem::path& path,
                                             double extrude_2d = 1.0);

// Fluent Common Fluids Format HDF5 mesh (.msh.h5, or the mesh inside .cas.h5).
PolyMesh read_fluent_cff(const std::filesystem::path& path,
                         double extrude_2d = 1.0);

// Formatted ASCII Plot3D multi-block grid, 2D or 3D (one k-plane is extruded).
// Gunzip .gz inputs before calling. Auto-loads a sibling <stem>.nmf when
// present, applying its BC ranges and one-to-one connections. Without a map,
// blockN_imin/imax/jmin/jmax[/kmin/kmax] remain independent boundaries.
// Extruded blocks share one front and one back Symmetry patch.
PolyMesh read_plot3d(const std::filesystem::path& path,
                     double extrude_2d = 1.0);
// An empty explicit map path disables auto-discovery.
PolyMesh read_plot3d(const std::filesystem::path& path,
                     const std::filesystem::path& nmf_path,
                     double extrude_2d = 1.0);

// Dispatch on file extension.
PolyMesh read_mesh(const std::filesystem::path& path);

struct CellField {
  std::string name;
  int components = 1;              // 1 (scalar) or 3 (vector, interleaved xyz)
  std::span<const double> values;  // n_cells * components
};

// VTK XML unstructured grid (polyhedra), readable by ParaView.
// Optional provenance is preserved as an XML header comment (no double hyphens).
void write_vtu(const std::filesystem::path& path, const PolyMesh& mesh,
               std::span<const CellField> fields = {}, std::string_view header_comment = {});

// Structured hex box [lo, hi] with nx*ny*nz cells. Patches, in order:
// xmin, xmax, ymin, ymax, zmin, zmax (all PatchKind::Wall; callers re-tag).
PolyMesh make_box_hex(int nx, int ny, int nz, Vec3d lo, Vec3d hi);

// Conforming Freudenthal split: six tets per hex, sharing the 000->111
// diagonal.
PolyMesh make_box_tet(int nx, int ny, int nz, Vec3d lo, Vec3d hi);

// Planar xy mesh. Directed edges leave owner on the left. neighbour and
// edge_patch have one entry per edge: -1 denotes boundary / internal
// respectively. patches carry metadata only; extrusion fills their start/size
// ranges.
struct PlanarMesh {
  std::vector<Vec3d> points;
  std::vector<std::array<Index, 2>> edges;
  std::vector<Index> owner, neighbour, edge_patch;
  Index n_cells = 0;
  std::vector<Patch> patches;
  std::vector<CellZone> cell_zones;
};
PolyMesh extrude_2d(const PlanarMesh& mesh, double thickness = 1.0);

// Move matching faces into a new patch immediately after the original. The
// original retains nonmatches. Metadata is inherited; returns the new patch id.
// Face ids change; rebuild any face-indexed fields after splitting.
Index split_patch(PolyMesh& mesh, Index patch, const std::string& name,
                  const std::function<bool(Vec3d)>& predicate);

}  // namespace kutta
