#pragma once

// Mesh readers, writers and generators. The solver is 3D only: 2D meshes are
// extruded one cell thick, and the two extrusion faces become Symmetry patches
// named "front" and "back".

#include <filesystem>
#include <span>
#include <string>
#include <vector>

#include "kutta/mesh/mesh.hpp"

namespace kutta {

// Fluent legacy mesh/case (.msh, .cas), ASCII or binary sections; 2D or 3D.
PolyMesh read_fluent_msh(const std::filesystem::path& path, double extrude_2d = 1.0);

// Fluent Common Fluids Format HDF5 mesh (.msh.h5, or the mesh inside .cas.h5).
PolyMesh read_fluent_cff(const std::filesystem::path& path, double extrude_2d = 1.0);

// Plot3D multi-block structured grid (NASA TMR format), whole or 2D (one k-plane).
PolyMesh read_plot3d(const std::filesystem::path& path, double extrude_2d = 1.0);

// Dispatch on file extension.
PolyMesh read_mesh(const std::filesystem::path& path);

struct CellField {
  std::string name;
  int components = 1;       // 1 (scalar) or 3 (vector, interleaved xyz)
  std::span<const double> values;  // n_cells * components
};

// VTK XML unstructured grid (polyhedra), readable by ParaView.
void write_vtu(const std::filesystem::path& path, const PolyMesh& mesh,
               std::span<const CellField> fields = {});

// Structured hex box [lo, hi] with nx*ny*nz cells. Patches, in order:
// xmin, xmax, ymin, ymax, zmin, zmax (all PatchKind::Wall; callers re-tag).
PolyMesh make_box_hex(int nx, int ny, int nz, Vec3d lo, Vec3d hi);

}  // namespace kutta
