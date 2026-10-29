#pragma once
#include <cstdint>
#include <vector>

#include "kutta/mesh/mesh.hpp"

namespace kutta {
struct CellFaces {
  std::vector<Index> offsets;
  std::vector<Index> faces;
  std::vector<std::int8_t>
      signs;  // +1 owner, -1 neighbour; outward area = sign*Sf
};
struct CellNeighbours {
  std::vector<Index> offsets;
  std::vector<Index>
      cells;  // sorted and unique within each row; no self entries
};
// Inputs must satisfy PolyMesh's invariants.
CellFaces build_cell_faces(const PolyMesh& mesh);
CellNeighbours build_cell_neighbours(const PolyMesh& mesh);
struct MeshPermutation {
  std::vector<Index> cells_new_to_old;
  std::vector<Index> faces_new_to_old;
};
// Deterministic reverse Cuthill-McKee, including disconnected components.
// Reorders internal faces by owner/neighbour and boundaries within their patch;
// reverses winding on owner swaps and updates cell zones. Point ids stay fixed.
// Use returned new->old maps to permute existing cell/face fields. An oriented
// face field also needs a sign flip wherever the old owner becomes neighbour.
MeshPermutation renumber_cells_rcm(PolyMesh& mesh);
}  // namespace kutta
