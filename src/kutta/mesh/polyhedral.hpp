#pragma once

#include "kutta/mesh/mesh.hpp"

namespace kutta {

struct PolyhedralConversion {
  PolyMesh mesh;
  Index tets = 0;         // input tetrahedra replaced
  Index polyhedra = 0;    // dual cells created, one per tet-mesh vertex
  Index kept = 0;         // input cells retained (prisms, hexes, pyramids, polyhedra)
  Index split_vertices = 0;  // vertices kept as hex quarters (see below)
};

// Fluent's tet-to-polyhedra conversion (mesh/polyhedra/convert-domain): the
// median dual of the tetrahedral region.
// Each vertex of a tetrahedron becomes one polyhedron made of the quarter of
// every tetrahedron around it, cut at edge midpoints, triangle centroids and
// tetrahedron centroids. One polygon face separates the cells of an edge's two
// vertices. Every other cell (prism layers, hexes) is kept; a triangle it
// shares with a tetrahedron, and every boundary triangle of a tetrahedron, is
// split into the three vertex quads. Volume, boundary patches, their kinds and
// the region's surface are preserved exactly. A dual cell that wraps around a
// sharp convex edge can have its centroid behind one of its faces, which
// validate() rejects; such a vertex instead keeps one hex quarter per
// tetrahedron (as a tet-to-hex split would). A mesh without tetrahedra is
// returned unchanged. Throws when tetrahedra of different cell zones share a
// vertex or touch an internal face zone.
PolyhedralConversion convert_tets_to_polyhedra(const PolyMesh& mesh);

}  // namespace kutta
