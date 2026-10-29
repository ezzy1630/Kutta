# Mesh import and topology

`mesh.hpp` is unchanged. The I/O additions and `mesh/topology.hpp` provide the
construction and gather operations without changing the solver's face-based
mesh contract. All geometry is computed in double precision on the host.

## Geometry and adjacency

Faces use a triangle fan about the arithmetic mean of their vertices. The area
vector is the sum of the triangle area vectors; the face centroid is weighted by
triangle area magnitude. Cell volume and first moments are integrated by signed
pyramids (tetrahedra) from a local estimated center to those same triangles.
Local coordinate differences avoid unnecessarily losing precision on translated
meshes. This is an original implementation of the geometric decomposition, not
copied from OpenFOAM or another GPL source.

`validate()` checks CSR sizes/ranges, finite coordinates, valid distinct face
nodes, owner/neighbour ordering, complete contiguous patch coverage, cell-zone
indices, positive volumes, outward winding, area closure, and paired oriented
cell edges. Closure tolerance is `1e-10 * sum(magnitude(face area))` per cell.
Faces and cells are expected to be simple and star-shaped about their estimated
centers; arbitrary self-intersecting geometry is not supported.

```cpp
auto faces = kutta::build_cell_faces(mesh);
for (auto i = faces.offsets[c]; i < faces.offsets[c + 1]; ++i) {
  auto f = faces.faces[i];
  auto sign = faces.signs[i]; // +1 owner, -1 neighbour
  // Gather using outward area sign * geometry.face_area[f].
}
auto neighbours = kutta::build_cell_neighbours(mesh); // sorted unique CSR
```

Call `renumber_cells_rcm(mesh)` before constructing solver fields. It supports
multiple disconnected components, updates cell zones, reverses faces on owner
swaps, sorts internal faces, and groups boundaries by their original patches.
Its returned maps are `cells_new_to_old` and `faces_new_to_old`. For existing
oriented face fields, also negate a value when
`cells_new_to_old[new_mesh.owner[new_f]] != old_mesh.owner[old_f]`.
The tests demonstrate graph bandwidth improvement on a scrambled chain; GPU
performance is not measured by these host tests.

`make_box_tet()` uses six conforming tetrahedra per structured hex. `extrude_2d()`
accepts a planar xy mesh with each directed edge leaving its owner on the left.
It adds one layer in positive z and appends `front` and `back` Symmetry patches.

## Plot3D

`read_plot3d()` reads formatted ASCII multi-block 2D `.p2dfmt` and 3D grids,
including Fortran D exponents and one-k-plane extrusion. Whole 3D grids and
single-block dimension headers are supported. Gzip, binary Plot3D, iblank and
function files are unsupported: callers must **gunzip** compressed grids first.
No zlib dependency is added.

Each block has `blockN_imin`, `blockN_imax`, `blockN_jmin`, `blockN_jmax` boundary
patches (and `kmin`/`kmax` for 3D). Extruded caps are named `front`/`back`.
Blocks are independent; no interfaces are inferred or stitched from coincident
coordinates. Case setup must apply block connectivity when needed.

NASA TMR boundaries can cover only part of a block side. The flat-plate wall,
for example, starts at x=0. Split the relevant patch before building any
face-indexed fields:

```cpp
auto plate = kutta::split_patch(mesh, jmin_patch, "plate",
                               [](kutta::Vec3d c) { return c.x >= 0; });
mesh.patches[plate].kind = kutta::PatchKind::Wall;
```

Matching faces move to a new patch immediately after the original, inheriting
source metadata; nonmatches remain in the original. Empty results are allowed.
The predicate receives the area-weighted face centroid.

## Fluent legacy

ASCII sections 10/12/13/18 and binary 2010/3010, 2012/3012, 2013/3013,
2018/3018 are supported, along with zone names/types in 39/45. The section-4
machine flag selects byte order (60 little endian; other flags big endian).
Binary payload lengths come from headers, never parenthesis scanning inside
raw bytes. Node coordinates may be floats or doubles; connectivity is 32-bit.

Grid ids are hexadecimal. Zone metadata ids in 39/45 are decimal. Mixed and
polygonal faces carry per-face node counts. File node winding points *toward*
c0 in 3D; the reader reverses it when c0 becomes the owner. A 2D edge leaves c0
on its left. Both zero-c0 boundaries and c0>c1 internal faces are tested.
Source zone ids, raw named boundary type strings, and cell-zone names are
retained. Numeric BC codes have canonical strings when named metadata is absent.
Unspecified cell material is retained as `fluent-cell-zone-type-N`, rather than
inferring fluid or solid. Named cell metadata, when present, is authoritative.

Use `read_fluent_msh_with_metadata()` to retain periodic face pairs. Pair ids
refer to the returned mesh after face sorting and extrusion. The simpler
`read_fluent_msh()` returns just the mesh. Neither reader installs a solver
periodic boundary condition. Adapted/hanging-node and nonconformal interfaces
are rejected explicitly; surface-only files are rejected.

## Fluent CFF / HDF5

CMake discovers HDF5 optionally and publishes `KUTTA_HAS_HDF5=1` or `0`. Disable
it explicitly with `-DKUTTA_ENABLE_HDF5=OFF`. The stub reader throws a clear
installation/configuration error when HDF5 is unavailable. HDF5's imported CMake
target carries the C library's platform-specific include/link definitions.

The reader uses the HDF5 C API and `/meshes/1` layout: dimension, node
zoneTopology/coords, cell zoneTopology, face zoneTopology, faces/nodes sections
with `nnodes`/`nodes`, and sectioned c0/c1 datasets with minId/maxId attributes.
It supports fixed semicolon-delimited and variable-length names, polygonal
faces, sparse boundary c1 sections, 2D extrusion, and mesh data embedded in
`.cas.h5`. Boundary kinds come from `zoneType`, not element `faceType`.
Cell material is kept as the numeric source type instead of being guessed.
Adapted, nonconformal and overset layouts are rejected. CFF periodic pairing
metadata is not imported; raw periodic boundary types are retained.

Layout and convention references (read for documentation; code is original):

- [Ansys legacy file-format appendix](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_app_form_case_data.html)
- [VTK vtkFLUENTCFFReader, BSD](https://github.com/Kitware/VTK/blob/master/IO/FLUENTCFF/vtkFLUENTCFFReader.cxx)
- [Ansys CFF mesh data model](https://developer.synopsys.com/docs/common-fluids-format-2024-r2/_mesh_d_m.xhtml)
- [Ansys CFF boundary-type enumerations](https://developer.synopsys.com/docs/common-fluids-format-2024-r2/group__enums.xhtml)
- [PyFluent case/mesh reader documentation](https://fluent.docs.pyansys.com/version/0.38/api/filereader/case_file.html)

**Real Fluent acceptance remains unverified.** CFF tests create small synthetic
HDF5 files independently through the C API. Legacy tests use small hand-made
ASCII and binary files. Before claiming Fluent compatibility on a production
case, re-verify a real Fluent `.msh` and `.msh.h5`/`.cas.h5`, including polygonal
cells, zone metadata, winding and total volumes.

## VTU and tests

`write_vtu()` writes ASCII VTK XML UnstructuredGrid, with tetrahedron (10),
hexahedron (12), and polyhedron (42) connectivity. Polyhedra have outward face
streams and cumulative faceoffsets; ordinary cells in a mixed file use -1 face
offsets. Scalar and interleaved xyz cell fields are supported. Tests load/check
XML nesting, connectivity sizes, positive tet/hex ordering, polyhedron streams,
field components and escaping. No ParaView or binary/base64 dependency is needed.

Run `cmake -S . -B build -G Ninja && cmake --build build && ./build/kutta-tests`.
Repeat with `-DKUTTA_ENABLE_HDF5=OFF` in another build directory to check the
stub. Optional NASA tests read the coarsest available formatted flat-plate and
bump grids directly from
`/Users/ezzy/.Others/Kutta/validation/data/nasa_tmr/*/grids/`; they report a skip
when that external directory is absent. No large validation files are copied.
