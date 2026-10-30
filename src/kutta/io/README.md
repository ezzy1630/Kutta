# Mesh import and topology

The I/O additions and `mesh/topology.hpp` provide the
construction and gather operations without changing the solver's face-based
mesh contract. All geometry is computed in double precision on the host.

## Geometry and adjacency

Faces use a triangle fan about the arithmetic mean of their vertices. The area
vector is the sum of the triangle area vectors; the face centroid is weighted by
signed triangle area projected onto the resultant face normal. Cell volume and first moments are integrated by signed
pyramids (tetrahedra) from a local estimated center to those same triangles.
Local coordinate differences avoid unnecessarily losing precision on translated
meshes. This is an original implementation of the geometric decomposition, not
copied from OpenFOAM or another GPL source.

`validate()` checks CSR sizes/ranges, finite coordinates, valid distinct face
nodes, owner/neighbour ordering, unique patch names, complete contiguous patch coverage, cell-zone
indices, positive volumes, outward winding, area closure, and paired oriented
cell edges. Closure tolerance is `1e-10 * sum(magnitude(face area))` per cell.
Planar simple faces may be concave; warped-face moments follow the center-fan
decomposition. Validation additionally requires each face normal to point away
from its cell centroid. Arbitrary self-intersecting geometry is not supported.

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
Plot3D adds no zlib dependency.

`read_plot3d(grid_path, thickness)` automatically loads a sibling `<stem>.nmf`
when present. `read_plot3d(grid_path, nmf_path, thickness)` selects an explicit
map; an empty explicit path disables discovery. A present but invalid map fails
instead of falling back to unclassified block boundaries.

The [NASA neutral map format](https://turbmodels.larc.nasa.gov/nmf_documentation.html)
uses one-based inclusive **point** ranges. Both range reversal and one-to-one
`Swap` orientation are supported. Declared interface points must coincide within
`1e-12 * max(1, |point1|, |point2|)`; only those points are merged. Paired faces
become internal faces with ordered owners/neighbours and outward winding.
Every boundary face must have exactly one BC or interface assignment. Invalid
indices, mismatched dimensions/ranges, overlaps, gaps and unsupported topology
(`patched`, collapsed/pole, `unprocessed`) fail with map path/line context.
Comments, quoted or bare case-insensitive type names, and the trailing backslashes
in TMR maps are accepted. Custom BC names retain `PatchKind::Other`.

Ordinary 3D maps use the grid's dimensions. Extruded maps may use either
`ni × nj × 2` or the TMR span-first layout `2 × ni × nj` (map i/j/k become
mesh z/x/y). If both shapes coincide, the ordinary layout takes precedence.
BC patches are named `<type>_block<B>_f<F>_<S1>-<E1>_<S2>-<E2>` with lowercase
and hyphens normalized to underscores; `source_type` retains the original type.
Wall, symmetry, farfield, inflow and back-pressure types carry `PatchKind` hints;
case setup still assigns the solver's actual boundary condition. All extruded
blocks share one `front` and one `back` patch, retaining the map's symmetry type.
Mapped caps must be symmetry boundaries with consistent source types.

Without a map, each block has independent `blockN_imin`, `blockN_imax`,
`blockN_jmin`, `blockN_jmax` patches (also `kmin`/`kmax` for 3D). Coincident
coordinates alone do not imply connectivity. Extruded caps still share one
`front`/`back` pair. `split_patch()` remains available for other partial boundaries:

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
raw bytes. Uniform binary cell sections may omit the payload and still carry
the binary trailer. Node flags (2041/3041) are skipped as one 32-bit word per
node in their declared range. Unknown binary sections are rejected rather than
scanned. ASCII settings and trailing zone lists are skipped with balanced
parentheses, respecting quoted strings and escaped quotes. Node coordinates may
be floats or doubles; connectivity is 32-bit.

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

CMake discovers HDF5 1.10.2 or newer optionally and publishes
`KUTTA_HAS_HDF5=1` or `0`. Disable it explicitly with `-DKUTTA_ENABLE_HDF5=OFF`. The stub reader throws a clear
installation/configuration error when HDF5 is unavailable. HDF5's imported CMake
target carries the C library's platform-specific include/link definitions.
CFF builds also require zlib 1.2.9 or newer for the legacy compressed-chunk
compatibility path; HDF5-disabled builds require neither library.

The reader uses the HDF5 C API and `/meshes/1` layout: dimension, node
zoneTopology/coords, cell zoneTopology, face zoneTopology, faces/nodes sections
with `nnodes`/`nodes`, and sectioned c0/c1 datasets with minId/maxId attributes.
Coordinate dataset names can be zone IDs or section ordinals: placement always
uses the ranges, independently of names or traversal order. Cell `ctype`
sections (uniform or mixed) are checked by range against cell-zone topology.
Overlaps, missing coverage, inconsistent element types and mismatched root
`nodeCount`/`cellCount`/`faceCount` attributes are rejected.
It supports fixed semicolon-delimited and variable-length names, polygonal
faces, sparse boundary c1 sections or full-length sections with zero boundary
neighbours, 2D extrusion, and mesh data embedded in `.cas.h5`. Boundary kinds
come from `zoneType`, not element `faceType`.
Cell material is kept as the numeric source type instead of being guessed.
Adapted, nonconformal and overset layouts are rejected. CFF periodic pairing
metadata is not imported; raw periodic boundary types are retained.

Layout and convention references (read for documentation; code is original):

- [Ansys legacy file-format appendix](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_app_form_case_data.html)
- [VTK vtkFLUENTCFFReader, BSD](https://github.com/Kitware/VTK/blob/master/IO/FLUENTCFF/vtkFLUENTCFFReader.cxx)
- [Ansys CFF mesh data model](https://developer.synopsys.com/docs/common-fluids-format-2024-r2/_mesh_d_m.xhtml)
- [Ansys CFF boundary-type enumerations](https://developer.synopsys.com/docs/common-fluids-format-2024-r2/group__enums.xhtml)
- [PyFluent case/mesh reader documentation](https://fluent.docs.pyansys.com/version/0.38/api/filereader/case_file.html)

Older Fluent writers can deflate only the valid rows of the final HDF5 chunk.
HDF5 2.2 rejects this layout ([upstream issue 6691](https://github.com/HDFGroup/hdf5/issues/6691)).
Only an ordinary numeric read failing with HDF5 2.2's specific
"chunk size is incorrect after being unfiltered" diagnostic permits the raw-chunk
fallback. It requires rank 1/2, complete rows, a single deflate filter, and a
partial final chunk. Each complete stream must expand to exactly a full chunk
or the valid final rows, with at least one unpadded final chunk. Decoded and
converted buffer sizes are bounded against the output. Other HDF5 read failures
propagate without compatibility recovery. HDF5 performs
numeric conversion, including byte order. Corrupt, truncated and unsupported
chunks fail; inputs are never rewritten.

### Verified real Fluent files

The optional `real Fluent external mesh acceptance` test reads these exact
[Ansys example-data](https://github.com/ansys/example-data/tree/master/pyfluent)
files (MIT license), in place under
`<checkout>/validation/fluent_examples/`. It skips cleanly when
that external directory is absent. No external fixture is copied into the repo.
Counts below are for the returned 3D mesh; cylinder and cavity are extruded one
unit thick from their source 2D grids. Other files are native 3D.

| File relative to the examples directory | Cells | Faces | Points | Patches | Total volume |
| --- | ---: | ---: | ---: | ---: | ---: |
| `flow_over_cylinder/flow_over_cylinder.msh` | 46,949 | 188,392 | 95,090 | 6 | 449.214654586013 |
| `lid_driven_cavity/lid_driven_cavity.msh` | 3,550 | 14,321 | 7,344 | 5 | 1.00000000000001 |
| `simple_duct/simple-duct.msh.h5` | 10 | 51 | 44 | 3 | 9.99999904 |
| `box/poly.msh.h5` | 86,803 | 591,932 | 505,132 | 6 | 0.500000000000006 |
| `aero/wing.msh.h5` | 114,219 | 654,343 | 478,848 | 5 | 6995.01687397708 |
| `file_session/elbow1.cas.h5` | 13,852 | 44,500 | 16,968 | 5 | 0.00263365409702418 |
| `mixing_elbow/mixing_elbow.cas.h5` | 17,822 | 91,581 | 66,417 | 6 | 0.00250065507937294 |
| `static_mixer/Static_Mixer_Parameters.cas.h5` | 22,771 | 119,140 | 82,247 | 4 | 3.95970656064711e-08 |
| `static_mixer/Static_Mixer_Parameters_legacy_text.cas` | 22,771 | 119,140 | 82,247 | 4 | 3.95970656064711e-08 |
| `static_mixer/Static_Mixer_Parameters_legacy_binary.cas` | 22,771 | 119,140 | 82,247 | 4 | 3.95970656064711e-08 |

Expectations in `tests/io/test_fluent_real.cpp` come from `h5dump` root count
attributes and zoneTopology ranges/names, or legacy declaration and zone
headers. Before extrusion, the cylinder has 94,494 edges and 47,545 points;
the cavity has 7,221 edges and 3,672 points. Both `.msh` files contain binary
sections despite the extension. Tests assert total/internal face counts, point
counts and every boundary/cell zone's ID, name and size. They run `validate()`
and `compute_geometry()`, require positive cell volumes, check every face's
area vector points away from its owner centroid, and require global boundary
area closure within `1e-10 * boundary area`. The measured minimum face-normal
cosine is positive in every file (lowest: wing, 0.21637); the largest relative
boundary closure residual is 2.0e-14 (cylinder).

The static mixer's three formats have identical connectivity/owners/neighbours,
matching point coordinates within `1e-12` of mesh scale, matching patch
names/IDs/kinds/ranges, and total volume equal within `1e-12` relative. These
checks establish import/topology/geometry compatibility for the listed files,
not solver boundary-condition setup or every Fluent layout.

The companion `file_session/elbow1.dat.h5` and
`mixing_elbow/mixing_elbow.dat.h5` are also checked: they are data-only files,
so mesh dispatch rejects their extension and the direct CFF reader reports no
`/meshes/1` mesh. With HDF5 disabled the four legacy meshes still run; the six
CFF mesh cases and two data-only files are explicitly skipped.

## VTU and tests

`write_vtu()` writes ASCII VTK XML UnstructuredGrid, with tetrahedron (10),
hexahedron (12), and polyhedron (42) connectivity. Polyhedra have outward face
streams and cumulative faceoffsets; ordinary cells in a mixed file use -1 face
offsets. Scalar and interleaved xyz cell fields are supported. Tests load/check
XML nesting, connectivity sizes, positive tet/hex ordering, polyhedron streams,
field components and escaping. No ParaView or binary/base64 dependency is needed.

Run `cmake -S . -B build -G Ninja && cmake --build build && ./build/kutta-tests`.
Repeat with `-DKUTTA_ENABLE_HDF5=OFF` in another build directory to check the
stub. Optional NASA tests read all seven available uncompressed flat-plate,
bump and NACA0012 grids directly from
`<checkout>/validation/data/nasa_tmr/*/grids/`; they report a skip
when that external directory is absent. They check every NMF patch name/type/size,
stitched point/internal-face counts, absence of coincident boundary faces, and
global boundary area closure. No large validation files are copied.
