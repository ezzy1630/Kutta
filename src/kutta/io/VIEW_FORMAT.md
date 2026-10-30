# Kutta view format, version 1

`kutta solve` writes `view-input.bin` alongside `result.vtu`. The completion
marker is `summary.json`, written last. `kutta export-view DIR` reads only that
snapshot; the original case/data paths need not exist. Runs from older versions
without a snapshot must be solved again. The snapshot contains the solver's
RCM mesh and float32 cell/boundary fields, including reference fields gathered
by exact Fluent global IDs before the source files are closed. It is private
and is not an app input.

`view/manifest.json` is the public entry point. Its `v` and `format.version` are
1. `zones` and `slices` each describe one binary file, its vertex/index counts,
and an ordered `fields` array. Every descriptor has `name`, `unit`, `components`
(1 or 3), `location: "vertex"`, `offset_bytes`, `min`, and `max`. Top-level
`fields` combines ranges across files. Scalar ranges are scalar values; vector
ranges are vector magnitudes, including delta vectors. Empty/all-missing ranges
are JSON null. Bounds and plane coordinates are in metres. Zone filenames use
patch ordinals, avoiding path characters and duplicate/unavailable source IDs;
`id` and `name` retain source identity. Slices have zero-based filenames.

## Public binary layout

All integers are **uint32 little endian**, all values **IEEE-754 float32 little
endian**. No padding, compression, strings or interleaving between arrays.
Offsets in the manifest are absolute byte offsets from the start of the file.
Writers are `view_geometry.cpp::write_geometry` and `binary.hpp`.

| Byte offset | Type | Meaning |
| --- | --- | --- |
| 0 | 4 bytes | ASCII `KUTV` |
| 4 | uint32 | version = 1 |
| 8 | uint32 | header bytes = 32 |
| 12 | uint32 | vertex count V |
| 16 | uint32 | index count I (multiple of 3) |
| 20 | uint32 | field count F |
| 24 | uint32 | reserved = 0 |
| 28 | uint32 | reserved = 0 |
| 32 | float32[3 V] | positions, interleaved xyz |
| 32 + 12 V | uint32[I] | triangle indices (zero-based) |
| 32 + 12 V + 4 I | float32[V C0] | first field; vector components interleaved xyz |
| subsequent | float32[V Cn] | remaining fields, in manifest order |

Triangles follow outward boundary winding or the requested slice normal.
Concave faces are ear-clipped in their dominant projection; warped faces retain
3D vertex positions. There are no normals in the file: the viewer can compute
triangle normals. Absolute positions must be representable as finite float32.

## Fields and comparison semantics

| Name | Unit | Placement |
| --- | --- | --- |
| `pressure` | Pa | all surfaces and slices; gauge pressure |
| `pressure_coefficient` | 1 | when saved reference pressure, positive density and speed exist; `(p - p_ref) / (0.5 rho_ref U_ref^2)` |
| `wall_shear` | Pa | walls; xyz stress **on the fluid** from the active closure (laminar: viscous stress) |
| `wall_shear_magnitude` | Pa | walls; face stress magnitude |
| `y_plus` | 1 | walls, SST only; u-star definition |
| `y_star` | 1 | walls, k-epsilon only; wall-function y* = rho Cmu^0.25 k^0.5 y / mu |
| `velocity` | m/s | non-wall surfaces and slices; xyz |
| `velocity_magnitude` | m/s | non-wall surfaces and slices |
| `k` | m^2/s^2 | slices; SST and k-epsilon |
| `omega` | 1/s | slices; SST |
| `epsilon` | m^2/s^3 | slices; k-epsilon |
| `viscosity_ratio` | 1 | slices; dynamic turbulent/molecular viscosity; turbulent models |

Surface vertex values are face-area-weighted averages **within that zone**.
Magnitude and Cp are derived before averaging. Kutta boundary pressure uses the
solver's reconstructed pressure; non-wall velocity uses the inlet value,
tangential symmetry value, or extrapolated outlet owner value. Wall stress and
y-plus/y-star use the solver's wall diagnostics, with no new wall law in the exporter.
Laminar runs export pressure, velocity and wall shear only.

Slices intersect cell faces with the requested plane, connect intersection
segments into closed polygons, and triangulate the polygons. Fields use the
intersected cell's constant value, with separate vertices on cell interfaces;
this intentionally preserves finite-volume discontinuities without inventing
nodal gradients. Coincident mesh faces are emitted once, with the owner's
values. Tangential point/edge intersections have zero area and are omitted.
An outside plane produces an empty file with null ranges. Planes are six finite
comma-separated numbers `ox,oy,oz,nx,ny,nz`; normals are normalized and must be
nonzero. The tolerance scales with mesh bounds and coordinate roundoff.

When `--data` seeded the run, `reference` is true. Available saved fields use
`<field>@fluent`; `<field>@delta` is **Kutta minus Fluent**, component by
component for vectors, after the same surface/slice interpolation. Missing
reference quantities are omitted; partial coverage uses quiet NaNs in the binary
and propagates to delta vertices. The viewer must mask nonfinite samples.
A missing quantity is never silently reported as zero. Saved boundary pressure
and velocity take priority; where absent, saved owner-cell values are used.
Fluent `SV_WALL_SHEAR` is N/face on the fluid: divide by face area before
comparison in Pa. `SV_MU_T` is divided by the adapted molecular viscosity.
Reference results are the supplied saved data, not an extra converged run.

The provenance string retains the approximation notice. A view is written to a
staging directory and replaces the previous view only after every file closes.

## Private snapshot layout

For reproducibility, the private snapshot is also explicitly little endian:
`KUTR`, uint32 version 1, uint32 reference flag, length-prefixed UTF-8 provenance;
six int32 counts (points, faces, internal faces, cells, face-node entries,
patches); float64 xyz points; int32 CSR face offsets, node IDs, owners and
internal neighbours. Each patch stores two UTF-8 strings (name, source type),
uint32 kind, int32 source zone/start/size. Cell fields then boundary fields each
start with uint32 count; each field stores UTF-8 name/unit, uint32 components,
and its float32 array of known row count. Strings use uint32 byte lengths.
The reader rejects unknown versions, trailing bytes, invalid sizes, ranges and
indices before exporting. Full solver geometry validation already ran before
the trusted snapshot was written.

## Polyline geometry, KUTL version 1

Streamlines, pathlines and probes use the same scalar representation as KUTV.
Their header is eight uint32 words: ASCII `KUTL`, version 1, header size 32,
vertex count V, index count I, field count F, line count L, reserved zero.
Positions (3 V float32) and indices (I uint32) follow. Next are **L+1 uint32
line offsets**, starting at zero and ending at I. A line uses the concatenated
index interval `[offset[j], offset[j+1])`; one-point lines and zero lines are
valid. Vertex fields then follow in manifest order, at their absolute byte
offsets. Geometry descriptors include `primitive:"polylines"` and `line_count`.
KUTV headers/payloads remain unchanged and add `primitive:"triangles"` metadata.
Polyline index counts need not be multiples of three. Never feed KUTL to a
triangle renderer. Missing probe/reference values are quiet NaNs and must be
masked; every stored geometry coordinate is finite.

## Face edges

Zone files triangulate polygon faces, so their triangle edges include
diagonals that are not mesh edges. Exports (`export-view`, `export-geometry`;
not live snapshots) also write `zone-N.edges.bin` per zone, a KUTL v1 file with
no fields. Its positions are byte-identical to, and in the same order as, the
zone file's, so a viewer may index the zone's vertex buffer with it. Its
indices hold every polygon edge of the zone exactly once (no diagonals, no
duplicates), chained greedily into edge-disjoint polylines, starting at
odd-degree vertices: I = E + L for E edges and L lines. Shared edges between
zones appear in each zone. The zone descriptor gains `edges`, the descriptor of
that file: `file`, `primitive:"polylines"`, `vertices`, `indices`,
`line_count` and `fields:[]`.

Manifest key order is not significant. A case-geometry manifest
(`export-geometry`) adds `geometry_only:true` and `case:{"sha256"}`, with
`reference:false`, empty `provenance`, no slices and no fields.

## Private snapshot v2

New solves write KUTR version 2. The v1 layout remains unchanged through the two
field lists. Seven optional float64 properties follow, each encoded as uint32
presence (0 or 1) then a float64 value when present, in this order: actual density,
actual molecular viscosity, reference density, reference velocity, reference
pressure, reference area, reference length. A `mass_flux` boundary field stores
the solver's integrated outward kg/s per boundary face. Readers accept both
versions. Public protocol and public KUTV/KUTL versions remain 1.


## Live surface publication

`solve --live-view DIR` uses the same KUTV v1 zone files and view manifest under
`DIR/live/<run_id>/snapshot-<sequence>/`. The small `DIR/live/manifest.json` pointer
has `v`, `run_id`, `sequence`, `iter`, `elapsed` and a relative `manifest` path.
`live_view.cpp` closes and renames the snapshot staging directory before replacing
that pointer by same-filesystem rename. Existing immutable snapshots remain
readable. A snapshot progress event is emitted after publication and uses an
absolute UTF-8 view-manifest path. The app never reads an in-progress payload or
loads the volume. See the app protocol for cadence, budget and retention.
