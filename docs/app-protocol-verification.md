# Sidecar protocol verification

Evidence from the implementation checkout on Apple M4 Pro, 2026-10-30. This
covers CLI/data behavior, not acceptance of the solver's physical accuracy.
Existing strict FV gates and the `wfgc` opt-in remain unchanged.

- Quick tier with `KUTTA_FLUENT_EXAMPLES` set: 170 cases passed, 12.29 s in a
  warmed build. The first baseline run before these changes took 17.56 s;
  startup/cache state affects that comparison, so it is not a speedup claim.
- Native `ctest -LE slow`: tests and portability lint passed. CPU and Metal
  compiled-GPU availability checks ran on the actual Apple device.
- CUDA: built in `nvidia/cuda:13.0.2-devel-ubuntu24.04`, amd64, with the CI
  configure flags `-DKUTTA_ENABLE_CUDA=ON
  -DCMAKE_CUDA_ARCHITECTURES="75-real;86-real;89-real;120"`. The CPU oracle suite
  passed with `compiled GPU backends*` excluded, as in CI. This Docker host has
  no NVIDIA driver; the GPU-presence test fails if explicitly included. No CUDA
  runtime or performance claim follows from compilation.
- Windows: the console control handler is limited to setting a flag for
  CTRL_C/CTRL_BREAK. Named-pipe EOF uses nonblocking `PeekNamedPipe`. Windows
  compilation and runtime could not be exercised on this machine.

Quick-tier coverage includes command schemas, setup/iteration/warning/done/error
ordering, real POSIX SIGINT after an iteration event, early and active stdin EOF,
final outputs after stopping, projected force coefficients, direct downforce
reports, strict override validation, and modified-output provenance. A generated
inlet/outlet case changes actual integrated inlet mass flow from -1 to -2 kg/s
when speed is overridden from 1 to 2 m/s. The source case hash stays unchanged.
The saved mixing-elbow case inventories 9,844 settings, including 26 editable
controls when seeded. Retained raw settings remain read-only.

The independent Python decoder verifies the view header, byte offsets, scalar
and vector fields/ranges, triangle areas on known box surfaces and slices,
area-weighted vertex interpolation, Fluent deltas, coincident/outside planes,
and replacement of stale exports. A concave prism checks nonconvex polygon and
slice triangulation. Generated saved face fields verify the N/face to Pa shear
conversion and face-field priority over owner-cell fallbacks. A real 17,822-cell
Fluent mixing-elbow seed/export produced all surface/slice comparison fields.

## Export measurements

`bench_export MESH OUTPUT_DIR` manufactures fields on imported meshes; it does
not solve them. `bench_export --10m OUTPUT_DIR` constructs a 250 x 250 x 160
structured hex mesh directly, without the general generator's face maps.

| Mesh | Cells | Faces | In-memory export, one mid-z slice | View bytes |
| --- | ---: | ---: | ---: | ---: |
| NASA NACA0012 897 x 257, extruded | 229,376 | 918,464 | 0.6840 s | 112,175,506 |
| Fluent wing | 114,219 | 654,343 | 0.0175 s | 1,767,414 |
| Manufactured structured volume | 10,000,000 | 30,142,500 | 0.3477 s | 29,700,865 |

The NACA mesh is the largest available external fixture by cell count; its
compressed grid was decompressed into the ignored build directory. It has far
more boundary data than the three-dimensional wing because extrusion caps cover
every cell. The synthetic 10M run wrote its snapshot in 1.405 s and used a peak
resident set of 1.948 GB. It is a scalability probe with manufactured fields.

The actual CLI command, including loading that 10M snapshot and exporting an
oblique slice, was:

```sh
./build/kutta export-view build/evidence/export-10m --slice 0.5,0.5,0.5,1,1,1
```

It reported **2.05762 s**, 2.09 s process wall time and 2.003 GB peak RSS.
These are single local measurements, not cross-platform bounds. Every zone is
processed independently; the slice path builds one cell-face adjacency and one
plane's distances/candidate mask, with output storage proportional to the cut.
No nodal field arrays for the full volume or per-cell polyhedron objects are
created. Array payload I/O is batched and explicitly little endian.

## Volume post-processing follow-up

Implementation source: `ff000f0738a3966e854e91da670d3e0c9c2285dc`, on the same
Apple M4 Pro. The follow-up changes app/export, I/O and the new `post/` module,
plus documentation, tests and the CMake source list. It changes no FV, linear
algebra, runtime or desktop UI source.

Twenty additional quick tests cover affine reconstruction and neighbour walking,
uniform-flow and solid-body-rotation streamlines, deterministic plane/surface
seeds, time-interpolated pathlines, saved-time knots, particles starting at rest,
physical flight time at length limits, analytic iso-plane area/winding/fields,
concave-cell kernel selection, line and ordered closed-section probes, distant
probe samples, analytic pressure/viscous forces and moments, constrained pressure
centres, absolute mass weighting, null undefined reports, reference overrides and
nonlinear Fluent deltas. The independent Python decoder exercises the actual CLI,
KUTL/KUTV layouts/ranges/indices, probe CSVs, reports, frame files and atomic
publication failures. The private snapshot reader accepts both v1 and v2.

The real saved Fluent mixing elbow is solved for **one CPU iteration**, seeded
from its `.dat.h5`, with the existing explicit `--approximate wfgc` opt-in. This
is post-processing smoke evidence, not a converged aerodynamic validation.
Its complete batch produced:

| Output | Observed result |
| --- | --- |
| Boundary zones | 6 |
| Surface-seeded traces | 2 curves, 31 vertices |
| Pressure iso-surface | 14,343 vertices, 27,985 triangles |
| Q=0 iso-surface | 51,441 vertices, 100,419 triangles |
| Volume line probe | 7 samples |
| Wall section probe | 80 vertices |
| Field catalog across artifacts | 59 fields, including Fluent/delta diagnostics |
| Aero reports | 6 zones plus their combined total |

The initial default diagonalization exposed warped Fluent faces whose triangles
invert about the Fluent face-pyramid centroid. Shared face diagonal selection
and a positive-kernel interior point resolve these without changing solver
geometry or discarding inverted tetrahedra. A narrow L-shaped manufactured prism
checks a centroid outside the cell kernel and an exact 0.9375 m² iso-section.
Tracing through its concavity stops at the first boundary even when the proposed
endpoint would lie inside another arm.

### Post-processing scalability

The existing structured generator constructs 10,000,000 cells, 30,142,500 faces
and 10,143,161 points. For this measurement, pressure is manufactured as `p=x`
Pa at cell/face centres, with constant `U=(1,1,1)` m/s. Density is 2 kg/m³,
viscosity 0.01 Pa s, reference density 2 kg/m³, reference speed 1 m/s, reference
pressure zero, reference area 1 m² and reference length 1 m. This is an analytic
export workload, not a CFD solve.

The actual CLI batch loads the saved volume, exports all boundary zones and an
oblique slice, adds Cp/Q fields, integrates two 0.2 m streamlines, extracts the
nonempty `pressure=0.5` iso-plane, samples a 51-point diagonal probe and computes
all-zone surface reports. It took **47.6557 s** (47.77 s process wall time), with
**4,634,722,304 bytes peak RSS**, while other verification was running. The view
occupies 73,572,145 bytes. Independent binary readback found 280,821 iso vertices,
315,382 triangles and area **1.00000000000004 m²**, with every x coordinate within
1e-5 m of 0.5. Probe pressures agree with analytic x to 1e-5 Pa. Both streamline
lengths are 0.2 m and six zone reports are finite.

The ordinary surface/single-slice in-memory path on that mesh took 0.422964 s,
with 1.950 GB process peak RSS including mesh generation and snapshot writing.
Volume reconstruction, gradients, point location and marching add the extra
post-processing cost; the app still reads only compact output geometry.

The v1 protocol defines steady and time-dependent pathlines separately, the
incompressible total-pressure normalization, absolute mass weighting, wall-stress
sign, pressure-only centre of pressure and residual couples. These definitions
avoid ambiguous Fluent comparisons. The public binary version stays 1; only
the private run snapshot advances to v2 to retain physical/reference properties
and solver face mass fluxes.

Final verification on this source:

- Native quick tier: **190 cases, 1,601,918 assertions passed in 12.70 s**;
  `KUTTA_FLUENT_EXAMPLES` points to the supplied external fixtures.
- Native `ctest -LE slow`: **2/2 passed in 13.59 s**, including portability lint.
  Real CPU/Metal backend availability checks are included.
- CUDA 13.0.2 amd64 container: the CI architecture configure line compiled
  without warnings; **189 CPU-oracle cases, 888,928 assertions passed**, with
  only the compiled-GPU-presence case excluded as in CI. The mixing-elbow fixture
  was mounted and its new post-processing smoke test ran. CUDA execution remains
  unverified because this host has no NVIDIA driver.
- Windows compilation and runtime remain unverified on this Mac. This follow-up
  adds no Windows-specific implementation code; the existing control-handler
  limitation above remains explicit.

Sample artifacts retained in the ignored build directory are
`build/evidence/post-elbow-final/{post.json,summary.json,view/manifest.json}` and
its compact binaries, CSV probes and `reports.json`. The summary records the
exact source commit above. The 10M readback is
`build/evidence/post-10m-readback.json`; its compact view, request, timings and
resource logs remain, while the large temporary volume snapshot was removed.
