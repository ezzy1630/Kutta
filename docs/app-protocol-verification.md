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
