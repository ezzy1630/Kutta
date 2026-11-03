# Peak memory of a full `kutta solve`

The undated sections describe the earlier host-memory round. The latest
[lossless device-compaction round](#2026-11-03--lossless-device-compaction-from-8011f3f)
is at the end; it did not reach the 2× target.

Peak process memory for a whole `kutta solve` (mesh read, case setup, solver setup,
iterations, VTU and loads output) is 37–41% lower on the Ahmed meshes. Every output
stays bitwise identical to main 96b30c9 on CPU and Metal. Device (GPU-buffer) bytes per
cell are unchanged. All of the saving is host memory that, on Apple silicon, shared
the unified pool with the device buffers.

| Case (Metal, SST, Coupled, 5 it) | Cells | Before | After | Change |
| --- | ---: | ---: | ---: | ---: |
| V6 mixing elbow | 17,822 | 0.31 GB (17,277 B/cell) | 0.22 GB (12,390 B/cell) | −28% |
| Ahmed L1 | 403,604 | 2.41 GB (5,972 B/cell) | 1.46 GB (3,622 B/cell) | −39% |
| Ahmed L2 | 690,104 | 3.75 GB (5,428 B/cell) | 2.40 GB (3,474 B/cell) | −36% |
| Ahmed L3 | 1,330,032 | 7.63 GB (5,735 B/cell) | 4.50 GB (3,380 B/cell) | −41% |
| Ahmed L4 (20 it) | 2,893,920 | not run (projected 16.6 GB) | 9.64 GB (3,329 B/cell) | |
| Ahmed "L4.25" (20 it) | 3,544,744 | not run (projected 20.3 GB) | **11.77 GB (3,321 B/cell)** | |

The elbow's per-cell figure is high because ~120 MB of fixed Metal and process
overhead is spread over only 18k cells.

The largest measured run on this Mac (M4 Pro, 24 GB unified memory) is **3.54M cells
at 11.77 GB**. Inside the 14 GB budget for one run, the after build fits about
**4.2M cells**, against about **2.4M** before. Both of those are projections from the
B/cell figures.

## Method

- **Measure.** `phys_footprint` is the number macOS uses for memory pressure and
  jetsam. On Apple silicon it includes every touched Metal buffer, because host and
  device share one pool. With `KUTTA_MEMORY_PHASES=1`, `kutta solve` now prints one
  line per phase:
  - the current footprint;
  - its peak since the previous phase, from `ri_interval_max_phys_footprint`, which
    is reset at each mark;
  - live malloc bytes;
  - live device-buffer bytes (each backend's buffer handle now tallies them).

  Without the variable, every mark is a no-op. The process peak is cross-checked
  with `/usr/bin/time -l` ("peak memory footprint").
- **Before and after.**
  - **Before** is main 96b30c9 plus only the probe commit (cb7ff05), built from
    `git archive` into a separate tree.
  - **Before-returned** is the same binary run with `MallocSpaceEfficient=1`. It
    separates the allocator effect from the structural changes.
  - **After** is this branch.

  Runs: `kutta solve CASE --backend metal --iterations 5 --approximate wfgc`, with
  2 runs per binary for L1–L3, taking the maximum. L4 and L4.25 used 20 iterations
  and one run each.
- **Ahmed cases.** The meshes are all-hex, from `tools/fsae_aero/ahmed.py` (35° slant).
  L4.25 used a scratch copy of that generator that accepts a fractional `--level`.
  `kutta solve` needs a case file, so `tools/memory_scale/ahmed_case.py` wraps each
  `.msh` in a legacy `.cas`, using the V6 elbow's settings and one zone thread per
  zone type. The meshes are not committed.
- **Box.** A 1.02M-cell synthetic hex box was also run with `bench_memory`
  (SIMPLEC, SST). See `raw/bench_box.log`.

## Where the bytes went

**Ahmed L3 before.** The peak was 5,735 B/cell, of which:

- **device buffers: 2,652 B/cell**, the solver's working set. This lane did not
  change it.
- **retained freed host memory: ~2,170 B/cell.**
  - macOS malloc keeps freed blocks of every size resident for reuse, and
    `phys_footprint` counts them.
  - Metal allocates its buffers outside malloc, so none of the freed setup memory
    (mesh reader maps, validation geometry, reconstruction scratch) was ever
    reused.
  - Each setup temporary therefore stacked on top of the device buffers allocated
    after it.
  - `malloc_zone_pressure_relief` returns nothing. The allocator reads its mode only
    at process start (`MallocSpaceEfficient` / `MallocLargeCache`).
  - Measured gap: 5,735 before vs 3,568 before-returned.
- **live host transients at the phase peaks:**
  - **read-mesh ~1 KB/cell.** The Fluent readers kept `std::map`s keyed by source
    id, one node vector per face and one `detail::Face` record per face.
  - **`validate()`: ~340 B/cell of face geometry.** It ran five times per solve.
  - **walk face reconstruction: ~870 B/cell.** It kept a 124-byte record per face,
    a full copy for upload and all four double moment levels at once.
  - **smaller items:** three float copies per Vec3d upload and two double copies of
    each VTU field.

**After.** The process peak is at `write-vtu`, at 3,380 B/cell on L3. It breaks down as:

- **device buffers: 2,652 B/cell.**
- **host heap: ~380–460 B/cell** in steady state.
- **fixed Metal overhead: ~117 MB.**
- **unattributed footprint: ~95 B/cell.** It scales with the mesh, and is probably
  driver-side.
- **export transient: 140 B/cell** above the iteration footprint (3,380 vs 3,240).

The host-side setup peak, estimated as phase peak minus device bytes, is now about
900 B/cell, at `stencil`. Before, even with the allocator retention removed, it was
about 1.4 KB/cell, at `face-reconstruction`.

### Ahmed L1, bytes per cell by phase

The peak column is the highest footprint since the previous mark. Device is the live
device-buffer bytes at the mark (after). The phases from `renumber` to `solver-ready`
are inside the solver's constructor.

| phase | before | before-returned | after | device (after) |
| --- | ---: | ---: | ---: | ---: |
| read-case | 214 | 206 | 207 | 0 |
| read-mesh | 1,443 | 1,085 | **431** | 0 |
| case-setup | 524 | 381 | **151** | 0 |
| renumber | 743 | 495 | 340 | 0 |
| geometry | 743 | 408 | 408 | 0 |
| connectivity | 746 | 416 | 416 | 0 |
| upload-mesh | 1,535 | 929 | 712 | 241 |
| stencil | 2,348 | 1,459 | 1,430 | 506 |
| fields | 2,764 | 1,339 | 1,335 | 784 |
| face-reconstruction | 4,272 | 3,041 | **2,158** | 1,471 |
| transport | 4,329 | 2,071 | 2,068 | 1,656 |
| lsq | 4,383 | 2,320 | 2,317 | 1,705 |
| turbulence-lsq | 4,636 | 2,563 | 2,560 | 1,948 |
| turbulence-transport | 4,730 | 2,462 | 2,459 | 2,114 |
| turbulence | 4,735 | 2,383 | 2,380 | 2,210 |
| pressure-solver | 4,834 | 2,611 | 2,612 | 2,287 |
| coupled-solver | 5,004 | 3,132 | 3,131 | 2,658 |
| solver-ready | 5,004 | 2,683 | 2,680 | 2,658 |
| initialize | 5,796 | 3,324 | 3,321 | 2,659 |
| first-iteration | 5,961 | 3,490 | 3,487 | 2,659 |
| iterations | 5,961 | 3,490 | 3,487 | 2,659 |
| write-vtu | 5,972 | 3,808 | **3,622** | 2,659 |
| forces-fluxes | 5,932 | 3,460 | 3,457 | 2,659 |
| loads-summary | 5,960 | 3,458 | 3,456 | 2,659 |
| **process peak** | **5,972** (2.41 GB) | **3,808** (1.54 GB) | **3,622** (1.46 GB) | |

### Ahmed L3, bytes per cell by phase

| phase | before | before-returned | after | device (after) |
| --- | ---: | ---: | ---: | ---: |
| read-case | 194 | 192 | 192 | 0 |
| read-mesh | 1,504 | 1,065 | **409** | 0 |
| case-setup | 762 | 358 | **126** | 0 |
| renumber | 990 | 480 | 318 | 0 |
| geometry | 990 | 378 | 375 | 0 |
| connectivity | 1,016 | 394 | 391 | 0 |
| upload-mesh | 1,572 | 908 | 688 | 241 |
| stencil | 2,250 | 1,442 | 1,416 | 506 |
| fields | 2,517 | 1,166 | 1,164 | 781 |
| face-reconstruction | 4,070 | 2,858 | **1,985** | 1,462 |
| transport | 4,126 | 1,903 | 1,902 | 1,646 |
| lsq | 4,179 | 2,063 | 2,062 | 1,696 |
| turbulence-lsq | 4,448 | 2,305 | 2,304 | 1,937 |
| turbulence-transport | 4,541 | 2,292 | 2,290 | 2,103 |
| turbulence | 4,545 | 2,200 | 2,199 | 2,197 |
| pressure-solver | 4,646 | 2,452 | 2,451 | 2,277 |
| coupled-solver | 4,838 | 2,993 | 2,992 | 2,652 |
| solver-ready | 4,838 | 2,514 | 2,514 | 2,652 |
| initialize | 5,563 | 3,080 | 3,080 | 2,652 |
| first-iteration | 5,731 | 3,241 | 3,240 | 2,652 |
| iterations | 5,731 | 3,241 | 3,240 | 2,652 |
| write-vtu | 5,735 | 3,568 | **3,380** | 2,652 |
| forces-fluxes | 5,726 | 3,241 | 3,240 | 2,652 |
| loads-summary | 5,685 | 3,242 | 3,240 | 2,652 |
| **process peak** | **5,735** (7.63 GB) | **3,568** (4.75 GB) | **3,380** (4.50 GB) | |

The elbow, L2, L4 and L4.25 tables can be regenerated from the raw records with
`python3 tools/memory_scale/tables.py validation/results/memory_scale/raw/phases.jsonl.gz`.
After the change, L4 and L4.25 follow L3 to within 3% per phase.

**Box via `bench_memory`** (1.02M cells, SIMPLEC SST, Metal). Device memory is
1,636 B/cell in all three runs.

| Run | Peak |
| --- | ---: |
| before | 4.47 GB |
| after, default allocator | 3.62 GB |
| after, space-efficient malloc | 2.53 GB |

`bench_memory` does not go through `kutta`'s `main()`, so it does not re-execute.
The third run sets the variable by hand.

## What changed

| Commit | Change | Effect |
| --- | --- | --- |
| cb7ff05 | Per-phase memory probe (`src/kutta/core/memory.*`), device-byte tally in the CPU, Metal and CUDA buffer handles, and phase marks in `solve.cpp` and the solver constructor. | Measurement only. |
| 69e3554 | On macOS, `kutta` re-executes itself once at start with `MallocSpaceEfficient=1`, so freed memory goes back to the system. An explicit `MallocSpaceEfficient` setting is respected. | L1: 5,972 → 3,808 B/cell. L3: 5,735 → 3,568 B/cell. Same solve time. |
| 899f02a | `validate()` recomputes each face's centroid and area where a check needs them, through a new `face_geometry()`. That is the function `compute_geometry` uses, so every check sees the same doubles. Previously `validate()` built whole-mesh face arrays, five times per solve. The case adapter's planar-inlet and mass-flux checks use `face_geometry()` too. The zone duplicate check uses a bitset. | −340 B/cell transient at every validate. |
| 8373ee4 | The Fluent `.msh` and `.cas.h5` readers store source ids in id-range tables instead of `std::map`s, and face nodes in one flat array. Assembly builds the CSR mesh directly, with the same numbering, stable face order and checks. The legacy case read limit is raised from 512 MiB to 2 GiB; a 2.9M-cell text case is 553 MB. | read-mesh: 1,085 → 431 B/cell (L1). Mesh read takes about half the time. |
| 81fe0d7 | Walk face reconstruction keeps only corrected rows per chunk and builds the moment levels one from the previous. `upload_mesh` and the stencil's `vectors()` convert one float component at a time. `result.vtu` fields are downloaded straight into the writer's arrays. | face-reconstruction: 3,041 → 2,158 B/cell. write-vtu: 3,808 → 3,622 B/cell (L1). |
| (final) | Report, scripts and evidence. `memory.cpp` guards `NOMINMAX` the way `stop.hpp` does (Windows only). | |

### Changes in the other lane's files

No arithmetic is changed, and no kernel or `*.kern.h` file is touched.

- `src/kutta/fv/steady_solver.cpp` has only `#include "kutta/core/memory.hpp"`
  and 16 `memory::phase("…")` calls between existing statements in the
  `Impl` constructor. They mark renumber, geometry, connectivity, upload-mesh,
  stencil, fields, face-reconstruction, transport, lsq, turbulence-lsq,
  turbulence-transport, turbulence, pressure-solver, coupled-solver and
  solver-ready. The diff is 16 lines and can be dropped, or merged line by line.
- `src/kutta/fv/operators.cpp` changes only the anonymous `vectors()` helper.
  It converts x, then y, then z into one reused float vector instead of three
  vectors. The float values and upload order are unchanged.
- `src/kutta/fv/face_projection.cpp` (walk reconstruction) and
  `src/kutta/fv/device_mesh.cpp` (`upload_mesh`) change only how
  intermediate rows are held and freed. The double expressions, fits, slot order
  and uploaded floats are the same.

## Larger meshes elsewhere (projected, not measured)

CUDA and Windows were not executed here: CUDA compiles only in CI, so neither
platform's memory has been measured. The RTX figures below are projections from
the Metal device tally, with the same 0.6 GB CUDA-context-and-display reserve
as `validation/results/performance_scale`.

- **Device limit.** On a discrete GPU the limit is device memory, and this lane did
  not change the solver's device arrays.
  - Ahmed-like external aero, SST and Coupled (2,653 B/cell): about **5.8M cells
    on 16 GB** and **8.8M cells on 24 GB**.
  - SIMPLEC on a box-like hex mesh (1,636 B/cell): about **9.4M / 14.3M cells**.
- **Host memory (64 GB).** It is not the limit. The host side peaks at about
  0.9 KB/cell during setup and holds about 0.4–0.7 KB/cell during the run and
  export, so 8.8M cells need roughly 8 GB of host memory.
  - The retention problem above comes from macOS malloc sharing one pool with
    Metal. On a discrete GPU, retained host memory does not take device memory.
    That is inferred from the platform design, not measured.
- **This Mac (24 GB unified).** A single run at the measured 3,321 B/cell reaches
  about 4.2M cells at 14 GB. That figure is projected; 3.54M cells is the largest
  measured run.

## Bitwise parity

`tools/memory_scale/parity.py` compares main 96b30c9 against the candidate. For each
case and backend it hashes:

- `fluxes.csv`, `forces.csv`, `loads.json` and `result.vtu`;
- `residuals.csv`, without its seconds-per-iteration column.

Every run used `kutta solve CASE --iterations 40 --out DIR`. The candidate was a
frozen copy of the 81fe0d7 code. The final commit changes only a preprocessor guard
inside the `_WIN32` branch of `memory.cpp` and adds no code to the macOS build. The
check was run twice and both passes were identical. Full hashes of all five files
for the second pass are in `raw/parity.jsonl`.

| Case | Backend | Identical | `result.vtu` SHA-256, first 16 hex (both binaries) |
| --- | --- | --- | --- |
| V6 mixing elbow, `--approximate wfgc` | Metal | yes | `1e9f77eaf776650d` |
| V6 mixing elbow | CPU | yes | `4f04ea3939b20f41` |
| elbow_param, `--approximate compressibility,ke-options` | Metal | yes | `d5dd146b8d27220a` |
| elbow_param | CPU | yes | `359a2b4519b5a051` |
| Ahmed L1 (403,604 cells), `--approximate wfgc` | Metal | yes | `50aff8d6755a386a` |
| Ahmed L1 | CPU | yes | `23418bf010a41e7b` |

## Quick tier and lint

- `./build/kutta-tests --test-suite-exclude=slow`, with fixtures:
  - 250 test cases passed, 29 skipped, 0 failed;
  - 7,417,230 assertions passed, 0 failed.
- `python3 tools/lint_portability.py` passes.

Details are in `raw/quick_tier.txt`.

## Open items

- **Device bytes.** These are now the limit everywhere.
  - Largest consumers on L1 (`raw/device_sites_ahmed_l1.log`):
    - walk-reconstruction coefficients: 343 B/cell;
    - solver per-cell arrays: 280 B/cell, across 70 arrays;
    - cell-face-entry vectors: 216 B/cell;
    - walk geometry: 160 B/cell;
    - `TransportEquation`: 233 B/cell;
    - Coupled block solver: 371 B/cell, including ten 4n Krylov vectors at
      160 B/cell.
  - A shared `face_correction`/`zero_` scratch across transport equations would save
    about 64 B/cell. It was left out because it changes the API the other lane is
    editing.
- **Stencil transient.** The stencil phase holds about 900 B/cell of host data
  (cell faces, double geometry, pattern) at the setup peak. This is next if the
  host side matters.
- **Export transient.** `write-vtu` adds 140 B/cell over the iteration footprint,
  mostly from its `validate()` and the field arrays. The VTU could be streamed field
  by field.
- **Unattributed footprint.** About 95 B/cell scales with the mesh but appears in
  neither the heap nor the device tally. It is probably Metal driver-side.
- **Legacy case reader.** `.cas` text files are read whole, so a 3.5M-cell case is
  ~0.7 GB, briefly, at `read-case`.
- **Re-exec on macOS.**
  - `kutta` re-executes once with the same pid, arguments and file descriptors, so
    the app's `QProcess` client is unaffected.
  - Child processes inherit `MallocSpaceEfficient=1`.
  - If `execv` fails, kutta runs on in the default mode.
  - Anything that links the solver in-process, such as `bench_memory` and the test
    binary, does not get the re-exec.
- **CUDA and Windows** memory is unmeasured. The Windows probe uses `PrivateUsage`
  and `PeakPagefileUsage` and compiles only in CI.

## Reproduce

```sh
# Ahmed case files (meshes are gitignored, not committed)
python3 tools/fsae_aero/ahmed.py --slant 35 --level 1 --out runs/mesh/ahmed35_L1.msh
python3 tools/memory_scale/ahmed_case.py --mesh runs/mesh/ahmed35_L1.msh \
    --template validation/fluent_examples/mixing_elbow/mixing_elbow.cas.h5 \
    --out runs/mesh/ahmed35_L1.cas
# per-phase memory, appended as JSON lines
python3 tools/memory_scale/phases.py --kutta build/kutta --label after \
    --case runs/mesh/ahmed35_L1.cas --cells 403604 --iterations 5 --jsonl phases.jsonl
python3 tools/memory_scale/tables.py phases.jsonl
# parity against a main build
python3 tools/memory_scale/parity.py --baseline BASE/build/kutta --candidate build/kutta \
    --out /tmp/parity --backends cpu,metal
```

## 2026-11-03 — device storage, baseline e014e36

Metal SST Coupled device storage on Ahmed L3 falls **2,865 → 2,710 B/cell
(−5.4%)**. Peak process memory falls **4.685 → 4.479 GB**. The 20–40M-cell car
still does not fit an 8–12 GB GPU.

Before is e014e36 plus allocation-site diagnostics only. Candidate solver source
is e63e157; the frozen binaries were built before that commit and retain e014e36
as their configure-time version string. Their solver source is identical to
e63e157. [manifest.json](device-2026-11-03/manifest.json) records source commits
and binary hashes. The older report's 2,652 B/cell is not this round's baseline:
e014e36 has additional solver state and caches. Every pair uses the same input.

GB below means 10^9 bytes. Device figures are requested buffer payload, rather
than process footprint or measured available RTX VRAM. CUDA was not executed.

### Measured process and device memory

All CLI pairs use Metal, SST Coupled and five iterations, including mesh read,
setup, stored-case initialization, iterations and output. One run per binary
per case; L3 ran once per binary and alone in this lane. The process peak is
macOS phys_footprint cross-checked by /usr/bin/time -l. Device is the live
iteration payload, including the backend reduction buffers in the CLI tally.

| Case | Cells | Peak process GB, before → after | Device B/cell, before → after | Device reduction |
| --- | ---: | ---: | ---: | ---: |
| V6 polyhedral elbow | 17,822 | 0.236 → 0.233 | 3,867.3 → 3,704.4 | 4.21% |
| Ahmed L1 | 403,560 | 1.555 → 1.493 | 2,875.5 → 2,723.2 | 5.30% |
| Ahmed L2 | 689,796 | 2.542 → 2.436 | 2,867.0 → 2,713.4 | 5.36% |
| Ahmed L3 | 1,329,048 | 4.685 → 4.479 | 2,865.2 → 2,710.4 | 5.40% |

The existing Ahmed L1 fixture has 403,560 cells; newly generated L2 and L3 have
689,796 and 1,329,048 cells. These differ from the older report. Input hashes
and the exact baseline/candidate binary hashes are in the accompanying records.
No wall-time speedup is claimed: all runs shared the machine with other lanes.

### Inventory totals and savings

The inventory benchmark runs one first iteration and one additional iteration,
uses default paired state, plain SST and fixed benchmark settings, and maps mesh
patch hints to velocity inlet, pressure outlet or wall. This is a memory probe;
its boundary setup differs from the CLI case adapter. Backend-owned allocations
are excluded from these totals. The CSV sums reconcile with MEMORY_TOTAL at its
printed precision; the table uses the exact sums.

| Mesh | Coupling | Device B/cell, before → after | Face scratch + zero saving, B/cell | LSQ saving, B/cell | Total saving, bytes |
| --- | --- | ---: | ---: | ---: | ---: |
| Ahmed L1 hex | SIMPLEC | 2,351.71 → 2,199.35 | 64.53 | 87.83 | 61,488,000 |
| Ahmed L1 hex | COUPLED | 2,875.34 → 2,722.98 | 64.53 | 87.83 | 61,488,000 |
| V6 polyhedra | SIMPLEC | 3,235.62 → 3,072.65 | 98.22 | 64.76 | 2,904,600 |
| V6 polyhedra | COUPLED | 3,862.73 → 3,699.75 | 98.22 | 64.76 | 2,904,600 |
| Ahmed native tet/prism | SIMPLEC | 2,036.73 → 1,897.38 | 52.76 | 86.58 | 16,147,488 |
| Ahmed native tet/prism | COUPLED | 2,512.33 → 2,372.99 | 52.76 | 86.58 | 16,147,488 |

Boundary-adjacent cell counts are 26,592 (hex), 4,490 (poly) and 8,800
(native). For **every pair**, the measured device-byte reduction equals
`16*(F+N) + 96*N - 124*C` exactly. All other allocation lengths are unchanged.
The native mesh contains 55,914 tetrahedra and 59,968 prisms; its two iterations
are a storage measurement, not a convergence or force-accuracy acceptance.

### Projected GPU capacity

These are projections, not RTX measurements. GB means 10^9 bytes. The budget
reserves 0.6 GB for context/display/driver use (an assumption inherited from the
previous report), then uses the source-derived hybrid-startup bound described
below. Real allocations and a different mesh/case can reduce capacity further.

| Mesh / coupling | Iteration Mcell/GB, before → after | Startup-bound B/cell after | 8 GB, Mcells | 12 GB, Mcells | 24 GB, Mcells |
| --- | ---: | ---: | ---: | ---: | ---: |
| Hex / SIMPLEC | 0.425 → 0.455 | 2,502.0 | 2.96 | 4.56 | 9.35 |
| Hex / COUPLED | 0.348 → 0.367 | 3,025.6 | 2.45 | 3.77 | 7.73 |
| Poly / SIMPLEC | 0.309 → 0.325 | 3,495.6 | 2.12 | 3.26 | 6.69 |
| Poly / COUPLED | 0.259 → 0.270 | 4,122.8 | 1.79 | 2.77 | 5.68 |
| Tet/prism / SIMPLEC | 0.491 → 0.527 | 2,159.0 | 3.43 | 5.28 | 10.84 |
| Tet/prism / COUPLED | 0.398 → 0.421 | 2,634.6 | 2.81 | 4.33 | 8.88 |

For the larger L3 hex case, the measured Coupled iteration payload is
2,710 B/cell after, about 0.369 million cells per GB before startup/reserve.
A 20–40M-cell car still needs roughly 54–108 GB for the Coupled iteration
payload alone (and more at startup). **The 8–12 GB full-car target remains unmet.**
SIMPLEC is smaller, but this round does not put that mesh range within a consumer
card either. The largest unchanged costs are walk reconstruction, redundant
geometry layouts, gradients, Coupled vectors and coarse block state.

### Byte-identical verification

All ten case/backend pairs completed **40 actual iterations** with identical
`fluxes.csv`, `forces.csv`, `loads.json` and `result.vtu` under `cmp`. Residuals
also pass `cmp` after the prescribed `grep -v '^#' | cut -d, -f1-7`. Comments and
the timing column are excluded. These compare before versus after **within each
backend**; they do not claim that CPU and Metal match one another.

[parity.jsonl](device-2026-11-03/parity.jsonl) contains the full hashes and counts.
The VTU hashes below are the first 16 hex digits, equal for both binaries.

| Case | CPU: all five outputs | Metal: all five outputs | CPU VTU hash | Metal VTU hash |
| --- | --- | --- | --- | --- |
| V6 Coupled, `wfgc` | identical | identical | `affcc3c8b86189f5` | `7d6a40de90a0406b` |
| elbow_param SIMPLE, `compressibility,ke-options` | identical | identical | `6da600a6ac7b28a1` | `85e1941771284b3c` |
| Ahmed L1 Coupled, V6-derived wrapper | identical | identical | `8a1fd6d93377bde9` | `3acd4286b6d55a67` |
| Ahmed L1 SIMPLEC, V6-derived wrapper | identical | identical | `1b95a4d50953e855` | `c55dabbe7cac7e9e` |
| Ahmed L1 SIMPLEC, 40 m/s air | identical | identical | `51544af0ca87d55f` | `106c301a507e0275` |

The scale wrappers inherit V6 settings and are memory/parity fixtures, not the
Ahmed wind-tunnel setup. The separate native
[Ahmed SIMPLEC case](../../../tools/memory_scale/ahmed_simplec.kutta.json) uses
40 m/s air, SST, 1% inlet turbulence, viscosity ratio 10, stationary body/ground
and symmetry on the remaining farfield patches. It uses standard initialization,
automatic pseudo time, momentum relaxation 0.7 and pressure relaxation 1. Its
settings inspection succeeds, and it needs no approximation acknowledgement.
Forty-iteration parity is not drag or convergence acceptance.

- Full default build passes; portability lint passes.
- Quick tier: **308 cases, 7,990,274 assertions, zero failures**. The 37 excluded
  cases are exactly the slow suite (`--test-suite=slow --count` confirms 37).
- The focused shared-inverse test passes 120 byte comparisons on CPU/Metal,
  including alternating FaceValue/Omitted/NormalFootValue fits on hex and tet
  meshes. It compares both inverse words and resulting gradients with `memcmp`.
- The phase runner now propagates a failed child exit; a controlled child-failure
  check confirms it. Parity rejects a run that stops before the requested count.

The brief allows the byte-identical exception to a new full slow-tier run. This
round changes buffer ownership and copies stored words, with no arithmetic,
model, default control, reference, threshold or iteration-budget change. The
specified parity procedure, native external-aero pair and quick tier pass. The
full slow tier/comparator was **not rerun**; slow regressions, recoveries and
changed assertion values are **not measured**, and no zero-count claim is made.
[checks.json](device-2026-11-03/checks.json) records the verification scope.

CUDA compilation/execution, CUDA/Windows host RSS, controlled timing of the extra
copy dispatches and full-scale car runs remain unverified. The allocations are
smaller; this is not a claim that every physical or performance gate was rerun.

### What changed

1. The three momentum and two SST transport equations share one face-correction
   array and one zero field. Each assembly writes every face correction before
   reading it. Paired residuals also rebuild the correction before using it.
   RHS, diagonals, limiters and gradients stay separate. Saving: `16*(F+N)` bytes
   for plain SST, independent of the coupling scheme.
2. Five boundary least-squares fits share one dense six-plane inverse. Each fit
   stores the original six float words only for cells touching a boundary. A
   copy kernel restores those rows immediately before a fit is consumed. Interior
   rows come from the same `cell_fit` expression with no boundary observations;
   they are identical across fits. Saving: `96*N - 124*C` bytes, where `C` is the
   number of distinct boundary-adjacent cells. This includes the shared cell-id
   map. No setup-only host inverse copy remains alive after construction.

The LSQ change adds a boundary-row copy dispatch per fit application. No numerical
threshold, solver control, iteration budget or reference value is changed. The
legacy standalone-fit path remains available to independent callers; sharing is
explicit and bound to one mesh and Device. Green-Gauss keeps its existing path.

### Exact buffer inventory

The accompanying `device-2026-11-03/*-buffers.csv.gz` files enumerate **each live
solver-owned allocation** before and after, for SST SIMPLEC and Coupled on an
Ahmed hex mesh, V6 polyhedra, and an Ahmed native tet/prism mesh. Every row has
allocation source, column, instance number, exact bytes, B/cell and B/face. These
are allocation records, not counts inferred by matching array lengths. Loop
allocations (AMG levels and vector components) each have their own row. Their exact byte
sums reconcile with `MEMORY_TOTAL` at its printed precision; exact totals are
in [inventory-summary.json](device-2026-11-03/inventory-summary.json).

For the formulas below, `N` is cells, `F` all faces, `I` internal faces,
`B=F-I` boundary faces, `E=F+I` cell-face incidences, `P=N+2I` fine CSR entries,
and `C` boundary-adjacent cells. Every device element is four bytes: float32 or
signed int32. A buffer with `L` elements costs exactly `4L/N` B/cell and `4L/F`
B/face. Vector rows mean three individually enumerated buffers, each of that
length. Neither global counts nor CSR offsets here fit in 16 bits.

| Owner | Individual buffers | Element count per buffer | Lifetime |
| --- | --- | --- | --- |
| DeviceMesh connectivity | owner; neighbour; cell_face_offsets; cell_faces; cell_face_sign | F; I; N+1; E; E | Whole solve |
| DeviceMesh geometry | area.{x,y,z}; owner_to_face.{x,y,z}; neighbour_to_face.{x,y,z}; owner_to_neighbour.{x,y,z} | F; F; I; F respectively | Whole solve |
| DeviceMesh scalar geometry | weight; volume; wall_distance; wall_normal_distance | I; N; N; B | Whole solve |
| Stencil connectivity | adjacent; matrix_entry | E each | Whole solve |
| Stencil geometry | lsq_delta.{x,y,z}; face_delta.{x,y,z}; boundary_lsq_delta.{x,y,z}; nonorthogonal.{x,y,z}; lsq_inverse; diffusion | E; E; B; F; 6N; F | Whole solve |
| Five boundary LSQ fits, before | boundary_delta.{x,y,z}; inverse | B; 6N per fit | Whole solve |
| Five boundary LSQ fits, after | boundary_delta.{x,y,z}; boundary_inverse | B; 6C per fit | Whole solve |
| Shared LSQ workspace, after | inverse; boundary_cells | 6N; C, once | Whole solve |
| Face jump LSQ | delta.{x,y,z}; inverse | E; 6N | Whole solve |
| Walk reconstruction | coefficients; geometry[0..3]; observed; slot; mask; moment[0..1] | 30*max(corrected faces,1); 10N each; B; F; F; 12N each | Whole solve; moment contents are scratch |
| Walk paired work (only compensated momentum) | paired_moment; paired_bounds[0..1] | 12N; 2N each | Lazy allocation, retained after first paired request |
| Shared fine matrix | row_ptr; col; diag; val | N+1; P; N; P | Whole solve, one structure/value array shared across scalar equations |
| Each transport equation | rhs; diagonal; consistent_inverse; grad.{x,y,z}; limiter | N each | Whole solve; overwritten each assembly |
| Transport scratch, before | face_correction; zero | F; N, five copies | Whole solve |
| Transport scratch, after | face_correction; zero | F; N, one copy | Whole solve, sequential contents |
| Omega projected diffusion scratch | diffusion_correction | F | Whole solve, scratch during assembly |
| Optional compensated transport | rhs_low; diagonal_low; scalar_delta; scalar_residual | N each per enabled equation | Whole solve |
| Optional scalar steady stopping | residual; absolute; scale | N each per enabled equation | Whole solve |
| Scalar Krylov workspace | R; Z; P; Q; Shadow; S; T; Zs; Best; Scratch | N each, shared across all scalar solvers | Whole solve, sequential scratch |
| Jacobi preconditioners | inverse; invalid | N each; invalid aliases scalar workspace Scratch | Whole solve, five separate inverses |
| Scalar solver controls | scalars; status; projection | 24; 3; 2 only for constant nullspace | Whole solve, per solver |
| Pressure AMG level l | l1; r; x; b; CSR row_ptr/col/diag/val on coarse levels | n_l each; fine x/b empty; CSR n_l+1, p_l, n_l, p_l | Whole solve, fine CSR comes from caller |
| Pressure AMG transfer l | aggregate; member_offsets; members; entry_offsets; entries | Exact level lengths in CSV | Whole solve; needed by refresh/restrict/prolong |
| Pressure AMG K-cycle only | p; q; z; t; scalars; status | n_l each on coarse levels; 12; 1 | Absent for these V-cycle runs |
| Coupled fine state | cell; fine_coefficients; coefficient; fine_r; fine_tmp | 27N; 5E; F; 4N; 4N | Whole solve; coefficients refreshed each outer solve |
| Coupled Krylov | x; b; r; shadow; p; q; s; t; z; zs | 4N each | Whole solve, scratch; x also used for update |
| Coupled controls | scalars; gauge; status | 24; 2; 3 | Whole solve |
| Coupled coarse level l | rows; cols; diag; a; c; g; h; d; x; b; r; tmp | n_l+1; p_l; n_l; p_l; p_l; 3p_l; 3p_l; 11n_l; 4n_l each last four | Whole solve |
| Coupled transfer l | aggregate; member_offsets; members | Exact level lengths in CSV | Whole solve |
| SST fields | mu_t; k_source; k_sink; omega_source; omega_sink; strain; vorticity; f1; f2; cdkw; diffusivity_k; diffusivity_omega; production_k; fixed_omega; wall_production | N each | Whole solve |
| SST cell masks | fixed_omega_mask; selected_wall_face | N each | Whole solve |
| SST face diffusion | face_diffusivity_k; face_diffusivity_omega | F each | Whole solve |
| SST wall fields | wall_mask; wall_tau; wall_u_star; wall_u_tau; wall_y_plus; wall_y_plus_utau; wall_shear.{x,y,z}; wall_shear_coefficient; wall_face_omega; wall_cell_omega; face_production | B each | Whole solve |

The steady-solver fields are also separately enumerated in the CSV:

- `N` floats each: pressure, pressure_rhs, correction, residual_u/v/w, continuity,
  provisional_continuity, scale, speed, wall_diagonal, wall_rhs, physical_diagonal,
  consistent_inverse, kinetic, omega, rho_cell, mu_cell.
- `N` vectors: velocity, body, source, limited_gradient[0..2], grad_p,
  grad_correction, absolute_residual, fluent_scale, strain_base[0..2],
  strain_gradient[0..2], and the five transport gradients.
- `F` floats each: mass, viscosity, coefficient, face_pressure, pressure_jump;
  `F` vectors: face_velocity and transpose_stress.
- `B` floats each: boundary_pressure, prescribed_pressure, correction_boundary,
  boundary_k, boundary_omega, turbulence_a, turbulence_b; `B` vectors:
  boundary_velocity, prescribed_velocity, inlet_direction.
- `B` ints each: boundary_kind, scalar_kind, total_pressure, open_boundary_kind,
  turbulence_spec, turbulence_kind; `N` ints each: zero_mask and poor when present.
- Paired-state low words and work: velocity_lo and pressure_gradient_lo (`N`
  vectors), pressure_lo/imbalance/momentum_delta (`N` each), mass_lo/previous_mass
  (`F` each), and pressure_boundary_lo (`B`). These remain allocated in both
  schemes. The inventory benchmark uses default paired state and does not enable
  compensated momentum or compensated scalar transport.
- Compensated momentum adds pressure_force_lo/canonical_velocity and six
  canonical gradients (`N` vectors), face_pressure_lo/paired_diffusion_work/
  paired_advection_work/canonical_mass (`F` floats), high_face_velocity (`F`
  vector), boundary_velocity_lo/canonical_boundary_velocity (`B` vectors), and
  scalar_kind_lo (`B` ints). They are absent unless that setting is enabled.
- The differentiable limiter adds six `N` scratch arrays; optional projected
  laminar gradients add three `F` vectors. They are absent in this SST benchmark.

Solver buffers are retained through output and destroyed with the solver. Host
setup vectors, coarsening graphs, double walk moments, the initial CSR values and
LSQ construction inverses are transient. Hybrid initialization additionally owns
a temporary transport equation/matrix, face zero/one fields, scalar iterates,
boundary data and potential flux; the CLI process measurements include that peak.
The pressure AMG hierarchy is used by hybrid initialization even with Coupled;
removing it unconditionally would break that path.

The CSV excludes backend-owned reduction/command buffers and allocator/driver
rounding. The CLI `DeviceBytes` tally additionally counts the Metal reduction
buffers: partial storage is 16,384 bytes; batched storage starts at 65,536 bytes
and grows with the largest simultaneous reduction batch. Metal sequence command
buffers are allocated outside this tally. Kernel threadgroup memory, Metal
pipelines, CUDA context/display use and hardware allocation granularity are not
solver buffer payload and must be budgeted separately.

### Host memory on a discrete GPU

CUDA's `CudaBuffer` uses `cudaMalloc`; upload/download use `cudaMemcpyAsync` and
complete before returning. It has no permanent host staging mirror for every
solver buffer. The shared LSQ dense buffer and boundary overrides therefore
remain device-only; their construction vectors are released after upload.

The solver deliberately retains host mesh/output data. Mirrored items are
owner/neighbour (4F/4I bytes each side), fine CSR pattern (row_ptr/col/diag:
`4*(2N+1+P)` bytes each side), and geometry in double on the host versus float
on the device. Host `MeshGeometry` is `48F+32N` bytes, excluding capacity and
allocator overhead. Host points, face-node CSR, zones, patches and permutations
also remain alive. Device area and volume mirror host quantities at float
precision; centroid differences and wall distances are derived data. Scalar
fields and Krylov vectors do not have permanent full host copies. Output
transiently downloads fields and constructs double VTU arrays.

At the L3 iteration mark, live host heap was 512,293,312 → 512,131,344 bytes
(385.46 → 385.34 B/cell): essentially unchanged.

The host heap/phase readings in the records are **macOS observations**, not CUDA
RSS measurements. Subtracting device bytes from a unified-memory footprint is
not a CUDA host-memory measurement. An RTX runtime and a Windows/Linux process
measurement remain unverified.

### Candidates left unchanged

- Recomputing face displacements from rounded coordinates can change bits.
  `face_delta` actually duplicates the already computed owner/neighbour face
  differences in a cell-incidence layout. Replacing those loads with face-side
  gathers could save `12E` bytes without new subtraction, but needs coordinated
  limiter/jump-kernel interfaces and its own parity measurement.
- Global indices are already 32-bit. Sixteen-bit global indices cannot represent
  any of the measured aero meshes. A mixed-width local format would need new
  kernel interfaces; no width or acceptance gate was weakened here.
- Fine matrix structure and scalar Krylov vectors already share storage.
- AMG maps and Coupled coarse block arrays remain intact. Restriction/prolongation
  and coefficient refresh consume them throughout the solve; setup ownership is
  not evidence that the arrays are setup-only.
- The Coupled fine-edge cache is `20E` bytes (about 120 B/cell on hex). Removing
  it trades repeated arithmetic for storage and overlaps the speed lane. It is
  unchanged here.
- Hybrid initialization creates a temporary transport linear solver but calls
  the separate pressure solver directly. An assemble-only transport path could
  avoid its unused Jacobi inverse, ten Krylov vectors and controls (`44N+108`
  requested bytes). This is a source-only opportunity, not a measured saving
  in this round; the temporary arrays remain in the startup capacity bound.
- Float16 storage was not introduced. No approximate-memory mode or unmeasured
  accuracy trade-off is hidden in these defaults.

Hybrid initialization's device transients, from the actual allocation sites in
`initialize_hybrid`, are unchanged by this round. Each row below again uses
`4L/N` B/cell and `4L/F` B/face for each listed buffer of length `L`:

| Temporary buffers | Length per buffer | Lifetime |
| --- | --- | --- |
| potential; old; boundary; zero; one | N; N; B; F; F | Hybrid initialization |
| potential_kind; prescribed_flux; potential_value; potential_flux; source_device | B; B; B; F; N | Latter part of hybrid initialization |
| solve_scalar bk; bv; bd; rhs | B; B; B; N | One scalar initialization call |
| solve_scalar LSQ delta.{x,y,z}; boundary_delta.{x,y,z}; inverse | E; B; 6N | One scalar initialization call |
| temporary transport row_ptr; col; diag; val | N+1; P; N; P | Hybrid initialization |
| temporary transport rhs; diagonal; consistent_inverse; grad.{x,y,z}; limiter; zero | N each | Hybrid initialization |
| temporary transport face_correction; Jacobi inverse | F; N | Hybrid initialization |
| temporary transport Krylov R/Z/P/Q/Shadow/S/T/Zs/Best/Scratch | N each | Hybrid initialization |
| temporary transport scalars; status | 24; 3 | Hybrid initialization |

For SST, the turbulence bootstrap reuses its persistent k/omega fields and does
not allocate another pair of cell fields. Zero-length placeholders request no
payload bytes; Metal/CUDA still create at least a four-byte backing allocation
for each such handle, outside the requested-byte totals.

`DeviceBytes::live()` is a live tally, not a high-water counter. The CLI phase
marks do not record the temporary hybrid arrays while `initialize_hybrid` is
inside a call. For capacity planning I therefore also use a source-derived
upper bound on concurrent hybrid payload: persistent bytes plus
`124N + 8P + 16F + 12E + 40B + 112` bytes. This conservatively includes the
second scalar initialization's LSQ and all later potential-flux fields together,
even when the stored case does not request that scalar solve. It applies to these
meshes (all have more than the scalar solver's 2,048-row speculation limit),
plain SST, the present V-cycle setup and the listed hybrid initializer. It is not
a measured CUDA peak. The GPU capacity projections use this bound, rather than
assuming that the iteration working set is also the startup peak.


### Reproduce this round

Use the input hashes in [input-hashes.json](device-2026-11-03/input-hashes.json).
The generated meshes and case wrappers stay outside git. The generator used
`gmsh==4.15.2` and `numpy==2.5.4`; the actual cell counts are recorded above.
Run each heavy solve serially, including both L3 runs.

```sh
uv run --with gmsh==4.15.2 --with numpy==2.5.4 python tools/fsae_aero/ahmed.py \
  --slant 35 --level 1 --out runs/mesh/ahmed35_L1.msh
python3 tools/memory_scale/ahmed_case.py --mesh runs/mesh/ahmed35_L1.msh \
  --template validation/fluent_examples/mixing_elbow/mixing_elbow.cas.h5 \
  --out runs/mesh/ahmed35_L1.cas
python3 - <<'PYCASE'
from pathlib import Path
text = Path("runs/mesh/ahmed35_L1.cas").read_text()
assert text.count("(flow/scheme 24)") == 1
Path("runs/mesh/ahmed35_L1_simplec.cas").write_text(
    text.replace("(flow/scheme 24)", "(flow/scheme 21)"))
PYCASE
# Allocation records: omit KUTTA_MEMORY_COUPLED for SIMPLEC.
MallocSpaceEfficient=1 KUTTA_MEMORY_BUFFERS=1 KUTTA_MEMORY_COUPLED=1 \
  build/bench_memory metal --mesh runs/mesh/ahmed35_L1.msh 1 > inventory.log
python3 tools/memory_scale/inventory.py inventory.log --revision REV > buffers.csv
# The default cases include both V6-derived Ahmed wrappers and native air/SST.
python3 tools/memory_scale/parity.py --baseline BASE/build/kutta \
  --candidate build/kutta --out runs/memory-parity --backends cpu,metal
```

The allocation-site provenance is permanent diagnostic support, not a retained
scratch probe. It does not change allocation sizes or numerical decisions.


## 2026-11-03 — lossless device compaction from 8011f3f

This records the original recompute-by-default trial, rejected during integration
review. The integration-review section below supersedes its default-mode
recommendation; these measurements are retained as trial evidence.

**The 2× target was not met.** Hex SIMPLEC goes from 2,199.35 to
2,067.36 device B/cell (6.00%), still well above 1,300. Hex Coupled saves 9.94%;
polyhedral Coupled saves 10.78%. The candidate is slower in the measured whole
CLI pairs. These results do not establish that a 5–20M-cell FSAE mesh fits an
8–12 GB RTX card.

The baseline is an unchanged build of `8011f3f`. The candidate solver is
`04b95595d8d54a64327f97e4f7e906408228179b`, built in the isolated checkout with
AppleClang 21, C++20 Release, HDF5 and Metal, on the M4 Pro with 24 GB unified
memory and macOS 26.6.2. The three local source commits are:

- `6b0e5ff`: canonical face-side gathers and sequential strain scratch.
- `b2db0a0`: Coupled edge recomputation and V-cycle residual scratch.
- `04b9559`: exact walk geometry row tables.

The later tools/report commit does not change this solver source. Full binary
hashes, the frozen intermediate checkpoints and measurement boundaries are in
[manifest.json](device-cut-2026-11-03/manifest.json). Input hashes are in
[input-hashes.json](device-cut-2026-11-03/input-hashes.json). Every comparison
uses the same input and the same iteration budget.

### Exact device inventory

The twelve `baseline-*-buffers.csv.gz` and `candidate-*-buffers.csv.gz` files in
[the evidence directory](device-cut-2026-11-03/) enumerate every live
solver-owned allocation at the benchmark iteration mark. Their exact integer
byte sums reconcile with `MEMORY_TOTAL`; the table uses those sums, not the
rounded log values. Backend reductions, command storage, allocation granularity
and thread-local storage remain outside this inventory. The benchmark uses
plain SST, default paired state, and its existing inlet/outlet/wall mapping;
it is a storage probe rather than the native aero boundary setup.

| Mesh | Coupling | Device B/cell, 8011f3f → candidate | Reduction |
| --- | --- | ---: | ---: |
| Ahmed L1 hex | SIMPLEC | 2,199.35 → 2,067.36 | 6.00% |
| Ahmed L1 hex | COUPLED | 2,722.98 → 2,452.42 | 9.94% |
| V6 polyhedra | SIMPLEC | 3,072.65 → 2,893.02 | 5.85% |
| V6 polyhedra | COUPLED | 3,699.75 → 3,301.01 | 10.78% |
| Ahmed native tet/prism | SIMPLEC | 1,897.38 → 1,783.19 | 6.02% |
| Ahmed native tet/prism | COUPLED | 2,372.99 → 2,149.83 | 9.40% |

[inventory-summary.json](device-cut-2026-11-03/inventory-summary.json) records
cells, faces, internal faces, buffer counts and exact totals. The meshes have
403,560 hex cells, 17,822 polyhedral cells and 115,882 native tet/prism cells.

### Changes, lifetimes and time cost

As above, `N` is cells, `F` faces, `I` internal faces and `E=F+I` incidences.

1. **Face displacement gathers: `12E` bytes saved.** The removed incidence
   vectors contained `float(face_centroid - cell_centroid)`. Those exact words
   already exist in `owner_to_face` and `neighbour_to_face`. Limiters and the
   jump correction select the face side using the existing incidence sign.
   No subtraction from rounded coordinates is introduced.
2. **Sequential scratch: `24N + 16(N + sum n_l)` bytes saved with Coupled.**
   The three base strain gradients become one vector, consumed by each
   component's jump correction before the next component overwrites it. The
   corrected gradients remain separate. The ordinary strain path adds two
   scalar gradient launches and two inverse restores versus the fused vector fit. In the Coupled
   V-cycle, the smoothing temporary holds the residual until restriction has
   consumed it; post-smoothing then overwrites it. Fine and coarse non-alias
   assertions check the actual solution/scratch handles. The base-gradient part
   saves `24N` with SIMPLEC; the Coupled residual part saves
   `16(N + sum n_l)`, where `n_l` are coarse-level cell counts.
3. **Coupled edge recomputation: `20E` bytes saved.** Five directed-edge
   coefficients are reconstructed from the mass, viscosity, weights and geometry
   frozen for the blocking solve. Their expressions and accumulation order stay
   the same. Volatile final float values preserve the former cache's rounding
   boundary before multiplication into the Krylov vector. The block inverse and
   coarse operators remain stored. This is the clear measured time penalty.
4. **Exact walk geometry rows: about `36N` bytes saved.** Rows are interned by
   all ten float words, including signed zero. Encounter order fixes row ids.
   An indexed level costs `40U + 4N` bytes for `U` unique rows; dense storage
   stays in use whenever it is smaller. On these meshes only hop zero compacts;
   the other three levels stay dense. The row dictionaries and construction
   maps are freed on the host after upload.

The per-change byte totals below reconcile exactly with the inventory delta;
[change-savings.json](device-cut-2026-11-03/change-savings.json) has integer
bytes, including the native mesh. Time attribution uses three alternating
pairs on V6 Coupled, 40 iterations, timing the final batch (21–40). Each row
compares successive frozen checkpoints. Intermediate binaries were built
before the source commits; their hashes and changes are recorded in the
manifest. The final whole-case measurements use the committed solver.

| Change | Hex B/cell saved, SIMPLEC / Coupled | Poly B/cell saved, SIMPLEC / Coupled | V6 paired median time change (range) |
| --- | ---: | ---: | ---: |
| Face-side gathers | 72.00 / 72.00 | 120.10 / 120.10 | -1.5% (-1.7% to -1.0%) |
| Sequential strain / V-cycle scratch | 24.00 / 42.57 | 24.00 / 42.95 | -6.2% (-53.5% to +4.2%) |
| Coupled edge cache removed | 0.00 / 120.00 | 0.00 / 200.16 | +22.8% (+14.2% to +22.9%) |
| Walk geometry row tables | 35.99 / 35.99 | 35.53 / 35.53 | -0.9% (-9.3% to +9.5%) |

The scratch and walk timing ranges cross zero and are inconclusive. Small
negative medians are not claimed as speedups. Coupled recomputation was slower
in every isolated pair. The exact data are in
[timing-summary.json](device-cut-2026-11-03/timing-summary.json) and the
`*-paired-runs.jsonl.gz` records. These are contended Mac timings, not RTX costs.

### Whole CLI process and iteration measurements

Each case uses Metal, SST Coupled and 40 actual iterations: three sequential
AB/BA pairs, with the last residual CSV batch (21–40) timed. Every run checks
actual mesh cells, completed iterations and all five output hashes. Peak process
memory is the maximum of the three runs per binary under `/usr/bin/time -l`,
including mesh read, setup, initialization, iterations and output. GB is decimal.
Device here is the live iteration payload including Metal reductions, so it
is slightly larger than the allocation inventory above. It is not a measured
CUDA allocation high-water mark.

| Case | Cells | Device B/cell, before → after | Peak process GB, before → after | Median s/iteration, before → after | Paired median time change |
| --- | ---: | ---: | ---: | ---: | ---: |
| V6 polyhedral elbow | 17,822 | 3,704.4 → 3,305.6 | 0.233 → 0.226 | 0.0711 → 0.0762 | +7.2% |
| Ahmed L1 | 403,560 | 2,723.2 → 2,452.6 | 1.492 → 1.382 | 0.2006 → 0.2373 | +17.3% |
| Ahmed L2 | 689,796 | 2,713.4 → 2,442.9 | 2.440 → 2.253 | 0.2846 → 0.3214 | +9.9% |
| Ahmed L3 | 1,329,048 | not measured | not measured | not measured | held |

**L3 remains unfinished.** A solo window was unavailable: the separate main
slow-tier process was still active with a 4.46 GB physical footprint at the
recorded gate. RSS understated that footprint. [l3-hold.json](device-cut-2026-11-03/l3-hold.json)
records the observation. The earlier section's L3 numbers are historical;
they are not substituted for a new before/after run. No owned solve remains
running.

The whole-case paired time-change ranges were +5.3% to +22.2% for V6,
+12.0% to +21.3% for L1, and +6.4% to +15.9% for L2. Paired medians need not
equal the ratio of the two separate timing medians.

### Projected RTX capacity

The same conservative hybrid-startup bound is retained:
persistent payload plus `124N + 8P + 16F + 12E + 40B + 112` bytes, with
`P=N+2I` and `B=F-I`. Pressure AMG still serves hybrid initialization under
Coupled, so it has not been removed. The table reserves 0.6 GB for context,
display and driver use, as before. Values use the measured mesh ratios and
are projections, not card measurements.

| Mesh / coupling | Startup-bound B/cell, before → after | 8 GB Mcells, before → after | 12 GB Mcells, before → after | 24 GB Mcells after |
| --- | ---: | ---: | ---: | ---: |
| Hex / SIMPLEC | 2,502.0 → 2,370.0 | 2.96 → 3.12 | 4.56 → 4.81 | 9.87 |
| Hex / COUPLED | 3,025.6 → 2,755.1 | 2.45 → 2.69 | 3.77 → 4.14 | 8.49 |
| Poly / SIMPLEC | 3,495.6 → 3,316.0 | 2.12 → 2.23 | 3.26 → 3.44 | 7.06 |
| Poly / COUPLED | 4,122.8 → 3,724.0 | 1.79 → 1.99 | 2.77 → 3.06 | 6.28 |
| Tet/prism / SIMPLEC | 2,159.0 → 2,044.8 | 3.43 → 3.62 | 5.28 → 5.58 | 11.44 |
| Tet/prism / COUPLED | 2,634.6 → 2,411.4 | 2.81 → 3.07 | 4.33 → 4.73 | 9.70 |

The 8 GB poly Coupled projection is still only about 2.0M cells. The
consumer-card full-car target remains open.

### Host mirrors and CUDA boundary

There are **zero new persistent host mirror bytes**. The current CUDA source
uses `cudaMalloc` and completed upload/download copies; it does not keep a
whole-solve staging mirror for every device allocation. Walk dictionaries are
construction temporaries. Coupled recomputation adds no iteration-time host
transfer. The removed stencil face-displacement construction also removes its
`24E`-byte host double-vector temporary.

Existing host owner/neighbour arrays, fine CSR pattern, mesh points, face-node
CSR, zones, patches and permutations remain. `MeshGeometry` remains
`48F + 32N` host bytes, excluding capacity/allocator overhead, against float
geometry on device. Scalar fields and Krylov vectors have no permanent full
host copies. In the measured CLI runs, maximum live heap at the iteration
mark was 408.730 → 408.729 B/cell on L1 and 397.029 → 397.047 on L2:
essentially unchanged. These are macOS heap observations, not CUDA host RSS.
CUDA compilation, actual VRAM high-water usage, Windows/Linux host RSS and
RTX iteration cost remain unverified.

### Parity and checks

All twelve final-source case/backend pairs completed 40 iterations and passed
`cmp` for `fluxes.csv`, `forces.csv`, `loads.json`, `result.vtu`, and residuals
after the literal `grep -v '^#' | cut -d, -f1-7` pipeline. Baseline outputs were
from actual 40-iteration runs of `8011f3f`; candidates were checked against the
final `04b9559` artifact. CPU and Metal are compared separately, not to one
another. [parity.jsonl](device-cut-2026-11-03/parity.jsonl) has full hashes.

| Case | CPU, 40 iterations | Metal, 40 iterations |
| --- | --- | --- |
| V6 Coupled, `wfgc` | identical | identical |
| elbow_param, `compressibility,ke-options` | identical | identical |
| V6 converted SIMPLEC, `wfgc` | identical | identical |
| Ahmed L1 Coupled, V6-derived wrapper | identical | identical |
| Ahmed L1 SIMPLEC, V6-derived wrapper | identical | identical |
| Native Ahmed L1 SIMPLEC, 40 m/s air | identical | identical |

- Full default Release build and portability lint pass.
- Quick tier: 308 cases and 7,990,274 assertions passed, zero failures.
  The 37 excluded cases match the separately counted slow suite.
- No test threshold, reference, control, iteration budget or skip was changed.
- Full slow tier and comparator were **not rerun** under the brief's
  byte-identical exception. Slow assertion changes and zero-regression counts
  are not measured. The separately running main slow baseline was not modified.
- Forty-iteration output parity is not converged force or Fluent acceptance.

[checks.json](device-cut-2026-11-03/checks.json) records the scope.

### Rejected work and remaining costs

Exact coefficient-row interning saved **zero bytes**: hex coefficient storage
stayed 138,654,120 bytes and poly storage stayed 10,268,640. It adds host table
construction, so it was removed. Its setup time was not measured in controlled
pairs. Geometry row interning remains because it does save payload.

The SST source/sink fields are consumed after `update()` by successive
assemblies and physical residuals. Its diagnostics are publicly exposed through
`SteadySolver::turbulence()` and remain valid after update; transition also reads
strain and vorticity. They cannot simply become overlapping scratch without
changing that contract or proving a replacement. Even hypothetical removal of
all fifteen cell float fields would save only 60 B/cell, versus the remaining
767 B/cell gap to 1,300 on hex SIMPLEC.

Walk coefficients, the three remaining dense geometry-moment levels, gradients,
fine state and Coupled Krylov/coarse block state still dominate. The coarse
operators and transfer maps remain live through refresh, restriction and
prolongation. A larger reconstruction/storage change is still needed to reach
2×; this round has not proved one. No float16/bfloat16 mode was implemented.
Its converged Ahmed/wing force effects and slow-tier results are unmeasured.

To reproduce the final paired measurement on the same generated inputs:

```sh
python3 tools/memory_scale/parity.py --baseline BASE/kutta \
  --candidate build/kutta --out runs/memory-parity --backends cpu,metal
python3 tools/memory_scale/paired.py --baseline BASE/kutta \
  --candidate build/kutta --case runs/mesh/ahmed35_L1.cas --cells 403560 \
  --iterations 40 --pairs 3 --out runs/memory-paired-l1
```

Use the existing fixture-link procedure and test environment from the brief.
The paired runner requires a fresh output directory, verifies cells and actual
iterations, and stops before another solve if a process peak exceeds its
3 GB default limit. L3 still requires a solo window and one run per binary.

## 2026-11-03 — integration review: cached default and opt-in memory mode

**Cached edges and the fused ordinary strain fit are restored as defaults.**
The final solver source is `95215f5`. `SolverSettings::coupled_edge_storage`
uses `CoupledEdgeStorage::Cached` by default; `Recompute` is the explicit
memory-limited alternative. The CLI selects it with
`--coupled-edge-storage recompute`, for either Fluent or Kutta case JSON inputs,
and records the choice in `summary.json`. No new case JSON property was needed.
The setting also selects sequential ordinary strain scratch with SIMPLEC.

The original combined `b2db0a0` is replaced in this branch by `2b12d9c`
(V-cycle residual scratch, default) and `ae31bee` (optional edge recomputation).
`6b0e5ff`'s canonical face gathers stay default. The exact walk tables from
`04b9559` were reapplied unchanged as `30aae18`. `95215f5` gates only the
ordinary sequential strain part of `6b0e5ff`. The compensated strain path was
already sequential and keeps its shared scratch without added launches.

The frozen `8011f3f` baseline and final binary hashes, source commits and
measurement methods are in [manifest.json](device-review-2026-11-03/manifest.json).
The later report/tools changes do not change the measured solver.

### Device payload and the remaining target

These are exact sums of the live solver-owned buffers in the benchmark,
excluding backend reduction/command storage and allocation granularity.
All twelve new default/memory snapshots reconcile with the formulas below.
The six before snapshots are the existing exact `8011f3f` inventory; all nine
current input hashes match that dataset. All eighteen compressed per-buffer
CSVs are retained in [the evidence directory](device-review-2026-11-03/).

| Mesh / coupling | 8011f3f B/cell | Default B/cell | Memory mode B/cell |
| --- | ---: | ---: | ---: |
| Ahmed L1 hex / SIMPLEC | 2,199.35 | 2,091.36 | 2,067.36 |
| Ahmed L1 hex / COUPLED | 2,722.98 | 2,596.42 | 2,452.42 |
| V6 polyhedra / SIMPLEC | 3,072.65 | 2,917.02 | 2,893.02 |
| V6 polyhedra / COUPLED | 3,699.75 | 3,525.17 | 3,301.01 |
| Ahmed tet/prism / SIMPLEC | 1,897.38 | 1,807.19 | 1,783.19 |
| Ahmed tet/prism / COUPLED | 2,372.99 | 2,264.18 | 2,149.83 |

The 2× target is still unmet. Hex SIMPLEC is 2,091 B/cell by default and 2,067
in memory mode, above the 1,300 target. This work does not establish that a
5–20M-cell FSAE mesh fits an 8–12 GB card.

For `E=F+I`, default savings are `12E` for canonical face gathers, the exact
walk-row saving, and `16(N + sum n_l)` for Coupled V-cycle residual reuse.
Restriction consumes the smoother temporary before post-smoothing overwrites
it; fine/coarse handle assertions guard solution/scratch non-aliasing.
Memory mode additionally saves `24N` for ordinary strain scratch and `20E`
for Coupled coefficients. The corrected strain gradients stay separate.

| Mesh | Face gathers B/cell | Walk rows B/cell | V-cycle residual B/cell | Optional strain B/cell | Optional edges B/cell |
| --- | ---: | ---: | ---: | ---: | ---: |
| Ahmed L1 hex | 72.000 | 35.987 | 18.565 | 24.000 | 120.000 |
| V6 polyhedra | 120.096 | 35.533 | 18.954 | 24.000 | 200.160 |
| Ahmed tet/prism | 54.210 | 35.985 | 18.615 | 24.000 | 90.350 |

[change-savings.json](device-review-2026-11-03/change-savings.json) retains
integer byte counts. Coupled Krylov vectors and coarse block coefficients are
unchanged. No reduced precision is introduced.

### Whole CLI paired timing and process peak

All runs use Metal and 40 actual iterations on the same contended M4 Pro.
Each Coupled case/mode has three alternating AB/BA pairs against `8011f3f`.
The harness checks cell counts, completed iterations and all five numerical
output hashes on every run. Whole-CLI wall time includes process startup,
mesh read, setup, initialization, iterations and export. The separate iteration
sample is the residual CSV batch for iterations 21–40. Peak process GB is the
maximum `/usr/bin/time -l` physical footprint, in decimal GB. Device here
includes Metal reductions, unlike the buffer inventory above.

| Case | Mode | Live device B/cell, before → after | Peak process GB, before → after |
| --- | --- | ---: | ---: |
| V6 Coupled | default | 3,704.4 → 3,529.8 | 0.223 → 0.220 |
| V6 Coupled | memory | 3,704.4 → 3,305.6 | 0.223 → 0.215 |
| Ahmed L1 Coupled | default | 2,723.2 → 2,596.6 | 1.484 → 1.434 |
| Ahmed L1 Coupled | memory | 2,723.2 → 2,452.6 | 1.485 → 1.376 |
| Ahmed L2 Coupled | default | 2,713.4 → 2,586.9 | 2.434 → 2.347 |
| Ahmed L2 Coupled | memory | 2,713.4 → 2,442.9 | 2.434 → 2.247 |
| Ahmed L1 converted SIMPLEC | default | 2,199.6 → 2,091.6 | 1.277 → 1.233 |
| Ahmed L1 converted SIMPLEC | memory | 2,199.6 → 2,067.6 | 1.277 → 1.223 |
| Ahmed L1 aero JSON SIMPLEC | default | 2,199.6 → 2,091.6 | 1.220 → 1.177 |
| Ahmed L1 aero JSON SIMPLEC | memory | 2,199.6 → 2,067.6 | 1.220 → 1.167 |

| Case | Mode / pairs | Whole-CLI paired median change (range) | Median s/iteration, before → after | Iteration paired median change (range) |
| --- | --- | ---: | ---: | ---: |
| V6 Coupled | default / 3 | -1.0% (-10.1% to +29.2%) | 0.0244 → 0.0235 | -2.4% (-3.5% to -1.5%) |
| V6 Coupled | memory / 3 | +1.9% (+0.7% to +2.4%) | 0.0394 → 0.0412 | +4.6% (+0.5% to +9.2%) |
| Ahmed L1 Coupled | default / 3 | -3.7% (-4.6% to +5.2%) | 0.1547 → 0.1495 | -3.4% (-12.2% to -0.2%) |
| Ahmed L1 Coupled | memory / 3 | +4.6% (-2.8% to +5.0%) | 0.1872 → 0.1946 | +0.9% (-14.5% to +11.7%) |
| Ahmed L2 Coupled | default / 3 | -3.1% (-13.6% to +7.5%) | 0.2721 → 0.2201 | -7.9% (-19.1% to +4.3%) |
| Ahmed L2 Coupled | memory / 3 | +4.4% (+3.9% to +7.8%) | 0.2280 → 0.2468 | +10.1% (+0.8% to +59.2%) |
| Ahmed L1 converted SIMPLEC | default / 1 | -0.3% (-0.3% to -0.3%) | 0.1956 → 0.1919 | -1.9% (-1.9% to -1.9%) |
| Ahmed L1 converted SIMPLEC | memory / 1 | -3.6% (-3.6% to -3.6%) | 0.2260 → 0.2020 | -10.6% (-10.6% to -10.6%) |
| Ahmed L1 aero JSON SIMPLEC | default / 2 | +1.8% (+0.2% to +3.4%) | 0.2397 → 0.2447 | +2.4% (-2.4% to +7.2%) |
| Ahmed L1 aero JSON SIMPLEC | memory / 1 | +0.6% (+0.6% to +0.6%) | 0.2405 → 0.2394 | -0.5% (-0.5% to -0.5%) |

All three Coupled defaults meet the requested whole-CLI criterion: absolute
paired median within about 3%, or a paired range crossing zero. Their iteration
medians are −2.4%, −3.4% and −7.9%; V6 is within the band and L2's range
crosses zero. L1's batch median is just outside a strict ±3% band on the faster
side. These contended measurements do not establish a
speedup. Memory mode's iteration medians are +4.6%, +0.9% and +10.1%; L1's
range crosses zero, while V6 and L2 are slower in every pair.

Converted SIMPLEC has one AB pair per mode. The native aero JSON default's
first AB pair was +3.4% whole CLI and +7.2% in the iteration batch. One additional
reverse-order BA pair was declared to investigate that result, and both samples
are retained above; the memory mode has one AB pair. The two native default
samples and their combination are in
[native-simplec-combined.json](device-review-2026-11-03/native-simplec-combined.json).
Both SIMPLEC default checks have paired medians within the requested band.
No sample was discarded. Paired medians need not equal ratios of the separate
before/after timing medians.

[timing-summary.json](device-review-2026-11-03/timing-summary.json) and the
compressed paired-run records contain every final sample, including setup and
heap measurements. The before peak in each row comes from that comparison's
baseline cohort, rather than from an unrelated earlier run.

### What failed as a default

The original isolated edge-cache removal was +22.8% on V6, slower in every
pair. Its recomputation path is now optional. The cached coefficient builder,
row, SpMV, residual and smoother sources are exactly those from `8011f3f`;
the host chooses cached or recomputed kernels outside their row loops.

Keeping sequential ordinary strain with cached edges then missed the L2
iteration gate: +23.1% median, range +14.9% to +39.7%. That triggered restoring
the fused vector fit. This is a whole-case observation under contention, not
an isolated causal estimate of the four extra launches. The pre-gate V6/L1
sets, the L2 miss and the first V6 set with cold/mixed Metal compilation are
retained in [pre-gate-timing-summary.json](device-review-2026-11-03/pre-gate-timing-summary.json).
The final-source tables above use a new complete set after the source change.

No new isolated V-cycle-only or walk-only time claim is made. Their combined
default cost is covered by the whole-case pairs. The earlier exact coefficient
row interning experiment still saves zero bytes and remains removed.

### Projected capacity and host mirrors

The conservative hybrid-startup bound remains persistent payload plus
`124N + 8(N+2I) + 16F + 12(F+I) + 40(F-I) + 112`, with a 0.6 GB context/display
reserve. Values are mesh-ratio projections, not measured RTX capacity.

| Mesh / coupling | Startup B/cell, before / default / memory | 8 GB Mcells, before / default / memory | 12 GB Mcells, before / default / memory | 24 GB Mcells, default / memory |
| --- | ---: | ---: | ---: | ---: |
| Ahmed L1 hex / SIMPLEC | 2,502.0 / 2,394.0 / 2,370.0 | 2.96 / 3.09 / 3.12 | 4.56 / 4.76 / 4.81 | 9.77 / 9.87 |
| Ahmed L1 hex / COUPLED | 3,025.6 / 2,899.1 / 2,755.1 | 2.45 / 2.55 / 2.69 | 3.77 / 3.93 / 4.14 | 8.07 / 8.49 |
| V6 polyhedra / SIMPLEC | 3,495.6 / 3,340.0 / 3,316.0 | 2.12 / 2.22 / 2.23 | 3.26 / 3.41 / 3.44 | 7.01 / 7.06 |
| V6 polyhedra / COUPLED | 4,122.8 / 3,948.2 / 3,724.0 | 1.79 / 1.87 / 1.99 | 2.77 / 2.89 / 3.06 | 5.93 / 6.28 |
| Ahmed tet/prism / SIMPLEC | 2,159.0 / 2,068.8 / 2,044.8 | 3.43 / 3.58 / 3.62 | 5.28 / 5.51 / 5.58 | 11.31 / 11.44 |
| Ahmed tet/prism / COUPLED | 2,634.6 / 2,525.8 / 2,411.4 | 2.81 / 2.93 / 3.07 | 4.33 / 4.51 / 4.73 | 9.26 / 9.70 |

Both modes add zero persistent host mirror bytes. The default edge cache is
written directly on the device; memory mode reconstructs it without host
transfers. Strain buffers are device allocations. Existing mesh/CSR mirrors
remain; host `MeshGeometry` payload is still `48F + 32N` bytes besides topology/CSR.
Walk dictionaries remain temporary host construction data. CUDA
uses a separate device pool, so the Mac process peak is not an RTX VRAM peak.
CUDA compilation, actual RTX memory high-water and throughput remain unverified.

L3 remains held while the independent main slow baseline runs; no new L3
numbers are substituted. [l3-hold.json](device-review-2026-11-03/l3-hold.json)
records that gate. Walk coefficients, remaining LSQ geometry, corrected gradients,
ten Coupled `4N` vectors, coarse block state and most SST cell/wall fields remain.
The 2× target remains open.

### Fresh parity and checks

All twelve default case/backend pairs against `8011f3f` completed 40 iterations
and passed literal `cmp` on fluxes.csv, forces.csv, loads.json and result.vtu,
plus residuals.csv after `grep -v '^#' | cut -d, -f1-7`. The six cases are V6
Coupled (`wfgc`), V6 SIMPLEC converted by the baseline, elbow_param
(`compressibility,ke-options`), Ahmed L1 Coupled, converted Ahmed L1 SIMPLEC,
and the native aero SIMPLEC JSON. Each is compared separately on CPU and Metal.

Twelve fresh memory-mode candidates also pass the same comparisons against
the just-produced final default outputs. These reuse the fresh default reference
outputs; they do not substitute the original rejected candidate.
[parity-default.jsonl](device-review-2026-11-03/parity-default.jsonl) and
[parity-memory.jsonl](device-review-2026-11-03/parity-memory.jsonl) retain hashes.

The full quick tier passes: 309 cases, 7,990,306 assertions, zero failures;
37 slow cases excluded by the unchanged suite boundary. The new SST field-byte
test passes 32 CPU/Metal assertions. The Release build and portability lint pass.
No thresholds, references, iteration budgets, skips or expected failures changed.
The slow tier was not rerun: this change uses the byte-identical procedure and
makes no new full-slow or Fluent-accuracy claim.
[checks.json](device-review-2026-11-03/checks.json) binds the checks to this artifact.
