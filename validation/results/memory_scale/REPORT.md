# Peak memory of a full `kutta solve`

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
