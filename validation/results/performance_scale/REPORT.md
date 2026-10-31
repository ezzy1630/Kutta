# Performance at scale: per-iteration overhead, memory and startup (2026-10-30)

Branch `t3/reduce-cfd-iteration-overhead`. Baseline is `main` at `63287b9`; the
candidate is `c21e59c` unless a row names another commit; `7eea93c` differs from it only in RCM's face renumbering. The machine is an
Apple M4 Pro with 14 CPU cores, 24 GB unified memory, a Metal
`recommendedMaxWorkingSetSize` of 19.07 GB, and Darwin 25.6.

Up to four other lanes' slow tiers shared the GPU and CPUs throughout. Every
timing comparison alternates baseline and candidate runs, and absolute times
include that contention. Raw logs are in `raw/`. The benchmarks are
`bench/bench_overhead.cpp` and `bench/bench_memory.cpp`. The baseline binaries
were built from a `63287b9` worktree with a copy of these bench sources;
`main`'s own sources were unchanged.

## What changed

| Commit | Change | Numerics |
| --- | --- | --- |
| `0174935` | `Device::reduce_many` queues several host reductions behind one wait. Krylov entry checks and the AMG refresh validation use it. Each Krylov batch poll also computes the reliable residual speculatively in dead work vectors, so a candidate exit costs one wait, not two. Solver limits come from a kernel instead of a synchronizing upload. Metal finishes one-threadgroup reductions in one dispatch by replaying `runtime_reduce_into`'s additions. | Each value bitwise identical |
| `a4e1f84` | Metal records runs of small dispatches (at most 2048 rows) and replays them in one `runtime_sequence` threadgroup, with a device barrier between steps. `kernelgen.py` generates the dispatcher. `KUTTA_METAL_SEQUENCE=0` turns replay off. | Same kernels in the same order |
| `1ae5929` | One `KrylovWorkspace` shares work vectors across the momentum, k, ω and pressure solves, which never interleave. | Unchanged |
| `d5d8610` | Boundary LSQ fits share the stencil's internal rows and store only their boundary-face rows. | Same floats, from the same double expression |
| `aec240e` | `validate()` builds labels only on failure, finds repeated nodes in place, and pairs edges with one reused sorted vector instead of a per-cell `std::map`. Wall distance runs on threads with per-cell results identical to the serial loop. | Unchanged |
| `c21e59c` | RCM permutes the face CSR arrays directly, giving the same order and patch ranges as `detail::pack`, instead of building one `detail::Face` per face twice. | Unchanged |
| `e17f7ea`, `e509246`, `7eea93c` | `bench_memory` gains a mesh-file mode, in-place transfers, a peak footprint figure, and a GCC include fix. | Bench only |

## Bitwise parity

Fields are hashed with FNV-1a after a fixed number of iterations
(`bench_overhead parity`, `raw/parity.log`). Every hash is identical between
`63287b9` and `7eea93c`. They were rechecked at `c21e59c` before commit.

| Case | CPU | Metal |
| --- | --- | --- |
| k-ε channel, 20 cells, 600 steps | `4263a6c55de206aa` | `72586688e49e9606` |
| SST tet channel, 16×10×4, 40 iterations | `1ed4c187c03fa89b` | `253bfa09c0471f56` |
| Laminar lid box, 12×12×6, 30 iterations | `51924bc0604c9eff` | `d34d0604f98f7a95` |
| Laminar lid box, 54×54×48 (140K cells), 3 iterations | `7132b2ddd3a62835` | `1dda4da81f0c54de` |

For the Fluent mixing elbow (`kutta solve --approximate wfgc`), every main and
candidate run produced the same residual history (timing column removed), forces,
fluxes and `result.vtu` (`raw/elbow_cli.log`):

| Backend | Iterations | Runs | residuals | forces | vtu |
| --- | --- | --- | --- | --- | --- |
| Metal | 150 | 3 main + 3 candidate | `b31cc5c309a5b305` | `62620087e01e830d` | `9272b5efcd6fcf37` |
| CPU | 30 | 3 main + 3 candidate | `66e236eb661851e9` | `c64a075e3fb91892` | `6e015e0383fd20e5` |

The CPU elbow hashes were rechecked at `c21e59c`. In the slow tier, the V8
channel's printed defect, change, τ and U error match main to every printed digit
on all backends and sizes (`raw/slow_*.log`).

## Per-iteration overhead

### V8 realizable k-ε channel on Metal

The test's own timing comes from `raw/slow_baseline.log` and
`raw/slow_candidate.log`. The bench rows are three alternating runs of
`bench_overhead ke-channel metal N counted` (`raw/ke_channel_metal.log`).
Steps and results are identical in every run.

| Cells | Steps | main, slow tier (s) | candidate, slow tier (s) | main, bench (s) | candidate, bench (s) |
| --- | --- | --- | --- | --- | --- |
| 20 | 600 | 13.6 | 6.7 | 13.8 / 13.4 / 14.4 | 6.6 / 6.9 / 7.3 |
| 40 | 1900 | 48.1 | 22.3 | | |
| 80 | 5700 | 180.7 | 78.6 | 173.2 / 140.9 / 158.6 | 102.7 / 71.9 / 172.8 |

V8's 180 s budget is unchanged. The candidate converged in all four 80-cell
runs. The slowest, 172.8 s, came while the machine was most contended: the main
run just before it took 158.6 s. In that run, ~50 waits at ~0.5 ms per wait set
the per-step cost.

Per-step counts at 80 cells, identical across runs:

| | main | candidate |
| --- | --- | --- |
| host waits | 100.96 | 49.75 |
| host reductions requested | 54.75 | 55.75 |
| status downloads | 43.2 | 43.2 |
| uploads | 3.01 | 0.012 |
| logical launches | 3915 | 3919 |

Replay reduces GPU dispatches. In an earlier 20-cell profile, about 3.0M
dispatches became about 28K.

### Dispatch and wait latency on Metal

Microseconds per operation from three alternating `bench_overhead latency metal`
runs (`raw/latency_metal.log`):

| Operation | main | candidate |
| --- | --- | --- |
| queued dispatch | 0.93–1.09 | 0.85–0.86 |
| queued device reduction | 4.37–4.97 | 1.62–1.76 |
| queued scalar kernel | 1.21–1.98 | 1.05–1.12 |
| queued SpMV | 1.66–1.84 | 1.26–1.31 |
| queued mixed dispatch | 2.09–2.51 | 1.30–1.37 |
| BiCGStab Jacobi update (one batch poll) | 79–91 | 40–43 |
| dispatch + host reduction / synchronize / upload | 84–98 | 82–98 |

The last row is the cost of one host wait. It is unchanged: removing waits is
what pays.

### Mixing elbow, Metal, 150 iterations

Three alternating runs (`raw/elbow_cli.log`):

| | main | candidate |
| --- | --- | --- |
| solve wall time (s) | 2.60 / 2.56 / 2.58 | 1.97 / 1.97 / 1.97 |
| setup to first iteration (s) | 1.81 / 1.81 / 1.80 | 1.57 / 1.57 / 1.56 |
| process real time (s) | 4.65 / 4.61 / 4.62 | 3.76 / 3.75 / 3.75 |

CPU solve times are unchanged within noise: 30 iterations took 2.75–2.85 s on
main and 2.79–2.92 s on the candidate. The CPU 20-cell V8 run took 40 ms on both
(`raw/ke_channel_cpu.log`).

At 1M and 2M cells, iteration time is unchanged (0.18–0.20 s and 0.39–0.41 s
on Metal): those meshes are bandwidth-bound, not overhead-bound.

## Memory

`bench_memory` runs an SST steady solve and wraps the device to classify every
allocation by length (`raw/memory_*.log`). Device bytes are identical on CPU and
Metal, and the device peak equals the steady value.

| Mesh | Cells | Faces/cell | main B/cell | candidate B/cell | Change |
| --- | --- | --- | --- | --- | --- |
| Hex box 160×128×100 | 2,048,000 | 3.02 | 2169.8 | 1633.3 | −24.7% |
| Hex box 80×128×100 | 1,024,000 | 3.03 | 2172.1 | 1636.5 | −24.7% |
| Fluent mixing elbow (polyhedral) | 17,822 | 5.14 | 3023.6 | 2262.5 | −25.2% |

Device bytes per cell by array length at 2M hex. On a hex mesh, the six LSQ
inverse planes have the same length as the cell-face entries.

| Array length | main arrays | main B/cell | candidate arrays | candidate B/cell |
| --- | --- | --- | --- | --- |
| cells | 191 | 764 | 146 | 584 |
| cell-face entries (and 6×cells) | 35 | 840 | 20 | 480 |
| faces | 32 | 387 | 32 | 387 |
| matrix entries | 3 | 83 | 3 | 83 |
| internal faces | 5 | 60 | 5 | 60 |
| boundary faces | 36 | 7 | 54 | 10 |
| other (AMG coarse levels and small) | 148 | 29 | 152 | 29 |

On the elbow, the candidate's entries are 520 B/cell (13 arrays), faces
658 B/cell (32), cells 584, LSQ inverse planes and AMG levels 208, and the matrix 129.

### Largest mesh on this Mac

On Apple unified memory the limit is the process footprint, not device bytes:
it adds the host mesh, geometry and setup temporaries to the device buffers.

| Hex box (Metal) | main peak footprint | candidate peak footprint |
| --- | --- | --- |
| 1,024,000 cells | 4.13 GB | 3.51 GB |
| 2,048,000 cells | 8.14 GB | 6.90 GB |
| 3,072,000 cells | not run | 10.83 GB (device 5.01 GB, 0.57 s/iteration) |

The largest mesh measured here is 3.07M hex cells, run while only about
9.3 GB was free on the shared machine (`raw/memory_3m_metal.log`; swap grew
0.2 GB). A 4M run would have needed about 14 GB, so it was not attempted while
other lanes were running.

On an otherwise idle 24 GB Mac with about 19 GB available to one process, the
measured peak footprint (about 3.5 KB/cell on the candidate, 4.0 KB/cell on main)
gives roughly 5.4M hex cells on the candidate and 4.8M on main. Polyhedral meshes
fit fewer. The Metal working-set limit of 19.07 GB covers device buffers alone
(11.7M cells at 1633 B/cell) and does not bind first.

Footprint readings at the same size vary by up to about 1 GB between sessions.
macOS keeps freed malloc regions resident: after RCM at 2M cells, malloc held
1.59 GB for 0.36 GB in use, and `malloc_zone_pressure_relief` returned none of
it. Metal buffers cannot reuse that memory. `c21e59c` removes the largest
source, RCM's per-face node vectors, which lowered the footprint after RCM from
2.27 to 1.26 GB at 2M.

### Extrapolation to a 16–24 GB RTX card (not measured)

The CUDA device allocates through the same `Device` API, so device bytes per
cell carry over. The host mesh and setup temporaries live in system RAM.
Assuming about 0.6 GB for the CUDA context and display, roughly 15 GB is usable
on a 16 GB card and 23 GB on a 24 GB card.

| Mesh type | B/cell (main → candidate) | 16 GB (main → candidate) | 24 GB (main → candidate) |
| --- | --- | --- | --- |
| Hex-like (3.0 faces/cell) | 2170 → 1633 | 6.9M → 9.2M cells | 10.6M → 14.1M cells |
| Polyhedral, elbow-like (5.1 faces/cell) | 3024 → 2263 | 5.0M → 6.6M cells | 7.6M → 10.2M cells |

These figures are for SST steady with the default AMG pressure solver. A
production FSAE mesh with prism layers sits between the two rows, depending on
its polyhedral fraction.

## Startup

### Mixing elbow, `kutta solve` to the first iteration

Phase timing on Metal, from temporary probes that are not committed:

| Phase | Time |
| --- | --- |
| hash, read case, setup | 0.03 s |
| `read_mesh` | 0.25 s |
| SteadySolver construction | 0.23 s |
| hybrid initialization | 1.09 s |

Hybrid initialization runs 10 AMG-PCG potential solves of 86–300 iterations at
about 0.44 ms per iteration. Those iteration counts are numerics and must not
change, so the speedup comes from cheaper iterations.

| Backend | main setup (s) | candidate setup (s) |
| --- | --- | --- |
| Metal | 1.81 / 1.81 / 1.80 | 1.57 / 1.57 / 1.56 |
| CPU | 6.74 / 6.75 / 7.01 | 6.70 / 6.78 / 6.96 |

On CPU, startup is the hybrid solve on the single-threaded CPU backend and is
unchanged.

### Large mesh (hex box, Metal)

At 1M cells, solver setup took about 6 s. These phases were profiled on `d5d8610`, whose mesh and setup code was still identical to main's:

| Phase | Time |
| --- | --- |
| RCM | 2.3 s |
| `upload_mesh` | 1.85 s |
| geometry | 0.43 s |
| AMG | 0.34 s |
| stencil | 0.32 s |
| LSQ fits | 0.35 s |

RCM and `upload_mesh` time was mostly `validate()` (0.82 s per call, called 4–5
times per solve) and the serial wall-distance search (0.86 s). Alternating runs
(`raw/memory_startup_final.log`):

| Cells | main mesh + setup (s) | candidate mesh + setup (s) |
| --- | --- | --- |
| 1,024,000 | 2.32 + 6.18, 2.15 + 5.61 | 1.55 + 2.71, 1.52 + 2.62 |
| 2,048,000 | 4.34 + 12.46, 4.77 + 12.82 | 3.31 + 6.11, 3.38 + 6.28 |

The mesh column is the bench's box generator, which also calls `validate()`.
Wall distance in this box has only two walls. On FSAE meshes, where wall faces
are a large share, the threaded search matters more.

## CUDA

At `7eea93c`, which carries every CUDA-relevant change, a compile-only build
ran in the local `nvidia/cuda:13.0.2-devel-ubuntu24.04` container (nvcc
13.0.88, GCC 13.3, architectures 75/86/89, `--fmad=false`). All 148 targets
built (`raw/cuda_build_7eea93c.log.gz`). The six warnings are in files this
branch does not touch.

That build found a missing `<cmath>` in `bench_memory.cpp`, fixed in
`7eea93c`. CUDA execution was not possible: there is no NVIDIA GPU here.
`CudaDevice::host_reductions` and the CPU==CUDA parity under it are
unexecuted.

## Tests

- Quick tier at `c21e59c`: `ctest -LE slow` passed (kutta-tests and
  lint-portability). Run directly with fixtures, it reported 225 cases and
  7,413,874 assertions, 0 failed (`raw/quick_*.log`). The same counts were
  recorded at `7eea93c`.
- `test_mesh.cpp` gains two `validate()` rejection subcases: the first repeated
  node is named, and a hanging node is caught by edge pairing alone.

- Full slow tier, using the same command, at `63287b9` (baseline) and
  `c21e59c` (candidate), with `KUTTA_FLUENT_EXAMPLES` set. Both runs reported
  24 cases (18 passed, 6 failed) and 10,816 passed / 68 failed assertions.
  `tools/forensics/compare_slow_assertions.py` found 0 regressions, 0
  recoveries, 0 missing and 0 new assertions (`no_new_failures: true`;
  `raw/slow_comparison.json`, XML in `raw/slow_*.xml.gz`).
- Both runs fail the same 6 cases: V6 mixing elbow (29), V12 fan import (1), V2
  tet refinement (6), V5 plate (6), V5 bump (9) and V7 SA (17).
- All 8 changed expanded values come from V7 SA bump level 2 on CPU. It stops
  on its 600 s wall-clock budget (`budget=1` in both logs), and the candidate
  reached 2,600 steps where main reached 2,500. Every line not limited by the
  budget matches.
- Wall time for the slow tier: 5,039 s on the baseline, 4,390 s on the
  candidate.

## Open limits

- **Host waits still set small-mesh time under contention.** About 50 waits
  per V8 step remain, mostly Krylov status polls after candidate exits (about
  14 restarts per solve). The restart decision compares a residual finished in
  double on the host, which Metal cannot reproduce on the device. Fewer waits
  need pipelined restarts: queue the next restart and batch behind the
  speculative residual, keep `x` for rollback, and finish host reductions
  asynchronously from retained partials. That is the next step if V8 must stay
  far below 180 s under heavy sharing.
- **AMG-PCG at 17.8K cells is near the dispatch floor.** About 0.44 ms per
  iteration, roughly 60 dependent dispatches, and the fine level exceeds the
  2048-row replay threshold. Going further needs kernel fusion.
- **Further memory savings are measured but not implemented**, because they
  touch transport and reconstruction kernels another lane owns:
  - Sharing interior-cell LSQ inverse planes across the five fits would save
    about 85–100 B/cell.
  - Deriving `face_delta` from `owner_to_face` and `neighbour_to_face` would
    save 72 B/cell on hex and more on polyhedral meshes.
  - Each `TransportEquation` keeps its own `zero_` field (4 B/cell each).
- **Host footprint on macOS.** Allocator retention makes the footprint noisy.
  `make_stencil` and `upload_mesh` still build `Vec3d` temporaries that could be
  converted to float as they are produced.
- **CUDA and Windows/RTX are not executed here.** The figures above are
  extrapolations.
