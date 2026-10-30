CSR operations, Krylov iteration, smoothing, transfer and coarse-operator refresh
use Device buffers and generated kernels. The host builds aggregation structure
once. Solve/refresh do not download or traverse field arrays. Iteration scalars,
nullspace sums and optional K-cycle corrections stay on the device; the host
polls a three-int status buffer between batches. Kernels gather, own their output
row/entry, use float arithmetic and require no atomics.

A normal caller uploads one matrix and reuses the solver:

```cpp
auto a = kutta::upload(dev, host_matrix);
auto solver = kutta::make_solver(kutta::SolverKind::AmgPcg, dev, host_matrix,
    {.nullspace = kutta::NullSpace::Constant}); // connected Neumann pressure
// Use NullSpace::None (also the convenience overload's default) for pinned SPD pressure.
auto stats = solver->solve(dev, a, x, b, {.rel_tol = 1e-6});
// Assemble/upload new a.val with the same pattern and nullspace.
solver->refresh(dev, a);
stats = solver->solve(dev, a, x, b, {.rel_tol = 1e-6});
```

The original three-argument make_solver signature remains available. SolverOptions
adds an explicit NullSpace::None/Constant choice; AmgOptions has the same field.
The equation's boundary conditions determine this choice. Float row sums do not:
random coefficients and a diagonal assembled as -sum(float offdiag) need not
sum to exactly zero. Debug builds check a requested constant nullspace against
symmetry and row sums within 8 float epsilons times the absolute row sum.
This is a consistency diagnostic, never an automatic mode selection.

`SolverOptions::structure` and `AmgOptions::structure` default to
`OperatorStructure::General`. FV pressure explicitly selects
`SymmetricMMatrix`: the caller guarantees symmetric nonpositive off-diagonals
and diagonal dominance of the stored float rows, including on refresh. Setup
validates this declaration. Nullspace selection remains independent and explicit.

Constant projection removes the mean from the initial guess, residuals and
preconditioned corrections, and from coarse RHS/corrections at every level.
Reported residuals remain norms of the full stored operator's b-A*x. An
incompatible nonzero-mean RHS is rejected. Rounding can still set a residual
floor. Multiple disconnected pressure nullspaces or a nonsymmetric nullspace
need problem-side constraints; one global constant projection is insufficient.

The solver belongs to its original Device, which must outlive it. Subsequent
matrices must retain the exact original pattern and declared nullspace.
Dimensions are checked at the API boundary; the caller owns the immutable-pattern
invariant and must refresh after value changes. Factories initialize diagonals
and coarse values from the host matrix without uploading another fine CSR.
The caller's initial device values must match those host values. Pattern-only
upload initializes values to zero for assembly. SpMV/residual outputs must not
alias inputs.

AMG defaults and aggregation
----------------------------

```cpp
kutta::AmgOptions options; // V-cycle, 5 pre + 5 post sweeps, 3 pairing passes
options.cycle = kutta::AmgCycle::K; // optional; useful on the serial CPU backend
options.pre_sweeps = options.post_sweeps = 2;
auto solver = kutta::make_amg_solver(kutta::SolverKind::AmgPcg, dev, host_matrix, options);
```

The default is a fixed V-cycle. Its five pre/post sweeps were selected with the
Metal measurements below and the <=15-iteration manufactured 3D Poisson gate.
Two-sweep V is faster at 128^3 for random forcing, but failed that precision gate:
32^3 took 17 iterations and 64^3 stagnated above 1e-8. Five sweeps pass both.
The smoother is damped L1-Jacobi, with damping 2/3 and
`D_l1 = max(|a_ii|, (|a_ii| + sum_{j != i}|a_ij|)/2)`.
Equal pre/post sweeps retain the SPD V-cycle property and are required for
AMG-PCG. AMG setup and refresh reject negative diagonals with
`AMG: requires positive diagonals`. Zero/invalid diagonals are also rejected,
except the scalar zero Neumann operator. Jacobi-BiCGStab supports negative
nonzero diagonals, including the regression that negates every seventh row.

Aggregation matches the strongest unassigned connection, with sorted index order
breaking ties and a row-relative threshold of 0.25. Strength is
`max(|a_ij|, |a_ji|)`, using CSR transpose lookup during setup. Incoming-only
structural edges are included. A large already-assigned upstream edge can no
longer hide all available downstream connections. Three matching passes produce
aggregates of up to eight cells and follow the strong direction on stretched
cells. Greedy matching still depends on order; symmetric strength prevents the
one-level failure, not every permutation-dependent difference in quality.

Isolated diagonal rows have aggregate -1, receive smoothing, and are omitted
from restriction/prolongation and the coarse grid. Only connected rows count
in the 0.85 coarsening-stall check. Defaults stop at <=24 rows or 24 levels;
the coarsest solve is 32 device sweeps, fused into one dispatch at <=24 rows.
A connected graph that stalls falls back to smoothing at its current size.

Interpolation is piecewise constant, restriction its transpose. Each coarse
entry owns a sorted list of contributing fine entries, implementing
`A_c = P^T A P` with deterministic compensated sums. Initial host setup uses
the same float summation order as device refresh. Intermediate pairing passes
build only the operator, without entry maps; each intermediate Coarsening is
released before the next is built. Refresh rebuilds only values and diagonals,
retaining aggregates/maps. Large coefficient changes can reduce hierarchy
quality; rebuild explicitly if needed.

For `SymmetricMMatrix`, opposite coarse edges gather the same fine edges in
canonical order, preserving exact stored symmetry. Each coarse diagonal is
reconstructed as the upward-rounded sum of its stored edge magnitudes and
aggregated fine-row surpluses. This is the graph form of `P^T A P`, rounded to
preserve the M-matrix property; it does not add an arbitrary diagonal shift.
Host setup/intermediate pairing and device refresh use the same upper-bound
float accumulator. General matrices retain compensated entrywise Galerkin sums.

A K-cycle uses up to two recursive inner Krylov corrections, with an adaptive
0.25 residual-reduction exit. Inner SPD corrections use conjugate directions;
nonsymmetric corrections use GCR. Outer PCG uses a one-direction flexible-CG
recurrence for this variable preconditioner. Both cycles remain selectable.

Reductions, exits and storage
----------------------------

SolveControl uses `max(abs_tol, rel_tol * initial_residual)`. Norms are absolute
L2 norms of b-A*x, including nonzero initial guesses. Device reductions return
double when requested by the host. Asynchronous `sum_into`, `dot_into` and
`dot2_into` write normalized `(hi, lo)` float pairs to a caller-owned device
buffer. CPU/CUDA accumulate in FP64; Metal uses compensated pairs including FMA
product residuals. Both stages have fixed reduction order and no atomics. Scalar
kernels retain both words through coefficient arithmetic before rounding the
coefficient to float. Pairs have about 48 significant bits and float's exponent
range; solver fields still need O(1) scaling on Metal. Residual kernels use
Neumaier summation and Dekker product correction in float; above 1e30 the product correction falls back to a rounded product to
avoid splitting overflow. Backends must preserve expression order and disable
fast-math reassociation.

All solvers honor `residual_check_interval` (default 8). PCG/BiCGStab evaluate
recursive norms and coefficient/breakdown guards on the device every update,
then poll status once per batch. A convergence candidate, divergence candidate,
or 50-update reliable replacement freezes updates until the host verifies the
full stored operator's residual. Frozen dispatches do not count as iterations;
large stencil/preconditioner sweeps skip their memory work. A failed candidate
within ten times the requested target switches recovery to one-update batches,
avoiding unused queued cycles near a float storage floor. This changes polling,
not the target or acceptance test. BiCGStab fuses t.t and t.s in one reduction.
The optional K-cycle also evaluates its adaptive exit and scalar guards on the
device. Jacobi still computes its residual field every sweep and verifies its
norm every batch and on the final sweep. The existing three-field aggregate
initialization remains valid. Iterations count attempted Krylov updates or
Jacobi sweeps; a satisfied initial tolerance takes
zero updates. A verified residual above 1000 times the initial residual, or a
nonfinite residual, exits with converged=false. Numerical exits restore the
best *verified* x, including the initial guess, rather than returning a worse
or nonfinite iterate. One recovery vector pays for this guarantee.

Negative curvature, BiCGStab breakdown/unrepresentable scalars, or twelve
stagnant reliable-residual replacements also return non-converged stats.
Replacement occurs at a convergence candidate or every 50 iterations and
restarts the recurrence. All numerical exits recompute the actual final
residual, including max_iter=0. Invalid CSR/control, nonfinite inputs,
incompatible Neumann forcing and invalid diagonals throw Error. The one-row
7*x=1 test prevents false convergence below float representability.

Device-scalar batching measurements (2026-10-29)
------------------------------------------------

Machine: Apple M4 Pro, 14-core CPU, 24 GB, macOS 26.6.2 (25G83), AppleClang 21,
Release CMake/Ninja build. Same deterministic benchmark RHS and matrix family.
Before: `bdbb18a77b8feebb75b8a1af54f6d1e71792c8d4`; after:
`0430a94d64921824833e96d12ea9233c53e7d363`, descended from that base.
Times below are **best of three** sequential runs per phase (setup/refresh
excluded). Raw records, including every sample, full residuals, source identities
and blocking/asynchronous reduction counts, are in
[benchmarks-device-scalars.jsonl](benchmarks-device-scalars.jsonl).
No convergence tolerance or test threshold changed. Commands for each sample:

```sh
cmake -S . -B build -G Ninja && cmake --build build
./build/bench_linalg
./build/bench_linalg --small 64 --large 64 --wind 20
```

Both commands run CPU and Metal. Default V5, rel_tol=1e-6, abs_tol=1e-12,
max_iter=2000, residual_check_interval=8. The default command also records the
anisotropic 64^3 case. Benchmark solve times include the final synchronization.
Baseline sample 1 overlapped the initial test-suite run; samples 2/3 and the
final three after samples had no other task-owned checks running.

| Backend | Case / solver | Before ms | After ms | Speedup | Iterations before / after | Final relative residual before / after |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| metal | 64^3, AmgPcg | 35.65 | 19.43 | 1.83x | 15 / 15 | 7.86701089e-07 / 7.86701089e-07 |
| metal | 128^3, AmgPcg | 229.11 | 215.35 | 1.06x | 21 / 21 | 8.61293301e-07 / 8.61293301e-07 |
| metal | 64^3, wind 20, BicgstabJacobi | 94.09 | 38.95 | 2.42x | 69 / 69 | 6.52857681e-07 / 6.52857681e-07 |
| metal | 64^3, wind 20, AmgBicgstab | 31.34 | 18.26 | 1.72x | 7 / 7 | 7.49268702e-07 / 7.49268736e-07 |
| metal | 64^3, PcgJacobi | 170.60 | 79.19 | 2.15x | 278 / 278 | 6.83041863e-07 / 6.83041863e-07 |
| metal | 128^3, PcgJacobi | 1016.03 | 830.00 | 1.22x | 640 / 640 | 7.16956407e-07 / 7.16956407e-07 |
| cpu | 64^3, AmgPcg | 653.82 | 717.33 | 0.91x | 15 / 15 | 7.86701089e-07 / 7.86701089e-07 |
| cpu | 128^3, AmgPcg | 7521.80 | 8099.05 | 0.93x | 21 / 21 | 8.61293301e-07 / 8.61293301e-07 |
| cpu | 64^3, wind 20, BicgstabJacobi | 196.26 | 223.68 | 0.88x | 69 / 69 | 6.52857681e-07 / 6.52857681e-07 |
| cpu | 64^3, wind 20, AmgBicgstab | 609.88 | 691.53 | 0.88x | 7 / 7 | 7.49268702e-07 / 7.49268736e-07 |
| cpu | 64^3, PcgJacobi | 386.57 | 452.23 | 0.85x | 278 / 278 | 6.83041863e-07 / 6.83041863e-07 |
| cpu | 128^3, PcgJacobi | 7258.99 | 7977.55 | 0.91x | 640 / 640 | 7.16956407e-07 / 7.16956407e-07 |

All these cases converge. Iteration counts are unchanged. Pressure and
Jacobi-BiCGStab final residuals match at the recorded precision. AMG-BiCGStab's
relative residual changes by 3.5e-14 absolute (4.6e-8 relative), consistent with
retaining about 48 rather than 53 coefficient bits in the device pair format.
The final norm is still freshly computed from b-A*x with the host double-return
reduction API. The exact CPU/Metal comparison regressions remain bitwise equal.

The serial CPU oracle is 8-17% slower on these cases: it executes the same scalar
state, guard kernels and frozen dispatch loops, without a GPU wait to remove.
Metal gains most on the scalar-heavy Jacobi Krylov solvers. At 128^3 the default
V-cycle already had no internal waits, so removing the 63 ordinary scalar waits
cannot account for most of its 229 ms baseline. The final implementation also
skips frozen memory sweeps, finishes reductions with a parallel tree, and fuses
coarse sweeps and smoothing. AMG smoothing uses alternating uniquely owned buffers and
the identical compensated residual arithmetic, removing one memory pass and
one dispatch per sweep. The <=24-row coarsest solve uses one owning thread and
local Jacobi arrays, with exactly the same 32 simultaneous sweep updates.
Larger coarsest levels retain the device sweep fallback.

The anisotropic forcing still does not converge at 1e-6: AMG-PCG takes 38
iterations and returns 1.51462927e-6 before and after. Metal improves from
69.63 to 57.92 ms; CPU changes from 1716.40 to 1808.91 ms. Recovery preserves
the verified best iterate and never reports that storage floor as convergence.

Synchronization inventory (Device call trace plus instrumented solve counts):

| Ordinary iteration or cycle | Before: blocking scalar reads | After: blocking scalar reads |
| --- | ---: | ---: |
| Fixed-preconditioner PCG iteration | 3 (p.Ap, r.r, next r.z); initial r.z once | 0; device coefficients and norms |
| Flexible PCG outer recurrence | 4, plus recursive K-cycle work | 0 |
| Full BiCGStab iteration | 6 (rho, rhat.v, s.s, t.t, t.s, r.r); s-only exit uses 3 | 0; t.t/t.s fused |
| Jacobi sweep | 1 every 8 sweeps, plus final check | same |
| V-cycle, NullSpace::None, including coarsest | 0 | 0 |
| V-cycle constant projection, L levels | 2(L-1), plus 1 fine correction sum at Krylov boundary | 0 |
| Optional K-cycle inner correction, excluding recursion/projection | 4 at early exit, up to 7 for both corrections | 0 |
| Coarsest 32-sweep solve | 0; 64 residual/update dispatches | 0; 1 dispatch at <=24 rows |

Alpha/beta/omega were host calculations fed by those reads; their ratios, signs,
zero/NaN/range guards now run on device pairs. Only the outer batch boundary
polls the three-int phase/count/short-step state. Convergence candidates freeze
updates, so polling at iteration 16 can truthfully report 15 attempted updates.
Reliable verification and best-x restoration retain their actual-residual
semantics, including max_iter=0, negative curvature and numerical breakdowns.

Whole-solve blocking calls below include input/final norm reads, status polls,
and the one boundary upload of scalar thresholds; they exclude setup, refresh
and the benchmark's final timing synchronize. CPU counts are the same API costs,
not GPU waits. Metal's bounded command-buffer backpressure may also wait after
8 pending 256-dispatch buffers; ordinary V-cycle batches stay below that bound.
Long K-cycle batches can still reach that resource-retention bound.

| Case | Before | After host norms | After status polls | After boundary uploads | After total |
| --- | ---: | ---: | ---: | ---: | ---: |
| 64^3 AmgPcg | 50 | 5 | 2 | 1 | 8 |
| 128^3 AmgPcg | 68 | 5 | 3 | 1 | 9 |
| 64^3 BicgstabJacobi wind 20 | 420 | 6 | 10 | 1 | 17 |
| 64^3 AmgBicgstab wind 20 | 44 | 5 | 1 | 1 | 7 |
| 64^3 PcgJacobi | 845 | 11 | 63 | 1 | 75 |
| 128^3 PcgJacobi | 1938 | 18 | 210 | 1 | 229 |

Verification on the implementation above:

- `cmake -S . -B build -G Ninja && cmake --build build && ./build/kutta-tests`:
  **127 tests / 1,439,199 assertions pass on CPU and Metal**.
- `ctest --test-dir build --output-on-failure`: both the full suite and
  `lint-portability` pass. Initially the untouched baseline had one missing-file
  inventory failure; all external Fluent examples were fetched from their
  documented source and verified against SOURCE.md hashes. No acceptance test
  was weakened or excluded from the native run.
- New tests cover device sum/dot/fused-dot precision, empty reductions,
  cancellation, scratch reuse and deterministic repeats; exact zero-omega and
  zero-shadow-product BiCGStab breakdowns at batch sizes 1/8/31; zero/NaN/Inf
  rho/omega guards; K/V-cycle and Neumann host-call budgets; and exact fused
  coarse/smoothing equivalence with ordinary Jacobi sweeps.
- Docker `nvidia/cuda:13.0.2-devel-ubuntu24.04` (local arm64 image) compiles the
  entire project with `-DKUTTA_ENABLE_CUDA=ON`,
  `-DCMAKE_CUDA_ARCHITECTURES="75-real;86-real;89-real;120"` and
  `-DKUTTA_ENABLE_HDF5=OFF`; the focused CPU-oracle selection passes 14 tests /
  155 assertions. The CI configure with HDF5 enabled stopped in pre-existing
  `case/hdf5.hpp`: Ubuntu's HDF5 lacks `H5L_info2_t`/`H5Lget_info2`. The optional
  reader was disabled for the CUDA compile check; its source was not changed.
  CUDA runtime correctness/performance remains unverified: no NVIDIA GPU here.
- `LinearSolver` public signatures, Buffer element types and convergence
  tolerances are unchanged. No files under `src/kutta/fv/**` were changed, and
  the other branch's FV flow benchmark was not available for measurement.

Historical measurements before device-scalar batching
----------------------------------------------------

The following measurements and verification counts predate the device-scalar
work described above. The current comparison is recorded separately above.

Measured reduction/allocation regressions in that earlier work, before -> after:

| Check | Before | After |
| --- | ---: | ---: |
| 64^3 random pressure, solve reductions (K2 -> default V5) | 726 / 10 iterations | 50 / 15 iterations |
| 128^3 random pressure, solve reductions (K2 -> default V5) | 1308 / 10 iterations | 68 / 21 iterations |
| Ordinary pinned-pressure PCG iteration, excluding setup/final checks | ~73 at 64^3 with K | 3 with V |
| BiCGStab, four full iterations | 32 reductions | 28 reductions |
| Jacobi, seventeen sweeps | 21 reductions | 7 reductions |
| Neumann AMG-PCG, one 8x8 iteration including all checks | 21 reductions | 15 reductions |
| Divergence regression, PCG final residual | 7.49e6 from initial 1.41 | 1.41; original x restored after one update |
| Upwind wind 5/20/100, increasing order | 1 level, complexity 1.000 | 6 levels at 64^3, both directions |

The three-reduction ordinary iteration count is for NullSpace::None. Constant
projection adds level-dependent sum reductions. Fine corrections are projected
only once at the Krylov boundary, and PCG does not build an unused direction
when max_iter has been reached. A device-resident compensated sum consumed by
the shift kernel would be needed to remove the remaining projection round-trips.

Device currently exposes only scalar sum/dot/max_abs reductions. Calling dot
twice cannot fuse host round-trips. Recommended additive API:
`dot2(x, y, z, w, n) -> std::array<double, 2>`, returning both dot products
from one ordered reduction/readback, with the same compensated precision.
This would fuse (t.t, t.s) and (rhat.r, r.r) in BiCGStab and reduce six ordinary
iteration reductions to four. A default implementation could call dot twice
for backend compatibility; a Metal override must actually fuse both channels.
No Device API or backend implementation was changed in this work.

At 128^3, AmgPcg solver-owned device payload is **132.134 -> 84.705 bytes per
fine row** (36% less); peak setup payload is **195.759 -> 84.705 bytes/row**
(57% less). These are counted allocation requests on both CPU and Metal, not
process RSS or driver-resident memory. They exclude the caller's fine CSR,
x/b, benchmark-only value backup, host arrays and backend reduction scratch.
Level 0 now owns only L1 and residual vectors, versus ten vectors. PCG owns
four working vectors plus recovery x, versus eight working vectors. Default
V coarse levels allocate only their smoother workspaces; K workspaces are
allocated only when K is requested. The duplicate fine CSR upload is gone.

bench_galerkin also measures a deterministic coarse-row gather that eliminates
the fine-nnz entry map and its offsets. On 128^3 eight-cell aggregates it saves
31.266 bytes per fine row at the first transfer, but is too slow to adopt:

| Galerkin first-transfer refresh | CPU ms | Metal ms |
| --- | ---: | ---: |
| Entry map (retained default) | 11.167 | 0.534 |
| Coarse-row gather (benchmark alternative) | 56.930 | 8.016 |

Outputs are bitwise identical. Each time is the median of five batches of
twenty queued launches plus synchronization, divided by twenty. This alternative
is 5.1x slower on CPU and 15.0x slower on Metal, so production retains the entry
map. The alternative kernel and benchmark remain available to reproduce this
memory/performance tradeoff.

Measured CPU and Metal evidence
------------------------------

Apple M4 Pro, 24 GiB, macOS 26.6.2 (25G83), AppleClang 21.0.0.21000101,
Release. Repeated timing/benchmark source:
`47a6987828369a98f6c1544b172a23595b085a51` on t3/aggregation-amg-solver.
The original baseline is `11f1533`; allocation baseline was instrumented at
`af06309` before workspace/factory changes. Raw records, sample numbers,
source identity and convergence flags are in benchmarks-2026-10-29.jsonl.
benchmarks-cpu.jsonl remains the historical CPU-only measurement.

```sh
cmake -S . -B build -G Ninja && cmake --build build
./build/kutta-tests
./build/bench_linalg
./build/bench_linalg --backend metal --cycle k --sweeps 2
./build/bench_linalg --backend metal --cycle v --sweeps 2
./build/bench_linalg --backend metal --small 64 --large 64 --wind 20
# --ordering 1 reverses rows/columns/RHS; --ordering 2 uses a deterministic permutation.
./build/bench_linalg --backend metal --small 64 --large 64 --wind 20 --ordering 2
./build/bench_linalg --tolerance 1e-5
./build/bench_galerkin
```

Tables use median wall times from three sequential runs, except the explicitly
single-sample 1e-5 precision check and the gather batches above. The deterministic
random RHS retains low modes and is permuted with the matrix for ordering tests.
rel_tol=1e-6, abs_tol=1e-12, max_iter=2000. Setup includes aggregation,
allocation/upload and initial diagonals/coarse values, excluding the caller's
fine matrix upload. Refresh excludes value upload and includes synchronization;
solve also synchronizes. Timings are subject to machine/driver variability;
raw samples preserve the spread. These are reusable solves, not setup+solve
speedups, and do not establish performance on a 20-million-cell mesh.

Metal cycle selection (AMG-PCG, solve milliseconds; iterations in parentheses):

| Case | K, 2+2 sweeps | V, 2+2 sweeps | Default V, 5+5 sweeps |
| --- | ---: | ---: | ---: |
| 64^3 | 157.85 (10) | 54.99 (20) | **38.22 (15)** |
| 128^3 | 306.50 (10) | **173.04 (28)** | 232.09 (21) |
| 64^3, 1:1000 | 286.15 (38) | 111.41 (43) | **70.42 (38)** |

V5 is 4.1x/1.3x faster than K2 on the two isotropic grids. V2 is 1.8x faster
than K2 at 128^3, but does not meet the manufactured 1e-8 gate. V5 is the
common default that passes that gate and wins against K2 on all measured
Metal cases; a hybrid K policy was therefore unnecessary. The anisotropic
1e-6 entries above are time to stagnation, **not converged solves**.

Final default comparison:

| Backend | Case | Solver | Setup ms | Refresh ms | Iterations | Solve ms | Relative residual | Converged |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| CPU | 64^3 | PcgJacobi | 1.86 | 1.49 | 278 | 408.53 | 6.83e-7 | yes |
| CPU | 64^3 | AmgPcg | 50.65 | 3.38 | 15 | 691.74 | 7.87e-7 | yes |
| CPU | 128^3 | PcgJacobi | 16.02 | 12.68 | 640 | 7360.32 | 7.17e-7 | yes |
| CPU | 128^3 | AmgPcg | 411.28 | 25.10 | 21 | 7606.88 | 8.61e-7 | yes |
| CPU | 64^3, 1:1000 | PcgJacobi | 1.58 | 1.40 | 354 | 560.36 | 1.55e-6 | no |
| CPU | 64^3, 1:1000 | AmgPcg | 51.09 | 2.86 | 38 | 1737.31 | 1.51e-6 | no |
| Metal | 64^3 | PcgJacobi | 1.65 | 0.47 | 278 | 278.77 | 6.83e-7 | yes |
| Metal | 64^3 | AmgPcg | 50.60 | 2.14 | 15 | 38.22 | 7.87e-7 | yes |
| Metal | 128^3 | PcgJacobi | 12.55 | 1.53 | 640 | 1310.16 | 7.17e-7 | yes |
| Metal | 128^3 | AmgPcg | 395.10 | 5.78 | 21 | 232.09 | 8.61e-7 | yes |
| Metal | 64^3, 1:1000 | PcgJacobi | 1.68 | 0.32 | 354 | 210.49 | 1.55e-6 | no |
| Metal | 64^3, 1:1000 | AmgPcg | 53.95 | 1.58 | 38 | 70.42 | 1.51e-6 | no |

Isotropic Metal solve speedups over PcgJacobi are 7.3x/5.6x. The serial CPU
backend favors explicit K2 for pressure: 10 iterations / 250.99 ms at 64^3 and
10 / 2074.14 ms at 128^3. The default follows measured GPU costs rather than
claiming to be optimal on CPU.

All 64^3 momentum cases below converge to 1e-6. Times are Metal solve medians;
J = BicgstabJacobi, AMG = default AmgBicgstab. Increasing/reversed order both
have six levels and complexities 1.137/1.136/1.134 for wind 5/20/100. Random
permutation yields 4/5/6 levels and complexities 1.373/1.340/1.324.

| Wind | Numbering | J iterations | J ms | AMG iterations | AMG ms |
| ---: | --- | ---: | ---: | ---: | ---: |
| 5 | increasing | 81 | 143.00 | 9 | 43.44 |
| 5 | reversed | 79 | 135.10 | 9 | 37.03 |
| 5 | random | 77 | 205.93 | 11 | 67.51 |
| 20 | increasing | 69 | 116.45 | 7 | 29.96 |
| 20 | reversed | 70 | 73.02 | 7 | 28.19 |
| 20 | random | 70 | 105.22 | 7 | 45.40 |
| 100 | increasing | 72 | 111.16 | 9 | 37.60 |
| 100 | reversed | 73 | 77.33 | 9 | 37.27 |
| 100 | random | 72 | 116.76 | 7 | 50.34 |

For increasing numbering, CPU J/AMG medians are 582.19/802.26 ms (wind 5),
494.05/620.23 ms (20), and 542.00/827.58 ms (100). Every Metal momentum case
favors AMG by 2.1x to 3.9x; CPU favors Jacobi-BiCGStab for these cases.

Recommended starting solvers for steady RANS:

| Equation/backend | Starting solver | Nullspace choice |
| --- | --- | --- |
| SPD pressure correction, Metal | AmgPcg, default V5 | Constant for connected all-velocity/wall Neumann; None when pinned/outlet |
| SPD pressure, serial CPU | AmgPcg with explicit K2 | Same equation-based choice |
| Upwind momentum/turbulence, Metal | AmgBicgstab V5 for the measured wind 5-100 family | None |
| Upwind momentum/turbulence, serial CPU | BicgstabJacobi for the measured family | None |
| Negative-diagonal nonsymmetric system | BicgstabJacobi; AMG rejects it | Problem-defined; no automatic projection |

Unstructured RANS coupling, different reaction terms and refreshed hierarchy
quality still need representative application measurements. BiCGStab may break
down on other nonsymmetric systems; non-converged stats must be handled.

The 1:1000 random-forcing case still cannot meet 1e-6: default AMG returns
1.51e-6, truthfully non-converged. At 1e-5 it converges in 19 iterations:
CPU 914.45 ms, Metal 30.78 ms, residual 8.89e-6 (single sample).
Previously, independent host-double CG reached 1.23e-13 on this forcing,
while rounding its solution to float gave 1.68e-6. That historical reference
supports a storage-rounding explanation, not a proof of the minimum residual
over all float vectors. The requested 1e-6 remains unmet on that forcing;
higher-precision solution storage would require a deliberate contract change.

Historical verification
-----------------------

All **34 tests / 1,254,455 assertions pass on CPU and Metal**. A separate
CPU Debug AddressSanitizer + UndefinedBehaviorSanitizer build passes all
34 tests / 627,205 assertions, including the debug nullspace consistency checks.
Metal compilation/runtime is verified here. CUDA is not available on this
machine and remains unverified.

New regressions failed before the fixes: directional strength/isolated-row
coarsening, reduction/allocation budgets, rounded Neumann projection, duplicate
projection/final-direction work, negative
AMG diagonals and divergence recovery. Refresh now changes coefficients by
non-power-of-two factors plus a cell-dependent reaction, comparing the reused
hierarchy with a fresh setup. Actual coarse matrices and actual P maps are
observed on the real solve path and checked against an independent double
host P^TAP, both before and after nonuniform refresh. Deliberately disabling
coarse refresh makes both strengthened tests fail (2,154 assertions); restoring
it passes the full suite. Test instrumentation is separate from solver code.

The representable manufactured 32^3/64^3 Poisson test at 1e-8 now has an
absolute <=15-iteration ceiling: default V5 takes 10/4 iterations, with
independent host-double residual checks. Iteration counts can be nonmonotonic
near a float rounding threshold. Random forcing is tested separately at 1e-6
and takes 11/15 iterations on those grids. Tests explicitly assert identical
CPU/Metal iteration counts, residuals and solutions for selected 2D/3D integer
stencils and upwind cases; generic random-coefficient reductions need not be
bitwise identical between double CPU accumulation and compensated Metal pairs.

Coverage also includes wind 5/20/100 in increasing, reversed and random order,
isolated cells, random-coefficient Neumann matrices with both diagonal assembly
orders, incompatible forcing, negative-row Jacobi solves, refreshed negative
AMG diagonals, controls/breakdown, empty matrices, preserved small reference
terms, anisotropy, dense double nonsymmetric references and double CG references.
