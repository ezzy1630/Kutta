CSR operations, Krylov iteration, smoothing, transfer and coarse-operator refresh
use Device buffers and generated kernels. The host builds aggregation structure
once; solve/refresh never download or traverse device arrays. Kernels gather,
write their own row/entry only, and use no atomics or double arithmetic.

A normal caller uploads one pattern and reuses the solver:

```cpp
auto a = kutta::upload(dev, host_matrix);
auto solver = kutta::make_solver(kutta::SolverKind::AmgPcg, dev, host_matrix);
// x and b are Device float buffers; a.val initially matches host_matrix.val.
auto stats = solver->solve(dev, a, x, b, {.rel_tol = 1e-6});
// Assemble or upload new values into a.val, keeping the pattern unchanged.
solver->refresh(dev, a);
stats = solver->solve(dev, a, x, b, {.rel_tol = 1e-6});
```

The solver belongs to its original Device, which must outlive it. All subsequent
matrices must have the exact original pattern and nullspace. Dimensions are
checked at the API boundary; the caller owns the immutable-pattern invariant.
Factories initialize value-dependent data from the host matrix. Call refresh
before solving after any value change. Pattern-only upload initializes values
to zero, ready for assembly. SpMV/residual outputs must not alias inputs.

`amg.hpp` adds optional configuration without changing csr.hpp or solver.hpp:

```cpp
kutta::AmgOptions options;
options.cycle = kutta::AmgCycle::V; // K is the default
options.aggregation_passes = 2;    // default 3; each pass matches pairs
auto solver = kutta::make_amg_solver(kutta::SolverKind::AmgPcg, dev, host_matrix, options);
```

Aggregation greedily matches the strongest unassigned connection, with index
order breaking ties and a default row-relative strength threshold of 0.25.
Repeated matching follows strong couplings first on stretched cells. Three
passes give aggregates of up to eight cells. Interpolation is piecewise
constant, restriction is its transpose. Each coarse entry owns a sorted CSR
list of contributing fine entries, implementing `A_c = P^T A P` by deterministic
sums. Refresh rebuilds only values and smoother diagonals on the device; it
neither rebuilds aggregates/maps nor reads values back to the host. Large
coefficient changes can reduce hierarchy quality; rebuild explicitly if needed.

The smoother is damped L1-Jacobi with
`D_l1 = max(|a_ii|, (|a_ii| + sum_{j != i}|a_ij|)/2)` and damping 2/3. This
bounds the symmetrically scaled SPD spectrum by two, even for non-diagonally
dominant SPD inputs. Equal pre/post sweeps preserve the SPD V-cycle property
and are required for AMG-PCG. Defaults are two sweeps on each side and 32
sweeps on the coarsest level (normally at most 24 rows). The coarsest solve
stays on the device; there is no dense host solve or transfer cost. A graph
that cannot coarsen falls back to smoothing at its current size.

A V-cycle is a fixed preconditioner. A K-cycle uses up to two recursive inner
Krylov corrections per coarse solve, with an adaptive 0.25 residual-reduction
exit. Inner SPD corrections use A-conjugate directions; nonsymmetric corrections
use GCR. Outer PCG uses the one-direction flexible-CG recurrence for K-cycles.
AmgBicgstab supports either cycle; BiCGStab can break down on difficult
nonsymmetric systems, which returns non-converged stats rather than success.
The nonsymmetric test exercises a three-level K-cycle and a V-cycle against an
independent double dense solve.

Pressure uses constant-nullspace projection, rather than pinning a row. Setup
detects symmetric operators whose stored-float row sums are exactly zero.
This supports connected Neumann Laplacians and leaves small nonzero reference
terms intact. Initial guesses, residuals and preconditioned corrections are
projected; coarse RHS and corrections are projected at every level. The result
has zero mean. A RHS with a nonzero mean beyond float rounding tolerance is
rejected. A general nonsymmetric nullspace or multiple disconnected-component
nullspaces requires an explicit problem-side constraint; automatic constant
projection does not solve those cases. Refresh must retain the nullspace.

SolveControl uses `max(abs_tol, rel_tol * initial_residual)`. Norms are absolute
L2 norms of `b - A*x`, reduced in double by Device, and include nonzero initial
guesses. Iterations count attempted Krylov updates or Jacobi sweeps; an already
satisfied tolerance performs zero updates. Negative curvature, zero/nonfinite
BiCGStab denominators, unrepresentable scalar coefficients, or 12 stagnant
reliable-residual replacements return `converged=false`. Invalid CSR/control,
nonfinite input vectors, incompatible Neumann forcing and invalid diagonals
throw Error. All numerical exits recompute the final residual, including
max_iter=0. Reliable replacement occurs at a convergence candidate or every
50 iterations, restarting the Krylov recurrence after replacement.

Residuals use Neumaier summation and Dekker product correction in float to
resolve stencil cancellation. Product correction falls back to a rounded
product for operands above 1e30 to avoid splitting overflow. Backends must
preserve floating-point expression order (disable fast-math reassociation),
as also required by compensated reductions. No Device API change is needed.
A float solution cannot satisfy arbitrarily small requested residuals: the
one-row `7*x=1` regression verifies that this never becomes false convergence.
Metal/CUDA compilation, backend comparisons and performance remain unverified
until those backends are integrated; the tests enumerate every available
backend automatically.

Measured on Apple M4 Pro, 24 GiB, macOS 26.6.2 (25G83), AppleClang
21.0.0.21000101, Release, serial CPU backend. Source commit:
`4f4619a25eb48bb2f4d46cfbb2bb228d32533227`. Commands:

```sh
cmake -S . -B build -G Ninja && cmake --build build && ./build/kutta-tests
./build/bench_linalg
```

All 15 tests / 100201 assertions pass. An AddressSanitizer + UndefinedBehaviorSanitizer
Debug build also passes the full suite. Coverage includes generated kernels,
2D/3D Poisson, 1:1000 and 1:10000 anisotropy, random-coefficient SPD diffusion,
upwind convection-diffusion against a partial-pivoted double dense solve,
Neumann projection against double CG, controls/breakdown, empty matrices,
small reference terms, and identical refreshed/fresh convergence after scaling.
No runtime, generator, existing contract or build-system files were modified.

At rel_tol=1e-8, the representable manufactured 32^3/64^3 Poisson solutions
take 11/11 AMG iterations, with independent host-double residual checks.
Operator complexity is 1.138/1.140 (5/6 levels). For deterministic random
forcing retaining low modes at rel_tol=1e-6, the same grids take 9/10 iterations.
These are distinct precision and scalability tests; the manufactured result
does not establish arbitrary-RHS convergence to 1e-8 in float.

On 64^3 isotropic random forcing, comparison runs of 1/2/3 pairing passes took
8/9/10 iterations, with operator complexities 1.989/1.329/1.140 and solve times
0.459/0.284/0.233 seconds. Three passes were selected for lower measured CPU
solve time and complexity. A two-pass V-cycle took 18 iterations / 0.415 s,
versus the two-pass K-cycle's 9 / 0.284 s. These single-sample comparisons are
CPU tuning evidence, not GPU performance claims.

Final benchmark uses identical deterministic random forcing for both solvers,
rel_tol=1e-6, max_iter=2000 and default three-pass K-cycles. Setup includes
host aggregation, device allocation/upload and initial refresh; the fine
matrix upload is outside the timer. Refresh times exclude value upload and
include synchronization; solve times also synchronize. These are single-run
wall times. Full machine-readable records are in benchmarks-cpu.jsonl.

| Case | Solver | Setup ms | Refresh ms | Iterations | Solve s | Relative residual | Converged |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 64^3 | PcgJacobi | 7.44 | 1.08 | 278 | 0.368 | 6.83e-7 | yes |
| 64^3 | AmgPcg | 45.36 | 2.47 | 10 | 0.233 | 4.31e-7 | yes |
| 128^3 | PcgJacobi | 72.92 | 9.35 | 640 | 6.858 | 7.17e-7 | yes |
| 128^3 | AmgPcg | 389.96 | 20.61 | 10 | 1.915 | 5.04e-7 | yes |
| 64^3, 1:1000 | PcgJacobi | 7.82 | 1.13 | 354 | 0.514 | 1.55e-6 | no |
| 64^3, 1:1000 | AmgPcg | 46.63 | 2.34 | 38 | 0.930 | 1.44e-6 | no |

AMG has 6/7 levels and operator complexities 1.140/1.142 at 64^3/128^3;
the 1:1000 case has 6 levels / 1.134. Isotropic solve speedup is 1.58x/3.58x.
The anisotropic random-forcing benchmark stagnates above its requested
tolerance for both solvers and is not an accepted 1e-6 convergence result.

As a separate precision check on the same anisotropic matrix and forcing,
independent host-double CG reached relative residual 1.23023e-13. Rounding
that solution to float raised the host-double residual to 1.67561e-6, before
any device arithmetic. The solver's 1.44e-6 stagnation is therefore consistent
with storage rounding, although the rounded reference is not a proof of the
minimum residual attainable over all float vectors. Achieving the requested
1e-6 on this forcing remains unmet; higher precision solution storage would
require a deliberate contract change, not a false convergence flag.
