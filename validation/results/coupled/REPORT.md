# Native pressure-based Coupled solver

Fluent's Coupled scheme ran in Kutta as SIMPLEC, with a convergence note. It now
runs natively. Fluent Coupled cases select it with no note, and SIMPLEC remains
selectable (`methods/pressure-velocity-coupling`, `--scheme` in bench).

- Branch `t3/native-coupled-solver`; implementation commit `9d0f838`, base `63287b9`.
- Machine: Mac M4 Pro (CPU and Metal). CUDA cannot run here (see [CUDA](#cuda)).
- Other agents' slow tiers and benches shared the CPU and GPU throughout. All
  timings are contended. Scheme order alternates per case
  (`run_comparison.sh`) so the contention does not favour one scheme.

## Method

Fluent Theory Guide 23.4.3.2 ([pressure-velocity coupling][pv], eqs. 23-58 to
23-62): momentum and continuity form one block system per outer iteration, with
an implicit pressure gradient (23-58) and implicit Rhie-Chow pressure
dissipation in the continuity rows (23-60). The block form (23-61/62) is a 4×4
block per cell in (u, v, w, p). Pseudo-time follows 23-101 ([pseudo-transient][pt]),
with Kutta's existing global time step (23.6.2) on the momentum diagonal.
Relaxation is Fluent's explicit form (the elbow saves
`pressure-coupled/{momentum,pressure} pseudo-explicit-relax 0.5`). Turbulence
stays segregated.

Kutta's choices:

- **Delta form.** The block RHS is the steady SIMPLEC residual itself: momentum
  `b − A u` with the reconstructed face pressure, and the Rhie-Chow mass
  imbalance with the physical (Majumdar) diagonal. The fixed point is therefore
  SIMPLEC's discretization exactly, whatever the Jacobian approximates. The
  Jacobian, pseudo-time and relaxation change only the convergence path.
- **Jacobian.** The momentum blocks are the assembled transport rows (first-order
  upwind part plus diffusion, pseudo-time on the diagonal, no implicit
  relaxation). The momentum-pressure block is the face-area gradient. It also
  carries the own-cell part of the least-squares face extrapolation,
  `Q_i = T_i M_i⁻¹`, which the residual's boundary and skew face pressures use.
  The continuity-velocity block is the face-area divergence, and the
  continuity-pressure block is the compact Rhie-Chow Laplacian
  `ρ V/a_P · |S|²/(S·d)`. The kernel header (`src/kutta/kernels/fv_coupled.kern.h`)
  gives every term.
- **Linear solver.** Right-preconditioned BiCGStab on the 4n block vector. It
  reuses the device scalar kernels of the existing Krylov solvers. The
  preconditioner is one aggregation AMG V-cycle: pairwise aggregation on the
  pressure-Laplacian strength graph, three passes per level, with Galerkin block
  coarse operators. The smoother is damped block Jacobi with an arrow-block
  Schur solve. Fluent uses coupled AMG with ILU smoothing ([AMG][amg]); block
  Jacobi is the gather-only, atomic-free GPU equivalent. Termination is Fluent's
  default AMG criterion, 0.1 relative reduction, with at most 20 iterations. A
  breakdown that grows the residual falls back to a single V-cycle, so a bad
  Krylov step never reaches the flow field.
- **Kernels.** Every kernel is a `kernelgen` kernel (CPU, Metal, CUDA). It is
  gather-only with no atomics, and the coarse operators are built by per-row
  gathers over aggregate members. CPU and Metal take the same iteration
  counts in `tests/fv/test_coupled.cpp` and agree within its 1e-5 field check.

[pv]: https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/flu_th_sec_uns_solve_pv.html
[pt]: https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/flu_th_sec_pseudo_auto_local.html
[amg]: https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/flu_th_sec_amg_multigrid.html

## Outcome against the acceptance criteria

| criterion | result |
|---|---|
| Same fixed point | **Yes where both schemes converge.** Hex channel and cavity agree to float floors in `tests/fv/test_coupled.cpp`. Elbow final drag agrees to 6e-6–8e-6 relative. Level-2 TMR plates agree to 5e-6 (SST) and 3e-5–8e-5 (SA) in drag and cf. On the bumps SIMPLEC is still unconverged after 20 000 iterations and is drifting toward Coupled's loads (0.4–0.6 % gap left); this is consistent with a shared fixed point but does not prove one. Delta form makes the shared fixed point structural (the residual is SIMPLEC's). |
| Metal wall-time speedup, elbow plus at least one TMR case | **Yes.** Elbow 1.59× (Metal) and 2.20× (CPU). SST plate 5.35× / 3.28×, SA plate 6.55× / 6.26× (Metal / CPU). Bumps: Coupled loads settle by 1300 iterations (12–27 s on Metal); SIMPLEC's do not settle in 3000 or in 20 000. |
| V6 velocity L2 ≤ 3.60 %, pressure ≤ 2.96 % | **Held.** The slow-tier V6 run with Coupled gives velocity 3.6010 %, pressure 2.9622 % (baseline SIMPLEC 3.6016 % / 2.9625 %), CPU and Metal. |
| Continuity gates retried | **Still failing, for both schemes.** See [Continuity gates](#continuity-gates). Coupled reaches lower raw continuity on every case but cannot pass a scaled 1e-5 gate whose target lies below the float32 floor. |
| No new failed slow assertion | **One changed outcome, not caused by this change.** The V5 bump SST CPU level-2 cf check moved 0.0992 → 0.1068 against its 0.1 gate. That test runs SIMPLEC under a 600 s wall budget: baseline stopped at 2300 iterations, candidate at 2200 under heavier GPU/CPU sharing. SIMPLEC is bitwise identical between the binaries on that case. 20 V6 per-iteration REQUIREs are "missing" because Coupled leaves the 50-iteration loop sooner. See [Slow tier](#slow-tier). |
| Quick tier | **Green.** `ctest -LE slow`: 2/2. `kutta-tests --test-suite-exclude=slow` with fixtures: 228 cases, 7 414 385 assertions. |
| CPU == CUDA parity | **Not executed (no CUDA here).** The kernels are kernelgen dialect only (no atomics, gather-only). A clang CUDA device-mode parse of `fv_coupled.kern.h` passes; it is not an nvcc build and `--fmad=false` was not exercised. CPU and Metal agree to within the test's 1e-5. |

## SIMPLEC vs Coupled

Commit `9d0f838`, `validation/results/coupled/run_comparison.sh build/bench_convergence <out>`,
summarized by `summarize.py` (criteria in its docstring; logs in `bench/`). Coupled
uses Fluent's pseudo-transient explicit relaxation 0.5/0.5 everywhere. The
elbow's case saves exactly that; the TMR and channel cases have no Fluent
case, and SIMPLEC uses each case's bench settings. Both schemes use the same
global pseudo-time step. "Its / s" is the first iteration meeting the case's criterion: V6's
criterion for the elbow (1e-6 transport, raw continuity below Fluent's saved
2.476e-6); transport defect ≤ 1e-5 plus loads changing ≤ 0.1 % per 100
iterations for TMR. Times are contended wall seconds including monitor
overhead.

| case | backend | cells | SIMPLEC its / s | Coupled its / s | speedup | loads settled S → C | Coupled vs SIMPLEC (final) | Coupled B/cell | peak RSS MiB S / C |
|---|---|---|---|---|---|---|---|---|---|
| bump-sa | cpu | 56320 | not reached in 3000 (723.1 s) | not reached in 3000 (241.8 s) |  | no → 1300 / 115.6 s | n/a: not both settled | 340 | 266 / 297 |
| bump-sa | metal | 56320 | not reached in 3000 (90.1 s) | not reached in 3000 (24.2 s) |  | no → 1300 / 12.2 s | n/a: not both settled | 340 | 249 / 257 |
| bump-sst | cpu | 56320 | not reached in 3000 (810.2 s) | not reached in 3000 (417.9 s) |  | no → 1300 / 194.6 s | n/a: not both settled | 340 | 296 / 299 |
| bump-sst | metal | 56320 | not reached in 3000 (137.5 s) | not reached in 3000 (55.8 s) |  | no → 1300 / 26.6 s | n/a: not both settled | 340 | 249 / 266 |
| channel-983k | metal | 983040 | not reached in 200 (27.7 s) | not reached in 200 (17.9 s) |  | no → no | n/a: not both settled | 343 | 2224 / 2448 |
| channel | cpu | 4608 | diverging (raw 4.2e-01 at 1000) | diverging (raw 9.2e+23 at 1000) |  | - → - | n/a: not both settled | 344 | 22 / 26 |
| channel | metal | 4608 | diverging (raw 3.9e-01 at 1000) | diverged at 73 (0.6 s) |  | - → - | n/a: not both settled | 344 | 33 / 33 |
| elbow | cpu | 17822 | 550 / 55.36 s (raw 1.45e-06) | 350 / 25.21 s (raw 5.73e-08) | 2.20x |  | drag 6.2e-06 | 389 | 154 / 167 |
| elbow | metal | 17822 | 550 / 13.34 s (raw 3.26e-07) | 350 / 8.40 s (raw 5.73e-08) | 1.59x |  | drag 7.9e-06 | 389 | 156 / 158 |
| plate-sa | cpu | 13056 | 1700 / 95.42 s (raw 1.40e-07) | 300 / 15.25 s (raw 6.84e-08) | 6.26x | 1700 / 95.4 s → 300 / 15.2 s | drag 3.7e-05, cf0 8.4e-05 | 340 | 72 / 79 |
| plate-sa | metal | 13056 | 1700 / 39.90 s (raw 1.08e-07) | 300 / 6.09 s (raw 7.18e-08) | 6.55x | 1700 / 39.9 s → 300 / 6.1 s | drag 3.4e-05, cf0 7.7e-05 | 340 | 73 / 77 |
| plate-sst | cpu | 13056 | 1500 / 83.26 s (raw 1.02e-07) | 400 / 25.37 s (raw 6.98e-08) | 3.28x | 1500 / 83.3 s → 400 / 25.4 s | drag 5.4e-06, cf0 4.9e-06 | 340 | 77 / 86 |
| plate-sst | metal | 13056 | 1500 / 54.99 s (raw 1.46e-07) | 400 / 10.27 s (raw 7.09e-08) | 5.35x | 1500 / 55.0 s → 400 / 10.3 s | drag 6.1e-06, cf0 5.2e-06 | 340 | 75 / 82 |

Per iteration Coupled costs about the same as SIMPLEC: 26 vs 37 ms on Metal,
63 vs 55 ms on CPU (level-2 SST plate), and 90 vs 138 ms on Metal for the
983 040-cell tet channel. The speedup comes from iteration count. The BiCGStab
solve stops at 0.1 relative reduction in 1–3 iterations.

Bumps after 20 000 SIMPLEC iterations (`bench/long_bump-*`), Coupled at 3000:

| | SIMPLEC drag | Coupled drag | SIMPLEC cf0 | Coupled cf0 | SIMPLEC raw continuity | Coupled |
|---|---|---|---|---|---|---|
| SST | 0.0036383 (still falling) | 0.0036229 | 0.0051743 | 0.0051705 | 7.0e-5 | 4.4e-6 |
| SA | 0.0036171 (still falling) | 0.0035962 | 0.0052720 | 0.0052678 | 4.2e-5 | 4.4e-6 |

V6 as run by the slow tier (contended; baseline and candidate were separate runs):

| mode / backend | SIMPLEC (baseline `63287b9`) | Coupled (`9d0f838`) |
|---|---|---|
| seeded CPU | 350 its, 38.6 s | 300 its, 22.1 s |
| cold hybrid CPU | 550 its, 55.5 s | 350 its, 26.0 s |
| seeded Metal | 350 its, 8.44 s | 300 its, 7.50 s |
| cold hybrid Metal | 550 its, 12.60 s | 350 its, 8.71 s |

### Tetrahedral channel (bench `--case channel`)

The bench SST channel is Freudenthal tetrahedra with stretched, skewed wall
layers at ν = 1e-4. **Both schemes diverge on it:** SIMPLEC slowly (raw
continuity 0.004 → 0.4 in 1000 iterations, CPU and Metal), Coupled by
iteration 73–80. At 983 040 cells, Coupled's raw continuity turns upward
after iteration 120 and SIMPLEC's is still falling at 200.

The cause is the discretization's Rhie-Chow jump, not Coupled. The jump is
the difference of the two reconstructed face pressures,
`p_P + ∇p_P·r_P − p_N − ∇p_N·r_N` (`fv_pressure_faces`). Diagnostic runs on a
laminar 24×8×4 stretched tet channel (ν = 1e-4) used a temporary,
uncommitted switch that drops the gradient terms from that jump only:

- with the switch, Coupled converges at Δτ = 0.01 and 0.001, and SIMPLEC stays bounded;
- without it, SIMPLEC also diverges (raw 0.008 → 9.2 by iteration 2000) and Coupled diverges sooner.

On tetrahedra the reconstruction gradients flip the sign of some cells'
Rhie-Chow self-coupling: the measured own-cell sensitivity ranges from −0.89
to 1.19 where the compact operator assumes 1. When the pseudo-time term
dominates the physical diagonal, pressure is controlled by that term alone.
Approaches tried in the Jacobian that did not fix it:

- the exact pressure linearization as the Krylov operator diverged faster; the system is near-singular in pressure;
- capping the momentum diagonal at 2–20 a_P only delays the blow-up;
- scaling the Rhie-Chow block by 0.1–4 does not help.

A fix belongs with the face-reconstruction work, not in the Coupled
Jacobian.

The own-cell pressure term Q was needed for tetrahedra regardless. Without
it a ν = 0.01 24×8×4 tet channel diverges at α = 0.5. With it, the channel
converges (u = 1e-7 by 2000 iterations), and SIMPLEC's trajectories are
unchanged. `FV Coupled converges on a tetrahedral channel` guards this case.

## Continuity gates

The TMR gates scale continuity by the largest raw value of the first five
iterations and require ≤ 1e-5. That scale depends on the start-up transient:

| level-2 run (Metal, 3000 its) | raw | first-five scale | scaled |
|---|---|---|---|
| SST plate, SIMPLEC | 1.6e-7 | 4.4e-6 | 3.6e-2 |
| SST plate, Coupled | 6.8e-8 | 5.3e-4 | 1.3e-4 |
| SA plate, SIMPLEC | 1.6e-7 | 5.3e-6 | 3.0e-2 |
| SA plate, Coupled | 6.9e-8 | 2.1e-4 | 3.3e-4 |
| SST bump, SIMPLEC | 2.7e-5 | 0.33 | 8.3e-5 |
| SST bump, Coupled | 4.4e-6 | 1.7e-3 | 2.6e-3 |

Coupled reaches lower raw continuity on every case, and its plate raw values
are flat from about iteration 300. They sit on the float32 floor of the
face-flux sums. Passing the gate on the plates needs raw ≈ 2e-9–5e-9, below
that floor. The gate is unchanged, so both schemes still fail it. V6's
`r0.continuity.raw<8e-4` concerns the initial state and fails identically in
baseline and candidate.

## Memory

`CoupledSolver::bytes()`, the device memory the Coupled solver owns beyond SIMPLEC's state:

| mesh | B/cell |
|---|---|
| hex / TMR quads | 340 |
| tet channel | 343–344 |
| elbow (polyhedral Fluent mesh) | 389 |

Composition per fine cell:

- packed diagonals, G, H and the 3×3 Q tensor: 76 B;
- ten 4-vector Krylov buffers: 160 B;
- two smoother buffers: 32 B;
- one Rhie-Chow coefficient per face: 8–16 B;
- aggregation levels and transfers: the remaining 50–60 B.

Peak process RSS (`/usr/bin/time -l`) for the 983 040-cell tet channel on
Metal: SIMPLEC 2224 MiB, Coupled 2448 MiB (+10 %, +239 B/cell). Scaled
linearly, that whole-process footprint (host mesh and setup included) is
about 2.4 KB/cell for SIMPLEC and 2.6 KB/cell for Coupled:

- 24 GB Mac: about 9 M tet cells either way;
- 30 M cells: about 75 GB.

The 10–30 M-cell FSAE target is therefore limited by Kutta's existing
footprint, not by Coupled. Coupled adds about 10 GB at 30 M cells.
Polyhedral cells have more faces per cell and will cost more than these tet
figures. Not measured: CUDA device memory on the 64 GB-RAM RTX PC.

Further reduction available: BiCGStab needs seven of the ten 4-vectors if `b`
aliases `r` and `z`/`zs` share storage with `p`/`s` (about −48 B/cell).

## Slow tier

Baseline `63287b9` (`slow/baseline.*`, one run, 80 min). Candidate `9d0f838`
in two chunks merged with `merge_doctest_xml.py` (`-tc="V5*,V7 SA*"` and
`-tce="V5*,V7 SA*"`, 55 + 25 min). Each chunk stays under the background-job
time limit; the merged report has the same 24 cases. Comparison:
`tools/forensics/compare_slow_assertions.py slow/baseline.xml.gz slow/candidate.xml.gz`
→ `slow/comparison.json`:

| | count | explanation |
|---|---|---|
| regressions | 1 | V5 bump SST, CPU, level 2, station 0.632: `cf_error<=.1` 0.0992 → 0.1068. SIMPLEC run stopped by its 600 s budget at 2300 vs 2200 iterations, unconverged; bench SIMPLEC on that grid is bitwise identical between the two binaries (STEP digests equal). |
| missing assertions | 20 | V6 `REQUIRE(std::isfinite(maximum))` and `REQUIRE(r.maximum_velocity<100)` inside the 50-iteration loop. Coupled converges in 300/350 iterations (seeded/cold) instead of 350/550, i.e. five fewer loop passes per backend. |
| changed values | 59 | V6 (Coupled's path) and the budget-bound V5 bump. |
| failed in both | V6 29, V12 fan import 1, V2 tet MMS 6, V5 plate 6, V5 bump 9 (10 in the candidate: the regression above), V7 SA 17 | Pre-existing. V6's are the strict-setup rejection without the WFGC opt-in, Fluent-accuracy thresholds, and the initial continuity check. |

SIMPLEC is bitwise unchanged by this branch. Bench STEP lines, wall time
excluded, are identical between the baseline and candidate binaries for:

- the level-0 plate (300 iterations);
- the level-2 bump (200 iterations);
- the elbow (100 iterations, CPU).

## CUDA

No CUDA toolchain on this machine; the CUDA backend was not built or run.
`clang++ -x cuda --cuda-device-only -nocudainc` with attribute stubs parses
`prelude.h`, `fv_flow.kern.h` and `fv_coupled.kern.h` without diagnostics.
That checks the dialect only, not nvcc or `--fmad=false` parity.

## Open limits

- Tetrahedral, low-viscosity channels diverge with either scheme (above);
  Coupled diverges sooner.
- The bump fixed-point agreement is approach, not proof: SIMPLEC does not
  converge within 20 000 iterations.
- Coupled's TMR turbulence residuals (SST k ≈ 5e-5 on the bump) stay above
  1e-5 while loads are settled. Turbulence is segregated and implicitly
  relaxed, as in Fluent.
- Scaled continuity gates remain unreachable for both schemes (float32 floor).
- Linear solver: block-Jacobi smoothing instead of Fluent's ILU, and a
  V-cycle instead of the F-cycle default. A diverged Krylov solve falls back
  to one V-cycle.
- The backflow total-pressure adjustment and FSAE interior pressure jumps
  are not in the Jacobian; they act through the residual only.
- `app/fixtures/cases/mixing_elbow/progress-*.jsonl` still contain the old
  "SIMPLEC replaces Coupled" warning; those fixtures belong to the app lane
  and were not regenerated.

## Commands

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build
export KUTTA_FLUENT_EXAMPLES=<fixtures>   # validation/fluent_examples
ctest --test-dir build -LE slow
build/kutta-tests --test-suite-exclude=slow
build/kutta-tests -tc="FV Coupled*"
validation/results/coupled/run_comparison.sh build/bench_convergence <out>
python3 validation/results/coupled/summarize.py <out>
build/bench_convergence --case bump --model sst --level 2 --backend metal --scheme simplec --iterations 20000 --max-seconds 3000
build/kutta-tests --test-suite=slow -tc="V5*,V7 SA*" --success --reporters=xml --out=part1.xml
build/kutta-tests --test-suite=slow -tce="V5*,V7 SA*" --success --reporters=xml --out=part2.xml
python3 validation/results/coupled/merge_doctest_xml.py candidate.xml.gz part1.xml part2.xml
python3 tools/forensics/compare_slow_assertions.py baseline.xml.gz candidate.xml.gz
```
