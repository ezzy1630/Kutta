# Float-pair state: continuity below the float32 floor

Branch based on main `1c4ff13`. Apple M4 Pro, CPU reference backend and Metal. CUDA was not run; it compiles in CI only. The Fluent corpus is read from `KUTTA_FLUENT_EXAMPLES` and is not committed. Logs are under [`runs/`](runs/).

## Result

SIMPLE and SIMPLEC now carry velocity, pressure and the face mass flux as float pairs hi + lo (about 48 bits) on every backend, Metal included. The discrete equations are unchanged. The [residual floor report](../residual_floor/REPORT.md) showed that the cylinder's 1e-5 continuity criterion could not be reached in float32; with the pair it is reached.

| Cylinder (Fluent 18 adjoint tutorial, `reverse-flow`), Metal | main (float32) | this branch |
| --- | ---: | ---: |
| scaled continuity at iteration 1000 | 4.6e-4 | 4.2e-4 |
| first iteration below the case's 1e-5 | never | 1700 |
| plateau, raw kg/s (iterations 2000–3000) | 1.0e-2 | 4.0e-5 – 5.0e-5 |
| plateau, scaled | 3.9e-4 | 1.5e-6 – 1.9e-6 |

- The plateau is about 200× below main's.
- The CPU and Metal histories are bitwise identical.
- Until iteration 1000 the history tracks main to two or three digits; the run is not slower to converge. Main's raw continuity at iterations 500, 1000 and 1500 is 0.0756, 0.0119 and 0.0101; this branch's is 0.0757, 0.0108 and 0.0030.

The pair is on by default for SIMPLE and SIMPLEC (`SolverSettings::paired_state`). It is off for Coupled and for cases with interior jumps (fans, porous jumps), where the state stays float32. Every kernel other than the few listed below, and every output, reads hi, which stays within one ulp of the velocity, pressure and flux scales of hi + lo.

Cost and checks, against main:

- **Memory:** +53 to +78 bytes per cell of device memory on SIMPLE and SIMPLEC (+2.5% to +5.1%). Coupled is unchanged to the byte.
- **Time per iteration:** +0.3% on Metal and +2% on CPU while the solution still moves. Once main reaches its float floor, the paired run is cheaper per iteration, because main's momentum solves then run 50 to 58 Krylov iterations against rounding noise.
- **Tests:** the quick tier passes (259 cases; main's 258 plus one new test). The slow tier has the same 65 failing assertions as main, with no regressions and no recoveries.
- **Converged Fluent references:** the cylinder and `elbow_param` gates pass on CPU and Metal, with field errors against Fluent's saved fields unchanged to three or four digits. The nozzle and V6 run Coupled and are bitwise unchanged.

Details are in "Tests and Fluent references" and "Cost" below.

## Which state sets the floor

A CPU build with every float replaced by double, which rounds chosen states back to float after every write, attributes the floor. Cylinder, raw continuity in kg/s; the criterion needs 2.6e-4 at Kutta's scale (see the residual floor report):

| Rounded to float after each write | Raw continuity at iteration 3000 |
| --- | ---: |
| nothing (all double) | 4.0e-8 |
| pressure correction p′ | 4.0e-8 |
| pressure | 2.4e-4 |
| velocity | 2.6e-3 |
| momentum matrix coefficients | 3.1e-3 (iteration 2000) |
| momentum right-hand side | 3.5e-4 (iteration 2000) |
| momentum matrix and right-hand side | 3.3e-3 |
| face velocities | 4.7e-3 |
| face mass flux | 5.8e-3 |
| all of the above | 8.2e-3 |

Logs: [`runs/fp64/`](runs/fp64/); the experiment is [`runs/fp64/experiment.patch`](runs/fp64/experiment.patch), not part of the branch.

Velocity, the mass flux, the face velocities and the momentum matrix each set a floor of their own. The pressure correction does not; the pressure alone almost does. Storing only velocity at higher precision cannot work, because the matrix and fluxes are rebuilt every iteration from float values. What does work is keeping each stored state as hi + lo and making every quantity derived from hi stop changing near convergence, so that lo converges against a fixed float operator.

## Design

Per cell the solver keeps lo words for velocity (3), pressure (1) and the pressure gradient of lo (3), plus two work arrays (imbalance and the momentum update). Per face it keeps the mass flux lo word and the previous mass flux; per boundary face, the pressure lo boundary value. Each SIMPLE iteration then does the following ([`fv_paired.kern.h`](../../../src/kutta/kernels/fv_paired.kern.h)):

1. **Momentum in delta form.** The residual r = b − A(hi + lo) + (a_relaxed − a_P)·lo − V∇p_lo is summed compensated (TwoSum and an fma product error). The linear solve finds δ, which is added to the pair.
2. **Rhie-Chow mass flux** F(hi) in float, plus its linearization in lo, which goes into the flux lo word.
3. **Continuity** summed compensated over hi and lo fluxes. The pressure-correction right-hand side is this pair imbalance.
4. **Flux and velocity correction** add into the pairs.

What each piece is for, measured on the cylinder:

- **Hysteresis.** hi moves only when lo exceeds a slack. Round-to-nearest flips hi whenever the true value crosses a rounding midpoint. Each flip re-rounds every float quantity built from hi and kicks the solution by about one ulp, which sustains a chaotic floor. With a plain renormalized pair (slack relative to each cell's own value) the floor was 7.6e-3 kg/s, barely below main.
- **Slack set by the scale hi feeds**, not hi's own value. A small velocity component or flux still enters sums at the scale of the largest. The scales are the maximum speed, max(ρU², max|p|) and ρU|S| per face, from the last monitor. Floor: 1.8e-3.
- **Flux hysteresis across recomputation.** F(hi) is recomputed every iteration through diagonal → interpolation weights → flux, and 10 to 25 faces cycled by one ulp on that loop. The flux kernel keeps the previous hi while the difference plus lo stays within the slack. Floor: 1.2e-3.
- **Pair-aware limiter.** With hi otherwise frozen, four cells near a limiter kink (x ≈ 16, y ≈ 3) still alternated with period 2. Running with hi frozen everywhere converged to 1.8e-8, and so did running without the limiter (1.7e-8), so the operators for lo were consistent and the limiter was the cause. Adding hysteresis to the limiter made it worse (5.0e-4 at 1 ulp of the face value, 7.4e-4 at 16). The Barth-Jespersen limiter now takes neighbour differences of the pair, (hi_j − hi_c) + (lo_j − lo_c), and is bitwise equal to the float limiter when lo is zero. Floor: 4e-5.

Tried and dropped: an exact linearization of the deferred correction in lo (a lo gradient, the skew term and the limiter's derivative) raised the floor to 1.9e-4 and barely changed the middle of the run.

### Slack

Slack trades how soon hi freezes against how much hi lags the true value while the solution still moves. Cylinder, Metal, raw continuity; slack in ulps of the velocity, pressure and flux scales ([`runs/slack/`](runs/slack/)). The sweep set the three slacks with a temporary environment variable that is not on the branch. The 1-ulp row matches the committed code bitwise: compare [`runs/floor/cylinder-metal-paired.log`](runs/floor/cylinder-metal-paired.log).

| Slack (velocity, pressure, flux) | Iteration 500 | 1000 | 1500 | First below 2.6e-4 | Iteration 2000 |
| --- | ---: | ---: | ---: | ---: | ---: |
| main, float32 | 0.0756 | 0.0119 | 0.0101 | never | 0.0106 |
| 16, 16, 16 | 0.128 | 0.0278 | 1.07e-4 | 1500 | 4.1e-5 |
| 16, 16, 4 | 0.0901 | 0.0213 | 1.8e-3 | 1700 | 4.6e-5 |
| 16, 16, 1 | 0.0861 | 0.0167 | 1.5e-3 | 1700 | 1.6e-5 |
| 4, 4, 4 | 0.0814 | 0.0157 | 1.6e-3 | 1700 | 2.7e-5 |
| 4, 4, 1 | 0.0759 | 0.0112 | 1.9e-3 | 1700 | 3.4e-5 |
| **1, 1, 1 (default)** | 0.0757 | 0.0108 | 3.0e-3 | 1700 | 5.0e-5 |

16 ulps reaches the criterion 200 iterations sooner but is 2.3× slower than main through the middle of the run, and leaves hi up to 16 ulps of scale from the pair. In the slow tier that took the 257² cavity from 900 to 3000 iterations before its steady test passed, and moved V4 cavity outputs by 2–4e-6, past that test's 2e-6 gate (the imported and generated grids froze hi differently). At 1 ulp both are as on main, and the float32 outputs that read hi are within one ulp of scale of the pair, as rounding would leave them.

## Initialization: why Fluent's iteration-1 continuity is 448 and Kutta's 13

Fluent's first continuity residual on the cylinder is 448.29 kg/s. Kutta's initialized state has a Rhie-Chow imbalance of 13.0 kg/s, and its first iteration reports 23.7 (its scale, the largest of the first five, is 26.0). Fluent scales by the largest of its first five too, so the scaled histories have denominators 17× apart. Findings ([`runs/init/`](runs/init/), CPU):

| | Rhie-Chow imbalance of the initial state, kg/s | Continuity at iteration 1, kg/s |
| --- | ---: | ---: |
| Kutta, least-squares ∇φ (current) | 13.0 | 23.7 |
| Kutta, Green-Gauss ∇φ (the case's gradient scheme) | 131 | 48.9 |
| Fluent | not saved | 448.29 |

- **Documented parts match.** Fluent's hybrid initialization (Fluent 2024 R2 Theory Guide, Hybrid Initialization, equations 23-149 to 23-158) solves a Laplace equation for a velocity potential φ with velocity ∇φ. Pressure comes from a Laplace solve only when at least one inlet and one outlet carry pressure information; otherwise it is the averaged boundary value. The cylinder has a velocity inlet, so its initial pressure is a constant, as in Kutta.
- **The gradient used for ∇φ is not documented.** Kutta uses a boundary least-squares gradient. Green-Gauss on this mesh gives a ten-times larger initial imbalance. Neither reproduces 448.
- **Fluent's iteration-1 residuals belong to the state after its first update, not to the initial state.** Fluent normalizes momentum by Σ|a_P u_P|. Kutta's value for the initial state is 2.66e7. After its first update it is 1.335e7, against Fluent's 1.339e7 at iteration 1. Fluent ran Coupled at Courant 200; Kutta runs this case under SIMPLEC (an existing named approximation), so its first update differs and its first continuity residual differs with it.

Nothing was changed. Kutta follows the documented behaviour; the remaining difference comes from an undocumented gradient choice and from Fluent's Coupled first step, and matching 448 would mean fitting a number, not reproducing a documented algorithm. The criterion is met at either denominator: at Fluent's scale of 448.29, this branch's plateau is 1e-7 scaled.

## Tests and Fluent references

Quick tier (`--test-suite-exclude=slow`, fixtures present): 259 of 259 test cases, 7,454,253 assertions, 0 failures. Main: 258 of 258, 7,453,864. The extra case is the new `FV SIMPLEC float-pair state falls below the float32 continuity floor` ([`tests/fv/test_steady.cpp`](../../../tests/fv/test_steady.cpp)). It runs a 24×8 channel for 200 iterations on both backends. Float32 stays above 1e-6 scaled continuity; the pair is below 1e-7 (1.95e-8 observed). The two agree in u to 1e-5 relative, and the paired u is bitwise equal on CPU and Metal.

Slow tier (`--test-suite=slow -r=xml --success`), compared with `tools/forensics/compare_slow_assertions.py`:

| | Assertions passed | Failed | Regressions | Recoveries | Missing |
| --- | ---: | ---: | ---: | ---: | ---: |
| main `1c4ff13` | 39,502 | 65 | | | |
| this branch | 39,502 | 65 | 0 | 0 | 0 |

The 65 are main's known failures: the V6 elbow's 25, V2 tet refinement 8, V5 bump 9, the backflow default 6, SA TMR 16 and the V12 fan 1. 916 assertions changed value but kept their outcome.

An earlier version with a 16-ulp slack failed 45 V4 cavity assertions; see "Slack" above. It also recovered 4 (3 backflow, 1 SA TMR), which the 1-ulp version does not. The gate was not changed.

Converged Fluent references, errors against Fluent's saved fields from the slow tier. Values are CPU; Metal agrees to the digits shown except `elbow_param`'s final residuals, which are 1.11e-5 → 6.8e-7 (continuity) and 4.4e-8 → 3.1e-9 (momentum) on Metal.

| Gate | main | this branch |
| --- | ---: | ---: |
| cylinder velocity L2 (< 9.7e-4) | 6.448e-4 | 6.448e-4 |
| cylinder pressure L2 (< 7.8e-3) | 5.175e-3 | 5.173e-3 |
| cylinder pressure drop (< 7.6e-3) | 5.086e-3 | 5.060e-3 |
| cylinder wall force (< 1.1e-3) | 7.051e-4 | 7.056e-4 |
| cylinder \|Cd/Cd_saved − 1\| (< 1.05e-3) | 6.944e-4 | 6.951e-4 |
| cylinder \|Cl − Cl_saved\|/Cd_saved (< 1.85e-4) | 1.226e-4 | 1.211e-4 |
| `elbow_param` velocity (< 0.016) | 0.010528 | 0.010528 |
| `elbow_param` pressure (< 0.087) | 0.05988 | 0.05988 |
| `elbow_param` pressure drop (< 0.017) | 0.01112 | 0.01114 |
| `elbow_param` final scaled continuity (< 1e-4) | 1.08e-5 | 9.2e-7 |
| `elbow_param` final momentum (< 1e-6) | 4.3e-8 | 3.4e-9 |

The nozzle and V6 run Coupled, so the pair is off. None of their assertions changed. A V6 Coupled `kutta solve` (`--approximate wfgc`, Metal, 300 iterations) writes byte-identical `result.vtu`, residuals, forces, fluxes and loads on both builds.

With convergence checking on, `elbow_param` stops at iteration 120 on both builds, scaled continuity 2.03e-4 (main) and 2.01e-4 (this branch). Its criteria are loose enough that the float floor never mattered there.

## Cost

Measured with `kutta solve` on an otherwise idle machine, main and this branch interleaved, two repeats each (the repeats agree within 0.5%). Device bytes are from `KUTTA_MEMORY_PHASES=1` after the iterations; RSS is from `/usr/bin/time -l`.

| Case, backend | Cells | Device MB, main → branch | Per cell | RSS MB, main → branch |
| --- | ---: | ---: | ---: | ---: |
| cylinder SIMPLEC, Metal | 25,740 | 36.47 → 38.33 (+5.1%) | +72 B | 97.1 → 97.1 |
| V6 forced to SIMPLEC, Metal | 17,822 | 56.35 → 57.75 (+2.5%) | +78 B | 138.2 → 138.4 |
| `elbow_param` SIMPLE, Metal | 45,637 | 87.71 → 90.13 (+2.8%) | +53 B | 147.8 → 148.5 |
| V6 Coupled (`wfgc`), Metal | 17,822 | 63.37 → 63.37 (0) | 0 | 149.0 → 149.5 |

The added memory is 9 floats per cell, 2 per face and 1 per boundary face, allocated only when the pair is active. The 3.54M-cell memory case runs Coupled, so it is unchanged; on V6 Coupled the device bytes are identical. A SIMPLE case of that size would need about 190 to 280 MB more at the measured 53 to 78 bytes per cell.

Mean time per iteration:

| Case, backend | Iterations 6–100 | 101–500 | 601–1000 |
| --- | ---: | ---: | ---: |
| cylinder, Metal, main | 12.61 ms | 12.48 ms | 26.99 ms |
| cylinder, Metal, branch | 12.65 ms | 12.38 ms | 10.42 ms |
| V6 SIMPLEC, Metal, main | 15.83 ms | 15.99 ms (to 300) | |
| V6 SIMPLEC, Metal, branch | 15.88 ms | 15.92 ms (to 300) | |
| `elbow_param`, Metal, main | 24.45 ms | 83.14 ms (to 300) | |
| `elbow_param`, Metal, branch | 24.54 ms | 48.28 ms (to 300) | |
| cylinder, CPU, main | 74.80 ms | 71.87 ms (to 300) | |
| cylinder, CPU, branch | 76.04 ms | 73.12 ms (to 300) | |
| V6 SIMPLEC, CPU, main | 106.3 ms (to 100) | | |
| V6 SIMPLEC, CPU, branch | 109.2 ms (to 100) | | |

While the solution still moves, the pair costs 0.3% on Metal and 1.7% to 2.7% on CPU. Wall times for the full runs: cylinder Metal, 1000 iterations, 18.6 s → 11.5 s; `elbow_param` Metal, 300 iterations, 19.1 s → 12.1 s; V6 SIMPLEC Metal, 300 iterations, 4.80 s → 4.80 s; cylinder CPU, 300 iterations, 21.9 s → 22.3 s.

The later speedup comes from the momentum solves. They stop at a 1e-2 relative reduction or 100 iterations. At main's float floor, the float residual b − Ax is mostly rounding noise and cannot be reduced 100×. The paired residual is compensated and is not noise there. Cylinder, Metal, last-step x-momentum Krylov iterations ([`runs/perf/`](runs/perf/)):

| Iteration | 300 | 500 | 700 | 900 |
| --- | ---: | ---: | ---: | ---: |
| main | 5 | 8 | 50 | 58 |
| this branch | 5 | 6 | 5 | 6 |

With `paired_state` off, the branch reproduces main's cylinder history bitwise ([`runs/perf/cylinder-metal-krylov-float32-option.log`](runs/perf/cylinder-metal-krylov-float32-option.log)).

### Default or opt-in

The pair is the default, for these reasons:

- It costs 2.5% to 5% of device memory and under 3% of time per iteration.
- It changes converged answers by at most 2.6e-5 in error (the cylinder pressure drop). That is the size of main's own CPU-to-Metal spread on `elbow_param` (2.1e-5).
- It turns a hard float32 floor into one about 200× lower.
- It is faster in wall time on any run that reaches main's floor.

Opt-in would leave every SIMPLE and SIMPLEC run of the app unable to meet Fluent's own criteria. Turning it off (`SolverSettings::paired_state = false`) restores main exactly.

## Open items

- **CUDA** is compiled in CI only and was not run. The new kernels use only the portable prelude. Their exact TwoSum steps have no multiplies to contract, the one product error uses an explicit `kfma`, and CUDA builds with `--fmad=false`. Metal compiles in safe math with fast math off.
- **Coupled is not paired.** V6, the nozzle and the 3.54M-cell memory case run Coupled. They are unchanged and keep the float32 floor. V6's floor, 5.0e-5 scaled, is already below its 1e-3 criterion. Pairing Coupled would need the same treatment of its block system.
- **Interior jumps** (fan curves and porous jumps) keep the state in float32. Their face flux is rebuilt from hi alone, so the pair would need a lo linearization of the jump first. MRF zones do run paired: their relative-flux term is applied to F(hi) before the pair's offset, and the V9 MRF gates pass unchanged.
- **Fluent's iteration-1 continuity of 448 kg/s is not reproduced** (see "Initialization"). The cause is undocumented in Fluent: its ∇φ gradient scheme and its Coupled first step. Kutta's smaller denominator makes its scaled residuals larger, so it is the conservative side of the criterion.
- **The slack is a trade-off measured on one case.** 1 ulp of scale was chosen over 16 because 16 slowed the middle of the run and failed the V4 cavity gate. A case whose velocity scale is far below its maximum speed in most of the domain could freeze hi less tightly than the cylinder does. No such case was measured.
