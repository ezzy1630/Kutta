# Coupled convergence: the external-aero continuity plateau, float pairs, and the flow Courant number

**The external-aero plateau is diagnosed but not fixed. Float pairs under
Coupled fix the internal cases whose plateau was the float32 floor. Coupled
now honours the case's flow Courant number.**

- **Plateau on external aero (wing in ground effect, Ahmed body).** Coupled
  stops at raw continuity about 1e-3. That is about 4 times the float floor,
  so the stall is algorithmic. The cause is a limit cycle in the
  Barth-Jespersen momentum limiter (about 130 leading-edge cells flip every
  iteration). The limiter enters the mass-flux face velocity, and the block
  Jacobian does not see it. Freezing the limiter drops continuity to the
  float floor. Starting Coupled from SIMPLEC's converged field, it leaves
  that fixed point and enters the cycle.
- **What I tried and reverted.** A limiter that falls at once and recovers
  slowly ("held") took the wing to the float floor, and with pairs below
  SIMPLEC's raw continuity. On the Ahmed body it only appeared to converge:
  - Under the held limiter its continuity read 6e-4 to 8.5e-4. The same state
    reads 1.76e-3 under the true (fresh) limiter.
  - Its Cd moved from +0.24-0.29% to +0.64-0.71% away from converged
    SIMPLEC.
  - Two Coupled unit tests failed.

  It is reverted (`fd53f8e`). The real fix changes the discretization: make
  the mass-flux face velocity independent of a non-smooth limiter. That would
  move SIMPLEC and every Fluent reference with it, so it needs its own
  validation lane (open item 1).
- **Float-pair state under Coupled** (`04c9e5f`). These cases were held at the
  float32 floor:

  | case | before | after |
  | --- | --- | --- |
  | SST nozzle (raw continuity) | 2.5e-4 | 6e-7 by 3000 iterations |
  | V6 mixing elbow | 7.4e-5 | 3.1e-5 |
  | tet channel (raw) | stalls at 2.8e-3 | 4.0e-4 at 11,000, flattening |

  External aero is unchanged within the plateau's own wander: Ahmed Cd
  within 3e-5 relative, lift 0.24%; wing lift 0.05%.
- **Flow Courant number** (`04c9e5f`, test `a8f5421`). The saved
  `pressure-coupled/courant-number` now sets the Coupled implicit relaxation.
  It used to be fixed at Fluent's default of 200, which every corpus case
  saves anyway. `controls/courant-number` is editable under pseudo time
  (docs/app-protocol.md). Explicit relaxation factors were already honoured.
- **SIMPLE and SIMPLEC are byte-identical to main** on CPU and Metal. The
  quick tier passes. Slow tier: see [Slow tier](#slow-tier).

## Cases and builds

- Main is `e308968`, the merge base, built in Release. Main has since gained
  `a88729c` and `789b535`: the app, test fixture paths and CLI output format,
  with no solver change. The branch merges cleanly. The candidate is this
  branch.
- Runs are on Metal (Apple M4 Pro) unless marked CPU, from a uniform start,
  under `bench_fsae_aero`'s FSAE settings: SST, Fluent's Coupled defaults
  0.5/0.5 and CFL 200.
- Wing: Zerihan/Tyrrell 026 at h/c 0.224. Level 1 has 56,463 cells and level 2
  has 195,581 (`--re 4.62e5 --downforce --moving-ground --tu 0.001 --ratio 1`).
- Ahmed body: 35 deg slant, level 1, 50 um first layer, 403,604 cells, 10 of
  them poor-quality (`--re 2.78e6 --length 1.044 --area 0.112 --walls body
  --lift-z --tu 0.01`).

"Scaled continuity" is Fluent's: the raw imbalance over the largest raw
imbalance of the first five iterations. Coupled's first iterations are far
better balanced than SIMPLEC's, so its scale is about 375 times smaller on the
Ahmed body (external_3d_robustness, open item 1). The same 1e-3 criterion is
therefore much harder for Coupled in raw terms. Fluent defines it the same
way, so I kept it.

## 1. Reproduction (main)

| case | scheme | iterations | scaled continuity | raw continuity | forces |
| --- | --- | ---: | ---: | ---: | --- |
| wing L1 | SIMPLEC | 4000 | 8.9e-7 (met) | 1.07e-4 | lift 1.38446, drag 0.032357 |
| wing L1 | Coupled | 4000 | 2.4e-3 | 1.03e-3 | lift 1.38339, drag 0.032436 |
| Ahmed 35 L1 | SIMPLEC (slow test) | 2000 | met | | Cd 0.3531 |
| Ahmed 35 L1 | Coupled | 2000 | 4.0e-3 | 1.28e-3 | Cd 0.35393, lift 0.21196 |

Both Coupled histories are flat from about iteration 800 (wing) and 500
(Ahmed). The velocity residuals sit at 2e-7 and 1e-6.

## 2. Diagnosis: algorithmic, not the float floor

**The float floor.** With the limiter frozen, Coupled in float32 stops at raw
continuity 2.5e-4 on wing L1. The variant runs that removed the cycle
(`flim`, `lr01`, `lr002`, `first_c`) all stop there too. With float pairs
added, the held-limiter run reached 4.75e-5 by 3000 iterations. So I take
2.5e-4 as the float32 floor of this iteration (inferred from the frozen runs).
The 1.03e-3 plateau is about 4 times that floor, so it is not rounding.

**Ablations** on wing L1, Coupled, with scratch hooks
(`evidence/diag-hooks.patch`, not committed). Rows marked "a different
discretization" also change SIMPLEC's answer, so they locate the cause but
are not fixes.

| change | scaled continuity | raw | lift | reading |
| --- | ---: | ---: | --- | --- |
| none | 2.4e-3 | 1.03e-3 | 1.38339 | plateau |
| turbulence frozen from 600 | 2.6e-3 | 1.09e-3 | 1.37917 | not turbulence |
| `a_P` frozen from 600 | 2.9e-3 | 1.25e-3 | 1.38259 | not the diagonal's lag |
| CFL 20 / CFL 5 | 2.1e-3 / 1.3e-3 | 8.9e-4 / 5.5e-4 | 1.3841 / 1.3844 | not the implicit relaxation |
| pseudo-time scale 0.1 | 2.5e-3 | 1.05e-3 | 1.38399 | not the pseudo step |
| explicit relaxation 0.25 | 1.6e-3 | 5.8e-4 | 1.38339 | scales with it, does not converge |
| block solve to 1e-3, 100 iterations | 3.8e-3 | 1.06e-3 | 1.38461 | not the linear solve |
| standard pressure interpolation | 3.2e-3 | 1.32e-3 | 1.38595 | not the pressure scheme |
| no limiter | non-finite before 50 | | | needs the limiter |
| first-order momentum | 8.1e-3 | 3.5e-3 | 1.4288 | a different discretization; still cycles |
| **limiter frozen from 600** | **5.8e-4** | **2.5e-4** | 1.38245 | **reaches the float floor** |
| face velocity without the limited gradient and the walk fit | 5.8e-4 | 2.3e-4 | 1.43804 | float floor; a different discretization |
| walk fit without the limiter blend | 8.0e-4 | 3.5e-4 | 1.36635 | a different discretization |

**The cycle.** Cell dumps over iterations 1000-1011 (wing L1, second order)
show:
- 136 cells hold half of the imbalance, all on the leading-edge suction peak
  (x/c 0.01-0.05).
- 111-164 cells change their momentum limiter by more than 0.05 between
  consecutive iterations. The worst flips 1, 0, 0, 1.
- The total imbalance cycles between 1.02e-3 and 1.34e-3 with a period of
  about 10 iterations.

The limiter reaches the mass flux twice: through the limited gradient in
`fv_momentum_faces` and through the θ blend of the walk-fit face velocity
(`WalkFaceReconstruction::momentum_values`). The block Jacobian treats both
as constant.

**The fixed point is unstable under Coupled.** SIMPLEC was converged to
2000 iterations (raw 1.07e-4, lift 1.38446), then switched to Coupled:
- With the fresh limiter, raw continuity rose to 9.9e-4 and lift drifted to
  1.38398.
- With the limiter relaxed toward the fresh value by 10% of the gap per
  iteration, both ways, it stayed at lift 1.38446 at the float floor.

SIMPLEC iterates the same discretization stably. Coupled's iteration has a
mode, carried by the limiter, that its Jacobian omits. This is the "mode the
Jacobian under-predicts" that the external-3D report left open.

**Ahmed body.**
- Removing the walk fit's limiter blend reaches 3.7e-4, but lift falls from
  0.212 to 0.157: a different discretization.
- Removing the velocity faces does not help (5.1e-3).
- Freezing the limiter from iteration 1, 100 or 200 diverges.
- Relaxing it symmetrically (it lags when falling as well) diverges.

## 3. The held limiter: tried, measured, reverted

**Variant.** The cell limiter may fall at once but recovers by a fraction r
of the gap per iteration. It was applied only under Coupled, to the momentum
convection and the face velocity alike. Optionally it snaps to the fresh
value after N iterations without a fall. Results at 2000 iterations:

| variant | Ahmed scaled | Ahmed Cd | wing L1 scaled | wing lift | unit tests |
| --- | ---: | ---: | ---: | ---: | --- |
| main (fresh limiter) | 4.0e-3 | 0.35393 | 2.5e-3 | 1.38231 | pass, 75 iterations |
| r 0.1 | 8.5e-4 | 0.35479 | 8.9e-4 | 1.38366 | pass, 200 iterations |
| r 0.02 | 6.4e-4 | 0.35545 | 6.0e-4 | 1.38377 | **fail** |
| r 0.02, snap after 10 | 8.5e-3 | 0.35855 | | | pass |
| r 0.02, snap after 30 | 1.3e-3 | 0.35676 | | | pass |
| r 0.02, snap after 50 | 8.5e-4 | 0.35559 | 1.1e-3 | 1.38331 | pass, 150 iterations |
| r 0.02, snap after 100 | 6.4e-4 | 0.35537 | 3.0e-4 | 1.38379 | pass, 250 iterations |

With float pairs, r 0.02 reached raw 4.75e-5 on the wing by 3000
iterations, below SIMPLEC's 1.07e-4 (lift 1.38377). On the Ahmed body it met
the 1e-3 criterion.

**Why it is reverted.** The monitor assembles with the held limiter, so it
reports the residual of the held discretization. These match only once the
cycle stops.
- **Wing.** I checked this only on the symmetric variant (relaxed 10% per
  iteration both ways from iteration 600, 2000 iterations). Evaluating its
  monitor with the fresh limiter gave the same continuity, 5.8e-4 against
  6.0e-4. On the wing the reduction looks real. The held run ended at lift
  1.38377 with raw continuity 4.75e-5, 0.05% from SIMPLEC's 1.38446. Started
  from SIMPLEC's field, the relaxed-limiter Coupled run stays at 1.38446. So
  the limited discretization appears to have more than one fixed point here.
  I did not check the held run's own state with the fresh limiter.
- **Ahmed body.** At iteration 1900 (snap after 50), the state reads 8.5e-4
  under the held limiter. With the fresh limiter from 1900 on, the same
  state reads 1.76e-3 and returns to main's plateau by 2000 (3.5e-3, Cd
  0.3541). The cycle never stopped. The held limiter lowered the limiter in
  the cycling cells, which converges a more limited, different
  discretization. Its Cd moved away from converged SIMPLEC (0.3531): main's
  plateau samples sit at +0.24% and +0.29%, the held states at +0.64% to
  +0.71%. Reporting "Ahmed meets its criterion" would have been false.
- **Unit tests.** With r 0.02 and no snap, two Coupled tests failed:
  - `FV Coupled reaches SIMPLEC's fixed point…`: the run stopped at 625
    iterations with v 4.7e-4 from SIMPLEC, against a 2e-4 gate. The held
    limiter was still recovering, so the per-iteration residual was small
    while the state still drifted.
  - `FV Coupled converges on a tetrahedral channel`: residual 9.6e-7 at
    iteration 400, against 2e-7.

  A snap from 50 iterations passes both but leaves the Ahmed problem.

## 4. Float pairs under Coupled

Coupled now carries velocity, pressure and face mass flux as float pairs
(`SolverSettings::paired_state`), as SIMPLE and SIMPLEC already did.
- The block RHS is the pair's momentum residual and the compensated
  continuity.
- The update adds into the pairs.
- The Jacobian is unchanged, so this changes only the fixed point's
  precision.

Results to 3000 iterations from the case's start (`kutta solve --approximate
wfgc --backend metal --iterations 3000`). Both cases save convergence
checking off, so they run to the limit.

| case | build | continuity at 500 / 1000 / 3000 | x-velocity at 3000 | k, omega at 3000 | wall | peak RSS |
| --- | --- | --- | ---: | --- | ---: | ---: |
| SST nozzle | main | 2.49e-4 / 2.47e-4 / 2.48e-4 | 3.6e-8 | 4.3e-8, 3.6e-8 | 362 s | 158.0 MB |
| SST nozzle | candidate | 1.38e-6 / 5.83e-7 / 6.18e-7 | 9.1e-10 | 3.2e-8, 2.9e-8 | 362 s | 159.3 MB |
| V6 mixing elbow | main | 8.59e-5 / 8.21e-5 / 7.42e-5 | 2.4e-8 | 1.2e-7, 1.3e-7 | 295 s | 150.7 MB |
| V6 mixing elbow | candidate | 5.22e-5 / 4.02e-5 / 3.10e-5 | 5.7e-9 | 9.1e-8, 1.0e-7 | 297 s | 149.9 MB |

The nozzle's main plateau was the float32 floor. The pairs remove it, and
continuity falls 400 times further. The elbow gains 2.4 times. Its remaining
floor is not continuity's float floor; I did not trace it. Agreement with
Fluent's saved fields is in [Slow tier](#slow-tier).

**Tet channel** (`bench_convergence --case channel --scheme coupled`, Metal,
294,912 cells). The earlier reports record it unconverged under both schemes
at 400 iterations. Run longer:

| build | raw continuity at 2000 / 4000 / 6000 | Cd at 6000 |
| --- | --- | ---: |
| main Coupled | 4.2e-3 / 2.9e-3 / 2.9e-3 (stalled) | 0.019800 |
| candidate Coupled | 2.9e-3 / 7.4e-4 / 4.8e-4 (falling) | 0.019765 |
| candidate SIMPLEC | 1.8e-4 / 7.2e-5 / 8.0e-5 (stalled) | 0.019786 |

Run on under `bench_convergence`'s default 600 s budget:
- Main Coupled stays at raw 2.8e-3 through 15,500 iterations.
- Candidate Coupled reaches 4.4e-4 at 8,000 and 4.0e-4 at 11,000. That is
  6-7 times lower than main, but flattening.
- SIMPLEC reaches 5.0e-5 at 7,000.

Drag at the end: 0.019798 (main Coupled), 0.019767 (candidate Coupled),
0.019793 (SIMPLEC). Neither scheme meets a 1e-3 scaled criterion here:
SIMPLEC's scaled continuity is 0.025 and Coupled's 0.19. The pairs remove
main's precision stall, but the remaining gap to SIMPLEC is the Rhie-Chow
jump on skewed tet wall layers, as the Coupled report found.

**Tet duct** (`duct.py tet`, 9,729 cells, Re 1e5, CPU, 3000 iterations). No
scheme converges. Raw continuity is 2.4e-2 for main Coupled, 1.9e-2 for
candidate Coupled and 2.5e-2 for SIMPLEC, with Cd 0.0165, 0.0165 and 0.0164.
SIMPLEC fails too, so it is not a Coupled defect; I did not pursue it.

**External aero** is unchanged, because its plateau is the limit cycle.
Main and candidate ran concurrently with the same load:

| case | build | scaled continuity | raw | forces | wall | peak RSS |
| --- | --- | ---: | ---: | --- | ---: | ---: |
| Ahmed 35 L1, 2000 iterations | main | 4.02e-3 | 1.28e-3 | Cd 0.35393, lift 0.21196 | 519 s | 1,519 MB |
| | candidate | 3.96e-3 | 1.26e-3 | Cd 0.35392, lift 0.21146 | 529 s | 1,517 MB |
| wing L1, 3000 iterations | main | 2.84e-3 | 1.21e-3 | lift 1.38266, drag 0.032344 | 236 s | 249 MB |
| | candidate | 2.70e-3 | 1.15e-3 | lift 1.38335, drag 0.032433 | 251 s | 252 MB |

The plateau states cycle, so their forces wander with the phase. Two
plateau samples of main's wing differ by as much as main and candidate do
(lift 1.38339 at 4000, 1.38266 at 3000).

**FSAE 2D sweep, level 2, Coupled** (`bench_fsae_aero --coupled
--iterations 3000`; it stops when the residuals are met). Wing settings as
above; S1223 at `--re 2e5 --tu 0.001 --ratio 1`. Iterations given are those
at which the criteria were met:

| case | main: iterations, lift, drag | candidate: iterations, lift, drag |
| --- | --- | --- |
| wing, free air | 850, 0.86123, 0.016139 | 700, 0.86123, 0.016139 |
| wing h/c 0.448 | not met at 3000, 1.14035, 0.021487 | not met at 3000, 1.14036, 0.021484 |
| wing h/c 0.224 | 800, 1.40151, 0.030516 | 750, 1.40149, 0.030517 |
| wing h/c 0.134 | 750, 1.46728, 0.045230 | 800, 1.46714, 0.045222 |
| wing h/c 0.090 | 1300, 1.41017, 0.059070 | 1800, 1.41056, 0.059070 |
| S1223 alpha 4.76 | 850, 1.58469, 0.030068 | 850, 1.58469, 0.030068 |
| S1223 alpha 10 | 1200, 1.93375, 0.047744 | 1200, 1.93387, 0.047738 |

The same 6 of 7 cases meet their criteria on both builds. h/c 0.448
plateaus at scaled continuity 2.0e-3 (main) and 1.8e-3 (candidate). The iteration
counts move with the plateau's phase, by up to 500 in either direction.
Forces change by at most 3e-4 relative. The two sweeps ran at different
times under different machine load, so their wall times are not compared.

## 5. Flow Courant number and the protocol

Fluent User's Guide 36.3 (v25.1) applies the Flow Courant Number with
explicit relaxation to Coupled with the pseudo-time method off. With
Fluent's global pseudo time step, it uses the pseudo-time explicit
relaxation factors instead. Kutta runs Coupled only under pseudo time, and
adds the CFL implicit relaxation a_P(1 + 1/CFL) on top of the pseudo step
(Theory Guide 18.4.4).
- Before, the CFL was fixed at 200. Now the case's saved value is used
  (`SolverSettings::coupled_courant`, default 200). A non-positive or
  non-finite saved value is rejected at setup.
- `controls/courant-number` is editable under pseudo time (min the smallest
  positive float, max 1e15) and active only under Coupled. Without pseudo
  time it stays inactive; its `reason` now says the case runs SIMPLEC.
  docs/app-protocol.md documents this.
- Every corpus Coupled case saves 200, so no corpus result changes.
- `FV Coupled flow Courant number changes the path, not the fixed point`
  checks that the value reaches the solver: CFL 5 and 200 differ by 2.6% in
  u after 5 iterations and converge to the same field.
- The explicit relaxation factors (`pressure-coupled/{mom,pressure}/
  pseudo-explicit-relax`) were already honoured and are unchanged.

**Coupled without pseudo time** (not shipped, measured): the Fluent 18
laminar cylinder and the Dash and PPTX SST elbows save Coupled with pseudo
time off, and Kutta runs them under SIMPLEC. Running them natively under
Coupled with their own CFL and explicit relaxation (a scratch hook, reverted)
gave:

| case | SIMPLEC | native Coupled |
| --- | --- | --- |
| cylinder, 2000 iterations | continuity 1.9e-6, Cd 1.2708606, 22.5 s | continuity 8.3e-7, Cd 1.2708618, 8.0 s |
| Dash elbow, 400 iterations | continuity 4.6e-6, 9.2 s | continuity 6.2e-5, 8.6 s |

Adopting it would change app expectations. The Inspector says a pseudo-off
Coupled case "runs SIMPLEC" (`app/src/ui/Inspector.cpp` couplingNote, about
lines 796-807 and 877-885), and `app/tests/tst_realcli.cpp` checks that text.
Another lane owns `app/`, so it stays a recommendation (open item 2).

## Slow tier

The tier was run in thirds (doctest `--first/--last` 1-11, 12-22, 23-33), main
and candidate concurrently, with `-r=xml --success`.
`tools/forensics/compare_slow_assertions.py` compared each third
(`slow-{main,cand}-*.xml.gz`, `slow-comparison-*.json.gz`):
- 33 cases, **0 regressions, 0 recoveries, 0 missing assertions, 0 new
  assertions**.
- Main and candidate both pass 39,544 assertions and fail 66.
- The failures are the same assertions on both sides:
  - the known 65: V6 25, V2 tet 8, V5 plate 6, V5 bump 9, V7 SA 16, V12
    polynomial fan 1;
  - one load artefact: V8's Metal 80-cell channel has a 180 s wall-clock
    budget. With both tiers sharing the GPU it stopped unconverged at step
    4600 (181.5 s), on both builds. Rerun alone, it converges at step 5700 in
    98 s on both builds, with identical values, and V8 passes 18 of 18.

40 asserted values moved, none across its bound, all in Coupled cases:

| case | values | largest change |
| --- | ---: | --- |
| Fluent k-epsilon elbow and SST nozzle vs converged Fluent fields (CPU, Metal) | 12 | 4.9e-5 relative (nozzle k error 0.0491931 -> 0.0491907, bound 0.08) |
| V6 Fluent mixing elbow fixed point and cold start | 26 | 5.3e-4 relative (a fixed-point field error 0.00116857 -> 0.00116795); the 25 failing assertions stay failing, moving by <= 2.2e-6 absolute |
| Ahmed 35 deg at y+ ~ 1, SIMPLEC and Coupled | 2 | Coupled/SIMPLEC lift 0.60% -> 1.45% (bound 3%), drag 0.089% -> 0.077% (bound 1%) |

The Ahmed change is the Coupled plateau's phase moving with the pairs: the
state is not converged (section 4). The SIMPLE and SIMPLEC cases, the
laminar cylinder (Coupled without pseudo time, run as SIMPLEC), and the FSAE
and TMR cases are unchanged.

Wall time per third: main 70.3, 74.8 and 3.6 min; candidate 58.9, 74.9 and 3.6
min. These are not a timing comparison: main's first third started 25
minutes earlier, under heavier load. Peak RSS: main 3.28, 1.18 and 0.13 GB;
candidate 3.22, 1.18 and 0.14 GB.

## Quick tier, parity, lint, cost

- **Quick tier** (`kutta-tests --test-suite-exclude=slow`, fixtures present,
  `KUTTA_FLUENT_EXAMPLES` set, binary at `fd53f8e`): 260 cases and 7,454,647
  assertions, all passed (234 s, peak RSS 1.92 GB). The final binary (`a8f5421`)
  adds the Courant test: 261 cases and 7,454,665 assertions, all passed (190
  s). Main: 259 cases and 7,454,253 assertions, all passed.
- **Parity** (`validation/results/external_3d_robustness/parity.sh`, 40
  iterations, CPU and Metal): `elbow_simplec` (SIMPLEC) and
  `elbow_param_simple` (SIMPLE) are byte-identical to main in fluxes, forces,
  loads, result.vtu and residuals. `elbow_coupled` and `nozzle_coupled`
  differ, as intended: the float pairs change Coupled's iterates.
- `tools/lint_portability.py` is clean. The new kernels (`coupled_rhs`'s
  `net` input, `coupled_paired_update`, `coupled_paired_gauge`,
  `coupled_paired_pick`) use only `KUTTA_KERNEL` primitives already used by
  the SIMPLEC pair kernels. CUDA compiles only in CI.
- **Memory.** Coupled allocates the pair buffers that SIMPLEC already had:
  velocity_lo and pressure_gradient_lo (6 floats per cell), pressure_lo,
  imbalance and momentum_delta (3 per cell), mass_lo and previous_mass (2 per
  face), and pressure_boundary_lo. Process peak RSS moves by -0.5% to +1.4%
  on the runs above. Wall time under the same load is +0% to +6%.

## Open items

1. **The external-aero plateau.** The limit cycle is in the discretization's
   non-smooth limiter feeding the mass flux. Candidate fixes all change
   SIMPLEC's fixed point as well, so each needs checking against every
   Fluent reference:
   - a differentiable limiter (Venkatakrishnan's, for example) in the face
     velocity;
   - an unlimited, linear face velocity in the Rhie-Chow flux, as in the
     textbook scheme;
   - a smooth θ in the walk-fit blend.

   I recommend that as a separate lane, starting from the wing L1 dumps
   above.
2. **Native Coupled without pseudo time.** It is faster on the cylinder,
   slower to converge on the Dash elbow, and needs the app's couplingNote and
   `tst_realcli` updated. This is an app-lane decision.
3. **The V6 elbow's remaining floor**, about 3e-5 under pairs, is not traced.
4. **Tet duct and tet channel are not fixed.** No scheme converges the duct.
   On the channel, Coupled with pairs flattens near raw 4e-4 by 11,000
   iterations, 8 times above SIMPLEC (section 4).
5. **Stale app text (app lane; I did not edit `app/`).** Two places still say
   Kutta fixes the flow Courant number at 200:
   - `app/src/ui/Inspector.cpp:875-885`, the note under
     `controls/courant-number`;
   - line 809, the Coupled coupling note ("Fluent's default flow Courant
     number").

   With this change the case's value applies and the setting is editable.
   Both should say so. The setting's `editable` flag already reaches the app
   through the protocol.
6. **Non-uniqueness (inferred).** On wing L1 the limited discretization
   appears to have at least two fixed points, 0.05% apart in lift (1.38446
   and 1.38377). Which one a run lands on depends on its path.

## Reproduce

```sh
# external aero, main vs candidate (scripts as in evidence/commands.md)
bench_fsae_aero --mesh zerihan_h0.224_L1.msh --re 4.62e5 --downforce --moving-ground --tu 0.001 --ratio 1 --settle -1 --coupled --iterations 3000 --out w
bench_fsae_aero --mesh ahmed35_L1.msh --re 2.78e6 --length 1.044 --area 0.112 --walls body --lift-z --tu 0.01 --window 1000 --settle -1 --coupled --iterations 2000 --out ah
# internal Coupled
kutta solve "$KUTTA_FLUENT_EXAMPLES/mixing_elbow/mixing_elbow.cas.h5" --approximate wfgc --backend metal --iterations 3000 --out elbow
kutta solve "$KUTTA_FLUENT_EXAMPLES/flprj_3891/nozzle_para_named_1.cffdb/nozzle_para_named_1-Solve/Base DP/nozzle_para_named_1.cas.h5" --approximate wfgc --backend metal --iterations 3000 --out nozzle
bench_convergence --case channel --backend metal --scheme coupled --iterations 20000 --interval 500
# parity, then the slow tier per build in thirds (doctest --first/--last), compared with
validation/results/external_3d_robustness/parity.sh <main kutta> <candidate kutta> runs/parity
tools/forensics/compare_slow_assertions.py slow-main.xml.gz slow-cand.xml.gz > slow-comparison.json
```

`evidence/` holds every history named here (`history.csv` / `residuals.csv`
and run logs) and the scratch diagnosis hooks as a patch against main. It
holds no meshes and no fixture data.
