# Tet robustness: second-order pressure on wall cells with one-plane neighbours

Branch `t3/a057564e`, based on main `742ef86`. Commits:

| Commit | Change | Status |
| --- | --- | --- |
| `c9179bb` | pressure-scheme override; Coupled Jacobian | kept |
| `2df4adf` | `p_c` rows at every unprescribed face | superseded |
| `b8390e7` | inlet rows only where the extrapolation amplifies | superseded |
| **`80535e5`** | rows at any unprescribed face only where the extrapolation amplifies | the fix |

`2df4adf` and `b8390e7` stay in the history because the slow tier is what showed why each was wrong ([Diagnosis](#diagnosis)). Evidence is for `80535e5` unless a row says otherwise.

Hardware: Apple M4 Pro, CPU reference backend and Metal; CUDA was not run. Fluent fixtures are read through `KUTTA_FLUENT_EXAMPLES` and never committed.

## Result

| | main `742ef86` | `80535e5` |
| --- | --- | --- |
| `elbow_param`, second-order pressure, CPU and Metal | non-finite within 30 iterations | continuity 1.1e-5 by about 200 iterations, flat after |
| `elbow-2`, second-order pressure, CPU and Metal | non-finite within 30 iterations | the same, and flat through 1,200 iterations (Metal) |
| SST tet channel, SIMPLEC, CPU, 1000 iterations | diverges: raw continuity 0.004 → 0.42 | bounded: 8.9e-4 at 1,000, momentum still falling |
| SST tet channel, Coupled, CPU, 400 iterations | stalls at raw 0.025 | stalls at raw 0.015 (open) |
| FSAE tet duct (FSAE aero lane), SIMPLEC, CPU | non-finite at iteration 60 | bounded through 2,000 iterations, continuity stalled at about 1.5 (open) |
| V6, nozzle, standard `elbow_param`: 40-iteration outputs, CPU and Metal | — | **bitwise identical** |
| V6, cold, 500 iterations, CPU | — | **bitwise identical** |
| Rest of the slow tier | — | [below](#tests): one regression, the V2 tet velocity order |

**The fix fails one existing gate.** The V2 tet finest-pair velocity order (8→16) falls from 1.902 to 1.815 against a 1.9 bound, on CPU and Metal. Yet the velocity error falls 18–23% at every level, and a 32-level run shows second order. That run used `2df4adf`, which held rows at every face; it is within 0.2% of `80535e5` at n=8 and 16.

Every other variant tested fails in some way:
- it diverges on the elbow;
- it converges to a wrong pressure field;
- it fails other tet gates; or
- it breaks well-conditioned cases (tables below).

The gate is unchanged. Whether to accept the fix with that failure, or keep main's behaviour on the second-order path, is an owner decision; see [Decision needed](#decision-needed).

## Reproduction

`methods/pressure-scheme` (new, `c9179bb`) switches a case's Fluent `pressure/scheme` between `standard` and `second-order`:

```sh
echo '{"v":1,"overrides":{"methods/pressure-scheme":"second-order"}}' > so.json
E=$KUTTA_FLUENT_EXAMPLES/mixing_elbow
build/kutta solve $E/elbow_param.cas.h5 --setup so.json --approximate compressibility,ke-options \
  --backend metal --iterations 400 --out out
```

Main, with only that setting added, stops before the first 20-iteration monitor on both elbows and both backends: `solve: nonfinite input vector` ([`reproduction/`](evidence.tar.gz)).

The k-epsilon elbows are the only cases in the corpus with Fluent data on this mesh. The coverage lane made them run through the standard scheme's wall rows ([case_coverage REPORT](../case_coverage/REPORT.md), open limit 1).

No other test in `tests/` is excluded for diverging on tets. The second problem is the bench SST tet channel (`bench_convergence --case channel`, CPU, 4,608 cells). The [Coupled report](../coupled/REPORT.md#tetrahedral-channel-bench---case-channel) found it diverging under both schemes and attributed this to the Rhie-Chow jump's gradient terms. On current main it still diverges under SIMPLEC and stalls under Coupled (table above).

## Diagnosis

### The extrapolation loop

On the second-order path, the pressure LSQ holds prescribed outlet rows only. Every other boundary face is omitted, because its value is extrapolated (`p_c + G·r_b`) and so "carries no information". The loop that follows:
- A tet wall cell fits its gradient from its three interior neighbours alone.
- Its own wall pressure is then an extrapolation outside that stencil.
- Momentum takes the wall pressure through its Green-Gauss source, and Rhie-Chow takes the gradient.
- The next pressure correction closes the loop.

`tools/forensics/pressure_boundary_rows.py` measures the fit ([`forensics/`](evidence.tar.gz)). The gain `sum_j |r_b · M⁻¹ w_j d_j|` is the factor by which the extrapolated boundary pressure can amplify neighbour differences. At gain ≤ 1 the extrapolation stays within `p_c ± max|Δp_j|`.

| Main's fit (outlet rows only) | λmin/λmax < 1e-5 | < 1e-2 | gain median | max | gain > 1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `elbow_param` (45,637 tets): 5,584 wall faces | 57 cells | 327 cells | 1.58 | 276 | 5,510 |
| `elbow_param`: 228 inlet faces | | | 1.56 | 3.62 | 172 |
| V6 `mixing_elbow` (poly): 2,436 wall faces | 0 | 0 | 0.47 | 0.52 | 0 |
| V6: 2,001 symmetry faces | | | 0.40 | 0.76 | 0 |
| V6: 207 inlet faces | | | 0.41 | 0.60 | 0 |

The blow-up starts at the worst cells. After one iteration the fastest cell is at 230 m/s, and after two at 1,726 m/s, with 10 and 25 m/s inlets; it is Fluent cell 2991, with ratio 4.6e-6 and gain 276.

**But conditioning is not the whole mechanism.** Closing only the badly conditioned cells does not stop it; wall cells with gains of 1–20 keep the wall layer unstable or wrong:

| Scratch candidate on `elbow_param`, second-order | Cells changed | Outcome |
| --- | ---: | --- |
| main (no rows) | — | non-finite |
| drop eigen-directions below 1e-3 / 1e-2 λmax | 33 / 272 | non-finite / diverged (continuity 5e8 by 100) |
| `p_c` rows where λmin/λmax < 1e-3 / 1e-2 | 33 / 272 | non-finite / diverged |
| `p_c` rows at walls where gain > 5 / > 3 | 279 / 379 | diverged |
| `p_c` rows at walls where gain > 2 | 1,246 | converges, but pressure L2 31.8% and drop 47.8% against Fluent: a wrong field |
| `p_c` rows at any unprescribed face where gain > 1.5 (Metal) | 3,554 faces | converges, pressure L2 10.1%, drop 13.7% |
| **`p_c` rows at any unprescribed face where gain > 1 (`80535e5`)** | 5,682 faces | **converges, pressure L2 5.99%, drop 4.02%** |
| `p_c` rows at any unprescribed face where gain > 0.75 (Metal) | 5,722 faces | converges, pressure L2 5.96%, drop 3.68% |
| `p_c` rows at every unprescribed face (`2df4adf`) | all | converges, pressure L2 5.78%, drop 2.76%, but biases well-conditioned faces (below) |
| Barth-Jespersen limit on `∇p`, no rows | all | non-finite |
| rows only in the p′ correction gradient, only for face pressures, or only for the Rhie-Chow jump | all | non-finite |
| `p_c` rows at walls only (normal foot) | all walls | diverged by 400 |
| walls and inlets, no symmetry (foot or face) | all | converges; fails two quick tet gates |

So the rows are needed in every consumer of the pressure fit: momentum face pressure, wall extrapolation, the Rhie-Chow jump and p′. Field accuracy plateaus once every face with gain above about 1 has its row, and degrades progressively as faces with gain above 1 are left without one. The scratch switches that produced these tables are not in the branch; their logs are in [`evidence.tar.gz`](evidence.tar.gz) under `variants/`.

The tet channel has no ill-conditioned cell (minimum ratio above 1e-2, maximum gain 2.3), yet the same rows remove its SIMPLEC divergence. Its instability is the same wall-layer extrapolation loop, not a singular fit.

Inferred, not measured: the Rhie-Chow self-coupling sign flips that the Coupled report measured on this mesh come from the same wall-cell gradients, since the jump on wall-adjacent faces is built from them.

### Why not rows everywhere (`2df4adf`, `b8390e7`)

A `p_c` row imposes a zero normal pressure derivative at its face. Where the face carries a real normal gradient, that biases the cell's gradient. The full slow tiers of `2df4adf` and `b8390e7` showed it:

| Slow gate (CPU and Metal agree) | main | `2df4adf` (rows everywhere) | `b8390e7` (rows at walls, symmetry; inlets by gain) | `80535e5` |
| --- | ---: | ---: | ---: | ---: |
| FSAE V8 Darcy-Forchheimer pressure Linf (gate 7e-3): uniform Darcy gradient at a velocity inlet | 1.9e-5 | **1.54e-2, fails** | 1.8e-5 | 1.9e-5 (main's value) |
| FSAE V7.2 coarse hex Taylor-Couette pressure L2 (gate 1e-2): radial gradient at a rotating wall | 0.31% | 0.92% | 0.92% | 0.31% (main's value) |
| FSAE V9 solid-body rotation pressure L2 (gate 1e-2) | 0.07% | 0.31% | 0.31% | 0.07% (main's value) |

The Taylor-Couette error tripled to 92% of its gate. This is the error Fluent's own Theory Guide names for its standard scheme's zero wall gradient: valid "for boundary layers, but not in the presence of body forces or curvature" (§18.4.1, [mirror](https://www.afs.enea.it/project/neptunius/docs/fluent/html/th/node371.htm)). FSAE rotating wheels are such walls.

Inlet treatments tried on `b8390e7`'s way to the gain rule (`variants/inlet/`, `elbow-2` on Metal unless noted):

| Inlet treatment | Outcome |
| --- | --- |
| no inlet rows | `elbow_param` continuity 52 at 400 |
| row holding momentum's lagged extrapolation `p_c + G_prev·r_b` (unbiased at a fixed point) | grows from iteration 160, non-finite by 1,200; fails three quick tests (Coupled tet channel, SIMPLE fixed point, monitor-cadence invariance) |
| `p_c` rows only in the p′ correction fit | diverged by 200 |
| `p_c` row where λmin/λmax < 1e-2 | diverged by 300 |
| `p_c` row where inlet gain > 0.5, 1, 1.5 or 2 | converges |

`80535e5` applies that gain rule to every unprescribed face, computed from main's own fit.

V6 has no face with gain above 1 (forensic above), and V6, the nozzle, V8 and the annuli reproduce main exactly. On the elbow, 98% of the unprescribed faces keep their rows. It converges to within 0.2 points of the all-rows field L2s; the pressure drop is 4.0% against 2.8%. The window that keeps V6 bitwise and the elbow accurate is narrow: V6's largest gain is 0.76, and the elbow degrades above about 1. The threshold 1 is the boundedness limit, not a fit to these cases.

## How Fluent handles it

- **Least squares cell-based gradient.** Fluent 12.0 Theory Guide §18.3.3 ([mirror](https://www.afs.enea.it/project/neptunius/docs/fluent/html/th/node368.htm)): the Gram-Schmidt decomposition produces weights "for each of the faces of cell c0". The section does not mention boundary faces, so the boundary rows are not documented explicitly.
- **What Fluent stores at those faces.** Fluent's saved V6 data has `SV_P` on every face. At walls, symmetry and inlets it equals the cell pressure exactly ([momentum forensics](../../../docs/research/fluent-momentum-forensics.md), item 8).
- **Wall gradient.** TG §18.4.1 says, of the standard scheme, that Fluent "assumes that the normal pressure gradient at the wall is zero", and that this "is valid for boundary layers, but not in the presence of body forces or curvature".
- **Whether Fluent's fit uses the stored values: likely on V6.** Fluent's saved V6 state is converged with second-order pressure. Take its momentum residual with momentum's boundary pressure extrapolated; it is lowest when the pressure LSQ holds the stored values as rows. `pressure_boundary_rows.py`, section 1:

  | Pressure LSQ boundary rows, Fluent's saved V6 state | Wall cells Σ\|R\| | Symmetry cells | Inlet cells | All |
  | --- | ---: | ---: | ---: | ---: |
  | every boundary face, `p_b = p_c` at the face centroid | **1.429e-4** | 1.476e-4 | 6.79e-6 | 9.745e-4 |
  | walls at the normal foot | 1.446e-4 | 1.477e-4 | 6.76e-6 | 9.762e-4 |
  | walls, symmetry, inlets at the normal foot | 1.445e-4 | 1.440e-4 | 7.01e-6 | 9.718e-4 |
  | walls and symmetry at the centroid, no inlet rows | 1.433e-4 | 1.478e-4 | 7.90e-6 | 9.758e-4 |
  | no wall rows | 2.351e-4 | 1.495e-4 | 8.01e-6 | 1.087e-3 |
  | outlet rows only (main's second-order path, and `80535e5` on V6) | 2.376e-4 | 1.493e-4 | 9.16e-6 | 1.089e-3 |

  **So on V6, `80535e5` keeps main's treatment, which Fluent's residual disfavours.** Yet Kutta's converged V6 fields were further from Fluent's with rows than without them (`2df4adf`, [below](#well-conditioned-meshes)), and rows bias the analytic cases above. The residual evidence and the field evidence disagree, and the slow gates decide. This is open limit 3.
- **Second-order pressure on bad meshes.** TG §18.4.1 says the second-order scheme "may have some trouble if it is used at the start of a calculation and/or with a bad mesh". Fluent's elbow tutorials use the standard scheme on this mesh. There is no Fluent second-order run on it to compare against.
- **Rhie-Chow.** Fluent's jump uses the momentum equation's own pressure gradient (Green-Gauss of second-order face pressures, forensics coefficient 0.998). Kutta keeps its LSQ jump; that choice was rejected earlier on V2 tet pressure order and is not changed here.

## The fix

`80535e5`, in `src/kutta/fv/steady_solver.cpp`, applies to second-order pressure only:

- Wall, symmetry and velocity-inlet faces start omitted, as on main.
- `boundary_extrapolation_gain()` (new, `operators.cpp`) computes each face's gain from that fit. Where it exceeds 1, the face becomes a face-centroid row holding the cell pressure. That value is set by `fv_cell_pressure_boundary`, formerly the standard scheme's `fv_standard_pressure_boundary`.
- Momentum's boundary pressure stays `p_c + G·r_b`.
- The standard scheme keeps its normal-foot rows.
- `boundary_extrapolation_gain()` and `boundary_least_squares()` share one per-cell helper with unchanged arithmetic.

`c9179bb` also stops the Coupled Jacobian adding unprescribed boundary rows to `d g_c/d p_c`, since their difference is identically zero; this is bitwise neutral on main. It also adds `methods/pressure-scheme`.

The rule is not documented by Fluent. Where it adds rows, the value is Fluent's stored value. Where it does not, the result is main's, bit for bit. It is an undocumented stabiliser that activates only on cells whose extrapolation amplifies, as the brief allowed.

## Previously diverging cases

Each run is a `kutta solve` cold start from the case's hybrid initialization: 400 iterations, `--approximate compressibility,ke-options`, plus the setup override. Scaled residuals at iteration 400 are flat from about 200 on.

Errors come from `bench_fluent_reference` against Fluent's saved fields:
- the converged standard-pressure `elbow_param` (Fluent iteration 89);
- the unconverged `elbow-2` snapshot at iteration 50.

The schemes differ, so these bound agreement rather than measure a same-scheme error.

| Case | Backend | Continuity | Momentum | k / epsilon | Velocity L2 | Pressure L2 | k L2 | epsilon L2 | Pressure drop | Wall force |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `elbow_param` | CPU | 1.08e-5 | 2.2e-8 | 3.1e-8 / 3.1e-8 | 1.05% | 5.99% | 8.71% | 16.9% | 4.02% | 1.20% |
| `elbow_param` | Metal | 1.11e-5 | 2.2e-8 | 3.1e-8 / 3.2e-8 | 1.05% | 5.99% | 8.71% | 16.9% | 4.03% | 1.20% |
| `elbow-2` (vs Fluent it 50) | CPU | 1.12e-5 | 2.2e-8 | 3.1e-8 / 3.3e-8 | 4.59% | 6.67% | 15.0% | 15.6% | 8.23% | 0.36% |
| `elbow-2` (vs Fluent it 50) | Metal | 1.11e-5 | 2.2e-8 | 3.1e-8 / 3.2e-8 | 4.59% | 6.67% | 15.0% | 15.6% | 8.23% | 0.36% |

For scale:
- Kutta's standard-scheme `elbow_param` run (Fluent's own scheme) gives, by the same bench: velocity 0.88%, pressure 5.82%, k 9.10%, epsilon 17.6%, pressure drop 1.01% and wall force 0.32%.
- With rows at every face, `2df4adf` gave `elbow_param` drop 2.76% and `elbow-2` 6.92%; the other columns agree to within 0.2 points.

Solve time for 400 iterations: CPU 166–169 s and Metal 59–60 s. All four solves ran at once, alongside this branch's slow and quick tiers and another lane's test run.

The slow test's own conventions (`solve_and_compare`, which takes the pressure drop from Kutta's boundary pressures) give, on both backends: velocity 1.05%, pressure 5.99%, k 8.71%, epsilon 16.9%, pressure drop 1.11% and wall force 0.41%.

New gates:

- `Second-order pressure keeps the k-epsilon elbow tet wall layer bounded` (quick, CPU and Metal): 30 iterations; finite; maximum speed below 50 m/s (main reaches 3,600 m/s); continuity below 0.1.
- `Fluent k-epsilon elbow converges with second-order pressure` (slow): 400 iterations, with continuity < 1e-4 and momentum < 1e-6. Its field, pressure-drop and wall-force bounds were set at about 1.5× the `2df4adf` values, and `80535e5` is within them: pressure 5.99% against 8.7%, drop 1.11% against 1.7%.

## Well-conditioned meshes

These are bitwise identical to main ([`parity/`](evidence.tar.gz)), comparing `fluxes.csv`, `forces.csv`, `loads.json`, `result.vtu`, and `residuals.csv` without its timing column:
- the 40-iteration outputs of V6 (Coupled, `--approximate wfgc`), the SST nozzle, and standard-pressure `elbow_param`, on CPU and Metal;
- V6 at 500 cold iterations on CPU.

In the slow tier ([Tests](#tests)), no second-order case outside V2 tet changes beyond the eighth digit, apart from the timing-limited TMR level-2 runs.

V6 with rows everywhere (`2df4adf`, cold, CPU, 500 iterations), for the record:

| V6 vs Fluent | main = `80535e5` | `2df4adf` |
| --- | ---: | ---: |
| velocity L2 | 3.37% | 3.61% |
| pressure L2 | 2.66% | 2.70% |
| k L2 | 14.7% | 15.0% |
| omega L2 | 2.80% | 2.91% |
| elbow wall force | 0.75% | 0.67% |
| pressure drop | 1.90% | 1.91% |
| continuity at 500 | 9.3e-5 | 5.2e-5 |

**V2** manufactured solution, tet: velocity and pressure L2 errors on CPU. Metal agrees to the printed digits. The n=32 values are from `2df4adf`, run on Metal. It also held rows at this case's faces with gain ≤ 1; at n=4–16 the two builds differ by at most 3% (n=4) and 0.2% (n=8, 16).

| Tet n | main u | `80535e5` u | main p | `80535e5` p |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 0.1388 | 0.1065 | 0.2050 | 0.1461 |
| 8 | 0.04437 | 0.03407 | 0.1131 | 0.08184 |
| 16 | 0.01187 | 0.00969 | 0.05589 | 0.04187 |
| 32 (`2df4adf`) | 0.002811 | 0.002480 | 0.02589 | 0.02106 |
| velocity order 8→16 (gate ≥ 1.9) | **1.902** | **1.815** | | |
| velocity order 16→32 | 2.078 | 1.965 (`2df4adf`) | | |
| pressure order 8→16 (gate ≥ momentum order − .1 = .914) | | | 1.017 | 0.967 |

The error falls 18–29% at every level, but the improvement shrinks with refinement, so the 8→16 order drops below the gate. The 16→32 order is second order.

The three wall-force order assertions fail on both builds (gate 1.8):

| | main | `80535e5` |
| --- | ---: | ---: |
| wall-force order 1 | 1.05 | 1.05 |
| wall-force order 2 | 0.81 | 0.84 |
| wall-force order 3 | 0.55 | 0.62 |

V2 hex is unchanged.

V2 tet variants, finest velocity order 8→16 (CPU). Every variant that keeps the elbow stable fails the 1.9 gate:

| Variant | Order 8→16 | Other outcome |
| --- | ---: | --- |
| wall rows only | 1.882 | |
| symmetry rows only | 1.851 | |
| every face at the normal foot | 1.819 | |
| every face at the centroid | 1.812 | |
| walls only at the normal foot | 2.05 | diverges on the elbow |
| walls at the normal foot, plus inlet rows | | fails quick `FV laminar tet manufactured solution` (u 0.227 > 0.18) and `FV Coupled converges on a tetrahedral channel` |

## Tests

**Quick tier** (`kutta-tests --test-suite-exclude=slow`, CPU and Metal):

| | test cases | assertions | failed |
| --- | ---: | ---: | ---: |
| main `742ef86` | 249 | 7,417,189 | 0 |
| `80535e5` | 250 | 7,417,198 | 0 |

The extra case is the bounded-wall-layer test: one setup REQUIRE plus three checks per backend. The other two assertions check the `methods/pressure-scheme` setting. `b8390e7` gave the same counts. The rejected lagged-inlet variant failed 9 (`variants/inlet/`).

**Slow tier** (`kutta-tests --test-suite=slow -r=xml --success`), compared assertion by assertion with `tools/forensics/compare_slow_assertions.py`. Main and `2df4adf` ran concurrently; `b8390e7` and `80535e5` each ran alongside other work.

| Slow case | main failed / total | `80535e5` | Note |
| --- | ---: | ---: | --- |
| FV V2 manufactured Navier Stokes tet refinement | 6 / 18 | **8 / 18** | the regression: `ou>=1.9`, 1.902 → 1.815, CPU and Metal |
| V7 SA TMR flat plate and bump | 17 / 70 | 16 / 70 | bump level-2 `cf_error` 0.105 → 0.097: timing (below) |
| FV V6 Fluent mixing elbow fixed point and cold start | 25 / 2,565 | 25 / 2,565 | the 25 known failures; no value changed |
| V5 TMR flat plate SST | 6 / 36 | 6 / 36 | no value changed |
| V5 TMR bump SST | 9 / 47 | 9 / 47 | level-2 values changed: timing (below) |
| FSAE V12 native Fluent fan import | 1 / 1 | 1 / 1 | unchanged; missing optional fixture |
| Fluent k-epsilon elbow converges with second-order pressure (new) | — | 0 / 11,220 | |
| FSAE V7.1 translating-wall Couette | 0 / 104 | 0 / 104 | poly: 4 values move in the eighth digit |
| 19 other cases (Fluent elbow, nozzle, cylinder; V1–V4, V8; FSAE V7.2–V11; research downloads) | 0 | 0 | no value changed |

Changed values against main:
- V2 tet: 14 values, the expected ones.
- V5 bump and V7 SA bump, level 2 only: 15 values.
- FSAE V7.1 poly Couette: 4 values, moving in the eighth digit (velocity error 1.22228e-5 → 1.22077e-5).

Elbow, nozzle and cylinder are identical. Their assertions show as missing/new in the script's output only because this branch adds lines to `test_fluent_coverage.cpp`, shifting the source-line identity; compared in order, all 17,044 and 207 values match.

**The level-2 TMR changes are timing, not arithmetic.**
- The new rule converts 9 boundary faces on bump levels 0 and 1, and none on level 2 (56,320 cells). It was measured with a scratch print, since removed.
- Each level runs independently, so level 2 executes main's arithmetic.
- Level 2 stops on its 600 s wall-clock budget (`budget=1`). Main stopped at 2,000 steps. `80535e5` got 2,100 steps on both cases, and `2df4adf`'s run got 2,000 on V5 but 2,100 on V7. The extra steps move every level-2 value slightly, mostly towards the reference.
- Levels 0 and 1 assert no values.

The brief flagged V7 SA bump level 2 as timing-sensitive; the V7 recovery should not be credited to this change.

Superseded commits, same comparison:

| | V2 tet `ou` | FSAE V8 Darcy | FSAE V7.2 / V9 / V9.3 Taylor-Couette `ep` | V6, nozzle, V4, other values |
| --- | --- | --- | --- | --- |
| `2df4adf` | fails (1.812) | **16 failures** (pressure Linf 1.54e-2) | pass, at 0.92% against a 1% gate (main 0.31%) | changed throughout |
| `b8390e7` | fails (1.812) | pass | same as `2df4adf` | changed throughout |
| `80535e5` | fails (1.815) | pass, main's values | main's values | unchanged |

## FSAE aero lane reproductions

Source: `t3/620d7176` at `96b30c9`, used read-only (its `validation/results/fsae_aero/REPORT.md`, "Solver problems found").
- Its `bench_fsae_aero` and `validation/cases/fsae_aero.hpp` were built unchanged in scratch copies of `742ef86`, `2df4adf`, `b8390e7` and `80535e5`.
- That branch changes no solver source, so the `742ef86` build reproduces its results.
- Meshes came from its `tools/fsae_aero/` scripts, with its gmsh 4.15.2 environment.
- Logs: `fsae_aero/` in the evidence archive.

| Case | main `742ef86` | `80535e5` | In scope? |
| --- | --- | --- | --- |
| 1. gmsh tet duct (`duct.py tet`, 9,729 cells, SST Re 1e5, CPU) | non-finite at iteration 60 | **no divergence through 2,000 iterations**: one excursion at iterations 300–500 (scaled continuity up to 79, max speed 8.1), then drag 0.01618–0.01625, momentum about 1.3e-4, continuity stalled at 1.2–1.6 | yes: divergence fixed, convergence open |
| same duct in hexahedra (`duct.py hex`, 2,000 cells) | max residual 2.1e-4 at 300, drag 0.014525 | identical | untouched |
| 2. Ahmed L1, 15 µm first layer (Metal, 50 iterations) | diverges: drag −2.3e8 at 50 | **still diverges**: drag −3.4e5 at 50 | not fixed (below) |
| 3. Zerihan 2D wing, h/c 0.224, L2, `--coupled` (Metal) | non-finite Coupled residual before iteration 20 | same | not fixed (below) |
| Zerihan 2D wing, SIMPLEC (regression check, Metal, 400 iterations) | lift 1.29318, drag 0.039030 | lift 1.29310, drag 0.039030 | unchanged in substance |

1. **Tet duct.** The divergence is the wall-layer extrapolation loop diagnosed above, and the rows remove it. The remaining continuity stall matches the SST tet channel (open limit 2): bounded, forces settled, not converged. Both are tet meshes under SIMPLEC with walls on four sides.

   Rows at every wall (`b8390e7`) were calmer here: continuity 0.4–0.7 and drag 0.01611–0.01615. This, and the elbow pressure drop (4.0% against 2.8%), are the measured costs of leaving faces with gain at most 1 without rows.
2. **Ahmed thin layer.** This is a hex mesh (tets and prisms split into hexes) whose wall cells already have interior neighbours on every side. No boundary-row variant stops its divergence; the fault lies elsewhere on the path through cells of aspect ratio about 2,700 at a 90° convex edge. The candidates are:
   - the Rhie-Chow jump's skew extrapolation (`fv_pressure_faces` fades it beyond one link length, a Kutta choice);
   - the turbulent walk-weighted face projection.

   Neither was reduced here, so this stays open rather than widening this lane.
3. **Coupled on 2D wing meshes.** Reproduced: the Coupled block solve goes non-finite within 20 iterations on every build. Not reduced. It is in the Coupled solver, not in the fits this change touches.

## Decision needed

Recommendation: **accept `80535e5`**, and record the V2 tet 8→16 velocity-order assertion as a known failure, with this report as the reason. The case for it:

- It is the only treatment found that is stable and accurate on the tet elbow while leaving well-conditioned cases bitwise unchanged.
- It lowers the V2 tet error at every level and is second order at 16→32.
- It keeps every other slow gate that passes on main (below).
- No threshold, budget or tolerance was changed.

The alternatives:
- Drop `2df4adf`–`80535e5` and keep `c9179bb`. That leaves the second-order path divergent on tet meshes like the elbow; real FSAE meshes from Fluent Meshing can have the same wall cells.
- Re-base the V2 tet order gate on 16→32. That would change a gate and is the owner's call.

## Open limits

1. **V2 tet 8→16 velocity order** (above).
2. **Tet channel and tet duct do not converge.** The channel under Coupled stalls at raw continuity 0.015 (main 0.025). Under SIMPLEC it no longer diverges but converges slowly (8.9e-4 at 1,000). The FSAE tet duct stalls at scaled continuity 1.2–1.6 (0.4–0.7 with rows at every wall). The Rhie-Chow jump with Fluent's momentum gradient is the documented candidate.
3. **Fluent's boundary rows on well-conditioned meshes.** Fluent's saved V6 residual favours `p_c` rows at every face, but Kutta's converged fields and the analytic gates favour none there. A Fluent second-order case with a rotating wall or a strong inlet pressure gradient would settle it.
4. **Fluent's limiter on the pressure gradient.** Forensics slightly favour a Barth-Jespersen-limited pressure reconstruction (wall cells 1.429e-4 → 1.386e-4). It is not implemented.
5. **No same-scheme reference.** No corpus case runs second-order pressure on a tet mesh in Fluent. The elbow comparison is against Fluent's standard scheme.
6. **FSAE thin wall layers at sharp edges, and Coupled on 2D wing meshes** (FSAE aero lane items 2 and 3): reproduced, unchanged by this fix, open.
7. CUDA was not run. The kernel change is a rename and a comment; the rest is host code.

## Evidence

[`evidence.tar.gz`](evidence.tar.gz), with absolute paths replaced by `<worktree>`, `<main-742ef86>`, `<scratch>` and similar:

- `reproduction/`: main with only the setting added, both elbows, CPU and Metal.
- `runs/`: `80535e5`'s elbow solves (logs, `bench_fluent_reference` output, residuals, `summary.json`, forces). `runs_2df4adf/` and `runs_b8390e7/` hold the same for the superseded commits.
- `variants/`: the scratch candidates in the diagnosis tables.
  - `all_faces/` holds the gain-threshold sweep.
  - `inlet/` holds the inlet treatments.
  - The eigenvalue-ratio-gated wall logs were overwritten on this case-insensitive filesystem, so their outcomes above come from the session's console output, as does the "no inlet rows" elbow run.
- `v2/`, `channel/`, `v6/`, `forensics/`, `parity/`, `fsae_aero/`: the runs quoted above.

Beside this report: `slow-*.xml.gz`, the full slow tiers of `742ef86`, `2df4adf`, `b8390e7` and `80535e5`; and `slow-compare-742ef86-*.json.gz`, their comparisons with main. Absolute paths in both are replaced the same way.

## Reproduce

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build
export KUTTA_FLUENT_EXAMPLES=$PWD/validation/fluent_examples
build/kutta-tests -tc='Second-order pressure keeps*'                 # quick
build/kutta-tests -tc='Fluent k-epsilon elbow converges with second-order pressure,FSAE V8 Darcy*,FSAE V7.2*'   # slow
python -I tools/forensics/pressure_boundary_rows.py validation/fluent_examples/mixing_elbow mixing_elbow /tmp/cache
build/bench_convergence --case channel --scheme simplec --backend cpu --iterations 1000
```
