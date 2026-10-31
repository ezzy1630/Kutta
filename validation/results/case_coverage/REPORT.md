# Case coverage: more Fluent pairs through `kutta solve`

Branch `t3/ad404a5c`, based on `c24da6d`. Evidence commit: `cd68d0c` (all solves, comparisons and tests below).
Apple M4 Pro, CPU reference backend and Metal; CUDA was not run. The corpus is the
43 case/data pairs of [`fluent-reference-cases.md`](../../../docs/research/fluent-reference-cases.md),
read through `KUTTA_FLUENT_EXAMPLES` and never committed.

## Result

| | main `c24da6d` | this branch |
| --- | ---: | ---: |
| Pairs admitted strictly | 0 / 43 | 0 / 43 |
| Pairs admitted with named approximations | 1 / 43 (V6, `wfgc`) | **15 / 43** |
| New pairs compared with Fluent's saved fields, CPU and Metal | — | 14 (5 distinct runs) |
| Of those, against a converged Fluent state | — | 2 (`elbow_param`, nozzle) |

Every new pair needs a named approximation: `wfgc` (SST elbows and nozzle) or the
two new ones below (k-epsilon elbows). None is a silent substitution; each is
opt-in on the command line, recorded in `summary.json`, and listed with its
measured effect by `inspect --json`.

## Cases

Cold starts from the case's own initialization, `kutta solve` to a float plateau,
then [`bench_fluent_reference`](../../../bench/bench_fluent_reference.cpp): unweighted
relative L2 per field, area-weighted inlet-minus-outlet pressure drop, and total
wall force against Fluent's saved `SV_P` and `SV_WALL_SHEAR` (V6 conventions).
"Saved" is Fluent's stored iteration and scaled continuity (raw / scale of the
last stored residual row).

| Case (pairs) | Model, numerics | Before | After | Fluent saved state | Velocity L2 | Pressure L2 | k L2 | omega / epsilon L2 | Pressure drop | Wall force | Backend |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `mixing_elbow/elbow_param` (1) | standard k-epsilon, std wall functions; SIMPLE, standard pressure, FOU, no pseudo time | rejected (20 issues) | runs, `compressibility,ke-options` | it 89, continuity 3.9e-4 | 0.93% | 6.37% | 9.08% | 17.59% | 0.91% | 0.34% | CPU and Metal; fields agree to 4 significant figures |
| nozzle `flprj_3891 … Base DP` (1) | SST; Coupled with pseudo time, distance flux | rejected (11) | runs, `wfgc` | it 10, continuity 8.4e-2 | 1.62% | 0.58% | 5.05% | 5.55% | 0.52% | 0.62% | CPU and Metal |
| `fluent-mixing_elbow_steady-state/elbow-2` (10 snapshots) | as `elbow_param` | rejected (32) | runs, `compressibility,ke-options` | it 5–50 of one run, continuity 2.4e-3 at 50 | 4.55% (it 50) | 6.90% | 15.2% | 16.2% | 2.83% | 0.69% | CPU and Metal |
| PPTX elbow (1) | SST; Coupled without pseudo time, distance flux | rejected (15) | runs, `wfgc` | it 20, continuity 0.20 | 8.89% | 28.0% | 27.0% | 18.3% | 50.9% | 49.0% elbow wall, 20.2% inlet wall | CPU and Metal |
| Dash elbow (1) | as PPTX, but air | rejected (15) | runs, `wfgc` | data are not this case's: see below | no valid comparison | | | | | | CPU and Metal |
| V6 `mixing_elbow` (1) | SST | runs, `wfgc` | unchanged | — | 60-iteration residuals and fields bit-identical to `c24da6d` (timing column aside) | | | | | | CPU |

Reading the table:

- **`elbow_param` and the nozzle are the real agreement evidence.** Fluent stopped
  `elbow_param` at 89 iterations with scaled continuity 3.9e-4. The nozzle is a
  10-iteration Coupled restart. Kutta's cold solutions converge (continuity 1e-5,
  momentum and turbulence 2e-8 to 1e-7) and differ from Fluent by about 1% in
  velocity and pressure drop and below 0.7% in wall force. These two pairs are
  now slow-tier gates (`tests/fv/test_fluent_coverage.cpp`) at about 1.5x the
  measured values.
- **elbow-2 is one unconverged Fluent run saved every 5 iterations.** Same mesh
  and setup as `elbow_param`. Fluent's snapshots approach Kutta's converged field
  monotonically, so the remaining difference is mostly Fluent's own iteration
  error:

  | Fluent iteration | 5 | 10 | 20 | 30 | 40 | 45 | 50 |
  | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
  | velocity L2 | 14.1% | 11.9% | 9.47% | 7.69% | 6.10% | 5.32% | 4.55% |
  | k L2 | 91.9% | 84.2% | 62.4% | 41.0% | 25.3% | 19.6% | 15.2% |
  | wall force | 33.2% | 5.22% | 1.26% | 1.02% | 1.04% | 0.84% | 0.69% |

  Not gated: the reference is unconverged.
- **PPTX: Fluent saved 20 iterations at scaled continuity 0.20.** A 51%
  pressure-drop difference against that state says nothing about agreement.
  Reported, not gated. The DOE, ACE and Parametric elbows with the same mesh
  have no data file.
- **Dash: the saved data belong to a different case.** The `.dat.h5` is
  byte-identical to PPTX's (SHA-256 `5706cd2c…`), a water run, while the Dash case
  sets the fluid to air. The bench now prints `REFERENCE_MISMATCH` (saved density
  814x the case material). Kutta's air solution (pressure drop 0.117 Pa) runs and
  converges but has no Fluent reference.
- Before/after counts are `kutta solve --iterations 1` setting issues:
  strict on `c24da6d` ([`preflight-c24da6d.txt`](preflight-c24da6d.txt)), and
  with `--approximate wfgc,compressibility,ke-options` on this branch
  ([`preflight-branch.txt`](preflight-branch.txt)). On `c24da6d`, only V6 was
  blocked by `wfgc` alone.

## Numerics added, and their sources

All are selected from the saved case. Nothing is a default chosen by Kutta.

| Fluent setting | Kutta | Source |
| --- | --- | --- |
| `pressure/scheme` 10, standard | `p_f = (p0/a0 + p1/a1)/(1/a0 + 1/a1)` with the unrelaxed momentum diagonal; zero normal gradient at walls | Theory Guide 23-42; "Ansys Fluent assumes that the normal pressure gradient at the wall is zero" |
| `flow/scheme` 20, SIMPLE | correction coefficient `V/a_P` | TG "SIMPLE vs. SIMPLEC" |
| `pbs/flux-index` 0, auto-select off | distance-weighted face velocity | TG 23-48 |
| flux keys absent (Fluent 21.2 cases) | momentum-weighted, with a setup note | Fluent 12.0 TG 18.4.1, the only documented scheme for that era |
| `pseudo-transient-formulation? #f` | classic steady iteration: implicit momentum relaxation, explicit pressure relaxation, no pseudo time | TG pressure-based solver |
| Coupled, no pseudo time | SIMPLEC with Fluent's explicit relaxation factors; reported as a convergence-path change | — |
| `ke-std-on?` standard k-epsilon | TG 4.3.1 equations; every coefficient (`kecmu`, `kec1`, `kec2`, `kesigk`, `kesige`) checked against Cmu 0.09, C1 1.44, C2 1.92, sigma_k 1.0, sigma_e 1.3 | TG 4.3.1 |
| `ke-std-wall?` | standard wall functions, log law above y* = 11.225, laminar below, equilibrium G_k and epsilon | TG "Standard Wall Functions" |
| laminar | no turbulence closure; viscous wall shear exported | — |

Classic SIMPLE, SIMPLEC and pseudo time only change the iteration path:
Rhie-Chow uses the unrelaxed (Majumdar) diagonal, so they share one fixed
point. `tests/fv/test_fluent_schemes.cpp` checks this on a tet channel on
every backend (fields within 1e-5).

### Standard pressure

<a id="standard-pressure"></a>The Theory Guide states the zero normal pressure
gradient for walls only. Kutta applies it at walls and symmetry (exact there)
and, as an **inference from measurement**, at velocity inlets:

| `elbow_param`, CPU, 300 iterations | Velocity L2 vs Fluent | Pressure drop |
| --- | ---: | ---: |
| inlet face pressure = cell pressure (shipped) | 0.93% | 0.91% |
| inlet face pressure extrapolated with the gradient | diverges (continuity 5e3, pressure drop 5e8 Pa) | — |

The shipped choice is the one that agrees with Fluent's converged fields.
It is not documented by Fluent. If a future case disagrees at an inlet, this
is the first place to look.

### Wall cells with one-plane neighbours (why the k-epsilon elbows first diverged)

On the 45,637-tet elbow mesh, some wall-adjacent tets have three interior
neighbours whose centroids lie almost in a plane parallel to the wall. The
least-squares pressure fit, which omitted wall faces, then has eigenvalues
2.6e-10, 4.6e-5 and 5.7e-5 (cell 36390; condition about 2e5). Its wall-normal
gradient amplifies rounding noise. Velocity in those cells reached
227 m/s after one iteration and 3,600 m/s after two, with 10 and 25 m/s
inlets. That happened on every path tried, including Kutta's default pseudo-time
SIMPLEC with second-order pressure.

Two faithful conditions close the fit:

- **Standard pressure:** zero normal derivative rows at unprescribed boundary
  faces (above).
- **Wall-function k:** `dk/dn = 0` at walls (TG, wall treatment for
  epsilon-based models). Without it, the same cell drove k to 2e7 m^2/s^2 by
  iteration 12.

With both, every k-epsilon pair converges to continuity 1e-5. The default
second-order-pressure path still omits wall pressure rows, so an SST case on
such a mesh would still be exposed. See open limits.

## Named approximations

<a id="compressibility-bound"></a>**`compressibility`.**
`turb-compress-mod? #t` adds Sarkar's dilatation dissipation
`Y_M = 2 rho epsilon M_t^2`, `M_t^2 = k/a^2` (TG 4.3.8). The elbows are
constant-density air, and how Fluent obtains `a` for a constant-density
material is undocumented, so the term is omitted. Relative to dissipation it is
`2k/a^2`. From Fluent's saved fields, `k_max` is 9.99 m^2/s^2 for `elbow_param`
and 9.62 for elbow-2 at iteration 50. The coldest saved temperature is 299.4 K,
which gives an ideal-gas air sound speed `a = sqrt(1.4 * 287 * 299.4) = 346.9 m/s`.
The omitted term is therefore at most `2 * 9.99 / 346.9^2 = 1.7e-4` of epsilon,
everywhere in these cases. Other cases need their own bound.

<a id="standard-k-epsilon-elbows"></a>**`ke-options`.**
`ke-source-mod` and `ke-time-marching?` are `#t` in these cases and appear in
no Fluent guide. Kutta solves the documented equations. Their effect cannot be
isolated. The combined discrepancy (with `compressibility`) on `elbow_param` is
in the case table: velocity L2 0.93%, pressure L2 6.37%, k L2 9.08%, epsilon
L2 17.59%, pressure drop 0.91%, wall force 0.34%. This is a combined
solver discrepancy, not an isolated error or a general bound.

## Tests

| Gate | Result |
| --- | --- |
| Quick tier `kutta-tests --test-suite-exclude=slow`, corpus present | 236 / 236 test cases, 7,414,417 assertions pass |
| `ctest -LE slow` (adds `lint-portability`) | 2 / 2 pass |
| New slow test, `elbow_param` and nozzle, CPU and Metal | 17,044 / 17,044 assertions |
| Slow tier against my own `c24da6d` baseline, `compare_slow_assertions.py` | **no new failures**: 0 regressions, 0 missing assertions; the same 68 failed assertions in the same 6 test cases on both (V6 and others that fail on main); candidate 28,098 passed (17,044 new) against baseline 11,054 |

Eight slow values changed, all in V7 SA TMR bump level 2 on CPU. Both runs stopped on that test's 600-second wall-clock budget, unconverged (2,500 steps on the baseline, 2,700 on the candidate, with other lanes sharing the machine), so the step count and its cf/cd errors track machine load. SA code paths are unchanged, and pass/fail status is identical. Files: [`slow-c24da6d.xml.gz`](slow-c24da6d.xml.gz), [`slow-cd68d0c.xml.gz`](slow-cd68d0c.xml.gz), [`slow-compare.json.gz`](slow-compare.json.gz).

One existing quick test changed meaning. `Unsupported k epsilon retains its
active equation schemes` asserted that standard k-epsilon was rejected. It now
asserts that standard k-epsilon is admitted, and that RNG k-epsilon is still
rejected. The scheme-retention checks are unchanged. No threshold, budget or
gate was loosened.

## Open limits, by impact

1. **The default (second-order pressure) path has the same singular wall-cell
   fit.** Pressure-LSQ wall rows are only added for the standard scheme. Kutta's
   default path diverged on the elbow tet mesh before this change, and still
   would on such meshes. Fluent's second-order wall pressure is extrapolated,
   so a zero-derivative row is not faithful there. A fit that stays bounded
   without inventing data (for example, rejecting near-null directions relative
   to the largest eigenvalue rather than at an absolute 1e-10) needs its own
   validation against V6 and the car cases.
2. **2D laminar cylinder (Re 40, converged).** This is the best remaining
   external-flow reference. It still needs 2D admission (`dimension`, front/back
   patch identity), Green-Gauss cell gradients, `prevent-reverse-flow?`,
   `p-profile-multiplier`, and a legacy `.dat` reader.
3. **Legacy binary meshes (section 2040) plus legacy `.dat`:** Borealis (3
   laminar tet pairs, 0.5M cells) and the VKI turbine.
4. **Converged references are scarce.** 12 of the 15 admitted pairs are saved
   short of convergence or mismatched. Only `elbow_param` and the nozzle can
   gate. Agreement on SST polyhedral elbows still rests on V6.
5. `file_session/elbow1` (hex, k-epsilon): `ke-realizability-on? #t` has no
   documented meaning for the standard model, and old cases store
   vector-valued neutral zone keys that the allowlist reads as active.
6. Out of scope here: transient (elbow-2 transient, 5 pairs), MRF/expressions
   (axial fan), DPM and gravity (oil separator), conjugate heat transfer
   (exhaust), multiphase/species readers, nonconformal interfaces, and ablation
   face winding.

## Reproduce

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build
export KUTTA_FLUENT_EXAMPLES=$PWD/validation/fluent_examples
E=$KUTTA_FLUENT_EXAMPLES/mixing_elbow
build/kutta solve $E/elbow_param.cas.h5 --backend metal --iterations 400 \
  --approximate compressibility,ke-options --out build/elbow_param
build/bench_fluent_reference $E/elbow_param.cas.h5 $E/elbow_param.dat.h5 \
  build/elbow_param/result.vtu metal compressibility,ke-options
build/kutta-tests -tc='Fluent*'                      # quick admission test
build/kutta-tests -tc='Fluent standard k-epsilon*'   # slow comparison
```

Per-run logs are in [`runs/`](runs/): `solve.log`, `residuals.csv`,
`forces.csv`, `summary.json`, the bench output `compare.log`, and the elbow-2
snapshot sweep `snapshots.log`. Cold solves ran for 800 iterations (SST) or
400 (k-epsilon). Wall time with other lanes sharing the machine, CPU / Metal: PPTX 203 / 61 s,
nozzle 135 / 50 s, `elbow_param` 147 / 35 s. The slow-tier XML files are in
this directory, gzipped.
