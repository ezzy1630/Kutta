# Continuity residual floor: diagnosis

> **Superseded for SIMPLE and SIMPLEC.** Velocity, pressure and face mass flux are now stored as float pairs, and the cylinder reaches its 1e-5 criterion at iteration 1700. Its plateau is about 200× lower: 1.5e-6 to 1.9e-6 scaled. See [the float-pair state report](../paired_state/REPORT.md). The diagnosis below still describes float32 storage, which Coupled uses.

Branch based on main `f545f26`. Apple M4 Pro, CPU reference backend and Metal. CUDA was not run; it compiles in CI only. The Fluent corpus is read from `KUTTA_FLUENT_EXAMPLES` and is not committed.

## Result

The cylinder's continuity plateau is the rounding floor of Kutta's float32 velocity storage, not a measurement error and not a stall. Fluent's 1e-5 criterion for this case cannot be reached in float32 under Fluent's own documented scaling. That is not an implementation gap. So this branch changes no solver or residual code. It adds a diagnostic tool and this report, and it corrects the cylinder report's wording.

The premise that "Fluent reaches 1e-5 in single precision" does not hold. **All four Fluent reference runs were double precision** (see "How Fluent ran").

Four measurements establish this:

1. **Kutta's residual definition already matches Fluent's.** Fluent's documented definition is the raw continuity residual over the largest raw continuity residual of the first five iterations. The case stores `residuals/scale? #t`, `normalize? #f` and `n-maximize-norms 5`. Kutta divides its raw residual, Σ over cells of |Σ_f F_f| in kg/s, by the largest of its first five values. That is the same quantity.
2. **The floor is not measurement error.** At the cylinder's plateau:
   - the corrected imbalance summed in float on the device is 4.91e-3 kg/s;
   - recomputed in double from the stored face fluxes, it is 5.37e-3 kg/s (CPU, iteration 3000).
   
   Double summation does not lower it, so the summation is not the problem.
3. **The solution has stopped changing.** After iteration 1500 the largest velocity change over a 20-iteration monitor interval is at most 2.4e-5 m/s on a 47.6 m/s maximum: a few ulps of 3.8e-6 each. Raw continuity has sat at 0.0100–0.0111 kg/s since iteration 1500.
4. **The plateau equals the effect of one ulp of velocity rounding.** At the converged state, the residuals were re-evaluated with every velocity, pressure, or both moved by one float32 ulp in a random direction:

   | Re-evaluated state (cylinder, CPU = Metal) | Continuity raw, kg/s | x-momentum raw |
   | --- | ---: | ---: |
   | converged state | 9.50e-3 | 0.782 |
   | velocities ±1 ulp | 1.376e-2 | 1.576 |
   | pressure ±1 ulp | 9.70e-3 | 0.787 |
   | both ±1 ulp (3 seeds) | 1.391e-2 – 1.410e-2 | 1.586 – 1.592 |

   - One ulp of velocity adds √(1.376² − 0.950²)·1e-2 = 9.95e-3 kg/s in quadrature, the size of the whole plateau.
   - Pressure adds almost nothing. The pressure offset is also small (0 to 990 Pa), so gauge-pressure storage would not help.
   - In units of ε·Σ|F| (ε = 2⁻²³, Σ|F| = 3.797e5 kg/s summed over cell faces), the plateau is **0.224**.

What the criterion would need, against what float32 achieves:

| Cylinder raw continuity | kg/s | in ε·Σ\|F\| |
| --- | ---: | ---: |
| Kutta plateau (iterations 1500–3000) | 0.0100 – 0.0111 | 0.22 – 0.25 |
| Velocity-rounding contribution of one ulp | 0.00995 | 0.22 |
| Needed for 1e-5 at Fluent's scale (448.29 kg/s) | 0.00448 | 0.099 |
| Needed for 1e-5 at Kutta's scale (26.01 kg/s) | 0.000260 | 0.0057 |
| Fluent's raw value at its stop, iteration 66 (double) | 0.00418 | 0.092 |

- Even with Fluent's denominator, the criterion sits 2.3× below Kutta's float plateau.
- With Kutta's own denominator, it sits 39× below. A float32 state can't reach 0.006 ε·Σ|F|: rounding to nearest alone leaves errors of up to ½ ulp in every velocity.

## How Fluent ran

All four reference runs were double precision; two independent checks show it:

1. **The section numbers.**
   - The cylinder `.dat` stores its fields in binary sections `3300`.
   - Fluent's file format prefixes a section index with 20 for single-precision binary and 30 for double-precision binary (Fluent User's Guide, Appendix B.2.2, Legacy Case and Data File Formats). `3300` is therefore a double-precision data field.
2. **The values themselves.** A single-precision run writes float values. In the cylinder data:
   - none of the 25,740 cell values of `SV_U`, `SV_V` or `SV_P` is exactly representable in float32;
   - the constants (`SV_DENSITY`, `SV_MU_LAM`) are.

   The same check over `SV_U/V/W/P` gives 0 float32-exact values in V6 (17,822 cells), `elbow_param` (45,637) and the nozzle (19,429).

Fluent's own guide describes the limit that applies here: "For single-precision computations … residuals can drop as many as six orders of magnitude before hitting round-off. Double-precision residuals can drop up to twelve orders of magnitude" (Fluent 6.3 User's Guide §26.13.1, Monitoring Residuals).

## Why Kutta's denominator is smaller

Fluent's continuity scale is "the largest absolute value of the continuity residual in the first five iterations" (Fluent 6.3 User's Guide §26.13.1, Eq. 26.13-6), and the cylinder case stores exactly that setting. The value therefore depends on how unbalanced the first five iterations are. Kutta's are better balanced than Fluent's in every case:

| Case | Fluent scale, kg/s | Kutta scale, kg/s | Fluent / Kutta |
| --- | ---: | ---: | ---: |
| Cylinder (Fluent Coupled; Kutta SIMPLEC) | 448.29 (iteration 1) | 26.011 (iteration 2) | 17.2 |
| V6 (both Coupled) | 2.614e-3 (iteration 2) | 9.768e-4 (iteration 1) | 2.7 |
| `elbow_param` (both SIMPLE) | 3.503 (iteration 2) | 0.4501 (iteration 1) | 7.8 |
| Nozzle (both Coupled) | 8.833e-3 (iteration 1) | 5.066e-3 (iteration 1) | 1.7 |

On the cylinder:

- Kutta's hybrid initialization leaves a potential flux whose imbalance is 0.021 kg/s.
- The Rhie-Chow imbalance of that initialized state is 13.0 kg/s.
- Kutta's first SIMPLEC iterations peak at 26.0 kg/s.
- Fluent's Coupled first iteration reports 448 kg/s. Kutta reproduces Fluent's first-iteration momentum scale (1.335e7 against 1.339e7), but not this imbalance.

Fluent's guide names this effect: "The scaling factor could be low if the initial continuity residual is small and thus the scaled residual could fail to meet the absolute convergence criterion" (§26.13.1, Choosing a Convergence Criterion). The case selects the absolute criterion (`residuals/convergence-criterion-type 0`).

Kutta could enlarge its denominator only by making its early iterations worse or by substituting Fluent's number. Both would game the criterion rather than match its definition, so neither was done.

## Every case: histories and stopping

Convergence checking was on (`monitors/check-convergence`), with each case's own criteria and the default monitor interval of 20. Iterations to stop:

| Case | Criteria | Fluent saved at | Kutta stop, CPU | Kutta stop, Metal | Kutta plateau, scaled continuity | Plateau raw in ε·Σ\|F\| |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| Cylinder | 1e-5, continuity and u, v | 66 (9.32e-6) | never (3000) | never (3000) | 3.9e-4 – 4.2e-4 | 0.22 |
| V6 mixing elbow | 1e-3, all | 294 (9.47e-4) | 200 | 200 | 5.0e-5 | 0.32 |
| `elbow_param` | 1e-3, all | 89 (3.87e-4) | 120 | 120 | 1.0e-5 | 0.28 |
| Nozzle | 1e-3, all | 10 (8.4e-2, saved unconverged) | 80 | 60 | 2.6e-4 | 2.70 |

- **3D cases.** V6, `elbow_param` and the nozzle stop on their Fluent criteria on both backends. Their continuity floors sit 4–100× below the criteria, and they reach them within 200–400 iterations. They would meet a floor only under a criterion tighter than their plateaus: 5e-5 (V6), 1e-5 (`elbow_param`) and 2.6e-4 (nozzle).
- **Nozzle backends.** It stops at 80 on CPU and 60 on Metal because its iteration-60 continuity is 1.0012e-3 on CPU and 9.499e-4 on Metal. That is a borderline crossing; the Coupled path is not bitwise identical across backends.
- **Momentum.** The cylinder's x-velocity residual, as `kutta solve` reports it, falls below 1e-5 at iteration 360 and plateaus at 2.5e-8 (raw 0.64). Continuity alone keeps the cylinder from stopping.

Kutta's continuity history on Metal (CPU agrees; the cylinder is bitwise identical across backends). Values are scaled continuity with raw kg/s in parentheses:

| Iteration | Cylinder | V6 | `elbow_param` | Nozzle |
| ---: | --- | --- | --- | --- |
| 1 | 1 (23.7)\* | 1 (9.77e-4) | 1 (0.450) | 1 (5.07e-3) |
| 100 | 6.1e-2 (1.59) | 1.6e-2 (1.6e-5) | 1.05e-3 (4.7e-4) | 2.8e-4 (1.40e-6) |
| 500 | 2.9e-3 (7.6e-2) | 5.1e-5 (4.9e-8) | 9.9e-6 (4.5e-6) | 2.5e-4 (1.29e-6) |
| 1000 | 4.6e-4 (1.19e-2) | 5.0e-5 | 1.02e-5 | 2.6e-4 |
| 1500 | 3.9e-4 (1.01e-2) | 5.0e-5 (4.9e-8) | 1.02e-5 (4.6e-6) | 2.6e-4 (1.30e-6) |
| 2000 | 4.1e-4 (1.06e-2) | 5.0e-5 | 1.0e-5 (4.5e-6) | 2.6e-4 (1.31e-6) |
| 3000 | 3.9e-4 (1.01e-2) | n/a | n/a | n/a |

\* At iteration 1 the scale is that iteration's own 23.7. From iteration 2 it is 26.0.

Fluent's double-precision histories for comparison: the cylinder falls from 448 to 4.18e-3 kg/s in 66 iterations, still dropping about 9% per iteration. V6 falls from 2.6e-3 to 2.48e-6 in 294 iterations.

Ulp probe on the 3D cases:

- **Nozzle.** It re-evaluates its iterated value exactly (1.31e-6). Its floor is pressure-sensitive: one ulp of pressure triples it (4.2e-6), velocity does nothing. Its pressure has no large offset (0.6 to 25 Pa), so this is the Rhie-Chow pressure term's own sensitivity, not gauge storage.
- **V6 and `elbow_param`.** Re-seeding through `initialize_fields` rebuilds the flux with the physical momentum coefficient, not the iteration's own, and gives 9.0e-7 and 7.7e-6 kg/s against iterated 4.9e-8 and 4.5e-6. The probe therefore doesn't isolate their floors. Their plateaus are 0.32 and 0.28 ε·Σ|F|, the same size as the cylinder's.

## What was considered and not done

- **Compensated or pairwise summation, or double recomputation of the per-cell sum.** The double recomputation is not lower than the float sum (5.37e-3 against 4.91e-3), so there is nothing to gain.
- **Gauge or offset pressure storage.** One ulp of pressure moves the cylinder's continuity by 2% and does not move V6's. No case has a large pressure offset.
- **Area-vector closure.** The stored float area vectors do not close exactly: Σ over cells of |Σ_f S_f| is 2.36e-4 m² on the cylinder. A converged solve absorbs that into the pressure field; it sets no floor in exact arithmetic, and it is of the same order as the other rounding terms.
- **Changing the normalization.** Kutta already uses Fluent's documented definition, and the case's settings select it. The relative criterion Fluent offers for this exact problem is not what the case selects.
- **Mixed precision (the only real fix).** Reaching 1e-5 at Kutta's own scale needs the cylinder's raw continuity about 39× below the float32 floor, roughly 5–6 more mantissa bits. Two routes:
  - float-float (hi + lo) storage for velocity and face flux, with compensated Rhie-Chow flux and continuity kernels on CPU, CUDA and Metal;
  - a double state on CPU and CUDA alone, which would leave Metal behind.

  Either changes solver arithmetic and memory traffic everywhere velocity is written: the momentum solves, the SIMPLEC velocity correction and the Coupled update. The Coupled update belongs to another lane. This is an architecture decision for the owner, not a residual fix, and it was not started.

## Parity, accuracy and tests

- **No solver or residual code changed.** `git diff main -- src tests` is empty. Every field, residual and force is therefore that of main `f545f26` by construction, so accuracy against Fluent is unchanged. The cylinder's numbers are in [its report](../cylinder_2d/REPORT.md): velocity relative L2 0.0645%, wall force 0.0705%, Cd 0.069% on CPU and Metal. The runs above used a main build.
- **Slow tier.** Not rerun: with no source or test change, it would compare main with itself.
- **Quick tier:** 246/246 test cases and 7,415,453 assertions pass, with fixtures present; 26 skipped are the slow suite.
- **Portability lint:** `tools/lint_portability.py` passes.

## Reproduce

`bench/bench_residual_floor` (new) prints Fluent's saved residual rows and the float32 check of its fields. It then prints Kutta's raw continuity, scale, float and double corrected imbalance, Σ|F| and field change along a cold-start run. Finally it re-evaluates the final state with one-ulp perturbations and the area closure. It only reads the solver.

```sh
./build/bench_restore_stripped_cr "$KUTTA_FLUENT_EXAMPLES/adjoint_cylinder_tutorial/cylinder_tutorial.dat" build/cylinder_restored.dat
./build/bench_residual_floor "$KUTTA_FLUENT_EXAMPLES/adjoint_cylinder_tutorial/cylinder_tutorial.cas" build/cylinder_restored.dat metal reverse-flow 3000
./build/bench_residual_floor "$KUTTA_FLUENT_EXAMPLES/mixing_elbow/mixing_elbow.cas.h5" "$KUTTA_FLUENT_EXAMPLES/mixing_elbow/mixing_elbow.dat.h5" metal wfgc 2000
./build/bench_residual_floor "$KUTTA_FLUENT_EXAMPLES/mixing_elbow/elbow_param.cas.h5" "$KUTTA_FLUENT_EXAMPLES/mixing_elbow/elbow_param.dat.h5" metal compressibility,ke-options 2000
# nozzle: flprj_3891/nozzle_para_named_1.cffdb/nozzle_para_named_1-Solve/Base DP/nozzle_para_named_1.{cas,dat}.h5, wfgc, 2000
```

Convergence-checked runs used `kutta solve <case> --approximate <list> --backend <b> --iterations 3000 --setup check.json`, where `check.json` is `{"v":1,"overrides":{"monitors/check-convergence":true}}`.

Evidence in [`runs/`](runs):

- `floor/<case>-<backend>.log`: the tool's full output;
- `checked/<case>-<backend>.csv`: convergence-checked residual histories (the cylinder's CPU history to 3000 iterations is in `floor/cylinder-cpu.log`, identical to Metal);
- `long/<case>-metal.csv`: 2000-iteration histories without checking.

## Sources

- Fluent User's Guide, Appendix B.2.2, Legacy Case and Data File Formats (binary section prefixes 20/30, sections 300 and 302): https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_app_form_case_data.html
- Fluent 6.3 User's Guide §26.13.1, Monitoring Residuals (Eq. 26.13-6, round-off orders, low scaling factor): https://www.afs.enea.it/project/neptunius/docs/fluent/html/ug/node812.htm

## Open items

- **Owner decision:** whether Kutta should offer mixed precision (float-float velocity and flux) so that double-precision-tight criteria such as the cylinder's 1e-5 can be met. Until then, a float32 run of this case stops on its iteration limit or on a report condition (`monitors/report/<report>/active`), not on its residual criteria.
- **Unexplained:** why Fluent's Coupled first iteration reports 448 kg/s on the cylinder, against Kutta's 13 kg/s for the Rhie-Chow imbalance of the initialized state. This would need Fluent's hybrid-initialized pressure or its iteration-1 fields, which the corpus doesn't hold.
- **CUDA** was not run.
