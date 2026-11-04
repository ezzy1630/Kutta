# NASA TMR convergence investigation

The rebased source completes all 37 slow cases and matches all 201,633 assertion
identities against main `8866cc0`. TMR failures fall from 28 to
17; total slow failures fall from 62 to 51.
There are 11 recoveries, zero regressions and zero missing or new
assertions. All four previously failing TMR accuracy assertions pass.
The remaining TMR failures are convergence gates; this matrix does not fully pass.

Quick validation passes 307 cases and 7,990,154 assertions. All six requested
40-iteration elbow case/backend pairs match main `271b31f` byte for byte across
five outputs. The standard-k-epsilon JumpCorrected opt-in also matches main on
CPU and Metal in exported cell fields and wall diagnostics.

## Protocol and provenance

The original baseline is main `8866cc0e16df5b0642b2243611099ddccd7f3cd3`.
The branch is rebased onto `271b31f7a96bbbc39eda3dd59721c9667d51ff46`.
The final numerical source and validation revision is
`96d8b34615f88b28fb5db2cefe45baad13cb4446`. It adds the paired poor-interface integration
correction after the derivative correction in `355e80b`. The concluding commit updates evidence
only. Binary and XML hashes below bind the completed run to this source.

The machine is an Apple M4 Pro with 24 GB unified memory. CPU is the serial
reference backend and Metal uses the actual GPU. One heavy solve ran at a time.
Shared-machine wall times are not performance evidence.

Every NASA grid/backend keeps its 3,000-step outer budget. Uniform initialization,
meshes, freestream data, turbulence formulations, first-order turbulence
advection, second-order momentum/pressure reconstruction, boundaries,
reference values, and all convergence/accuracy gates are retained. No
reference field seeds a run. No assertion is removed, weakened or skipped.

The completed coordinator baseline on exact main supplies all 24 target
solves, each at 3,000 steps. This lane also ran a partial exact-main target
reproduction before that full baseline became available. None converged. The
full baseline has 28 failing target assertions: 24 convergence and four
aft-bump Cf accuracy assertions. Every failure is enumerated below.
All `failure.empty()`, finite-value and load-availability assertions pass;
there is no reproduced coarse-grid divergence on starting main `8866cc0`.

| Model/case | Backend | Level | Continuity | Largest transport defect | Load change | Classification of failed convergence assertion |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| SST plate | cpu | 0 | 0.0001487815 | 3.478709e-07 | 2.61297e-07 | Continuity floor; loads settled |
| SST plate | cpu | 1 | 0.0004876112 | 1.48377e-07 | 7.114497e-08 | Continuity floor; loads settled |
| SST plate | cpu | 2 | 0.007514194 | 8.368414e-08 | 7.325466e-07 | Continuity floor; loads settled |
| SST plate | metal | 0 | 0.0001835793 | 3.384738e-07 | 4.203027e-07 | Continuity floor; loads settled |
| SST plate | metal | 1 | 0.0004031023 | 1.404783e-07 | 0 | Continuity floor; loads settled |
| SST plate | metal | 2 | 0.005696822 | 8.21644e-08 | 1.662186e-06 | Continuity floor; loads settled |
| SST bump | cpu | 0 | 2.269245e-05 | 0.000208247 (progress maximum) | 0.0002074413 | Continuity and turbulence plateau |
| SST bump | cpu | 1 | 1.549432e-05 | 3.567575e-05 (progress maximum) | 0.0001271153 | Continuity and turbulence plateau |
| SST bump | cpu | 2 | 1.654531e-05 | 0.0001197847 (progress maximum) | 0.0123655 | Slow convergence: loads and transport unresolved |
| SST bump | metal | 0 | 3.513322e-05 | 0.0002082479 (progress maximum) | 9.510323e-05 | Continuity and turbulence plateau |
| SST bump | metal | 1 | 2.633376e-05 | 3.56816e-05 (progress maximum) | 0.0001055961 | Continuity and turbulence plateau |
| SST bump | metal | 2 | 1.397719e-05 | 0.0001197803 (progress maximum) | 0.01236376 | Slow convergence: loads and transport unresolved |
| SA plate | cpu | 0 | 0.0002559519 | 8.825898e-06 | 4.60888e-07 | Continuity floor; loads settled |
| SA plate | cpu | 1 | 0.0006237517 | 9.858658e-06 | 3.406633e-08 | Continuity floor; loads settled |
| SA plate | cpu | 2 | 0.006870226 | 6.943274e-06 | 4.694171e-06 | Continuity floor; loads settled |
| SA plate | metal | 0 | 0.0003023776 | 7.153247e-06 | 1.539584e-07 | Continuity floor; loads settled |
| SA plate | metal | 1 | 0.0004093038 | 9.564456e-06 | 1.105921e-07 | Continuity floor; loads settled |
| SA plate | metal | 2 | 0.005286517 | 6.710996e-06 | 5.631781e-06 | Continuity floor; loads settled |
| SA bump | cpu | 0 | 3.560498e-05 | 0.0007388555 | 0.0003563521 | Continuity and turbulence plateau |
| SA bump | cpu | 1 | 1.393743e-05 | 0.0002499332 | 0.0001048261 | Continuity and turbulence plateau |
| SA bump | cpu | 2 | 2.786607e-05 | 0.0002351986 | 0.01228497 | Slow convergence: loads and transport unresolved |
| SA bump | metal | 0 | 3.013505e-05 | 0.0007388047 | 0.0007059305 | Continuity and turbulence plateau |
| SA bump | metal | 1 | 2.966741e-05 | 0.0002499212 | 0.00013187 | Continuity and turbulence plateau |
| SA bump | metal | 2 | 1.204483e-05 | 0.0002351449 | 0.01228709 | Slow convergence: loads and transport unresolved |

For SST bump, the CPU maximum transport/progress defects at levels 0/1/2
are approximately 2.08247e-4 / 3.56757e-5 / 1.19785e-4. The fine-grid
load change is about 1.24%, exceeding 0.1%. These runs combine a residual
plateau with incomplete load convergence; the mechanism of each plateau
has not been isolated. A continuity floor is observed, but its entire cause
is not established.

| Additional failing accuracy assertion | Before relative Cf error | Frozen maximum | Classification |
| --- | ---: | ---: | --- |
| SST bump, CPU, level 2, x=0.8678025 | 0.145326 | 0.10 | Accuracy on an unconverged field |
| SST bump, Metal, level 2, x=0.8678025 | 0.145348 | 0.10 | Accuracy on an unconverged field |
| SA bump, CPU, level 2, x=0.8678025 | 0.149930 | 0.10 | Accuracy on an unconverged field |
| SA bump, Metal, level 2, x=0.8678025 | 0.149882 | 0.10 | Accuracy on an unconverged field |

All other target Cf, CD, Cp and three-grid-trend assertions pass on main.
Passing accuracy at an unconverged endpoint is not verification acceptance.

## Retained changes

`SolverSettings::compensated_momentum` preserves the momentum RHS sum/product
remainders, the under-relaxation subtraction remainder, and the difference
between the physical diagonal sum and its positive upward-rounded matrix
coefficient. The paired defect uses these low words. The stored positive
linear operator stays unchanged. The same option computes the existing
bounded Rhie-Chow pressure jump after cancelling the common cell pressure,
rather than subtracting two already-rounded reconstructed traces.

The option defaults off. `paired_state_slack_ulps` defaults to the previous
one feeding-scale epsilon. The SST and SA TMR setups opt into compensation,
16-epsilon high-word hysteresis, and momentum/pressure inner absolute
tolerances of `1e-20`. Relative tolerances and linear iteration limits stay
unchanged. The hysteresis retains low-word corrections while allowing the
high-word nonlinear coefficients to settle. It is a numerical setting, not
a change to the outer residual gate. Other case setups keep their prior
settings. The bump setup uses fixed flow pseudo-time .1 and momentum
relaxation .7, compared with .01 and .5 on main; pressure relaxation remains
.3. These change the iteration, not the steady spatial equation.

The precision path also evaluates the same low-word pressure, velocity,
projected/secondary diffusion, deferred advection, transpose-stress and wall
operators as the high word. Flux normalization retains products between
low flux and low velocity. Physical low-state reaction is evaluated directly,
instead of subtracting two large rounded pseudo-time products. Green–Gauss
force accumulation carries the high-word product/sum errors as well as low
face pressures. Prescribed low boundary values are zero; symmetry projects
the low velocity. Normal backflow selects its direction from the physical flux sign
and its incremental speed from the low flux. The compensated flux high word
is rebased when necessary to preserve that sign.

Two later checks corrected the nonlinear limiter treatment. An intermediate
full-low-word prototype bounded the high and low projected face fields
independently. With high reconstruction 1.5 inside [1,2], and a cloud low
perturbation 2^-20, it returned a low increment 2^-20 instead of the correct
4*2^-20. The complete reconstruction stays inside [1,2]; clipping its
increment loses 2.86102295e-6. `momentum_increment` now bounds the total pair
and subtracts the high reconstruction. Walk-cloud extrema retain point
identities, including prescribed boundary points, so an active bound carries
the low word of the actual extremum.

The cell limiter also used paired extrema with only the high gradient.
A local extremum carried entirely in the low gradient incorrectly received
limiter 1 instead of 0. The precision path now uses both gradients and the
actual low boundary observations. This includes projected symmetry data,
which need not equal the owner's low component. Both fixes leave the
legacy/default path unchanged. The new focused checks pass on CPU and Metal;
the walk predictor differs from reconstructing a representable total field
by at most 3.05532e-7 in the checked mesh.

The pressure least-squares reduction now carries both state words and the
sum/product errors through the same geometry rows. Independent float fits
of the two words lose the high fit's error. A cancellation check with a true
gradient 2^-40 fails that old arithmetic; the paired fit retains 2^-40 under
two equivalent state partitions. Canonical full face pressures retain curvature
2^-31 below a common pressure of one. Neumann boundary low traces also carry
the high trace's rounding remainder. This is an operator consistency fix;
fixed-velocity boundary fluxes ignore that pressure jump, so the boundary
trace correction did not change the measured TMR runs by itself.

Rhie–Chow accumulation now combines both face velocity words with carried
products and sums, and cancels full pressure pairs before rounding the small
pressure curvature. Previously, separately rounded high and low fluxes omitted
the high dot-product remainder. The native counterexample 1 + 2^-25 - 1
loses 2^-25 in the old float reduction and retains it in the corrected flux.
The coefficient and boundary flux policies stay the same.

Nonlinear turbulence rates, scalar mass fluxes and viscous forces in this
option read the physical velocity/flux pair independently of its held high
word. Velocity values round once for wall data. Velocity derivatives now
apply the existing least-squares and SST face-jump operators to both state
words before rounding the final derivative to float32. The constitutive
kernels and their outputs remain float32.

An intermediate version rounded the velocity pair before differentiating it.
That erased a physical derivative 2^-27 in the native affine-field check.
The corrected derivative retains 2^-27 under two equivalent state partitions
on CPU and Metal. The SST face-jump correction uses the same full traces;
its geometric rows and reconstruction order are unchanged. SIMPLEC's
compensated pressure update normalizes overlapping words before applying
local pressure hysteresis. An isolated pressure 2^-43 with a correction
2^-60 is retained instead of lost under a global feeding-scale slack.
Velocity and flux keep feeding-scale hysteresis. Coupled retains its prior
pressure-update policy; the final TMR driver uses SIMPLEC.

Scalar transport can opt into a compensated delta solve. TMR enables it,
omits scalar pseudo-time damping, lowers scalar absolute tolerance to
1e-20, and checks an inner physical L1 defect tolerance of 1e-6. This is ten
times below the unchanged outer 1e-5 gate. The source is re-evaluated at
every solve, so a previously settled component resumes when coupling changes
its equation. Fixed wall omega values are imposed before and after the
increment solve. Outer monitors use the carried physical equation.

Independent arithmetic checks in
[`test_paired_precision.cpp`](../../../tests/fv/test_paired_precision.cpp)
cover the following failures. The focused final run passes 23 cases / 470
assertions across CPU and Metal, including the independent SA correction and
the unchanged legacy SA equation oracle. These checks establish local
arithmetic behavior; the full matrix determines whether a flow passes.

| Reproduction | Before | After |
| --- | --- | --- |
| Steady source 0.125-q, pseudo-time diagonal eight orders larger, q=1 | Physical defect -0.875 erased | Error below 1e-8 at all three checked q values |
| Constant diffusion field, face coefficients 1, 1e-8, 1e-12 | False defect 4.36837e-7 | Defect below 1e-12; stored matrix identical |
| Pressure curvature with gauge offset 0 or 1 | True -2^-31 erased at offset 1 | Both jumps within 1e-14 of -2^-31 |
| Constant advected/diffused boundary field | Boundary factor/product carry lost | Defect below 1e-12 |
| Opposing integrated source and volume product | Remainder erased or counted twice in an intermediate candidate | Source and equal source/sink cancellation retained |
| Flux increment 2^-10, scalar increment 2^-20 | Missing physical product 9.31323e-10 | Expected total-flux advection within 1e-13 |
| Low-state reaction +/-2^-19 under strong pseudo-time | Physical reaction rounded to zero | Both defects within 1e-13 of the reaction |
| Uniform-pressure face sum with areas 1, 2^-25, -1, -2^-25 | High reduction leaves false -2^-25 | Carried force cancels it |
| Normal backflow, high flux -1, low flux 2^-20, area 2 | High speed -0.5 incorrectly enters low data | Incremental speed +2^-21 |
| Low projected perturbation inside the total bounds | Increment clipped from 4*2^-20 to 2^-20 | Full increment retained; active total bounds retain their own low word |
| Extremum carried only by the low gradient | Limiter 1 permits an undershoot | Limiter 0; actual low boundary observations retained |
| Scalar already below inner tolerance; then change its source | Repeated roundoff chasing | Initial hold, followed by a real solve when source changes |

## Integrated TMR measurements

These are the complete final integrated results, compared with the original
`8866cc0` baseline. Main's poor-cell correction is included in the final source,
so bump differences combine that change with this lane's retained changes.
No final field shift is attributed to an individual fix without an isolating
experiment.

| Model/case | Backend | Level | Converged before -> after | Steps after | Continuity before -> after | Maximum defect before -> after |
| --- | --- | ---: | --- | ---: | --- | --- |
| SST plate | cpu | 0 | False -> True | 300 | 0.0001487815 -> 1.086581e-06 | 0.0001487815 -> 1.086581e-06 |
| SST plate | cpu | 1 | False -> True | 800 | 0.0004876112 -> 4.278579e-06 | 0.0004876112 -> 4.278579e-06 |
| SST plate | cpu | 2 | False -> False | 3000 | 0.007514194 -> 1.90981e-05 | 0.007514194 -> 1.90981e-05 |
| SST plate | metal | 0 | False -> True | 300 | 0.0001835793 -> 1.86489e-06 | 0.0001835793 -> 1.86489e-06 |
| SST plate | metal | 1 | False -> False | 3000 | 0.0004031023 -> 3.045154e-05 | 0.0004031023 -> 3.045154e-05 |
| SST plate | metal | 2 | False -> True | 3000 | 0.005696822 -> 7.407654e-06 | 0.005696822 -> 7.407654e-06 |
| SST bump | cpu | 0 | False -> False | 3000 | 2.269245e-05 -> 3.143857e-05 | 0.000208247 -> 3.143857e-05 |
| SST bump | cpu | 1 | False -> False | 3000 | 1.549432e-05 -> 2.308443e-05 | 3.567575e-05 -> 2.308443e-05 |
| SST bump | cpu | 2 | False -> False | 3000 | 1.654531e-05 -> 1.210993e-05 | 0.0001197847 -> 2.573423e-05 |
| SST bump | metal | 0 | False -> False | 3000 | 3.513322e-05 -> 1.51518e-05 | 0.0002082479 -> 1.991881e-05 |
| SST bump | metal | 1 | False -> False | 3000 | 2.633376e-05 -> 2.81366e-05 | 3.56816e-05 -> 2.81366e-05 |
| SST bump | metal | 2 | False -> False | 3000 | 1.397719e-05 -> 1.08788e-05 | 0.0001197803 -> 2.572706e-05 |
| SA plate | cpu | 0 | False -> True | 400 | 0.0002559519 -> 4.356717e-08 | 0.0002559519 -> 9.231233e-07 |
| SA plate | cpu | 1 | False -> True | 900 | 0.0006237517 -> 4.993192e-06 | 0.0006237517 -> 4.993192e-06 |
| SA plate | cpu | 2 | False -> False | 3000 | 0.006870226 -> 4.282008e-05 | 0.006870226 -> 4.282008e-05 |
| SA plate | metal | 0 | False -> True | 300 | 0.0003023776 -> 2.887325e-06 | 0.0003023776 -> 2.887325e-06 |
| SA plate | metal | 1 | False -> False | 3000 | 0.0004093038 -> 2.417789e-05 | 0.0004093038 -> 2.417789e-05 |
| SA plate | metal | 2 | False -> False | 3000 | 0.005286517 -> 1.397759e-05 | 0.005286517 -> 1.397759e-05 |
| SA bump | cpu | 0 | False -> False | 3000 | 3.560498e-05 -> 3.439757e-05 | 0.0007388555 -> 3.439757e-05 |
| SA bump | cpu | 1 | False -> False | 3000 | 1.393743e-05 -> 2.368411e-05 | 0.0002499332 -> 4.124983e-05 |
| SA bump | cpu | 2 | False -> False | 3000 | 2.786607e-05 -> 7.859322e-06 | 0.0002351986 -> 6.128235e-05 |
| SA bump | metal | 0 | False -> False | 3000 | 3.013505e-05 -> 5.78428e-05 | 0.0007388047 -> 5.78428e-05 |
| SA bump | metal | 1 | False -> False | 3000 | 2.966741e-05 -> 2.251525e-05 | 0.0002499212 -> 4.124282e-05 |
| SA bump | metal | 2 | False -> False | 3000 | 1.204483e-05 -> 1.033504e-05 | 0.0002351449 -> 6.129666e-05 |

| SST plate raw continuity | Level | Before | After |
| --- | ---: | ---: | ---: |
| cpu | 0 | 2.726211e-08 | 1.989631e-10 |
| cpu | 1 | 1.33167e-08 | 1.16349e-10 |
| cpu | 2 | 3.258206e-08 | 8.092247e-11 |
| metal | 0 | 3.363881e-08 | 3.414859e-10 |
| metal | 1 | 1.100876e-08 | 8.280803e-10 |
| metal | 2 | 2.470181e-08 | 3.138772e-11 |

The full slow driver prints raw continuity for SST plates only. Unprinted raw
values are not inferred from separate representative probes. Complete public
station values and printed transport defects are in [`measurements.json`](measurements.json).

| Plate Cf | Backend | Before | After |
| --- | --- | ---: | ---: |
| SST | cpu | 0.002713813 | 0.002713727 |
| SST | metal | 0.002713817 | 0.002713699 |
| SA | cpu | 0.002726861 | 0.002727054 |
| SA | metal | 0.002726862 | 0.002727054 |

| Plate CD | Backend | Before | After |
| --- | --- | ---: | ---: |
| SST | cpu | 0.002876909 | 0.002876801 |
| SST | metal | 0.002876906 | 0.002876764 |
| SA | cpu | 0.002880514 | 0.002880676 |
| SA | metal | 0.00288051 | 0.002880676 |

| Bump aft Cf | Backend | Before | After |
| --- | --- | ---: | ---: |
| SST | cpu | 0.00298066 | 0.002553902 |
| SST | metal | 0.002980718 | 0.002553883 |
| SA | cpu | 0.003079214 | 0.002638858 |
| SA | metal | 0.003079086 | 0.002638848 |

| Finest bump aft Cf relative error | Backend | Before | After | Frozen limit |
| --- | --- | ---: | ---: | ---: |
| SST | cpu | 0.1453259 | 0.01865692 | 0.10 |
| SST | metal | 0.1453479 | 0.01866434 | 0.10 |
| SA | cpu | 0.1499296 | 0.01452083 | 0.10 |
| SA | metal | 0.149882 | 0.01452457 | 0.10 |

## Rebase and compatibility

The conflict in `update_turbulence()` was resolved around main's shared
`reconstruct_strain()`. For compensated flow, velocity derivatives and the
jump correction are formed from both state words before rounding. SST and
standard k-epsilon JumpCorrected use that same helper. Other standard k-epsilon
and SA paths consume `model_gradient()`; constitutive velocities, boundary
velocities and mass fluxes use their corresponding model accessors. When
compensated momentum is disabled, those accessors return main's original
buffers and reconstruction follows main's original operators.

The normal-foot pressure option and its gradient configuration guards from
main are retained; the TMR setup keeps the default extrapolated mode.
Main's harmonic-k transition setting and SST molecular omega-production
parameter are also preserved through the existing settings copies. TMR keeps
harmonic k disabled and the molecular ratio zero; no new model option is
enabled by this lane.
The standard k-epsilon opt-in branch, its variant/LSQ guards, the Péclet option,
poor-pressure-interface correction and main's performance caches are retained.
The duplicate comparator commit `6a74dc5` was dropped. The comparison uses
main's `tools/forensics/compare_slow_assertions.py` unchanged, including its
existing CLI-context normalization from `971fbb3`.

| Original SHA | On d8e26c7 | On 05ef218 | Final SHA on 271b31f | Change |
| --- | --- | --- | --- | --- |
| `deb6781` | `948993f` | `42b0137` | `7738484` | Preserve momentum assembly remainders in paired solves |
| `6c7c9e9` | `6eaa7f5` | `c26a3fd` | `a6bbb99` | Use compensated paired precision in NASA TMR cases |
| `4bfb039` | `75a597f` | `2d7d9ee` | `6344a69` | Carry boundary and source products consistently in paired assembly |
| `9479942` | `8faea41` | `aa9f521` | `2599c39` | Keep TMR SA production separate from the strain safeguard |
| `180c717` | `60afafb` | `921a0ce` | `3497c70` | Carry consistent low-word flow operators and scalar defects |
| `a789dc2` | `b2aef1d` | `d5809d7` | `60dd9f4` | Use physical scalar stopping in NASA TMR solves |
| `a5e8067` | `d594855` | `06899e8` | `dbf8dde` | Bound paired face predictors using the complete field |
| `3cac7f0` | `4a9067c` | `09a1464` | `577fdce` | Limit paired gradients with their low-word boundary observations |
| `b62295c` | `72984a2` | `fcb8f65` | `7f4a55e` | Carry complete paired pressure traces into momentum forces |
| `86a665a` | `c4d1461` | `a55bda8` | `1645a31` | Keep compensated flow operators consistent with the physical pair |
| `6a74dc5` | dropped | dropped | main supplies normalization | Normalize transient CLI context in slow assertion comparisons |
| `d823d6e` | `9bc07c2` | `6506725` | `355e80b` | Take constitutive velocity derivatives before rounding paired state |
| `8bbebfb` | `a50d195` | `681b728` | `8ea334f` | Record TMR convergence results and interrupted final validation |
| new integration fix | `4cf62d9` | `4a0050a` | `96d8b34` | Preserve poor-interface pressure correction in paired reconstruction |

The integration fix, `96d8b34` (previously `4cf62d9` and `4a0050a`), carries main's poor-interface rule into
the compensated pressure-face operator. The paired reconstruction previously
rebuilt face pressure after the default correction and retained the good
cell's extrapolation at a mixed interface. An isolated check produced 0.75
or 1.625 instead of the required 1.5 on CPU and Metal: eight assertions
failed. The corrected operator omits both extrapolations and preserves the
low-word cell mean; all eight assertions pass without changing an existing
expected value or threshold. The first rebased slow attempt was stopped
after this gap was found, and the final complete run was restarted on the
corrected source. The intermediate `4cf62d9` attempt was then stopped when
the requested main moved to `05ef218`. Both partial XML files are excluded
from acceptance evidence. The subsequent `4a0050a` run was also stopped
when the requested main moved to `271b31f`. A subsequent turn interruption
terminated the first `96d8b34` attempt; all 37 cases were restarted on its
unchanged source and binary. Only the complete final run is compared.

The requested CLI checks used the same 40-step budget, initialization,
approximations and allocator for both binaries. They compared `fluxes.csv`,
`forces.csv`, `loads.json` and `result.vtu` in full. For `residuals.csv`, only
`#` comment lines and the `seconds_per_iteration` column were excluded; the
last recorded iteration was checked to be 40.

| CLI case | Approximation | CPU | Metal |
| --- | --- | --- | --- |
| mixing_elbow | `wfgc` | all five outputs identical | all five outputs identical |
| elbow_param | `compressibility,ke-options` | all five outputs identical | all five outputs identical |
| mixing_elbow converted to SIMPLEC | `wfgc` | all five outputs identical | all five outputs identical |

The SIMPLEC conversion used `kutta convert --setup` with
`{"v":1,"overrides":{"methods/pressure-velocity-coupling":"simplec"}}`.
The separate opt-in check used `bench_elbow_turbulence`, cold initialization,
40 iterations and `jump` on elbow_param. `cells.csv` and `walls.csv` match main
on each backend, including strain, production and wall quantities.
All 30 final output hashes also match the archived 40-step `d8e26c7`
reference outputs, using the same residual normalization. This preserves the
original reference comparison alongside the fresh `271b31f` comparison.
Binary and output hashes are in [`elbow-parity.json`](elbow-parity.json).

## SA definition audit

The frozen [plate](https://tmbwg.github.io/turbmodels/flatplate_sa.html)
and [bump](https://tmbwg.github.io/turbmodels/bump_sa.html) datasets specify
standard SA with ft2, as recorded in the pinned
[fixture manifest](../../data/nasa_tmr/SOURCE.md). [TMR](https://tmbwg.github.io/turbmodels/spalart.html)
recommends SA-neg for robustness. Kutta clips nu-tilde at zero and lacks
SA-neg's negative branch. SA-neg normally uses Note 1(c); Kutta's reported
TMR safeguard is Note 1(b). The production safeguard error described below was corrected explicitly; ft2
and the reference variant remain unchanged.

| Term or condition | Current implementation and source |
| --- | --- |
| Eddy viscosity | rho nu-tilde fv1; chi=nu-tilde/nu, cv1=7.1 (`sa_coefficients`) |
| Production | cb1 (1-ft2) S-tilde nu-tilde; vorticity base with the fv2 correction |
| Destruction | [cw1 fw - cb1 ft2/kappa^2] (nu-tilde/d)^2; cw1 derived from cb1, cb2, sigma and kappa |
| Diffusion | div[(mu+rho nu-tilde)/sigma grad(nu-tilde)], plus rho cb2/sigma times the squared gradient |
| Damping | TMR ft2=1.2 exp(-0.5 chi^2); Fluent formulation omits ft2 |
| Safeguards | TMR Note 1(b), S-tilde at least 0.3 Omega, and r at most 10; the published production rate remains unclipped in the TMR setup |
| Wall and inflow | nu-tilde=0 at resolved walls; freestream nu-tilde=3 nu; geometric nearest-wall distance |
| Source linearization | Full local negative derivative is implicit with a compensating explicit term; net reaction source is unchanged |

An isolated positive-nu-tilde state with negative raw modified strain
exposed a model error: the Note 1(b) denominator safeguard was also applied
to production. The [TMR definition](https://tmbwg.github.io/turbmodels/spalart.html)
and [NASA CFL3D source at 776d5169](https://github.com/nasa/CFL3D/blob/776d5169f34d26ae65701c87038e63b125b7d86a/source/cfl3d/libs/spalart.F#L1380)
retain raw modified strain in production while safeguarding r. Kutta now
separates those quantities and their local derivatives. This is not evidence
that the downloaded historical CFL3D reference used the same optional r
safeguard; the current source's standard branch uses Note 1(a).

The lane first kept the clipped production behind a compatibility option.
At integration the option was removed: every `TmrStandard` SA solve now uses
the published form, and the legacy equation oracle in `tests/turbulence/test_sa.cpp`
was corrected to match. Fluent's existing form is unaffected. The independent test uses
nu-tilde=.003, nu=.001, wall distance=.1 and Omega=1. Before the correction,
the production-rate error was 1.9385 and the net reaction error 7.77496e-4.
Four assertions failed across CPU/Metal. The corrected form passes all eight
assertions, including nonnegative explicit production and implicit sink.
This proves the equation correction. The final field comparison combines
precision and iteration changes, so its load shifts are not attributed to
this production correction alone. The legacy equation oracle remains unchanged.

Fluent Theory Guide 2026 R1 sections 4.2.2 (4-15), 4.2.3 (4-16 to 4-18),
4.2.4 (4-19 to 4-25), 4.2.5 (4-26 to 4-29), 4.2.6 and 4.2.7
were checked against the implementation. The transport, production and destruction
SVG equations were downloaded and rendered because the text view omits their
mathematics. Fluent production and destruction omit ft2; the TMR branch
adds it explicitly. The automatic-wall blend described in section 4.2.7
is distinct from these resolved-wall TMR runs.

Resolved-wall, inlet and backflow data were checked against the setup. The
low-word boundary operator used zero prescribed low inlet/wall
values, homogeneous pressure-correction data and projected symmetry velocity.
It did not establish boundary treatment as the cause of the remaining plateau.

The existing independent double-equation test covers both legacy configured forms,
source splitting, boundary data and CPU/Metal parity. Passing that constitutive
test alone does not establish converged flow verification. Fluent automatic
walls remain unavailable; resolved molecular walls are used here.

Sources: [TMR SA equations and variants](https://tmbwg.github.io/turbmodels/spalart.html),
[TMR plate SA reference](https://tmbwg.github.io/turbmodels/flatplate_sa.html),
[TMR bump SA reference](https://tmbwg.github.io/turbmodels/bump_sa.html),
[Fluent Theory Guide 2026 R1, section 4.2.4, production](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_th/flu_th_sec_spal_turb_prod.html).
The Fluent formula contract and its other section citations are in
[`car_models.md`](../../../src/kutta/fv/turbulence/car_models.md).

## What did not work

These isolating experiments predate the rebase. Their numbers belong to the
configurations named here; the integrated matrix above is the final evidence.

Each exploratory solve used at most 3,000 outer steps. Intermediate variants
below are not the final solver. The retained operators and TMR controls are described above. Earlier
prototype numbers below are diagnostics, not final acceptance.

| Isolating experiment | Observed result |
| --- | --- |
| Lower momentum or pressure absolute tolerances alone | CPU coarse SST plate still failed; continuity 6.8803e-5 or 1.48553e-4. Lowering both gave 2.81011e-4. |
| Wider high-word hysteresis alone | 16 epsilons passed coarse CPU plate at 400 steps, continuity 3.49234e-6; middle/fine still failed at 5.74774e-5 / 2.5566e-3. 256 epsilons still failed the middle grid at 3.30734e-5. |
| Freeze turbulence after step 1000 | Middle CPU plate still failed, continuity 7.99189e-5. |
| Freeze all flow high words after step 1000 | Middle plate reached 2.92726e-6 at 1100. Freezing velocity, pressure or flux alone still failed. This was diagnostic instrumentation, removed. It is not an accepted solution. |
| Consistent low pressure, face-velocity, diffusion, advection, transpose-stress and wall terms | Fine CPU plate continuity improved from about 3.1e-3 with the earlier compensated-RHS prototype to 4.38298e-5; still above 1e-5. Adding the physical-diagonal remainder gave 3.01106e-5. That early prototype was discarded; the corrected low-word operators above were rebuilt and retained. |
| Tighter relative momentum, pressure and scalar solves with that full reconstruction | Fine CPU plate remained at 3.01642e-5. Coarse bump retained k near 2.09e-4. |
| Increase bump pseudo-time step from .01 to .05 alone | Coarse k improved to 6.38385e-5 but continuity was 8.688e-5. Fine Metal bump still had continuity 1.59196e-5, load change .0118142 and aft Cf .002941025. |
| Existing differentiable limiter | Coarse CPU plate continuity 1.65420e-4. With bump dt=.05 and momentum relaxation .7, continuity 2.32530e-5, k 1.45272e-5 and omega 1.31402e-5; still failed. Coupled bump also failed. No limiter change retained. |
| Coupled plus full paired reconstruction | Fine CPU plate passed at 600 with continuity 2.56369e-6, but coarse plate failed at 3.44269e-5. |
| Component steady stopping at 1e-7 with Coupled/full reconstruction | Coarse CPU plate passed at 300 with continuity 3.05581e-7, but middle plate remained at 1.33269e-5. Coarse bump developed a bounded limit cycle, continuity 2.4946 and load change 14.247 at 3000. Removed. |
| Tighter Coupled block solve | Middle plate continuity worsened to about 8.74e-5; tighter linear termination did not remove the nonlinear plateau. |
| Compensated scalar increments plus component stopping | Coarse SA plate passed in 200 (continuity 2.81086e-6); middle SST plate still failed at 1.08169e-5. Fine Metal SST plate passed in 600 (1.59899e-6), but this was not a passing common configuration. The final configuration uses a different inner stopping tolerance. |
| Disable scalar pseudo-time, retaining the full prototype | Coarse CPU bump at 3000 still failed: continuity 2.84891e-5, k 8.75590e-6, omega 2.08179e-5. This alone was insufficient. |

The final arithmetic path was also tested separately from those earlier
prototypes:

| Later controlled check, CPU | Result |
| --- | --- |
| Fine SST plate, full low operators, carried physical reaction, scalar inner stopping 1e-7 | Continuity 3.31229e-5 at 3000; failed |
| Same with compensated Green–Gauss force | 2.82884e-5, raw 1.19910e-10 at 3000; failed |
| Same binary, scalar inner stopping 1e-6 | Passed at 2900: continuity 4.66468e-6, raw 1.97728e-11, k 4.85604e-7, omega 9.17301e-7, zero measured load change |
| Wider flow hysteresis, 256 epsilons | Fine plate 4.30620e-5 at 3000; failed |
| Coarse SST bump, dt=.1, momentum relaxation .7, scalar stopping 1e-6 | Continuity 1.70541e-5, k 5.79433e-6, omega 1.21983e-5; failed |
| Same, both scalar relaxation factors .3 | Continuity 4.90779e-6, k 1.20001e-4, omega 7.61034e-6; failed |
| Same, only omega relaxation .3 | Continuity 4.08207e-5, k 2.34429e-5, omega 1.51188e-5; failed |
| Corrected SA coarse bump, dt=.1, scalar stopping 1e-7 | Continuity 3.18097e-5 and nu-tilde 2.67610e-5; failed |
| Corrected SA coarse bump, scalar relaxation .3 and stopping 1e-6 | Continuity 2.76500e-5 and nu-tilde 1.70707e-4; failed |

The scalar-stopping pair is an isolating experiment: the binary, uniform
start and all other inputs match. It shows that allowing an already small
scalar defect to settle can reduce the pressure–velocity iteration floor.
It does not show that all remaining plateaus have this cause. The final
full-suite run rechecks the chosen common controls on all grids/backends.
No source, load, iteration-budget or acceptance gate was frozen or bypassed.

After the limiter corrections, the following checks still did not produce a
passing common configuration:

| Check | Observed result |
| --- | --- |
| Coarse SST bump, retained scalar relative tolerance 1e-3 | Continuity 1.98756e-5, k 5.90183e-6, omega 2.27048e-5; scalar inner solves converged in 6 / 5 iterations |
| Same, scalar relative tolerance 1e-6 | Continuity 1.24490e-5, k 5.58984e-6, omega 1.16941e-5; inner solves converged in 12 / 9 iterations |
| Same tighter solves, scalar stopping 8e-6 | Continuity 1.32611e-5, k 7.90939e-6, omega 1.16934e-5; still failed |
| Fine Metal SA plate, scalar stopping 1e-6 / 5e-6 / 8e-6 | Continuity 5.72649e-5 / 1.52547e-5 / 1.84223e-5; all failed at 3000 |

A read-only CPU probe captured every scalar update in the last 100 of the
3,000 coarse-bump steps. With the tighter relative solve tolerance, the median
changed-cell count was 2,455 of 3,520 for k and 2,022 for omega. Among changed
cells, the median jump was 41.5 / 22 ulps; median per-step p95 jumps were
5,559.5 / 2,873 ulps. Median relative L1 field movement was 8.08e-6 / 1.84e-6.
This is continued nonlinear motion, not a demonstrated all-cell one-ulp
storage stall. The inner linear solves converge, and tightening them does
not remove the failure. These observations do not identify its sole cause.
No extra outer steps or extra monitor calls were used for this probe.

The final pressure/flux consistency work was measured separately on the
fine Metal SA plate and coarse CPU SST bump. These are representative
experiments, not substitutes for the final suite.

| Isolated later configuration | Fine SA continuity at 3000 | Coarse SST bump continuity / omega at 3000 |
| --- | ---: | ---: |
| Full paired pressure fit and canonical constitutive inputs, feeding-scale pressure hysteresis | 1.29405e-4 | 1.24805e-5 / 1.15592e-5 |
| Same, compensated full Rhie–Chow flux | 1.71804e-5 | 1.47958e-5 / 1.16993e-5 |
| Same, local compensated SIMPLEC pressure update (retained) | 1.64536e-5 | 1.22855e-5 / 2.17146e-5 |
| Same, classic steady flow without pseudo-time (rejected) | 4.07637e-5 | 2.66522e-5 / 1.17152e-5 |

The flux pair is an isolating check: it reduces the fine-plate raw
continuity defect from 6.67230e-10 to 8.85415e-11 with the same flow settings.
Local pressure updating then gives 8.47958e-11. Neither clears the gate.
An earlier local-pressure prototype, before completing flux accumulation,
had continuity about 9.32e-3. Applying local hysteresis to velocity and flux
as well was also rejected (about 8.93e-3). Gradient and boundary trace
corrections alone did not cure that behavior. These experiments establish
specific lost arithmetic terms and their measured consequences; they do not
establish that one term explains every remaining plateau.

Correcting the order of velocity rounding and differentiation was then
measured with the same uniform starts and settings:

| Representative solve | Before derivative correction | After |
| --- | --- | --- |
| Coarse CPU SST plate | Failed at 3000, continuity 2.79749e-5 | Passed at 300, continuity 1.08658e-6 |
| Fine Metal SA plate | Failed at 3000, continuity 1.64536e-5 | Failed at 3000, continuity 1.39776e-5 |
| Coarse CPU SST bump | Continuity 1.22855e-5, omega 2.17146e-5 | Continuity 3.19929e-5, omega 1.42903e-5, still failed |

The native affine-field check isolates the lost derivative. These flow
probes measure its effect without changing another setting. They recover
one convergence gate and leave the other two failing. The final full run
checks every grid/backend together.

In the classic-steady experiment, the raw fine-plate defect falls slightly
to 7.81911e-11, but its frozen initial continuity normalization also differs,
so normalized continuity is worse. The gate is always reported as defined;
a smaller raw residual is not substituted for acceptance.

Lower inner tolerances, larger time steps and the alternative limiter each
left failures. The paired operators recover specific arithmetic terms, but
the remaining nonlinear convergence failures still need a cause established
by a controlled experiment.

## Full validation

The complete final process ran `build/kutta-tests -ts=slow -r=xml --success`
in the foreground and was polled to completion. All 37 cases executed in one
process. The nonzero test exit reflects retained failing numerical gates;
the exact assertion comparison exits zero.

- 201,633 assertions in both baseline and final XML.
- 62 -> 51 failed assertions; 11 recoveries, zero regressions.
- Zero missing/new cases or assertions.
- 52 changed expanded values, all in TMR or already moved by main's poor-cell fix.
- All 34 off-target failures remain: V6 25, V2 tet eight, V12 one.

| Slow case | Assertions before / after | Failures before -> after | Changed expanded values |
| --- | ---: | ---: | ---: |
| Kutta cases converted from Fluent cases solve to the same bytes on every backend | 68 / 68 | 0 -> 0 | 0 |
| The NACA ground-effect wing from a bare mesh and a Kutta case solves like its donor Fluent case | 1 / 1 | 0 -> 0 | 0 |
| Native FSAE Fluent physics convert and solve to byte-identical Kutta outputs | 477 / 477 | 0 -> 0 | 0 |
| Corpus FSAE physics preserve strict solver holds or convert with identical outputs | 0 / 0 | 0 -> 0 | 0 |
| FV V3 Re1000 cavity refinement against Erturk | 16 / 16 | 0 -> 0 | 0 |
| FV V3 Re100 and Re400 cavity reported deviations | 4 / 4 | 0 -> 0 | 0 |
| FV V4 Fluent cavity matches its actual grid | 763 / 763 | 0 -> 0 | 0 |
| Fluent standard k-epsilon elbow and SST nozzle match converged Fluent fields | 17044 / 17044 | 0 -> 0 | 0 |
| Fluent k-epsilon elbow converges with second-order pressure | 11220 / 11220 | 0 -> 0 | 0 |
| Fluent laminar 2D cylinder matches converged Fluent fields | 207 / 207 | 0 -> 0 | 0 |
| FV V6 Fluent mixing elbow fixed point and cold start | 2565 / 2565 | 25 -> 25 | 0 |
| FSAE aero: Tyrrell 026 sectional Cl in ground effect vs thesis CFL3D k-omega | 28 / 28 | 0 -> 0 | 0 |
| FSAE aero: Ahmed body 35 deg drag, level 1 | 4 / 4 | 0 -> 0 | 1 |
| FSAE aero: Selig S1223 at Re 2e5, alpha 4.76 deg vs UIUC | 9 / 9 | 0 -> 0 | 0 |
| FSAE aero: Ahmed body 35 deg at y+ ~ 1, level 1, SIMPLEC and Coupled | 8 / 8 | 0 -> 0 | 4 |
| FSAE V12 native Fluent polynomial fan import acceptance | 1 / 1 | 1 -> 1 | 0 |
| V8 realizable KE scalable channel MKM DNS profile and momentum balance | 18 / 18 | 0 -> 0 | 0 |
| FV V2 manufactured Navier Stokes hex refinement | 16 / 16 | 0 -> 0 | 0 |
| FV V2 manufactured Navier Stokes tet refinement | 18 / 18 | 8 -> 8 | 0 |
| FV V1 manufactured scalar diffusion and upwind refinement | 48 / 48 | 0 -> 0 | 0 |
| FV SIMPLEC develops planar Couette flow with moving and slip walls | 14 / 14 | 0 -> 0 | 0 |
| V5 TMR flat plate SST accuracy and three-grid trend | 36 / 36 | 6 -> 2 | 11 |
| V5 TMR bump SST station accuracy and three-grid trend | 47 / 47 | 8 -> 6 | 12 |
| V7 SA TMR flat plate and bump accuracy convergence and three-grid trend | 70 / 70 | 14 -> 9 | 24 |
| Transition: ERCOFTAC T3A skin friction vs experiment (REGRESSION GUARD) and vs the documented equations (VERIFICATION) | 161470 / 161470 | 0 -> 0 | 0 |
| Transition: Selig S1223 at Re 2e5, alpha 4.76 deg vs UIUC | 11 / 11 | 0 -> 0 | 0 |
| FSAE V7.1 translating-wall Couette hex and poly | 104 / 104 | 0 -> 0 | 0 |
| FSAE V7.2 rotating-wall Taylor Couette hex and poly | 6668 / 6668 | 0 -> 0 | 0 |
| FSAE V8 Darcy Forchheimer pressure drop and oblique directions | 24 / 24 | 0 -> 0 | 0 |
| FSAE V9 absolute-velocity MRF analytic solutions and frame invariance | 48 / 48 | 0 -> 0 | 0 |
| FSAE V10 incompressible pressure and mass-flow inlets | 12 / 12 | 0 -> 0 | 0 |
| FSAE V11 internal fan polynomial and porous jump pressure discontinuity | 312 / 312 | 0 -> 0 | 0 |
| FSAE V7.3 SST moving-wall and translating-frame invariance | 16 / 16 | 0 -> 0 | 0 |
| FSAE V9.2 conformal MRF interface preserves absolute Taylor Couette | 4 / 4 | 0 -> 0 | 0 |
| FSAE V10.2 SST pressure and mass inlets match velocity-inlet transport | 12 / 12 | 0 -> 0 | 0 |
| FSAE V9.3 rotation and translation preserve Taylor Couette pressure | 6 / 6 | 0 -> 0 | 0 |
| Fluent research downloads retain published byte identities | 264 / 264 | 0 -> 0 | 0 |

### Changed-value audit

There are 5 changed values outside TMR. Each matches
both the identity and the exact before/after expanded strings in the
31-value comparison committed with main's `22e65c7`. They are the standard
Ahmed drag guard and four thin-layer Ahmed guards: SIMPLEC mean y+, SIMPLEC
drag error, Coupled/SIMPLEC drag agreement and lift agreement. There are no
other changed-value findings. Counts, identities and hashes are in
[`changed-value-audit.json`](changed-value-audit.json); restricted expanded
values remain in local build artifacts. Public TMR assertion changes are in
[`tmr-assertion-changes.json`](tmr-assertion-changes.json).

| Changed-value case | Count |
| --- | ---: |
| FSAE aero: Ahmed body 35 deg drag, level 1 | 1 |
| FSAE aero: Ahmed body 35 deg at y+ ~ 1, level 1, SIMPLEC and Coupled | 4 |
| V5 TMR flat plate SST accuracy and three-grid trend | 11 |
| V5 TMR bump SST station accuracy and three-grid trend | 12 |
| V7 SA TMR flat plate and bump accuracy convergence and three-grid trend | 24 |

Relative to the complete `86a665a` candidate, final TMR convergence outcomes
have 2 gains and 2 losses. This is separate from the required comparison
against main `8866cc0`; zero regressions against main does not mean every
historical candidate pass is retained. Compensated momentum remains opt-in;
the table records this trade-off within its TMR setup.

| TMR case | Backend | Level | Converged at 86a665a -> final | Steps at 86a665a -> final |
| --- | --- | ---: | --- | --- |
| SST plate | cpu | 0 | False -> True | 3000 -> 300 |
| SST plate | metal | 1 | True -> False | 800 -> 3000 |
| SST plate | metal | 2 | False -> True | 3000 -> 3000 |
| SA plate | metal | 1 | True -> False | 900 -> 3000 |


The historical complete `86a665a` run had 51 total failures, 17 TMR failures,
11 recoveries and 48 changed values, all in TMR. Its public measurements remain
in [`measurements-86a665a.json`](measurements-86a665a.json). The interrupted
`d823d6e` run is superseded by this complete rebased run; its partial historical
measurements remain explicitly labeled in
[`measurements-d823d6e-partial.json`](measurements-d823d6e-partial.json).
No partial XML is used for the final acceptance comparison.

| Evidence | Result |
| --- | --- |
| Quick tier | 307/307 cases and 7,990,154/7,990,154 assertions pass |
| Portability lint | passes |
| Elbow default CLI byte parity | 30/30 file comparisons pass, six CPU/Metal case pairs |
| Standard k-epsilon JumpCorrected | CPU/Metal cell fields and wall diagnostics identical to main |
| Full slow comparison | zero regressions and zero missing/new assertions |

| Artifact | SHA-256 |
| --- | --- |
| Main baseline XML | `9e60b5b5c8795f5aaacd4465e94cc67ee7bcccce12f0bf1f848bc70e164ce26a` |
| Final test binary | `133aed593c8ffe03df9c05b2bb4dd7b296f27355e39b0534cdfc20c90de816b0` |
| Final slow XML | `5162e44126c63274c64901f87230c9d7d36f26604eccd691c0e06c7969bc378d` |

The final process used `MallocSpaceEfficient=1`, matching the existing CLI
allocator setup. It recorded 1,253,163,008 bytes maximum RSS and 5,542,384,392 bytes
peak memory footprint. This exceeds the approximate 3 GB guidance; one heavy process ran at a time.
Baseline memory was not recorded in the supplied artifact. The earlier
complete run, before enabling allocator relief, had a 6,075,110,392-byte peak.
These are resource measurements, not performance claims.

Reproduction from the checkout root, with read-only fixtures linked from the
shared checkout and the fixture environment set:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
python3 tools/lint_portability.py
export MallocSpaceEfficient=1
build/kutta-tests --test-suite-exclude=slow
build/kutta-tests -ts=slow -r=xml --success > branch.xml
python3 tools/forensics/compare_slow_assertions.py "$BASELINE_XML" branch.xml
```

`BASELINE_XML` is the completed coordinator baseline for main `8866cc0`.
`KUTTA_FLUENT_EXAMPLES` and `KUTTA_RUNS` select its same read-only fixture
sources. Full XML, restricted outputs and off-target expanded values are local
ignored build artifacts. No home-directory paths or restricted data are added
to the committed evidence.

## Remaining work

There are 17 remaining TMR convergence assertions.
Every solve completed within its unchanged 3,000-step budget. All TMR accuracy
gates now pass, but the remaining nonlinear convergence mechanism is not
isolated. The arithmetic checks prove specific lost terms and derivatives;
they do not prove a sufficient fix for every floor. The scalar-update probe
also rules out claiming every bump failure is an all-cell one-ulp storage stall.

SA-neg's negative branch and Fluent automatic SA wall treatment remain
unsupported. They are model-coverage gaps; this investigation did not prove
them to cause the remaining convergence failures. The pinned reference variant,
ft2, numerical gates and budgets remain unchanged.

CUDA compilation/runtime is unverified. Default elbow parity and the full
slow comparison validate integration with main `271b31f` on CPU and Metal.
After validation, this lane's fixture links and temporary main-source worktree
were removed. No lane processes remain running. No push or PR was made.

## 2026-11-03–04 UTC: defaults and remaining convergence study

This follow-up is in progress on base `e014e36`. Five of twelve default
studies are complete. Momentum compensation and slack 16 regress gates;
tighter momentum and pressure tolerances preserve every gate. Scalar
compensation also preserves every gate outcome, with passing CLI INFO
identity changes retained in the comparison. No new default or
convergence fix has been accepted. The 3,000-step budget, all reference values and every
outer gate remain unchanged. The shared exact-main baseline and the first
complete default-candidate comparison are recorded below. Other candidates
and the causal convergence experiments remain unfinished.

The first candidate changes only the default of `compensated_momentum` to
true. Its quick tier passes 307 cases and 7,990,154 assertions. The complete
slow comparison rejects this default: the Fluent standard k-epsilon elbow
pressure-drop relative error becomes 2.00753% on CPU and Metal, failing the
unchanged 1.5% gate that main passes. It stays opt-in.
An earlier attempt was stopped in its first case because
memory profiling appeared in captured CLI assertion context; that interrupted
XML is excluded. The replacement runs without profiling. The comparator is
main's unchanged `tools/forensics/compare_slow_assertions.py`.

### Reboot and OS control, 2026-11-04 UTC

The owner paused this work. The Mac subsequently rebooted and changed from
Darwin 25.6.0 (recorded macOS 26.6.2) to macOS 27.2, build 26B123,
Darwin 27.2.0. The matrix controllers were gone on resumption; this agent
did not kill those CFD processes. The scalar-steady study had 22 complete,
receipted cases. Case 23 was interrupted and has no completed receipt.
The remaining six option/group studies had not run. The uncommitted
scalar-steady candidate header remained intact.

The frozen exact-main binary still loads and matches its archived hash.
Before reusing any completed case, it replayed V6's CPU/Metal seeded and
cold-start solves and all six SST flat-plate solves on the new OS. A third
control replayed the frozen scalar-steady candidate's V6 case against its
own archived result. Main's unchanged canonical comparator reports zero
regressions, recoveries, changed values, missing/new assertions or
missing/new cases for **all three controls**. Their iteration counts also
reproduce. The known convergence/accuracy failures remain; these controls
prove reproduction, not recovery. No binary rebuild was needed.

The scalar-steady candidate's quick tier also passes again on the new OS:
307 cases and 7,990,154 assertions. Its measured new-OS process peak is
1,331,791,768 bytes, below the guard.

The resume validator separately verifies all 22 complete cases against their
frozen test/CLI hashes, commands, resource receipts, manifest and original
assertion fragments. The interrupted case-23 XML and stderr are preserved
in a reboot archive and excluded. Case 23 restarts in full on the new OS.
The resumed study records old-OS cases 1–22 and new-OS cases 23–37 as
separate measurement segments. The controls support this numerical reuse;
they do not establish universal OS parity. A later failure to reproduce an
otherwise unchanged case requires a fresh comparison, not mixed values.

All cost pairs recorded so far (options 1–6) belong entirely to the old OS.
They remain historical measurements. Future timing pairs run **both sides
on the new OS**, using a fresh OS-specific baseline cache while preserving
the old cache. No timing ratio crosses the OS change. Process-memory
samples also retain their OS labels; an OS-related footprint change is not
claimed as solver memory savings. Control hashes, complete comparison
counts, CPU/Metal iteration samples and new-OS peaks are in
[os27-numerical-control.json](os27-numerical-control.json).

### Option decisions and allocation cost

Momentum compensation and slack 16 are rejected by their completed full slow tiers;
slack 16 also fails the quick tier. Momentum and pressure tolerances preserve
all slow gates and pass the quick tier. Scalar compensation also preserves
all outcomes and passes quick; its canonical comparison retains 76 passing
missing/new INFO identities. Promotion awaits the remaining
studies and combined-default proof. Other decisions remain pending. Allocation measurements
come from independent one-step CPU plate probes, changing one setting at a
time from ordinary settings, on 816 cells. The compensation row instead uses
the measured fine-bump delta after the cleanup below. These allocations are
not convergence or timing measurements. The exact data and probe hashes are
in [default-study-memory.json](default-study-memory.json) and
[lean-sst-parity.json](lean-sst-parity.json).

| Candidate | Slow-tier effect | Iteration cost | Extra device B/cell, SST / SA | Decision |
| --- | --- | --- | ---: | --- |
| Compensated momentum | Full: 2 new elbow failures, 0 recoveries; 1,173 changed values, 84 missing/new INFO identities | Startup time/step: CPU 2.91x, Metal 1.22x; ordinary SST plate: unconverged at 3,000 → converged at 300; SA still fails at 3,000 | 316.76 / 356.76, fine bump | Keep opt-in |
| Paired-state slack 16 | Full: 45 regressions, 0 recoveries, 1,237 changed values; 84 missing / 88 new identities; quick adds two failures | Paired elbow time: CPU 1.002x, Metal 1.559x; finest cavity steps: 900 → 2,800 CPU / 4,400 Metal | 0 / 0 | Keep opt-in |
| Momentum absolute tolerance 1e-20 | Full: 0 regressions/recoveries, 6 changed cavity values, 0 missing/new identities; quick passes | Ordinary SIMPLEC plate momentum inner work: 3.19–3.47x; all four solves still fail at 3,000. Late step time, SST / SA: CPU 0.988x / 0.936x, Metal 1.046x / 1.260x; Metal SA total 1.187x; larger T3A whole case 1.027x (1.006–1.048) | 0 / 0 | Numerically eligible; promotion pending remaining studies and combined proof |
| Pressure absolute tolerance 1e-20 | Full: 0 regressions/recoveries/changed values, 0 missing/new identities; quick passes | All four ordinary plates still fail at 3,000. Late time, SST / SA: CPU 1.021x / 1.012x, Metal 1.006x / 0.980x; SA inner work unchanged | 0 / 0 | Numerically eligible; promotion pending remaining studies and combined proof |
| Compensated scalar transport | Full: 0 regressions/recoveries, 231 changed values, 76 missing/new passing INFO identities; quick passes | All four ordinary plates still fail at 3,000. SST scalar inner work: 72,239→10,677 CPU / 73,502→10,910 Metal; SA unchanged. Late time, SST / SA: CPU 0.845x / 0.992x, Metal 0.634x / 1.079x | 32 / 16 | Gate outcomes eligible; promotion pending remaining studies and combined proof, with identity gaps disclosed |
| Scalar pseudo-time off | Pending | Pending | 0 / 0 | Pending |
| Scalar steady tolerance 1e-6 | Full study running; quick passes | All four ordinary plates still fail at 3,000. SST scalar inner work: 72,239→4,677 CPU / 73,502→4,530 Metal; SA unchanged. Late time, SST / SA: CPU 0.820x / 1.000x, Metal 0.629x / 1.014x | 24 / 12 | Pending |
| Scalar absolute tolerance 1e-20 | Pending | Pending | 0 / 0 | Pending |
| Fixed flow pseudo-time .1 | Pending | Pending | 0 / 0 | Pending |
| Three inner absolute tolerances together | Pending | Pending | Pending group measurement | Pending |
| Four scalar controls together | Pending | Pending | Pending group measurement | Pending |
| Combined precision controls | Pending | Pending | Pending group measurement | Pending |

Momentum relaxation .7 is already the ordinary solver default. The bump
harness's fixed time step is a separate candidate from its relaxation.
The planned groups are the three inner absolute tolerances, the four scalar
controls, and the combined precision controls. A group does not substitute
for the singleton checks. Any regression leaves the setting opt-in.

After the retained cleanup below, compensated momentum adds 17,840,128
device bytes on the 56,320-cell SST bump (316.76 B/cell). SA adds 20,092,928
bytes (356.76 B/cell). These deltas keep the TMR scalar controls fixed. With
walk reconstruction active, the increments are
`148*Ncells + 28*Nfaces + 28*Nboundary` for SST and
`188*Ncells + 28*Nfaces + 28*Nboundary` for SA. They include the 64 B/cell
paired walk scratch allocated on first use. The original allocation receipt
records SST before that cleanup: 19,867,648 bytes (352.76 B/cell); the linked
parity receipt records its exact 2,027,520-byte reduction. These are device
allocation deltas, not process-footprint savings.

### Remove unused SST fits without changing its solution

The compensated SST path allocated and evaluated three general velocity
gradient fits. SST instead consumes its dedicated strain fit, which overwrites
the shared low scratch before using it. The general fits have no consumer on
this path. Omitting their nine cell arrays saves exactly 36 B/cell and removes
three unused gradient calls per model-velocity refresh. Other models retain
their general fits. No time saving is claimed without paired timing runs.

On the fine bump, SST device bytes at solver-ready fall from 179,067,872 to
177,040,352, a 2,027,520-byte reduction. Compensation's extra device storage
therefore falls from 352.76 to 316.76 B/cell. SA remains at 163,938,028 bytes
with compensation active.

The following checks compare exact-main `e014e36` with this cleanup, retaining
the ordinary defaults. All checks pass on CPU and Metal. Hashes and complete
comparison results are in [lean-sst-parity.json](lean-sst-parity.json).

| 40-step path | Compared outputs | Result |
| --- | --- | --- |
| V6 mixing elbow, `wfgc` | fluxes, forces, loads, VTU, residual columns 1–7 | Identical on both backends |
| elbow_param, `compressibility,ke-options` | Same five outputs | Identical on both backends |
| Finest SST bump, compensation active | Cell fields, wall Cf, full history | Identical on both backends |
| Finest SST and SA bumps, compensation active | Binary downloads of public velocity, pressure, mass flux, boundary pressure, scalars and eddy viscosity | Identical for all four model/backend pairs |

The cleanup's full quick tier passes 307 cases and 7,990,154 assertions, and
portability lint passes. The completed isolated default-candidate experiment
includes this cleanup; its source and runtime qualifications are recorded
below. No default was promoted by the byte-parity checks.

### Ordinary-solve cost of compensation

Four alternating baseline/candidate pairs per backend (`AB, BA, AB, BA`)
measure the V6 SST mixing elbow with ordinary settings, a cold start and
40 outer steps. Only compensated momentum is enabled in the candidate,
which includes the unused-fit cleanup. Timing uses steps 21–40. The machine
is contended, so these are short-startup measurements, not an end-to-end
throughput claim. Raw measurements and binary hashes are in
[compensation-cost.json](compensation-cost.json).

| Backend | Median baseline ms/step | Median compensation ms/step | Median paired ratio | Paired ratio range |
| --- | ---: | ---: | ---: | ---: |
| CPU | 99.40 | 291.89 | 2.906 | 2.575–4.029 |
| Metal | 63.29 | 77.38 | 1.224 | 1.128–1.294 |

Both binaries use four Coupled block iterations at step 40 in every repeat.
The option materially increases this ordinary solve's measured startup cost,
in addition to its pressure-drop regression on the separate k-epsilon case.

Two resource-only experiments did not explain the previous 5.54 GB full-suite
peak. A single Metal device executing 102,400 launches remains near 111 MB
with and without an outer autorelease pool. Six create/use/destroy cycles
remain near 107 MB after destruction in both modes. These probes rule out
those specific simple lifetime paths; they do not establish that every
runtime path is free of retention.

A separate replay-buffer allocation prototype rounded each list to native
pages and reused the smallest adequate completed buffer. Its 40-step public
states and the required CLI output comparisons were byte-identical. Fine SST
bump peak footprint fell from 333,562,672 to 324,174,592 bytes; SA fell from
319,521,536 to 311,722,752. Four alternating ordinary-SST pairs measured a
Metal paired time ratio of 1.235 (range 0.942–1.432); the CPU control ratio
was 1.029 (0.891–1.132). These contended measurements do not justify shipping
the allocation change for a roughly 9 MB saving. The prototype was reverted.
See [sequence-buffer-experiment.json](sequence-buffer-experiment.json).

The original single-process compensated-default slow run was stopped after
14 complete cases when its footprint reached 3,852,733,848 bytes. That partial
run is not a full-tier comparison. The replacement executes all 37 slow test
cases in fresh processes, with their original assertions and solve budgets.
It preserves each actual test-case XML element and uses the unchanged main
comparator. There are no expression or INFO rewrites. This measures individual
case peaks without allocator/driver retention across unrelated cases. A
3,000,000,000-byte guard stops the run rather than hiding a resource excess.
The isolated thin-wall Ahmed case completes both SIMPLEC and Coupled with
eight passing assertions and a 2,307,148,032-byte peak. The complete 37-case
isolated run finishes below the guard; its heaviest case is this thin-wall
Ahmed case. This is a different execution protocol
from the original single-process suite, not a like-for-like solver memory
reduction.

The replacement compensation binary includes the byte-identical SST cleanup
and the since-reverted replay-buffer prototype. Its exact hashes are retained
in the experiment metadata. Its numerical gate effects are compared, while
its memory/time results are kept separate from the final source. Future
option binaries use the retained source only.

A turn interruption stopped the momentum-tolerance study's controllers after
25 fully receipted cases. The in-flight transitional S1223 process survived;
it was monitored to completion below 3 GB. Its output is preserved separately
and excluded from the combined tier because the original waiting parent's
exit receipt was unavailable. Only that unreceipted case is repeated under
the continuous guard. The first 25 cases are reused after checking their
frozen test/CLI hashes, manifest order, commands, resource receipts, canonical
assertion completeness, and exact agreement between each original XML case
and its saved fragment. Resume validation accepts that prefix and rejects a
different CLI binary before starting a solve. Gates and budgets are unchanged.

### Additional scratch reuse, validated privately

An isolated source variant reuses the low-word gradient scratch across
sequential component fits, and the base scratch across paired strain fits.
The fitted output gradients remain separate. It keeps the same kernel calls
and their order. On the 56,320-cell finest bump it saves another 48 device
B/cell for SST and 24 B/cell for SA, on both CPU and Metal:

| Model | Device bytes before → after | Additional saving |
|---|---:|---:|
| SST | 177,040,352 → 174,336,992 | 2,703,360 bytes / 48 B per cell |
| SA | 163,938,028 → 162,586,348 | 1,351,680 bytes / 24 B per cell |

The required 40-step CPU/Metal common-case comparisons keep all 20 file
checks identical against main. Six finest-SST-bump CSV comparisons and four
SST/SA public-state comparisons also remain byte-identical. These checks are
recorded in [gradient-scratch-parity.json](gradient-scratch-parity.json).
No timing improvement is claimed. This change remains private while the
per-option defaults matrix uses a frozen source; it is not yet a retained
solver saving. The currently retained change still saves 36 B/cell for SST.

Four subsequent alternating AB, BA pairs per model and backend measure the
same fine-bump configuration for 40 steps, timing only steps 21–40. All 32
runs reproduce the same physical reports and cumulative inner work within
each model/backend. The measured additional device saving remains exactly 48
B/cell for SST and 24 B/cell for SA in every pair.

| Model / backend | Paired time ratio, median (range) |
|---|---:|
| SST CPU | 1.0213 (0.9890–1.0303) |
| SST Metal | 0.9930 (0.9596–1.0907) |
| SA CPU | 0.9994 (0.9068–1.1303) |
| SA Metal | 0.9862 (0.9527–1.0054) |

These contended samples establish no speedup. They show a small SST CPU
median cost, rather than a consistently large penalty in this window.
Whole-process peaks must remain separate from exact buffer bytes. SST CPU
peaks span 154.6–198.6 MB before and 154.8–209.8 MB after; SA CPU spans
155.0–171.2 MB before and 153.5–200.0 MB after. CPU process-peak improvement is
unproven, and no cause is assigned to that variation. The unchanged baseline
binary itself has variable peaks. Metal process-peak medians fall from
323.8 to 320.9 MB for SST and 314.0 to 312.5 MB for SA; each of their four
pairs falls. The change remains unmerged until the matrix releases the root
source. Data, controls and frozen input hashes are in
[gradient-scratch-cost.json](gradient-scratch-cost.json).

### Complete compensation-default comparison

All 37 slow cases execute with 201,633 assertions. Main has 51 failures;
the compensated-default candidate has 53. The unchanged comparator reports
two regressions, zero recoveries and 1,173 changed expanded values. It also
reports 84 missing and 84 new identities because physical CLI residual logs
inside INFO changed. There are no missing or added test cases. These identity
changes are retained; this is not a clean canonical comparison.

The two new failures are k-epsilon elbow pressure-drop relative errors:
CPU `0.0101201 -> 0.0200753`, Metal `0.0101167 -> 0.0200753`, against the
unchanged `0.015` limit. All 153 TMR assertions retain their exact identities,
outcomes and expanded values. Compensated momentum stays opt-in.

[default-study-compensation.json](default-study-compensation.json) contains
the complete case-by-case assertion counts, failed-gate counts and changed-value
counts, the original binary/source hashes and the explicit runtime-prototype
qualification. The full XML and raw comparison stay outside version control.
The full quick tier remains green: 307 cases, 7,990,154 assertions.

### Slack 16: full slow-tier rejection and cost

Changing only the ordinary slack default from 1 to 16 adds two failures in
`Normal-foot pressure boundary gradients recover uniform flow on native tetrahedra`,
for the moving-wall variant. The errors are `1.48279e-5` on CPU and
`1.48278e-5` on Metal against the unchanged `1e-5` gate. The run executes all
307 quick cases: 306 pass, with 7,990,154 passing and two failing assertions.
The candidate test binary SHA-256 is
`19b03b060494386aa685a2ec911d51357ad7ef9c953d69fe7444f15141cd0591`.
Slack 16 stays opt-in. Its full slow tier completes all 37 cases with 201,637
assertions and 96 failures, against main's 201,633 assertions and 51 failures.
The unchanged comparator reports 45 regressions, zero recoveries, 1,237 changed
expanded values, 84 missing identities and 88 new identities. No test cases
are missing or added. The canonical comparison is not clean.

The 84 missing/new pairs are CLI-context identities in conversion cases.
The other four new assertions are passing V6 finite-value and maximum-speed
safety checks. The unchanged loop executes two extra 50-step checkpoints:
saved-state CPU convergence takes 200 → 300 steps. Cold-start CPU stays 250;
Metal saved/cold stay 200/250. These are measured iteration costs, not
contended wall-clock comparisons. All identity changes are retained.

All 45 regressions are in the Fluent cavity case, comparing imported-mesh
and equivalent generated-mesh velocity profiles. That case has 715 changed
expanded values and no missing/new identities. Main passes all 763 assertions;
the candidate passes 718 and fails 45. The failures split into u/v counts of
12/11 on CPU and 15/7 on Metal.

The corrected cost run checks every emitted setting against an exact-main
solver baseline. Four alternating elbow pairs measure median time ratios of
1.0017 on CPU (range 0.9439–1.0311) and 1.5586 on Metal (1.1153–2.4918).
The Metal slowdown is observed under contention; no cause is attributed.
At 40 steps, the coupled block uses 162 inner iterations for both settings
on both backends. This is stored in `IterationReport.pressure`, which the
benchmark labels `pressure`; it includes the joint momentum/pressure solve. CPU scalar work is 439 → 438, Metal 438 → 438. Peak footprint is
63,832,520 → 89,965,000 bytes on CPU and 205,537,976 → 205,538,000 on Metal.
These process peaks do not establish a setting-specific allocation increase;
the independent one-setting device measurement is unchanged at zero bytes.

The completed Re1000 cavity case also shows a material iteration cost. All
16 assertions pass, but the same stopping rule takes:

| Grid | CPU outer steps, before → after | Metal outer steps, before → after |
|---|---:|---:|
| 65² | 600 → 1,000 | 600 → 1,000 |
| 129² | 600 → 1,100 | 600 → 1,000 |
| 257² | 900 → 2,800 | 900 → 4,400 |

These are iteration counts on identical cases and gates, not contended
wall-clock ratios. The finest-grid cost increases by 211% on CPU and 389%
on Metal. The complete slow comparison confirms the rejection.

On the ordinary-control coarse plate, SST remains unconverged at 3,000
baseline steps but converges at 400 with slack 16 on both backends. SA remains
unconverged on CPU at 3,000, while Metal reaches the gate at 500. These are
separate ordinary-control probes, not the TMR harness or its 17 failures.
[default-study-slack.json](default-study-slack.json) preserves the hashes,
paired measurements, cumulative work and quick-tier rejection.

The completed thin-wall Ahmed case passes all eight assertions and peaks
at 2,215,168,184 bytes, below the 3,000,000,000-byte guard. This slack candidate
uses the retained original Metal runtime. It is a different default candidate
from the earlier compensation run, so their peak difference is not attributed
to a memory optimization. All 37 case processes remain below the 3 GB guard.

A private field-motion diagnostic observes the true high-plus-low flow pairs
and stored scalar fields without extra monitor or force calls. With downloads
enabled every step, four 40-step finest-bump SST/SA CPU/Metal comparisons keep
all public state bytes identical. The full slack run uses this instrumented
binary to gather per-component residuals and lag-one/lag-two field motion.
The diagnostic data cannot establish a cause until examined and tested with
an isolating experiment.

### Corrected ordinary-control iteration cost of compensation

These probes use the coarse TMR plate geometry with ordinary solver controls,
independently of the TMR harness opt-ins. The baseline is rebuilt from exact
main `e014e36` solver sources, with only the independent diagnostic benchmark
added; all 165 solver source files match main. Each sample's emitted settings
are checked with
[`check_solver_controls.py`](../../../tools/forensics/check_solver_controls.py).
That check accepts the corrected baseline log and rejects the actual mislabeled
log. Enabling only compensated momentum gives:

| Model/backend | Outer steps, before → after | Momentum inner steps | Pressure inner steps | Scalar inner steps |
|---|---:|---:|---:|---:|
| SST cpu | 3,000 → 300 | 8,974 → 1,093 | 60,635 → 5,329 | 72,239 → 4,453 |
| SST metal | 3,000 → 300 | 9,303 → 1,080 | 58,694 → 4,952 | 73,502 → 4,475 |
| SA cpu | 3,000 → 3,000 | 9,415 → 4,076 | 34,727 → 24,674 | 417 → 416 |
| SA metal | 3,000 → 3,000 | 9,750 → 9,442 | 34,514 → 24,677 | 417 → 417 |

The SST baseline remains unconverged at 3,000 steps on both backends;
compensation converges in 300. Both SA variants remain unconverged at 3,000.
SA's terminal maximum defect drops from `2.04057944461e-4` to
`1.31528557856e-5` on CPU, and from `1.83064395142e-4` to
`3.40673138094e-5` on Metal. No unpaired timing comparison is claimed.
The earlier alternating elbow timings remain a separate experiment.

This is a useful opt-in trade-off, not a universal default: the full slow
tier adds two elbow failures. Both iteration-cost binaries use the original
Metal runtime; the candidate includes the byte-identical SST memory cleanup.
[compensation-iteration-cost.json](compensation-iteration-cost.json) records
the hashes, controls and endpoints.

The first supplemental comparison mislabeled a compensation-enabled binary
as baseline. Its equal-work conclusion is withdrawn. The interrupted slack
cost run used that same incorrect baseline and is excluded. Both attempts
are retained privately, and the public receipt identifies the rejected
baseline. Neither invalid cost run changes the completed slow-tier comparison
or slack's quick-tier regression.

### Momentum tolerance: measured inner work

Changing only momentum absolute tolerance from `1e-10` to `1e-20` passes the
quick tier (307 cases, 7,990,154 assertions). The full slow tier now passes its
comparison with main, as recorded below.
The ordinary-control coarse plate uses SIMPLEC, where this tolerance is
active. Its cumulative momentum work and final maximum defect are:

| Model / backend | Momentum inner iterations, before → after | Maximum defect, before → after |
|---|---:|---:|
| SST CPU | 8,974 → 30,084 | 1.487815e-4 → 6.880296e-5 |
| SST Metal | 9,303 → 32,275 | 1.835793e-4 → 8.711917e-5 |
| SA CPU | 9,415 → 32,631 | 2.040579e-4 → 1.940796e-4 |
| SA Metal | 9,750 → 31,136 | 1.830644e-4 → 1.914472e-4 |

All four before/after solves remain unconverged at 3,000 steps. There is no
direct device allocation increase. These are iteration measurements, not
paired timing results. Data and binary identities are in
[momentum-absolute-cost.json](momentum-absolute-cost.json).

The existing four alternating elbow pairs use its saved Coupled scheme.
That path solves a joint block and does not use the separate momentum
absolute tolerance. Its CPU time ratio of 0.9984 and Metal ratio of 1.0098
therefore do not measure this option's active-path cost. They must not hide
the extra SIMPLEC work above.

`tools/forensics/measure_default_plate_cost.py` adds alternating AB, BA runs
on this ordinary-control plate, with both turbulence models and backends.
It retains early convergence and the original 3,000-step budget. It reports
time to the actual stop and late per-step time only when both runs reach
3,000; a failed baseline and successful candidate do not imply equal
accuracy. Every repeated numerical endpoint and inner-work total must
reproduce the checked receipt exactly before its timing is accepted.
Binary, receipt and mesh hashes are frozen, and actual controls are checked
against independent ordinary-setting definitions. Validation accepts the
three existing corrected receipts and rejects an already-compensated
baseline despite its claimed checked marker. The paired plate measurements below are complete for these first three
options.

### Complete momentum-tolerance comparison

Only the ordinary momentum absolute tolerance changes from `1e-10` to
`1e-20`. The full tier completes all 37 cases and 201,633 assertions, retaining
main's 51 failed assertions. Main's unchanged comparator exits zero:
zero regressions or recoveries, six changed expanded values, and zero
missing/new assertions or cases. The quick tier also passes all 307 cases
and 7,990,154 assertions. This candidate qualifies numerically for a default;
promotion awaits the remaining option studies and final combined-default
proof. Its active-path paired cost is recorded below. The matrix continues with ordinary defaults restored between candidates.

All six changed values are passing Re=1,000 cavity checks. Their changed
quantities are:

| Backend | Quantity | Before → after |
|---|---|---:|
| CPU | Middle-grid v deviation | 0.00711187 → 0.00711191 |
| CPU | v refinement order | 2.71436 → 2.71437 |
| Metal | Middle-grid u deviation | 0.00435289 → 0.004353 |
| Metal | Middle-grid v deviation | 0.00711203 → 0.00711196 |
| Metal | u refinement order | 2.87464 → 2.87468 |
| Metal | v refinement order | 2.71455 → 2.71454 |

All 153 TMR assertion identities, outcomes and expanded values match main
exactly: plate still has two convergence failures, bump six, and SA nine.
ERCOFTAC T3A uses 6,000 outer iterations in both runs. The full-tier elapsed
times are unpaired and are not attributed to this option. The measured
ordinary-plate inner-work increase above remains material.

The heaviest completed case is thin-wall Ahmed at 2,214,447,264 bytes; the
heaviest TMR case is 389,645,248 bytes. Both remain below the 3 GB guard.
After the controller interruption, the repeated S1223 case records its exit
and resource receipt and reproduces all 11 preserved assertion payloads
exactly. The excluded orphan contributes no assertion to the combined tier.
Complete per-case counts, source/binary/XML hashes, changed values and the
cost scope are in
[default-study-momentum-absolute.json](default-study-momentum-absolute.json).

### Paired ordinary SIMPLEC cost of the first three options

Each row uses two alternating AB, BA pairs on the 816-cell plate, with
ordinary controls and one changed option. Every endpoint and cumulative
inner-work total reproduces its earlier checked receipt. Binary, control,
mesh and receipt hashes are verified before and after the runs. The exact
rows, process peaks and hashes are in
[paired-ordinary-plate-cost.json](paired-ordinary-plate-cost.json).

Ratios are candidate / baseline. Late-step ratios cover steps 1,501–3,000
and exist only when both solves reach 3,000. The machine is shared; ranges
show both pairs. These small ordinary cases do not bound costs on other
meshes or models. A failed baseline and passing candidate achieve different
accuracy, so their time-to-stop ratio is not an equal-accuracy speedup.

| Option | Model / backend | Outer steps, before → after | Converged, before → after | Total solver time ratio | Late-step ratio (range) |
|---|---|---:|---|---:|---:|
| Compensated momentum | SST cpu | 3,000 → 300 | no → yes | 0.197 | — |
| Compensated momentum | SST metal | 3,000 → 300 | no → yes | 0.101 | — |
| Compensated momentum | SA cpu | 3,000 → 3,000 | no → no | 1.611 | 1.581 (1.580–1.582) |
| Compensated momentum | SA metal | 3,000 → 3,000 | no → no | 0.895 | 0.853 (0.828–0.879) |
| Slack 16 | SST cpu | 3,000 → 400 | no → yes | 0.125 | — |
| Slack 16 | SST metal | 3,000 → 400 | no → yes | 0.123 | — |
| Slack 16 | SA cpu | 3,000 → 3,000 | no → no | 0.762 | 0.711 (0.656–0.766) |
| Slack 16 | SA metal | 3,000 → 500 | no → yes | 0.161 | — |
| Momentum abs 1e-20 | SST cpu | 3,000 → 3,000 | no → no | 1.001 | 0.988 (0.985–0.990) |
| Momentum abs 1e-20 | SST metal | 3,000 → 3,000 | no → no | 1.087 | 1.046 (0.793–1.299) |
| Momentum abs 1e-20 | SA cpu | 3,000 → 3,000 | no → no | 0.975 | 0.936 (0.906–0.967) |
| Momentum abs 1e-20 | SA metal | 3,000 → 3,000 | no → no | 1.187 | 1.260 (0.984–1.536) |

Compensated momentum increases ordinary SA CPU late-step time by 58.1%
while decreasing it by 14.7% on Metal. It still fails both ordinary SA
convergence gates at 3,000. Slack 16 changes stopping behavior substantially,
but its full-tier and quick-tier regressions already reject a default.

The tighter momentum tolerance leaves all four ordinary solves unconverged
at 3,000, with no direct device allocation increase. In these pairs it adds
18.7% to SA Metal total time and 26.0% to its late-step time; the late-step
pair range is 0.984–1.536. SST Metal's total ratio is 1.087 and late ratio
1.046, also with a wide range. CPU totals are 1.001 for SST and 0.975 for SA.
The 3.19–3.47x momentum inner-work increase is deterministic. This option
must not be described as free, nor can these contended timings establish a
uniform slowdown across ordinary solves. No default has yet been changed.

### Momentum tolerance on a larger ordinary solve

The existing T3A transition case uses SIMPLEC, whose momentum solves read
this tolerance. Two alternating AB, BA pairs compare the frozen exact
main test binary with the archived momentum-tolerance candidate. All other
settings and runtime inputs match. This is the original T3A budget of 6,000;
the TMR budget remains 3,000. All four runs use Metal, reach 6,000 steps,
and reproduce all 161,470 original assertion payloads through the unchanged
canonical comparator.

| Pair and order | Baseline elapsed seconds | Candidate elapsed seconds | Candidate / baseline |
|---|---:|---:|---:|
| 0, AB | 367.03 | 369.39 | 1.0064 |
| 1, BA | 355.93 | 372.91 | 1.0477 |

The median paired increase is 2.7%, with a 0.6–4.8% pair range. These are
whole-test times including setup, output and assertions, on the shared
machine. They do not establish a cost bound on other cases. The candidate
peaks at 340,706,024 process bytes; baseline peaks vary from 340,443,856 to
355,762,920. No direct device allocation changes, and that process-peak
variation is not attributed to the tolerance. T3A's documented continuity
stall is not gated here; passing these unchanged assertions does not claim
physical convergence.

The donor checkout has advanced since `e014e36`. The baseline here is its
frozen copied binary, with its embedded runtime paths checked and its data,
mesh, ledger and independent equation reference matched to the candidate's
inputs. Its current C++ source and current build are not substituted. This
case calls the solver API and does not invoke the current root CLI. Protocol,
input and binary hashes, all four receipts and the scope limits are in
[momentum-t3a-cost.json](momentum-t3a-cost.json).

The ordinary SA plate's 18.7% total and 26.0% late-step Metal increases above
remain part of the decision. The larger SST transition result does not make
the option free across ordinary solves.

### Pressure tolerance: complete comparison and cost

Changing only pressure absolute tolerance from `1e-12` to `1e-20` passes the
quick tier: 307 cases and 7,990,154 assertions. The full slow tier completes
all 37 cases and 201,633 assertions, retaining main's 51 failed assertions.
The unchanged canonical comparator exits zero: no regressions, recoveries,
changed expanded values, missing/new assertion identities or missing/new
cases. Every recorded assertion payload matches main, including all 153 TMR
assertions and the 17 remaining convergence failures. This qualifies
numerically; promotion awaits the remaining studies and combined-default
proof. The matrix continues to restore ordinary defaults between candidates.

The heaviest completed case peaks at 2,201,946,248 bytes; the heaviest TMR
case peaks at 388,137,872 bytes. Both stay below the 3 GB guard. Full per-case
counts and source, binary, XML and comparator hashes are in
[default-study-pressure-absolute.json](default-study-pressure-absolute.json).

All four ordinary-control plate probes still fail convergence at 3,000.
SST CPU pressure inner work changes from 60,635 to 59,608, with final maximum
defect `1.487815e-4` → `1.485530e-4`. SST Metal changes from 58,694 to 59,039,
with defect `1.835793e-4` → `1.276474e-4`. Both SA endpoints and all three
inner-work counters reproduce baseline exactly; all 30 monitor and 30 cumulative-work samples
are identical. Direct device storage
is unchanged.

Two alternating AB, BA pairs per model/backend give the following late-step
time ratios over steps 1,501–3,000. They measure a fixed budget, with both
solves unconverged, and do not establish costs on other cases.

| Model / backend | Candidate / baseline, median (range) | Total solver time ratio |
|---|---:|---:|
| SST CPU | 1.021 (0.982–1.060) | 0.999 |
| SST Metal | 1.006 (0.929–1.083) | 1.025 |
| SA CPU | 1.012 (0.994–1.029) | 1.019 |
| SA Metal | 0.980 (0.962–0.998) | 0.982 |

These shared-machine timing variations accompany unchanged SA numerical work;
they are not evidence of a pressure-solver speedup. The paired probes peak
at 152,584,792 process bytes. The earlier elbow timing uses Coupled and its
separate pressure tolerance is inactive. Checked controls, endpoints, work,
binary hashes and paired ranges are in
[pressure-absolute-cost.json](pressure-absolute-cost.json).

### Scalar compensation: complete comparison and cost

The fifth singleton enables compensated scalar transport for SST and SA.
Its quick tier passes 307 cases and 7,990,154 assertions. The full slow tier
completes all 37 cases and 201,633 assertions, retaining main's 51 failures.
The unchanged canonical comparator reports zero regressions and recoveries,
231 changed expanded values, 76 missing and 76 new assertion identities,
and no missing/new cases. It exits **one**, because the 76 passing CLI
assertions carry changed physical solver histories in INFO. All 76 original
missing identities pass in baseline; all 76 new identities pass. No identity
is normalized or paired with a different one to turn that result green.

Every gate retains its outcome. All 153 TMR assertion identities and values
match exactly, so the 17 convergence failures remain. Scalar compensation
qualifies on gate outcomes, with the identity gaps disclosed; promotion still
awaits the remaining studies and combined-default proof. Full source,
binary and comparator identities are in
[default-study-scalar-compensation.json](default-study-scalar-compensation.json).

The 231 changes comprise 119 serialized-output equality assertions, 48
interpolation-bound assertions and 64 derived metrics. The first two groups
remain passing. The equality checks still confirm identical donor/converted
outputs within each candidate run; the output contents change across the
default study. Their full fields and interpolation reference endpoints stay
in the ignored local canonical comparison. All 64 remaining changed metric
pairs are in [scalar-compensation-effects.json](scalar-compensation-effects.json).

| Case | Changed expanded values |
|---|---:|
| Native Fluent physics conversion | 103 |
| T3A transition | 55 |
| V6 mixing elbow | 24 |
| General Fluent/Kutta conversion | 16 |
| Standard k-epsilon elbow / SST nozzle | 11 |
| SST moving-wall/frame invariance | 5 |
| Thin-wall Ahmed SIMPLEC/Coupled | 4 |
| SST pressure/mass inlets | 4 |
| Tyrrell lift | 3 |
| S1223 transition | 3 |
| S1223 without transition | 2 |
| Ahmed drag | 1 |

Eighteen existing failed V6 metrics change: ten improve and eight worsen.
Their largest absolute movement in the printed metrics is `1e-6`; none
recovers or introduces a failed gate. The maximum field L2 error is
`0.144782` → `0.144783` on both backends and still fails. T3A's verified
laminar/turbulent Cf errors change from `0.0169599` / `0.0134799` to
`0.0169701` / `0.0134984`, retaining their `0.03` gates. These are changed
values, not new accuracy acceptance. V6 still takes 200 seeded and 250
cold-start steps on both backends; T3A still takes its original 6,000 steps.

The heaviest completed case peaks at 2,216,233,072 bytes; the heaviest TMR
case peaks at 396,608,424 bytes. These fresh case processes stay below the
3 GB guard. They do not prove a monolithic-suite retention fix.

All four ordinary-control plate solves still miss convergence at 3,000.
SST scalar inner work falls from 72,239 to 10,677 on CPU and from 73,502
to 10,910 on Metal. SA scalar work remains 417 on both backends. Device
storage increases by 32 B/cell for SST and 16 B/cell for SA. This is scalar
storage alone, separate from compensated momentum's larger allocation.

Two alternating AB, BA pairs per model/backend measure late-step time over
steps 1,501–3,000. Candidate/baseline median ratios are:

| Model / backend | Late time ratio, median (range) | Total solver time ratio |
|---|---:|---:|
| SST CPU | 0.845 (0.841–0.849) | 0.855 |
| SST Metal | 0.634 (0.633–0.634) | 0.626 |
| SA CPU | 0.992 (0.982–1.002) | 0.996 |
| SA Metal | 1.079 (0.967–1.190) | 1.048 |

SST CPU's final maximum defect worsens from `1.487815e-4` to `1.649299e-4`;
SST Metal improves from `1.835793e-4` to `7.874356e-5`. Neither meets the
unchanged gate. Reduced inner work and faster failed-budget runs do not
establish a convergence or accuracy improvement. SA's paired Metal spread
includes both faster and slower steps, so its measured median is not a
uniform slowdown claim. The paired runs peak at 152,126,040 process bytes.
Actual settings, numerical endpoints, inner work, binary/receipt hashes and
ranges are in [scalar-compensation-cost.json](scalar-compensation-cost.json).

### Scalar steady tolerance: cost complete, full comparison running

The sixth singleton changes SST and SA scalar steady tolerance from zero to
`1e-6`. Quick passes 307 cases and 7,990,154 assertions. Its full slow tier
is running, so no numerical default decision is made yet.

On the ordinary-control plate, all four solves still miss convergence at
3,000. SST scalar inner work falls from 72,239 to 4,677 on CPU and from
73,502 to 4,530 on Metal. SA remains at 417; its terminal printed physical
monitors and all three inner-work totals match baseline on both backends.
This is printed-observation agreement, not full-state byte proof. Extra
device storage is 24 B/cell for SST and 12 B/cell for SA.

Two alternating AB, BA pairs per model/backend measure:

| Model / backend | Late time ratio, median (range) | Total solver time ratio |
|---|---:|---:|
| SST CPU | 0.820 (0.812–0.828) | 0.836 |
| SST Metal | 0.629 (0.619–0.638) | 0.650 |
| SA CPU | 1.000 (0.999–1.002) | 1.002 |
| SA Metal | 1.014 (1.011–1.017) | 1.014 |

Late means steps 1,501–3,000. Both baseline and candidate remain unconverged,
so reduced work and time are not convergence acceptance. In the separate
40-step Coupled elbow startup, scalar work increases from 439 / 438 to 464
on CPU / Metal. Four alternating pairs give late-step ratios of 1.012
(0.999–1.017) on CPU and 1.032 (1.023–1.045) on Metal. These shared-machine
measurements do not establish costs on other cases.

The paired plate runs peak at 149,570,088 process bytes. Actual controls,
reproduced endpoints/work, allocation deltas, hashes and all ranges are in
[scalar-steady-cost.json](scalar-steady-cost.json).

### Fresh remaining-failure baseline

The supplied exact-main run has completed all 24 TMR solves and reproduces
all 17 remaining convergence failures: SST plate 2, SST bump 6, SA 9. Their
maximum-defect and load-change endpoints match the preceding report. Main's
whole baseline subsequently completes all 37 slow cases with 51 failed
assertions in 16,993 seconds. The complete candidate comparison is above.
[default-study-baseline.json](default-study-baseline.json) records every
100-step observation and the supplied binary hash. The absent first load
comparison is null, not zero.

The completed TMR subset matches all 153 baseline assertion identities,
outcomes and expanded values exactly: SST plate 36 assertions with two
failures, SST bump 47 with six, SA 70 with nine. There are no changed values,
regressions or recoveries in that subset. This confirms that the existing TMR
opt-ins make the global compensation-default experiment numerically unchanged
on these cases. All 17 convergence failures remain.

### Initial classification of all 17 failures

The exact-main 100-step histories show 15 declining envelopes and two
plateaus. The ratio below is the maximum defect at step 3,000 divided by
its value at 2,000. This is an initial classification, not a causal result.
A declining envelope can contain a superposed two-cycle; a flat envelope
does not establish an arithmetic residual floor. Per-step diagnostics and
controlled experiments remain required for every row.

| Failing solve | Terminal maximum defect | Defect ratio, 3,000 / 2,000 | Initial observed shape |
|---|---:|---:|---|
| SST plate cpu L2 | 1.90981e-05 | 0.084034 | Slow approach in envelope |
| SST plate metal L1 | 3.045154e-05 | 1.0019 | Plateau |
| SST bump cpu L0 | 3.143857e-05 | 0.45292 | Slow approach in envelope |
| SST bump cpu L1 | 2.308443e-05 | 0.5394 | Slow approach in envelope |
| SST bump cpu L2 | 2.573423e-05 | 0.34281 | Slow approach in envelope |
| SST bump metal L0 | 1.991881e-05 | 0.30103 | Slow approach in envelope |
| SST bump metal L1 | 2.81366e-05 | 0.67352 | Slow approach in envelope |
| SST bump metal L2 | 2.572706e-05 | 0.34277 | Slow approach in envelope |
| SA plate cpu L2 | 4.282008e-05 | 0.071949 | Slow approach in envelope |
| SA plate metal L1 | 2.417789e-05 | 1.0005 | Plateau |
| SA plate metal L2 | 1.397759e-05 | 0.024291 | Slow approach in envelope |
| SA bump cpu L0 | 3.439757e-05 | 0.44899 | Slow approach in envelope |
| SA bump cpu L1 | 4.124983e-05 | 0.65014 | Slow approach in envelope |
| SA bump cpu L2 | 6.128235e-05 | 0.49282 | Slow approach in envelope |
| SA bump metal L0 | 5.78428e-05 | 0.70232 | Slow approach in envelope |
| SA bump metal L1 | 4.124282e-05 | 0.64986 | Slow approach in envelope |
| SA bump metal L2 | 6.129666e-05 | 0.49293 | Slow approach in envelope |

The cause of each row remains unestablished. No row is declared a stationary
floor or a two-cycle from even-spaced samples. The complete endpoints,
including load changes, are in
[initial-classification.json](initial-classification.json).

![Exact-main TMR maximum-defect histories, sampled every 100 steps](baseline-convergence.png)

The plot shows all 24 recorded solves, including the seven that converge.
A terminal cross means the unchanged convergence gates are not met by step
3,000; the dashed line shows the maximum-defect gate. Load change also applies.

### Per-step classification of the eight SST failures

All SST plate and bump field-motion captures have completed. The unchanged
canonical comparator confirms all 83 SST assertions retain their exact
identities, outcomes and expanded values: zero regressions, recoveries,
changed values or missing/new identities. The same eight convergence gates
still fail. This subset result is not the unfinished whole-tier comparison.

| Failing solve | k defect | Omega defect | Continuity defect | Observed failure mode |
|---|---:|---:|---:|---|
| SST plate cpu L2 | 9.906767e-07 | 4.136047e-07 | 1.90981e-05 | Slow approach with fluctuations |
| SST plate metal L1 | 9.306819e-07 | 5.835587e-07 | 3.045154e-05 | Pressure–velocity two-cycle |
| SST bump cpu L0 | 5.89877e-06 | 1.4441e-05 | 3.143857e-05 | Omega plateau + continuing flow motion |
| SST bump cpu L1 | 2.040069e-05 | 3.23637e-06 | 2.308443e-05 | Slow approach with fluctuations |
| SST bump cpu L2 | 2.573423e-05 | 8.517128e-07 | 1.210993e-05 | Slow approach with fluctuations |
| SST bump metal L0 | 6.84309e-06 | 1.991881e-05 | 1.51518e-05 | Omega plateau + continuing flow motion |
| SST bump metal L1 | 2.04067e-05 | 3.234025e-06 | 2.81366e-05 | Slow approach with fluctuations |
| SST bump metal L2 | 2.572706e-05 | 8.419657e-07 | 1.08788e-05 | Slow approach with fluctuations |

The fine CPU plate is a slow approach. Over steps 2,800–3,000, median relative
streamwise-velocity movement is `4.2498073e-9` per step; lag-two movement is
`1.9991` times lag-one. The pressure and mass-flux ratios are `1.3217` and
`1.7371`. Continuity remains the failing component.

The middle Metal plate has a strong pressure–velocity two-cycle with small
drift. Its median lag-two / lag-one movement ratios are `0.006692` for u,
`0.004219` for v, `0.001978` for pressure and `0.002364` for mass flux.
Pressure moves `2.00757324e-6` in relative L1 per step, almost returning after
two steps. k, omega and eddy viscosity are bitwise unchanged through all
201 sampled steps. This is observed alternation, not a stationary floor.

The coarse bumps need a more specific classification than their declining
maximum-defect envelopes suggested. Omega ends at `1.44410010558e-5` on CPU
and `1.99188050291e-5` on Metal. Its 3,000 / 2,000 ratios are `1.1753` and
`1.0494`; it plateaus while k and continuity decline. Omega still moves by
median relative L1 `1.97561e-6` / `2.27687e-6` per step. These are plateaus
with continuing nonlinear field motion; an arithmetic floor is not established.

On the middle and fine bumps, k continues moving with median lag-two /
lag-one ratios between `1.9947` and `1.9999`. The fine-grid k contraction over
steps 2,000–3,000 is about `0.99893` per outer step on both backends.
Fine Metal continuity also has a near-flat noisy tail and remains above its
gate. Accelerating k alone has not been shown to cure that component.

[field-motion-classification.json](field-motion-classification.json) records
all eight component profiles and motion metrics. The causes are still
unestablished. The SA classifications follow below. Controlled experiments remain;
no iteration budget or acceptance criterion has changed.

### SA plate classification

All three failing SA plate captures have completed. The fine CPU and Metal
plates show slow continuity approach: over the final 1,000 steps their
continuity ratios are `0.07195` and `0.02429`. Streamwise velocity lag-two /
lag-one movement is `1.9998` and `2.0001`; mass flux gives `1.9778` and
`1.9780`. Their final continuity defects remain `4.28200779547e-05` on CPU
and `1.39775922557e-05` on Metal. The transported scalar already meets its gate.

The middle Metal SA plate is a plateau with irregular paired-flow motion,
not the nearly repeating two-cycle seen in SST. Its pressure lag-two /
lag-one ratio is `0.9148`, with u `1.2461`, v `1.1657` and mass flux `2.2690`.
Transported viscosity and eddy viscosity are unchanged through all final
201 sampled steps. Continuity remains `2.41778890637e-05`. These observations
do not establish its cause or a stationary arithmetic floor.

The six SA bump classifications follow below. Controlled causal results
remain unfinished.

### All 17 classifications and complete TMR subset proof

The six SA bumps show slow transported-viscosity approach with continuing
field motion. The final 1,000-step scalar defect ratios and final 201-step
field-motion ratios are:

| Failing solve | Nu-tilde defect | Continuity defect | Scalar defect, 3,000 / 2,000 | Scalar field lag-two / lag-one |
|---|---:|---:|---:|---:|
| SA bump cpu L0 | 2.6725204e-05 | 3.4397567e-05 | 0.44044 | 1.9827 |
| SA bump cpu L1 | 4.1249825e-05 | 2.3684114e-05 | 0.65014 | 2.0001 |
| SA bump cpu L2 | 6.1282351e-05 | 7.8593215e-06 | 0.49282 | 2.0002 |
| SA bump metal L0 | 2.6721403e-05 | 5.78428e-05 | 0.44044 | 1.9826 |
| SA bump metal L1 | 4.1242816e-05 | 2.2515254e-05 | 0.64986 | 2 |
| SA bump metal L2 | 6.1296664e-05 | 1.0335044e-05 | 0.49293 | 2.0002 |

The fine CPU SA bump already meets continuity; nu-tilde alone prevents its
convergence. Fine Metal still misses both components. The coarser cases have
superposed flow fluctuations, while the transported scalar continues its
slow approach. No stationary scalar storage stall is demonstrated.

All 17 failing solves now have a per-step temporal classification: 13 slow
approaches, three plateaus and one near-two-cycle. This does not establish
their causes or exclude arithmetic contributions. No positive causal
experiment has run, and no convergence fix is accepted.

The unchanged canonical comparator confirms all 153 TMR assertions retain
their exact identities, outcomes and expanded values, with zero regressions,
recoveries, changed values or missing/new identities. All 17 convergence
failures remain. The complete 37-case slack comparison is reported above.

The heaviest completed TMR case peaks at 392,692,648 bytes with field-motion
diagnostics enabled, below 3 GB. These are fresh test-case processes; the
previous 5.54 GB figure measured the monolithic suite. This observation does
not claim that monolithic-process retention has been fixed. The byte-identical
36 B/cell SST cleanup remains the only retained solver memory change.

### Exact fine-bump reproduction

The finest SST bump on Metal completes 3,000 steps on exact base `e014e36`.
Its printed physical endpoint matches the preceding report: k defect
`2.57270639969e-5`, omega `8.41965737184e-7`, continuity
`1.08788010059e-5`, load change `0.000803104893041`, and aft Cf
`0.00255388300098`. It remains unconverged.

The k defect decreases from `2.61183114633e-4` at step 1,000 to
`5.95866563549e-5` at 2,200 and `2.57270639969e-5` at 3,000. Its median
per-step contraction inferred from the last ten 100-step intervals is
0.998928. This component is a slow approach. Its cause has not yet been
isolated. It must not be grouped with the plate continuity floors solely
because both fail the same convergence assertion.

This isolated complete solve peaks at 350,667,544 bytes of process footprint,
including shared Metal memory, below 3 GB. The preceding report's
5,542,384,392-byte figure measured the entire slow-suite process. It did not
isolate one TMR solve. No reduction of that full-suite peak is claimed here.

### Reproducible diagnostics

`bench_convergence` now exposes the compensation, slack, scalar stopping and
individual absolute-tolerance controls independently. It prints their actual
values, cumulative inner-iteration work and process/device bytes. Optional
`--motion-start` samples lag-one and lag-two field movement on every step
without inserting extra monitor or force calls. Scalar samples are their
actual storage; flow samples are high words and cannot prove the full pair
is stationary.

Against the preceding benchmark, coarse-plate CPU and Metal runs preserve
every printed physical value at 40 steps with motion sampling enabled for
the last two steps. This checks the diagnostic path; it does not earn a
solver default. Portability lint passes.

The following command selects the retained SST TMR controls explicitly;
changing one option makes the corresponding controlled probe. It retains
the 3,000-step limit and removes the benchmark's exploratory wall-time stop.

```sh
build/bench_convergence --case bump --model sst --backend metal --level 2 \
  --iterations 3000 --max-seconds 1000000000 --dt .1 \
  --momentum-relaxation .7 --pressure-relaxation .3 \
  --compensated-momentum 1 --paired-state-slack 16 \
  --absolute-tolerance 1e-20 --scalar-compensation 1 \
  --scalar-pseudo-time 0 --scalar-steady-tolerance 1e-6 \
  --motion-start 2700
```

The full per-option comparisons, cost measurements, classification and
isolating experiments for all 17 remaining failures, accepted fixes and
final verification remain unfinished. CUDA runtime remains unavailable.
