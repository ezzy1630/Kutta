# Elbow turbulence discretization

An opt-in reduces a measured standard-k-epsilon production deficit. The unchanged
default in this elbow uses plain cell-LSQ velocity gradients for production.
Reusing SST's derivative reconstruction (cell LSQ with wall
rows omitted, plus half the face-jump fit) reduces the cold elbow's k error
from 9.10% to 4.15% and epsilon from 17.61% to 10.37%. Every measured field
and pressure drop in that configuration improves. Select `RealizableKeSettings::standard_strain_reconstruction =
KeStrainReconstruction::JumpCorrected`. It applies to standard k-epsilon with
LSQ gradients; realizable k-epsilon and Green-Gauss retain their current fits.

V6 SST is **not fixed**. Its cold k error remains 14.48%. Holding the flow
at Fluent's saved state reduces it to 6.50%; that controlled turbulence solve
still disagrees with Fluent.
The wall-face omega offset contributes only 0.003 percentage points to k.
No SST formula or coefficient was changed.

## Scope and evidence

Measured 2026-11-01–02. Base: `8866cc0e16df5b0642b2243611099ddccd7f3cd3`.
Apple M4 Pro, Release,
float fields, double reductions. CPU and Metal run; CUDA runtime is unverified.
Solves ran one at a time in this lane; the full-tier memory overrun is recorded
below. Wall times are not
performance evidence on this shared machine.

Fixtures are the exact corpus `mixing_elbow/mixing_elbow` and
`mixing_elbow/elbow_param`, read without modification. Case/data SHA-256:

| pair | case | data |
| --- | --- | --- |
| V6 | `8145813f50f989a521ed8d1b02c891c7ca22ccc5fafa35a4eeba011949f19cbc` | `863a813227cf4274a9ca24fdad43e6c1f9f132dd62fcabc39a3be8cf493cff42` |
| standard k-epsilon | `d5795b18a7b8b3f17701609486cd1f288d9e2badb9d33b571f642b68afed8a96` | `b779bf1c946b33e7ee469bb05d4f52a8f7628e5a921a921c6c5737ec96b3372c` |

V6 uses `--approximate wfgc`; k-epsilon uses
`--approximate compressibility,ke-options`. These remain approximations.
The new derivative reconstruction does not establish Fluent's undocumented
`ke-source-mod` or `ke-time-marching?` semantics. WFGC is off in elbow_param.

`bench_elbow_turbulence` evaluates the real generated model at the saved
fields and exports cells in original mesh order, sources, sinks, strain,
cell diffusivities, wall diagnostics and regional comparisons. The source
change is confined to the gradients supplied to standard-k-epsilon production.
Constitutive coefficients, wall treatment, transport schemes, iteration
budgets, reference values and test thresholds are unchanged.

## Cold-start comparison

Metal, 300 iterations from each case's own initialization, cell-gradient
default versus the jump-corrected opt-in. V6 has no setting change.
Unweighted relative L2, percent. Wall means a cell owning a wall face;
interior means every other cell. The k-epsilon mesh has 5,584 wall and 40,053
interior cells; V6 has 2,436 wall and 15,386 interior cells. Pressure drop is area-weighted inlet minus
outlet pressure, compared as an absolute relative error. Velocity L2 combines
all three components. Energy is decoupled and not solved here.

| case | quantity | global before | global after | wall before | wall after | interior before | interior after |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| k-epsilon | velocity | 0.880755 | 0.869458 | 1.241857 | 1.238156 | 0.818342 | 0.805261 |
| | pressure | 5.820325 | 5.721179 | 7.768666 | 7.728186 | 5.128628 | 5.001067 |
| | k | 9.103990 | 4.148658 | 4.115781 | 1.420887 | 9.889197 | 4.543926 |
| | epsilon | 17.605960 | 10.371950 | 4.740352 | 1.686338 | 22.252059 | 13.230213 |
| | mu_t | 2.500773 | 2.450293 | 2.243617 | 0.819280 | 2.502888 | 2.459713 |
| | pressure-drop error | 1.011669 | 0.589807 | — | — | — | — |
| V6 SST | velocity | 3.361590 | unchanged | 8.923693 | unchanged | 3.037647 | unchanged |
| | pressure | 2.650766 | unchanged | 2.799907 | unchanged | 2.624367 | unchanged |
| | k | 14.476602 | unchanged | 12.203492 | unchanged | 15.000551 | unchanged |
| | omega | 2.774865 | unchanged | 1.524458 | unchanged | 8.445638 | unchanged |
| | mu_t | 10.448965 | unchanged | 11.212807 | unchanged | 10.448696 | unchanged |
| | pressure-drop error | 1.965305 | unchanged | — | — | — | — |

The fresh base differs from the older case-coverage report: its k-epsilon
velocity/pressure/drop errors are 0.88/5.82/1.01%, rather than
0.93/6.37/0.91%. The table uses freshly measured `8866cc0`, not those older
numbers. V6's slow test stops at 200 seeded / 250 hybrid iterations; the
fixed-300 diagnostic above is separate from its unchanged stopping criteria.

Cold k-epsilon transport defects after 300 iterations are about 3e-8;
momentum is below 7e-9 on both versions. V6's fixed-300 k/omega defects
are below 2e-7. Neither numerical convergence nor improved errors establishes
Fluent equivalence.

## Saved-state isolation

At Fluent's saved k-epsilon fields, mu_t relative L2 is 7.31e-8. Its
constitutive relation is already correct. Wall-cell epsilon from the wall
function has relative L2 7.02e-5 (0.00702%). All 5,584 wall cells have one
wall face; y-star ranges from 56.62 to 501.25. Thus neither a multi-wall
selection rule nor the laminar/log switch at 11.225 is active.

At V6's saved state, the current base's mu_t L2 is **1.0326%**, not the
3.4% in the initial brief. Wall-cell mu_t agrees to 5.25e-8; all of this
mismatch is interior. Fixed wall-cell omega relative L2 is 0.01237%.
The largest saved-state production/cap ratio is 0.3781, with zero active
cells. Changing the production limiter cannot repair this saved-state error.

Saved-state closure errors, percent relative L2 (same on/off for the opt-in):

| quantity | k-epsilon global | wall | interior | V6 global | wall | interior |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| mu_t | 7.30602e-6 | 7.30783e-6 | 7.30601e-6 | 1.03263 | 5.25268e-6 | 1.03280 |
| wall shear magnitude | — | 0.044079 | — | — | 0.007226 | — |
| y-star / saved SV_WALL_YPLUS | — | 0.041324 | — | — | 0.006739 | — |
| wall-cell epsilon / omega target | — | 0.007017 | — | — | 0.012370 | — |

Wall shear is compared with saved force magnitude divided by face area, using
the wall cell as identity (one wall face per cell on both meshes). The saved
`SV_WALL_YPLUS` is the u-star quantity here; `SV_WALL_YPLUS_UTAU` is different
and was not substituted for it.

Fluent's saved states have finite convergence. Its last recorded residual
rows give the following context; the physical normalization is the same
volume-integrated destruction used below. These are recorded iteration
residuals, not a saved term-by-term oracle at the imported fields.

| case / equation | Fluent recorded scaled residual | recorded sum abs(R) / sum Y |
| --- | ---: | ---: |
| k-epsilon k | 2.36294e-4 | 0.012998 |
| k-epsilon epsilon | 1.73502e-4 | 0.003350 |
| V6 k | 5.16072e-5 | 0.000522 |
| V6 omega | 2.10807e-5 | 0.000228 |

Kutta's assembled defects are substantially larger, but finite Fluent
convergence remains part of the comparison uncertainty, especially for k in
the standard elbow. No reference field or residual was adjusted.

The following are actual assembled Kutta residuals at the same saved fields,
including Kutta's face mass flux and walk diffusion. For these diagnostics,
temporary exports sampled the assembled CSR residual after monitor assembly;
they were removed before committing. Rows are mapped back through the solver's
cell permutation. Define R = convection + diffusion - production + destruction,
and Y as volume-integrated physical destruction. Wall epsilon/omega identity
rows are excluded. For SST omega, Y is rho*beta*omega^2, not the doubled Newton
diagonal. Values are ratios, not percent.

| model / equation / region | sum R / sum Y before | after | sum abs(R) / sum Y before | after |
| --- | ---: | ---: | ---: | ---: |
| k-epsilon k, all | +0.062805 | -0.003040 | 0.312935 | 0.293718 |
| k-epsilon k, wall | -0.009732 | unchanged | 0.243373 | unchanged |
| k-epsilon k, interior | +0.084277 | -0.001058 | 0.333526 | 0.308621 |
| k-epsilon epsilon, interior | +0.125049 | +0.023435 | 0.281773 | 0.240366 |
| SST k, all | -0.006280 | unchanged | 0.079018 | unchanged |
| SST k, wall | -0.005910 | unchanged | 0.021719 | unchanged |
| SST k, interior | -0.006628 | unchanged | 0.132763 | unchanged |
| SST omega, interior | -0.024835 | unchanged | 0.222462 | unchanged |

At these saved k-epsilon fields the reconstruction changes production alone:
mu_t, wall production, destruction and transport coefficients do not change.
The interior k net sink excess falls from 8.43% of destruction to -0.11%;
epsilon from 12.50% to 2.34%. This isolates an interior production deficit in
the former discretization. It does **not** identify Fluent's exact velocity
gradient algorithm or explain all cellwise residuals. The near-zero corrected
net k balance is also within the scale of Fluent's recorded finite convergence. The saved data cannot
distinguish a matching gradient reconstruction from compensation for an
undocumented source switch. That uncertainty is another reason to keep the
change opt-in.

An independent transport budget, `tools/forensics/ke_budget.py`, uses Fluent's
saved face flux to remove Kutta's momentum-flux disagreement. It reads the
actual generated sources, sinks and cell diffusivities from the bench.
Convection is the case's first-order upwind; diffusion uses a compact
nonorthogonal split rather than Kutta's walk fit. It is consequently a
separate diagnostic, not a replacement for the CSR residual above.

| k region / variant | C/Y | D/Y | P/Y | Y/Y | R/Y | abs(R)/Y |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| all, default | -0.331877 | -0.000612 | 0.604659 | 1 | +0.062851 | 0.135808 |
| wall, both | +0.244650 | +0.010560 | 1.259122 | 1 | -0.003912 | 0.006610 |
| interior, default | -0.502537 | -0.003919 | 0.410930 | 1 | +0.082614 | 0.174052 |
| all, opt-in | -0.331877 | -0.000612 | 0.670504 | 1 | -0.002993 | 0.101213 |
| interior, opt-in | -0.502537 | -0.003919 | 0.496265 | 1 | -0.002722 | 0.129216 |

For the interior epsilon equation, the independent budget is:

| variant | C/Y | D/Y | P/Y | Y/Y | R/Y | abs(R)/Y |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| default | -0.440269 | -0.006740 | 0.431045 | 1 | +0.121945 | 0.194323 |
| opt-in | -0.440269 | -0.006740 | 0.532659 | 1 | +0.020331 | 0.135002 |

The wall k budget nearly closes with the current wall production and
dissipation. The interior imbalance persists with Fluent's flux. Half the
face-jump correction reduces this independent interior abs(R)/Y to 0.129216
using the actual exported kernel source. The independent Python strain sweep
below gives 0.129212.

The documented equations are Fluent Theory Guide 2024 R2,
[standard k-epsilon, section 4.3.1](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/flu_th_sec_turb_ke.html),
and [wall treatment, section 4.18.2.1.4, equations 4-349–351](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/ch04s18s02.html):
zero wall-normal k gradient, log-law production, and algebraic wall-cell
epsilon. The new gradient reconstruction is an inference supported by the
isolating experiment, not an additional formula attributed to the guide.

## Rejected and non-fixing experiments

| experiment | measured result | disposition |
| --- | --- | --- |
| k-epsilon: add scalar symmetry normal rows | all cold-300 fields and pressure drop bit-identical | no change retained for this experiment |
| V6: divide wall-face omega reconstruction by 1.004187, holding the wall-cell target unchanged | k 14.476602% → 14.473258%; omega 2.774865% → 2.775797%; drop 1.965305% → 1.964580% | diagnostic sensitivity only; reverted, not a fitted correction |
| V6: freeze velocity, pressure and Kutta's initially reconstructed flux at Fluent's saved fields; solve k/omega for 600 steps | k global/wall/interior 6.4981/3.4437/7.0681%; omega 2.2445/0.0955/8.0346%; global mu_t 9.2128%; k/omega defects below 3.3e-8 | isolates flow feedback; not a cold-flow accuracy result |
| k-epsilon: use saved wall-face epsilon (median face/owner ratio 2.000000) as a normal-foot or face-centroid gradient observation | independent interior epsilon abs(R)/Y 0.135002 → 0.138685 / 0.138561; normal-foot net R/Y +0.020331 → +0.024992 | worse; no solver change retained |
| k-epsilon: raise the relative LSQ eigenvalue cutoff from 1e-10 through 1e-2 in the independent epsilon gradient | interior epsilon abs(R)/Y remains 0.135002 to printed precision | no gain in this diagnostic; shared LSQ untouched |
| k-epsilon: plain LSQ weighting 1/d^0 instead of 1/d^2 | independent interior abs(k residual)/Y 0.174627 vs 0.173906 using the reconstructed plain LSQ2 rate | no material gain; no solve retained |
| k-epsilon: add full Boussinesq mean-divergence terms using the numerical LSQ trace (not Sarkar dilatation dissipation) | independent interior abs(k residual)/Y: plain 0.173906 → 0.263519; half-jump 0.129212 → 0.270045; 9,484 / 10,039 interior cells have negative production | no improvement; no solver change retained |
| k-epsilon: full instead of half face-jump correction | independent interior abs(k residual)/Y 0.189720 vs 0.129212 with half; net R/Y -0.130994 vs -0.002713 | overproduction; rejected before coupled solve |

The independent gradient sweep used the same saved fields and scalar
transport terms. Interior abs(k residual)/Y, with the strain recomputed by the
Python operator (rather than imported from the generated kernel):

| LSQ weight | no jump | half jump | full jump |
| --- | ---: | ---: | ---: |
| 1 | 0.174627 | 0.131670 | 0.200794 |
| 1/d | 0.173510 | 0.129812 | 0.193695 |
| 1/d^2 | 0.173906 | 0.129212 | 0.189720 |

The V6 frozen-flow run retains Kutta's reconstructed saved-state mass flux,
not Fluent's saved SV_FLUX. Its change in L2 is a controlled sensitivity,
not an additive decomposition of the 14.48% error. No result from another
momentum lane was merged or assumed.

## Additional vendor-reported omega regularization

[Menter et al. (2022), equation 6, author-provided full text](https://www.researchgate.net/publication/362198958_An_Algebraic_LCTM_Model_for_Laminar-Turbulent_Transition_Prediction)
reports Fluent omega production as `gamma*rho*Pk/(mu_t + 0.001*mu)`.
This is a vendor-authored implementation statement, distinct from the printed
Theory Guide formula and not independently established for this saved 22.2
case. Kutta currently divides by mu_t without that added molecular term.

At the saved V6 fields, with fixed omega wall rows excluded, that denominator
would reduce the integrated production by 0.1641%, or 0.122995% of physical
omega destruction. The largest local production reduction is 86.68% in a
very-low-mu_t cell; the minimum interior mu_t/mu is 0.0001537. Applying only
this source difference to the exported assembled residual changes net R/Y
from -0.024835 to -0.023605, and abs(R)/Y from 0.222462 to **0.222553**.
Thus it does not close the aggregate saved-state defect. No source change is
retained here, and its cold-solve effect is unmeasured. This remains a concrete
constitutive difference to check, particularly for transition cases, rather
than an explanation established for the elbow's k error.

## Pressure-scheme trade-off

The second-order-pressure regression guard uses the same Fluent saved
**standard-pressure** reference; it is not an independently converged Fluent
second-order solution. At 400 iterations, CPU measurements are:

| quantity, relative error % | cell-gradient default | derivative opt-in |
| --- | ---: | ---: |
| velocity L2 | 1.052814 | 1.044297 |
| pressure L2 | 5.987973 | 5.929358 |
| k L2 | 8.712619 | 4.210473 |
| epsilon L2 | 16.935521 | 10.073011 |
| pressure drop | 1.113939 | 1.518024 |
| total wall force | 0.407685 | 0.406903 |

Both variants pass every original guard, including the 1.7% pressure-drop
bound. Nevertheless, the pressure-drop error grows by 0.404 percentage points.
The implementation therefore remains opt-in. There is no default selection
from Fluent's undocumented k-epsilon switches.

## Full verification

The numerical solver source is `24e6c5ad2afadcd74737f5459bd93064ebb56471`.
Later commits through `3984a6dcc448165b4eacb9680bf1042c06cbf36d` only add
diagnostic provenance and strict iteration parsing. The full-tier test binary
SHA-256 is `f502894e6c16315df2d4d7c1424bf7523478a365738c722a0c08eb0403ec9b38`;
it was unchanged throughout the run.

| check | observed result |
| --- | --- |
| full quick tier, default settings | 277/277 selected cases, 7,457,592 assertions passed |
| focused k-epsilon/nozzle and second-order-pressure slow cases, derivative candidate | 2/2 selected cases, 28,264 assertions passed on CPU and Metal |
| final explicit opt-in, cold Metal 300 | exported fields byte-identical to the measured derivative candidate |
| final V6, cold Metal 300 | exported fields byte-identical to the original default |
| portability lint | passed |
| full slow tier, default settings | 37 selected cases: 31 passed, 6 failed; 201,571 assertions passed, 62 failed |
| assertion comparer metadata controls | 9 Python tests passed |

The focused run preceded the explicit opt-in setting; final saved-state and
cold-field exports verify that the explicit selection retains that candidate's
numerics. The full slow tier exercises the unchanged defaults. CUDA runtime
is unverified.

The supplied completed main `8866cc0` baseline and this branch each contain
201,633 assertions. The full comparison has **0 regressions, 0 recoveries,
0 changed expanded values, 0 missing/new assertions and 0 missing/new cases**.
The unchanged failure counts are:

| slow case | main failures | branch failures |
| --- | ---: | ---: |
| V6 mixing elbow | 25 | 25 |
| V2 manufactured Navier–Stokes tet | 8 | 8 |
| V5 SST flat plate | 6 | 6 |
| V5 SST bump | 8 | 8 |
| V7 SA flat plate and bump | 14 | 14 |
| V12 native polynomial fan | 1 | 1 |
| all other selected cases | 0 | 0 |
| total | 62 | 62 |

The first comparison returned exit 1 with 266 missing/new identities: 68 in
Fluent/Kutta CLI parity and 198 in native-FSAE CLI parity. INFO included bare
process addresses, random `MeshTemp` directory ids and CLI wall times. Commit
`bf6902753137010898e58a5884ef84644ce0d3fa` fixes that diagnostic metadata
normalization. It preserves model labels, physical residual columns, iteration
counts, expressions, expanded values and outcomes. Its controls verify that
changes to those physical values or outcomes remain visible. Rerunning the
same comparer on the original, unmodified XML files returns exit 0 and matches
all 201,633 assertions. Both raw and corrected comparison JSON are retained
locally. No CFD test, threshold, reference or iteration limit changed.

Selected-case totals include two coverage limits present in both runs:

- Bare-wing starter parity checked the donor case hash, then returned because
  `runs/wing/naca4412_ground.msh` and `naca4412_ground.kutta.json` are absent.
  That parity solve is unverified.
- The corpus-FSAE case reports 21 existing setup holds, zero equivalences and
  zero assertions. Its hold output is identical between runs. These are held
  cases, not accepted native-physics solves. The separate synthetic native-FSAE
  CLI parity case did execute.

The full-suite process reported maximum RSS **2,415,460,352 bytes (2.42 GB)**
and macOS peak memory footprint **5,763,240,784 bytes (5.76 GB)**. The reported
footprint exceeds the requested approximately 3 GB limit. RSS-only monitoring
missed that overrun during execution. These are process metrics, not a separate
GPU-allocation measurement; the responsible case/cache contribution was not
profiled. No concurrent heavy solve was launched by this lane.

## Remaining discrepancies, ranked by measured contribution

Source ratios and field sensitivities are different measures; these effects
are not an additive decomposition of the total error.

1. **V6 flow feedback.** Freezing flow lowers global k error by 7.98 percentage
   points (14.48% → 6.50%) and wall k error by 8.76 points (12.20% → 3.44%).
   This is the largest measured effect. It does not isolate which momentum
   operator is responsible; that belongs to the parallel momentum lane.
2. **Remaining V6 turbulence transport/reconstruction.** Even with frozen
   flow, k is 6.50% off and interior omega 8.03% off. The saved-state interior
   omega abs(R)/Y is 0.2225. The frozen-flow run also retains a 9.21% global
   mu_t error and Kutta's initial flux, so it does not isolate source/diffusion
   from that remaining mass-flux difference. Gradient reconstruction, cross diffusion and
   diffusion are not separated by these measurements. Saved-state wall k
   abs(R)/Y is 0.0217. Fluent does not save wall production here, so the wall
   production blend remains unverified despite the matching wall targets,
   shear and y-star. The previous
   `k_wall_cells/REPORT.md` localizes an additional wall/symmetry-corner issue.
3. **Remaining k-epsilon interior transport/production.** After the change,
   cold interior epsilon is 13.23% off, k 4.54%; the saved-state interior
   epsilon net R/Y is still +0.0234 and abs(R)/Y 0.2404. Production improvement
   is isolated, but the remaining terms are not. `ke-source-mod` and
   `ke-time-marching?` remain undocumented and unisolated.
4. **Vendor-reported omega denominator regularization.** Its saved-state
   source contribution is 0.123% of omega destruction, and it slightly
   worsens the aggregate absolute residual. Local effects can be large at very
   low mu_t. Its presence in this case and cold-field effect remain unverified.
5. **V6 wall-face omega offset.** Removing its 0.42% bias moves global k by
   only 0.00334 percentage points and worsens omega by 0.00093 points. It
   cannot account for the target discrepancy.
6. **Inactive or already matched terms.** The V6 saved-state production
   limiter never activates. Standard k-epsilon wall epsilon and saved-state
   mu_t already match closely. The omitted compressibility term's existing
   case-specific bound is 1.7e-4 of epsilon (`case_coverage/REPORT.md`); it was
   not enabled or independently rederived in this lane.

The supplied cell-field datasets contain velocity, turbulence variables,
viscosities and pressure, but no velocity-gradient or turbulence-production
field. A direct comparison of those fields, or a controlled Fluent run with
one k-epsilon switch changed at a time, is needed to distinguish the remaining
reconstruction and source-option hypotheses. This lane does not identify
those undocumented switches from their names.

## Reproduction

```sh
export KUTTA_FLUENT_EXAMPLES=<fixture-checkout>/validation/fluent_examples
export KUTTA_RUNS=<fixture-checkout>/runs
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
build/bench_elbow_turbulence \
  "$KUTTA_FLUENT_EXAMPLES/mixing_elbow/elbow_param.cas.h5" \
  "$KUTTA_FLUENT_EXAMPLES/mixing_elbow/elbow_param.dat.h5" \
  build/ke-saved metal saved 0 compressibility,ke-options jump
python3 -I tools/forensics/ke_budget.py \
  "$KUTTA_FLUENT_EXAMPLES/mixing_elbow" elbow_param build/forensics-cache build/ke-saved
# Replace `saved 0` by `cold 300` for the cold comparison; use `cell` for
# the unchanged default. There is no Fluent rpvar mapping for this opt-in.
# For V6 use stem mixing_elbow and approximation wfgc.
build/kutta-tests --test-suite-exclude=slow
build/kutta-tests -ts=slow -r=xml --success > branch.xml
python3 tools/forensics/compare_slow_assertions.py <baseline>/slow.xml branch.xml
python3 tools/forensics/test_compare_slow_assertions.py
python3 tools/lint_portability.py
```

Raw logs, CSV diagnostics, scratch isolation patch, branch XML and full JSON
comparison remain uncommitted in `build/elbow-turbulence/`. They contain
reference-derived fields. No restricted data, mesh, slow XML or fixture link
is committed. Temporary fixture symlinks are removed at completion.
