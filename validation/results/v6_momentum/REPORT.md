# V6 momentum and face flux

The requested Fluent reproduction is **not achieved**. This branch adds an
explicit momentum-face reconstruction option and saved-state isolation tools.
It does not change the default equations or silently reinterpret a Fluent
setting.

With `FaceVelocityReconstruction::PecletSide` and the existing differentiable
momentum limiter, seeded V6 velocity error falls from 3.362% to 1.388%, pressure-drop
error from 1.965% to 0.106%, and elbow shear error from 1.641% to 0.568%.
The pressure-drop and both total-wall-force gates recover. Velocity, elbow
shear, the seeded all-field gate, and original WFGC compatibility still fail.

This makes the previously scratch-only Péclet experiment available through
`SolverSettings`. It is a partial improvement, not a newly identified exact
Fluent formula. It stays opt-in because k-epsilon elbow and cylinder errors,
and some MMS orders, get worse, including values inside passing gates.

## Identity and scope

- Numerical opt-in: `2ba69af`; probe identity: `e450904`; probe input and
  convergence guards: `f1e431b` and `4f017b8`; comparison identity repair:
  `87f9197`.
  These are the measured pre-integration commits. On main they are `a3fda32`,
  `d310f03`, `f165b48` and `d04b452`; the comparison identity repair was
  superseded by main's equivalent `971fbb3`, which normalizes the same CLI
  timing and scratch-path context for every source.
- Base: `8866cc0`, branch `t3/v6-momentum-discretization`, Release AppleClang,
  Apple M4 Pro, CPU and Metal. Solve timings are not used as evidence on the
  shared machine. Each solve runs alone in this lane.
- Fluent 22.2, `mixing_elbow`, 17,822 cells and 91,581 faces. Case SHA-256:
  `8145813f50f989a521ed8d1b02c891c7ca22ccc5fafa35a4eeba011949f19cbc`;
  data SHA-256:
  `863a813227cf4274a9ca24fdad43e6c1f9f132dd62fcabc39a3be8cf493cff42`.
- Both starts use the unchanged V6 stopping rule and 4,000-iteration ceiling.
  The cold start uses the stored hybrid initialization. All results explicitly
  permit the existing `wfgc` approximation; none represents the complete
  original Fluent settings.
- Device fields retain the base solver's float storage. Frozen saved values
  use that storage precision; the Python budgets use the saved doubles.
- No turbulence-model equation, constant, acceptance threshold, reference
  value, or iteration budget changes. There is no case-name/model/mesh switch.
- The opt-in uses device gather kernels generated for CPU, Metal and CUDA.
  CUDA compilation and runtime were not available or tested here.

## Retained change

For each internal face, the option uses

```text
outP = sum_faces max(outward_mass_flux, 0)
beta = min(1, (outP + outN) / (aP + aN))
side = w (uP + GP.rP) + (1-w) (uN + GN.rN)
u_face = existing_point_or_projected_value + beta (side - existing_value)
```

`aP` is the physical momentum diagonal, without pseudo-time or relaxation;
`w` follows the selected momentum/distance weights. The side gradients retain
the selected limiter. Boundary fluxes are unchanged. The option allocates one
float per cell and adds one outflow gather plus three scalar face kernels.
The default allocates none of that storage and dispatches none of those kernels.

At fixed Reynolds number in the viscous limit, beta is O(h). The O(h²)
side-minus-point term therefore becomes O(h³). This is Kutta's consistency
construction, not a Fluent-documented coefficient. The previous unweighted
side mean's hex wall-layer order loss is discussed in
[`fluent_face_velocity/REPORT.md`](../fluent_face_velocity/REPORT.md).

Caller usage:

```cpp
settings.face_velocity_reconstruction = fv::FaceVelocityReconstruction::PecletSide;
settings.limiter = fv::Limiter::Differentiable;
```

There is no automatic mapping to `pressure/ho-rc-scheme 2`. The setting's
meaning is undocumented and the corpus does not establish that it selects a
limiter. The default remains Point and Standard.

## Every V6 gate, before and after

Numbers below are CPU; Metal differs by less than 0.0001 percentage points
on the listed field/force/drop errors. S = saved-field start, C = cold hybrid.
Percentages are relative errors, not relative improvements. The new option is
combined with the existing differentiable limiter in the after columns.

| Gate | Bound | Before S / C | Opt-in S / C | After outcome |
| --- | --- | ---: | ---: | --- |
| Original settings accepted | no rejected active settings | WFGC rejected | WFGC rejected | fail |
| Stored initialization accepted | no rejection | pass | pass | pass |
| Initial mu_t L2 | <1.1% | 1.0326% / n/a | 1.0326% / n/a | pass |
| Saved-state continuity L1 | <8e-4 kg/s | 2.06914e-3 / n/a | 1.16938e-3 / n/a | fail |
| All steady defects and continuity | <1e-6; <2.475884e-6 kg/s | pass / pass | pass / pass | pass |
| Cell velocity L2 | <0.5% | 3.36162 / 3.36165 | 1.38844 / 1.38858 | fail |
| Seeded max(u,p,k,omega) L2 | <0.1% | 14.4782 / n/a | 8.11728 / n/a | fail |
| Elbow integrated shear | <0.5% | 1.64074 / 1.63789 | 0.56764 / 0.56739 | fail |
| Inlet integrated shear | <0.5% | 0.15238 / 0.15067 | 0.12198 / 0.12196 | pass |
| Elbow total wall force | <0.5% | 0.60239 / 0.60006 | 0.07961 / 0.07788 | recover |
| Inlet total wall force | <0.5% | 0.53384 / 0.52918 | 0.21057 / 0.21162 | recover |
| Inlet-to-outlet pressure drop | <0.5% | 1.96534 / 1.96604 | 0.10618 / 0.10823 | recover |
| Finite defects; maximum speed | finite; <100 m/s | pass / pass | pass / pass | pass |

The test's field/zone-size and source-identity checks also pass. The strict V6
test retains 25 failed assertions by default and 13 with the option, counting
both starts and both backends: 12 recoveries, no new failures.

| Additional metric | Before S / C | Opt-in S / C |
| --- | ---: | ---: |
| Iterations, both backends | 200 / 250 | 150 / 200 |
| Pressure L2 (%) | 2.65576 / 2.65111 | 0.93596 / 0.93701 |
| k L2 (%) | 14.4782 / 14.4777 | 8.11728 / 8.11832 |
| omega L2 (%) | 2.77599 / 2.77487 | 2.29492 / 2.29490 |
| Elbow pressure force error (%) | 0.65793 / 0.65558 | 0.05967 / 0.05793 |
| Inlet pressure force error (%) | 0.53392 / 0.52925 | 0.21009 / 0.21115 |

## What the saved fields isolate

`tools/forensics/v6_momentum_budget.py` holds all saved fields fixed. It uses
Fluent's saved conservative flux, viscosity and wall shear, and builds the
LSQ/SOU momentum imbalance. The Python assembly omits Kutta's boundary
extrapolation-gain guard; its residuals are operator probes, not exact
`SteadySolver::monitor()` results or claims of Fluent normalization parity.

The flux CSV from `bench_v6_operators` identifies a face by its input node IDs
and both adjacent cell IDs. The Python reader verifies every face once,
including its orientation. It verifies the case/data SHA-256 headers and never merges faces by cell pair or estimates
identity by nearest-neighbour coordinates. Geometry caches are keyed by case
SHA-256.

| Substitution at identical saved fields | Sum of absolute x/y/z momentum imbalances (N) |
| --- | ---: |
| Fluent's own final raw residual sum | 1.06141e-5 |
| Saved F, mu_t and wall shear; Kutta-like spatial operator | 1.02092e-3 |
| Replace only mu_t by Kutta's saved-state reconstruction | 1.02650e-3 |
| Replace only F by Kutta's default constitutive face flux | 2.81666e-3 |
| Replace only F by the opt-in constitutive face flux | 1.85496e-3 |

The last row keeps the budget's Standard-limited momentum transport unchanged;
it substitutes only the flux from the opt-in plus differentiable-limiter run.
These contributions are not additive error estimates.

Using the saved conservative flux reduces this default constitutive-flux
probe's imbalance by 63.8%. A large residual remains even with saved flux,
viscosity and wall shear: about 96 times Fluent's reported raw residual sum.
The saved-flux and reconstructed-flux defects are therefore distinct targets.

Direct `SteadySolver` face-flux measurements at the saved fields:

| Metric | Default | Opt-in plus differentiable limiter |
| --- | ---: | ---: |
| Relative L2 vs SV_FLUX, all 86,782 internal faces | 0.96294% | 0.33831% |
| Median absolute face error / absolute saved flux | 8.37279e-3 | 1.35629e-3 |
| Sum absolute divergence (kg/s) | 2.06914e-3 | 1.16938e-3 |

The native solver's saved-state scaled momentum defects (x/y/z) are
1.12783e-3 / 1.98653e-3 / 5.25224e-4 by default and
8.00413e-4 / 1.18141e-3 / 4.63558e-4 with the option and differentiable
limiter. These are direct monitor values, distinct from the dimensional
Python budget above. Neither formulation makes the saved fields a fixed
point of Kutta's momentum operator.

The converged retained-runtime probe puts 35.6% of squared velocity error in
the 2,436 wall-adjacent cells, whose local velocity L2 error is 4.8275%.
The three interior skew groups have local errors 1.3071%, 1.1269% and
1.0393%, respectively. This locates the remaining error; it does not identify
a particular wall or gradient formula as its cause.

The saved conservative flux itself has divergence 7.65770e-8 kg/s. That is not
Fluent's final constitutive continuity residual, 2.475884e-6 kg/s.

## Rejected trials and limits of attribution

**Isotropic stress.** Adding the explicit gradient of `2 rho k/3` to the
saved pressure gives imbalance 2.33891e-3 N; the opposite sign gives
2.35424e-3 N, against 1.02092e-3 N without it. The unconstrained least-squares
coefficient of that term is 0.00645, versus the proposed coefficient 1.
This rejects adding that term directly to this stored-pressure formulation.
It does not establish how Fluent defines or stores pressure in every mode.
No pressure-convention change is retained.

**Viscosity interpolation.** Arithmetic gives 1.03713e-3 N and harmonic
9.93331e-4 N, against distance-weighted 1.02092e-3 N. The 2.7% local residual
reduction does not close the mismatch. Neither alternative was adopted or
claimed to improve the converged solution.

**Two geometric hypotheses, not Fluent algorithms.** Integrating the linear
pressure over the actual triangle-fan area vectors, rather than using one
flattened face, gives 1.00179e-3 N: a 1.9% reduction. Integrating the covariance
of affine mass-flux and transported-velocity variations gives 1.08591e-3 N:
a 6.4% increase. Neither collapses the saved-state floor; neither is retained.
These are not implementations or isolations of Fluent WFGC.

**Green-Gauss pressure gradient in Rhie-Chow.** With the Péclet option and
differentiable limiter, replacing the RC reconstruction gradient by the
momentum pressure-force gradient gives, seeded Metal: velocity 1.37754%,
pressure 1.27933%, drop 0.53991%, inlet total force 0.89452%. Cold gives
velocity 1.37721%, pressure 1.27792%, drop 0.53842%, and inlet force 0.89357%. Both converge in 200 iterations.
It loses the recovered drop and inlet-force gates, so it was removed.

**Component momentum diagonals.** A separate offline face-flux proxy rebuilt
`aP_i` from saved convection, viscosity, wall shear divided by tangential
wall-cell speed, and the normal symmetry coefficient. With pure side
reconstruction, the differentiable limiter and production-like LSQ stencils,
using each component's diagonal changes flux relative L2 from 0.33408% to
0.33109%. Projecting the component inverses onto the face normal for RC instead
gives 0.33834%; both changes give 0.33546%. Across both boundary-row policies,
both limiters and LSQ/momentum-GG pressure jumps, component velocity weighting
removes at most 1.31% of the flux L2 error; the normal-mobility variant never
improves it. These are coefficient proxies, not native-solver runs or an
identification of Fluent's mobility. Neither change was adopted.

**Frozen turbulence, flow solved.** A scratch isolation held k and omega at
Fluent's saved values, skipped their transport solves and supplied Fluent's
exact saved mu_t to momentum viscosity. The wall closure still used the
current velocity and the fixed saved k/omega. These are momentum-only fixed
points, not accepted coupled SST solves:

| Momentum formulation, Metal seeded | Iterations | Velocity L2 | Pressure L2 | Drop error |
| --- | ---: | ---: | ---: | ---: |
| Default, saved k/omega/mu_t | 150 | 3.43929% | 3.22983% | 2.62058% |
| Péclet + differentiable, saved k/omega/mu_t | 100 | 1.34621% | 0.86098% | 0.12890% |

The scratch diagnostic used the same numeric flow thresholds but omitted the
deliberately frozen turbulence residuals from its stopping check. It is not a
V6 acceptance run. No such change remains in tests. All freeze and RC
experiment hooks were removed before the final build. The retained option
has no environment-variable hook in the solver.

The current base already has saved-state mu_t L2 1.0326%, not the historical
3.4%, owing to the retained half-jump strain operator. Supplying exact mu_t and
holding both turbulence quantities does not recover the velocity gate, even
with the improved face velocity. The saved cell mu_t mismatch therefore does
not explain the whole error in this formulation. Wall-closure differences
remain active in this experiment; it does not eliminate them. The
approximately 0.42% wall-face omega offset was not perturbed in this momentum
lane.

## Other cases and full verification

The quick tier on the clean numerical implementation passed **278 cases and
7,457,602 assertions**, with no missing-fixture skip messages. Its 37 excluded
cases are the slow tier. Peak process RSS was 1.53 GB. The dedicated kernel
case passes 10 assertions across CPU and Metal; portability lint is clean.
The forensic reader also rejected a mismatched input hash and a mismatched
face adjacency in negative CLI probes. It now rejects nonfinite saved inputs
and checks the original wall-shear rows before the existing reader's
non-wall zero-fill. Removing one wall-shear component in memory reproduced
the silent zero-fill and was rejected by this guard; source data was unchanged.
JSON output also rejects nonfinite numbers.

The legacy pair-indexed text flux input had a separate truncation bug:
removing its final row returned success and compared 86,780 faces instead
of 86,782. The probe now checks every row, declared counts, cell ranges,
finite values and source-mesh adjacency coverage. The truncated input and
a wrong-adjacency input are rejected before solver construction; the complete
input still compares all 86,782 faces at unchanged 0.96294% flux L2 error.
The hash-bound node-identity CSV remains the budget reader's input format.

The option was forced only for the focused comparison below, then the hook
was removed and the committed tree rebuilt. Against the completed shared
baseline, the canonical comparator filtered to these five cases reports:
28,505 assertions, eight known V2 tet failures, zero added failures, zero
recoveries, zero missing/new assertions or cases, and 72 changed expanded
values. The changes comprise 24 elbow/nozzle assertions, 16 second-order
pressure elbow assertions, 12 cylinder assertions and 20 MMS order assertions.
They include worse metrics inside passing gates and worse already-failing
tet orders; passing gate counts do not justify promoting the option.

The MMS orders below are the finest-grid orders, equal on CPU and Metal to
the shown precision. All hex gates pass. Tet velocity and three wall-force
gates remain failing; its pressure order retains the existing operator-order
comparison.

| MMS metric | Baseline order | Opt-in order | Existing bound |
| --- | ---: | ---: | --- |
| Hex velocity | 1.91640 | 1.91725 | >=1.9 |
| Hex pressure | 1.82115 | 1.88139 | >=1.8 |
| Hex pressure force | 1.91472 | 1.97990 | >=1.8 |
| Hex shear force | 1.95917 | 1.99665 | >=1.8 |
| Hex total force | 2.01196 | 2.00955 | >=1.8 |
| Tet velocity | 1.81502 | 1.80510 | >=1.9 |
| Tet pressure | 0.967086 | 0.940530 | >=operator order - 0.1 (0.914063 here) |
| Tet pressure force | 1.05200 | 1.07362 | >=1.8 |
| Tet shear force | 0.842617 | 0.852198 | >=1.8 |
| Tet total force | 0.622970 | 0.612269 | >=1.8 |

The four additional changed second-order-elbow assertions are its final
continuity/momentum defects: CPU 9.15578e-7 / 3.40308e-9 to
4.34845e-7 / 2.58047e-9, and Metal 6.80308e-7 / 3.08176e-9 to
4.22521e-7 / 2.68611e-9. All retain their 1e-4 / 1e-6 bounds.

The following CPU values compare this branch's default against the explicit
option. Bounds are the existing tests' unchanged bounds. Every listed gate
passes in both runs, while the k-epsilon fields and several cylinder metrics
get worse. The second-order-pressure elbow keeps its existing comparison
against the standard-pressure Fluent reference; this is a bounded method
variation, not a same-pressure-scheme equivalence claim.

| Case | Error (%) | Default | Opt-in | Bound |
| --- | --- | ---: | ---: | ---: |
| k-epsilon elbow | velocity | 0.88076 | 1.1588 | 1.5 |
| k-epsilon elbow | pressure | 5.8204 | 8.95 | 10 |
| k-epsilon elbow | k | 9.104 | 10.84 | 14 |
| k-epsilon elbow | epsilon | 17.606 | 20.403 | 27 |
| k-epsilon elbow | drop | 1.012 | 1.1193 | 1.5 |
| k-epsilon elbow | wall force | 0.31571 | 0.31007 | 0.6 |
| k-epsilon, second-order pressure | velocity | 1.0528 | 1.332 | 1.6 |
| k-epsilon, second-order pressure | pressure | 5.988 | 8.3137 | 8.7 |
| k-epsilon, second-order pressure | k | 8.7126 | 10.534 | 13.1 |
| k-epsilon, second-order pressure | epsilon | 16.936 | 20.02 | 25.4 |
| k-epsilon, second-order pressure | drop | 1.1139 | 0.22199 | 1.7 |
| k-epsilon, second-order pressure | wall force | 0.40769 | 0.11786 | 0.61 |
| SST nozzle | velocity | 1.6034 | 1.4298 | 2.5 |
| SST nozzle | pressure | 0.527 | 0.16499 | 0.9 |
| SST nozzle | k | 4.9191 | 4.4722 | 8 |
| SST nozzle | omega | 5.5184 | 5.4868 | 8.5 |
| SST nozzle | drop | 0.46207 | 0.070341 | 0.8 |
| SST nozzle | wall force | 0.56554 | 0.155 | 1 |
| cylinder | velocity | 0.064484 | 0.059651 | 0.097 |
| cylinder | pressure | 0.51733 | 0.44106 | 0.78 |
| cylinder | drop | 0.50597 | 0.74997 | 0.76 |
| cylinder | wall force | 0.070559 | 0.075006 | 0.11 |
| cylinder | Cd | 0.069512 | 0.075005 | 0.105 |
| cylinder | Cl / saved Cd | 0.012111 | 0.00030196 | 0.0185 |

The option was not forced through the full slow tier. Its global external-aero,
transition and TMR effects remain unverified here. The complete branch tier
uses the default Point/Standard settings.

The requested full default slow command completed all 37 selected cases:
201,571 passing assertions and the same 62 failing assertions in six cases
(31 passed cases). The canonical comparison reports **zero regressions,
zero recoveries, zero changed expanded values, and zero missing/new
assertions or cases**. The 62 failures remain V6 25, V2 tet 8, V5 plate 6,
V5 bump 8, V7 SA 14 and V12 1.

| Full slow case | Assertions main / branch | Failures main / branch |
| --- | ---: | ---: |
| Kutta cases converted from Fluent cases solve to the same bytes on every backend | 68 / 68 | 0 / 0 |
| The NACA ground-effect wing from a bare mesh and a Kutta case solves like its donor Fluent case | 1 / 1 | 0 / 0 |
| Native FSAE Fluent physics convert and solve to byte-identical Kutta outputs | 477 / 477 | 0 / 0 |
| Corpus FSAE physics preserve strict solver holds or convert with identical outputs | 0 / 0 | 0 / 0 |
| FV V3 Re1000 cavity refinement against Erturk | 16 / 16 | 0 / 0 |
| FV V3 Re100 and Re400 cavity reported deviations | 4 / 4 | 0 / 0 |
| FV V4 Fluent cavity matches its actual grid | 763 / 763 | 0 / 0 |
| Fluent standard k-epsilon elbow and SST nozzle match converged Fluent fields | 17044 / 17044 | 0 / 0 |
| Fluent k-epsilon elbow converges with second-order pressure | 11220 / 11220 | 0 / 0 |
| Fluent laminar 2D cylinder matches converged Fluent fields | 207 / 207 | 0 / 0 |
| FV V6 Fluent mixing elbow fixed point and cold start | 2565 / 2565 | 25 / 25 |
| FSAE aero: Tyrrell 026 sectional Cl in ground effect vs thesis CFL3D k-omega | 28 / 28 | 0 / 0 |
| FSAE aero: Ahmed body 35 deg drag, level 1 | 4 / 4 | 0 / 0 |
| FSAE aero: Selig S1223 at Re 2e5, alpha 4.76 deg vs UIUC | 9 / 9 | 0 / 0 |
| FSAE aero: Ahmed body 35 deg at y+ ~ 1, level 1, SIMPLEC and Coupled | 8 / 8 | 0 / 0 |
| FSAE V12 native Fluent polynomial fan import acceptance | 1 / 1 | 1 / 1 |
| V8 realizable KE scalable channel MKM DNS profile and momentum balance | 18 / 18 | 0 / 0 |
| FV V2 manufactured Navier Stokes hex refinement | 16 / 16 | 0 / 0 |
| FV V2 manufactured Navier Stokes tet refinement | 18 / 18 | 8 / 8 |
| FV V1 manufactured scalar diffusion and upwind refinement | 48 / 48 | 0 / 0 |
| FV SIMPLEC develops planar Couette flow with moving and slip walls | 14 / 14 | 0 / 0 |
| V5 TMR flat plate SST accuracy and three-grid trend | 36 / 36 | 6 / 6 |
| V5 TMR bump SST station accuracy and three-grid trend | 47 / 47 | 8 / 8 |
| V7 SA TMR flat plate and bump accuracy convergence and three-grid trend | 70 / 70 | 14 / 14 |
| Transition: ERCOFTAC T3A skin friction vs experiment (REGRESSION GUARD) and vs the documented equations (VERIFICATION) | 161470 / 161470 | 0 / 0 |
| Transition: Selig S1223 at Re 2e5, alpha 4.76 deg vs UIUC | 11 / 11 | 0 / 0 |
| FSAE V7.1 translating-wall Couette hex and poly | 104 / 104 | 0 / 0 |
| FSAE V7.2 rotating-wall Taylor Couette hex and poly | 6668 / 6668 | 0 / 0 |
| FSAE V8 Darcy Forchheimer pressure drop and oblique directions | 24 / 24 | 0 / 0 |
| FSAE V9 absolute-velocity MRF analytic solutions and frame invariance | 48 / 48 | 0 / 0 |
| FSAE V10 incompressible pressure and mass-flow inlets | 12 / 12 | 0 / 0 |
| FSAE V11 internal fan polynomial and porous jump pressure discontinuity | 312 / 312 | 0 / 0 |
| FSAE V7.3 SST moving-wall and translating-frame invariance | 16 / 16 | 0 / 0 |
| FSAE V9.2 conformal MRF interface preserves absolute Taylor Couette | 4 / 4 | 0 / 0 |
| FSAE V10.2 SST pressure and mass inlets match velocity-inlet transport | 12 / 12 | 0 / 0 |
| FSAE V9.3 rotation and translation preserve Taylor Couette pressure | 6 / 6 | 0 / 0 |
| Fluent research downloads retain published byte identities | 264 / 264 | 0 / 0 |

The test executable was built at `e450904`; its SHA-256 is
`b8e63fcfcb75661bd19d025598095c1c76ba8a2526694999d1af310f46e853f6`.
Subsequent commits change forensic tools and this report, without changing
the solver, kernels or C++ acceptance tests.

**Comparison-tool defect.** The initial canonical invocation returned 266
missing and 266 new assertion identities in two CLI tests, despite equal
counts and no changed matched values. Their INFO contexts contained bare
string-pointer addresses, random MeshTemp paths and timing columns. The
comparator now numbers those CLI addresses by first appearance and removes
only the named timing column, solve duration and generated temporary prefix.
It retains numerical columns, artifact names, original expressions and all
expanded values. Six focused regression checks pass; two failed before the
fix. Changed residuals, changed outcomes, artifact names and missing
assertions remain detectable. Both full XML files are unchanged; rerunning
the requested command with the corrected comparator gives the zero-delta
result above. No numerical acceptance criterion was changed.

**Fixture coverage.** Both original full runs logged missing
`wing/naca4412_ground.msh` and `wing/naca4412_ground.kutta.json`, so their
wing CLI case checked only the donor file identity. These were generated
in this worktree with the pinned project generator and matched the existing
ledger exactly. The omitted case was then run separately on the frozen
`8866cc0` binary and branch binary: **19/19 assertions pass on both**, with
zero changed values or unmatched identities. The shared baseline and runs
root were not modified. There are no other missing-fixture messages in
the full branch run, and none in the quick run.
The corpus-FSAE test produced zero assertions in both full runs; those runs
do not establish additional coverage from that corpus filter.

**Guarded runtime.** The final retained probe (`4f017b8`) converges in 150
Metal iterations with the Péclet/differentiable settings: maximum scaled
transport defect 9.149605e-7, raw continuity 6.255915e-8 kg/s and velocity
L2 1.388443%. Péclet side velocity with the Standard limiter exhausts the
unchanged 4,000-iteration ceiling: maximum scaled defect 2.370901e-4 and
raw continuity 7.309929e-5 kg/s. It returns status 2 and emits no converged
field-error attribution. The differentiable limiter is therefore part of
the measured successful option, not an automatically selected default.

**Memory limit.** Quick-tier peak RSS/physical footprint were 1.53/1.79 GB.
The full-tier process reported **1.89 GB RSS and 5.85 GB lifetime peak
physical footprint**, exceeding the requested approximately 3 GB process
budget. No per-case footprint trace was captured, so the responsible case
or allocation is unisolated. The V6 guarded probes were 0.31 GB peak
footprint each; supplemental wing runs were below 0.31 GB. Solves ran
serially in this lane. The footprint is the process lifetime metric read
by [Apple time](https://github.com/apple-oss-distributions/shell_cmds/blob/main/time/time.c#L237-L299),
distinct from its RSS statistic. No performance claim is made from these
shared-machine runs.

## What remains, ordered by the measured budget effects

The ordering uses the fixed-field momentum budget; its substitutions are not
additive estimates of converged field error.

1. **Face-flux/high-order velocity reconstruction.** The option removes much
   of the default flux-induced momentum defect, but its saved-state
   divergence is still 472 times Fluent's final raw continuity target. The
   independent all-face error and both cold/seeded solves establish a real
   mismatch; they do not identify Fluent's unpublished high-order formula.
2. **Momentum spatial reconstruction at fixed saved flux.** The remaining
   1.02092e-3 N imbalance survives exact saved viscosity and wall shear. The
   budget does not separate pressure, convection and diffusion reconstruction
   errors from a discrepancy between the saved fields and reported final
   residual. WFGC and saved-state drift have not been separated.
   Two consecutive Fluent saves, per-cell residuals, or a converged WFGC-off
   save would distinguish them. None was available here.
3. **Smaller spatial candidates.** Harmonic viscosity and the warped-face
   pressure moment reduce the fixed-state budget by 2.76e-5 N (2.7%) and
   1.91e-5 N (1.9%), respectively. Their converged effect and cross-case safety
   are unverified; neither was adopted.
4. **Saved-state constitutive substitution.** Replacing saved mu_t by
   Kutta's reconstruction changes that budget by 5.58e-6 N (0.55%). Conversely,
   exact-mu_t/frozen-k-omega flow solves still miss the velocity gate. This
   measures these substitutions, not the nonlinear sensitivity to every SST
   change. Wall-face omega remains unisolated here and belongs with the
   turbulence lane.

The explicit isotropic-pressure addition is contradicted by this experiment;
it is not retained as a leading explanation of this reference's error.

## Reproduction and documentation

The full tiers need the documented restricted datasets under `validation/data`
and generated meshes under `KUTTA_RUNS`, in addition to the Fluent corpus.
This run supplied matching fixtures through temporary read-only links.

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
mkdir -p build/v6-momentum
export KUTTA_FLUENT_EXAMPLES="<checkout>/validation/fluent_examples"
export KUTTA_RUNS="<checkout>/runs"
build/kutta-tests --test-suite-exclude=slow
build/kutta-tests -ts=slow -r=xml --success > build/v6-momentum/branch.xml
python3 tools/forensics/compare_slow_assertions.py BASELINE.xml \
  build/v6-momentum/branch.xml > build/v6-momentum/slow-comparison.json
python3 tools/forensics/test_compare_slow_assertions.py

# Saved-state comparison of the actual retained runtime option:
build/bench_v6_operators "$KUTTA_FLUENT_EXAMPLES/mixing_elbow" - metal \
  peclet-side differentiable flux=build/v6-momentum/pd-flux.csv
# In a Python environment with numpy, scipy and h5py:
python -I tools/forensics/v6_momentum_budget.py \
  "$KUTTA_FLUENT_EXAMPLES/mixing_elbow" --cache build/v6-momentum/cache \
  --kutta-flux build/v6-momentum/pd-flux.csv
```

For a converged field-error split, add `converge` to `bench_v6_operators`.
The supplemental wing proof uses the recorded donor case plus these two
locally generated artifacts, with `KUTTA_RUNS` pointed at that local root:

```sh
uv run tools/generate_fixtures.py --runs build/v6-momentum/complete-runs \
  --only wing/naca4412_ground.msh --only wing/naca4412_ground.kutta.json
```

The four full V6 option runs used the unchanged existing test with only the
two settings above forced in a temporary constructor hook. The focused
cross-case option runs used the same hook. It is absent from committed code.

Local numerical evidence stays under `build/v6-momentum`: baseline and option
V6 logs and CSVs, operator/flux dumps, JSON budgets, rejected-trial logs,
focused option-reference XML, full branch XML and comparison, quick/lint logs.
Restricted data, generated meshes, XML and fixture symlinks are not committed.
All 22 temporary fixture links were removed after verification. The captured
T3A scratch directory was already absent at cleanup; no other agent's
temporary directory was touched.

Fluent-documented statements are limited to the following. The pressure-based
face flux uses a momentum-weighted velocity and an unspecified high-order
correction ([Theory Guide 23.4.2](https://ansyshelp.ansys.com/public/views/secured/corp/v252/en/flu_th/flu_th_sec_solver_pb_cont_eq.html)).
WFGC is intended to improve gradients on warped, high-aspect-ratio and deformed
cells; its reproducible algorithm is not supplied in that description
([User's Guide 36.21.6](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_sec_solve_conv_stab.html)).
The Péclet construction and the rejected geometric terms above are our
hypotheses and experiments.
