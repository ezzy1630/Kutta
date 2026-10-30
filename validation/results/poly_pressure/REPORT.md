# Polyhedral face reconstruction: Taylor–Couette and regression evidence

The V7.2 pressure mode comes chiefly from skewed momentum interpolation in
the mass predictor. Controlled ablations distinguish this from the real viscous
flux defect: correcting diffusion alone leaves 1.678% pressure error, correcting
the predictor alone gives .0686%, and correcting both gives .0927%. The retained
change improves quadratic tangential reconstruction in both operators. Second
Order pressure interpolation, pressure LSQ, the momentum pressure source and the
Rhie–Chow pressure term are unchanged. The Taylor–Couette analytic solution is
the reference; no Fluent installation or same-mesh annulus run exists or is required.

## Converged solution and unchanged acceptance gates

CPU and Metal use the same original tests and limits. Pressure removes a volume-weighted constant gauge and is normalized by the analytic radial pressure range.

| Poly cells | Before pressure L2 (%) | After pressure L2 (%) |
| ---: | ---: | ---: |
| 384 | 9.076386 | 2.927974 |
| 1,536 | 3.607781 | 0.521156 |
| 6,144 | 1.691271 | 0.092681 |

The 6,144-cell gate passes below 1%; the coarse diagnostic remains above 1%. The three-level observed poly pressure order is about 2.49, versus 1.21 before. Hex pressure error is 0.311599%. The native rotating-wall adapter hold is removed; its quick test checks the imported angular velocity.

| V6 backend / start | Before velocity / pressure L2 (%) | After velocity / pressure L2 (%) |
| --- | ---: | ---: |
| cpu / seeded | 3.601602 / 2.962540 | 3.601605 / 2.962466 |
| cpu / cold-hybrid | 3.601460 / 2.961665 | 3.601465 / 2.961828 |
| metal / seeded | 3.601603 / 2.962296 | 3.601605 / 2.962460 |
| metal / cold-hybrid | 3.601460 / 2.961921 | 3.601464 / 2.961693 |

V6 remains approximately 3.60% velocity / 2.96% pressure error. Its same 29 accuracy assertions fail; this laminar reconstruction does not change its SST operators. Small integrated-source/rounding differences are retained in the JSON.

## Cause and formulation

All ablations use the same 6,144-cell poly mesh, analytic seed, viscosity one,
1,200 iterations and CPU solver, with only the indicated operator disabled.
The production source is unchanged by these diagnostics; isolated patches are
preserved as [diffusion-only](ablation_diffusion_only.patch) and
[predictor-only](ablation_momentum_only.patch).

| Active correction | Pressure relative L2 (%) | Velocity L2 |
| --- | ---: | ---: |
| Original operators | 1.691271 | .0003241313 |
| Projected diffusion only | 1.677778 | .0003017230 |
| Projected momentum predictor only | .068611 | .0003261246 |
| Both, retained | .092681 | .0003025337 |

The predictor-only pressure error is lower than the retained result. The
combined method is retained for its independent quadratic viscous-flux accuracy
and improved V2 finest-grid field/traction errors, rather than selecting the
smallest annulus pressure number. The pressure error's viscosity scaling is not
unique evidence against diffusion: in this diffusion-dominated flow, Rhie–Chow
mobility is proportional to `V/aP`, approximately `1/mu`. Compensating a spurious
face-normal predictor velocity therefore requires a pressure mode proportional
to viscosity. This is an inference from the published coupling formula and the
controlled ablations, not an independent Fluent observation.

The original kernel implements Fluent's published primary/secondary diffusion
split correctly. An independent double oracle agrees within float roundoff,
including with exact cell gradients. The defect is therefore not a missing
secondary term or the wrong primary normal distance. Its linearly reconstructed
secondary term does not resolve quadratic tangential variation on skewed faces.
The coarse/fine radial faces expose this most clearly: their exact-state viscous
normal-flux error converges at order 0.99 before the change and 1.95 afterward.

Let F be the face centroid, n its unit normal, P/N the cell centroids,
dp = n·(F−P), dn = n·(N−F), and alpha = A/(dp+dn). Define the projected points
P* = F−n dp and N* = F+n dn. A geometry-only quadratic least-squares fit qhat
supplies the deferred correction

```
C = alpha [qhat(N*) − qhat(N) − qhat(P*) + qhat(P)]
Fdiff / mu = alpha (qN − qP) + C.
```

At a Dirichlet wall N=F and dn=0. The implicit primary coefficient is unchanged.
For a linear field, C is exactly `grad(q)·(S−alpha*(N−P))`, the documented
secondary diffusion term. For a quadratic field the tangential shifts are exact.
The normal two-point curvature error is deliberately retained, including at
orthogonal walls; replacing it by a wide fitted derivative caused cavity
regressions in rejected trials.

Momentum face velocity uses the same normal-projected states with geometric
normal interpolation, while retaining the compact momentum-diagonal-weighted
cell coupling. If w=aP/(aP+aN), its deferred correction is

```
qstar = [dn*qhat(P*) + dp*qhat(N*)] / (dp+dn)
qface = w*qP + (1−w)*qN
      + [qstar−qhat(P)] − (1−w)*[qhat(N)−qhat(P)].
```

This corrects the skew curvature in the mass predictor. It does not replace the
Rhie–Chow pressure term. With limiting enabled the existing cell Barth limiter
bounds/blends this correction from the original face predictor; a zero limiter
recovers that predictor exactly.

Fits are centered at F, include both adjacent cells, and use symmetric weights
`1/(|sample−P|²+|sample−N|²)`. Three overlapping neighbour rings supply the
cloud, expanded to five only when the requested functionals are rank deficient.
Only independently prescribed boundary face values enter the fit. Roundoff-
orthogonal faces retain the original stencil without allocating reconstruction
weights. Rank-deficient faces and clouds containing discontinuous prescribed
wall traces retain the original method. This trace check uses rigid wall motion
before projection onto faceted normals, so smooth rotation is not misclassified.

Wall forces retain the fitted tangential derivatives and replace their normal
component with the exact flux assembled by the compact momentum stencil. The
full Newtonian traction therefore uses the same normal viscous flux as the
momentum equation, even for non-polynomial fields.

This is a standard quadratic geometric reconstruction extending the documented
linear split. It is **not a claim that Fluent uses these polynomial weights**:
Fluent's public documentation does not specify that reconstruction. The latest
task instruction permits a standard consistent treatment where that formula is
not published. The [Fluent Customization Manual, §3.2.5, Eq. 3.3](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/pdf/Ansys_Fluent_UDF_Manual.pdf)
specifies the primary/secondary split; the [Theory Guide, §23.4.1](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/flu_th_sec_solver_pb_mom_eq.html)
specifies central Second Order pressure and the separate body-force-weighted
option. No body-force-weighted pressure scheme was needed to pass this annulus.

## Operator evidence

The probe evaluates production kernels on exact Taylor–Couette cell fields at
8, 16 and 32 radial rings, for hex and poly meshes. It separates azimuthal,
regular radial, coarse/fine radial, inner/outer wall and cap faces, and separates
interior/inner-wall/outer-wall cells. The poly sequence has 384, 1,536 and 6,144
cells. Both analytic wall observations and the actual faceted rotating-wall
values are measured, with LSQ and exact velocity/pressure gradients, on CPU and
Metal. Every norm, maximum and adjacent refinement order is preserved in
[before](operators-baseline.json) and [after](operators-after.json); the readable
[operator table](operator_orders.md) contains every CPU analytic-wall/LSQ row.
The raw probe logs are [before](operators-baseline.log) and [after](operators-after.log).

Face RMS norms are area weighted; cell RMS norms are volume weighted. The probe
reports both centroid face fluxes and signed seven-point triangle-integrated
fluxes. Triangle orientation is retained for concave poly cap polygons. Pressure
source truncation is also compared against integrated analytic cell forcing,
not only the pressure gradient sampled at a centroid.

The first-order-operator heuristic needs qualification. The inner-wall viscous
normal derivative is first order on **both** hex and poly meshes; the unchanged
one-sided normal wall stencil has that local truncation error. Its cell residual
need not vanish pointwise with refinement. The faceted physical rotating-wall
trace adds another O(h²) velocity mismatch divided by an O(h) wall distance.
The fixed-thickness axial caps are not axially refined. None of these observations
can by itself diagnose the poly-only pressure mode. Local flux, local divergence
and converged solution orders are distinct; [Diskin and Thomas, AIAA 2012-0609](https://ntrs.nasa.gov/citations/20120001451)
explicitly separates these error measures on irregular meshes.

The retained reconstruction improves coarse/fine radial viscous flux and the
coupled mass predictor. Bulk continuity at the exact state improves from order
0.94 to 2.34. Some local errors remain first order, and some magnitudes increase:
for example bulk viscous divergence RMS at 32 rings rises from .05972 to .08409,
although its observed order improves from approximately zero to .95. Azimuthal
viscous face error also increases. These are recorded rather than hidden behind
the improved pressure solution. The complete before/after records are needed to
judge the remaining irregular-mesh truncation error.

The quadratic pressure checks also narrow the original hypotheses. Linear
pressure and boundary-cell gradients reproduce the manufactured field within
float tolerance. Quadratic cell LSQ gradients are only first order on this
irregular sequence, while quadratic face pressures approach second order. The
recent remote-extrapolation bound is inactive: the largest normalized skew is
.3101, below its activation threshold of one. Thus removing that bound cannot
repair V7.2. Pressure Gauss gradients retain first-order local errors, but the
same unchanged pressure source now supports a passing pressure solution. The
radial Rhie–Chow normal-velocity contribution is about 1.02e-6 at 32 rings,
compared with the original mass-predictor error of 1.12e-4; its observed order
is 3.35 and it is unchanged by the repair.

All 224 pressure-operator records and all 64 isolated Rhie–Chow pressure-term records are identical
before and after; so are the 80 isolated convective-flux records. The converged
coupled convective flux changes through the corrected mass predictor, as it
should. See [the comparison](unchanged_operators.json).

## Manufactured proofs

A paired manufactured velocity `(u,v,w)=(x²,−2xy,0)` uses identical
observations before and after, including exact cell gradients. The finest-poly
area-weighted RMS derivative errors and 16→32-ring orders are:

| Component / face | Before RMS | Before order | After RMS | After order |
| --- | ---: | ---: | ---: | ---: |
| u / azimuthal | 2.81935e-4 | 2.011 | 1.12880e-6 | float noise |
| u / coarse-fine radial | .00683831 | .995 | 8.56878e-5 | 1.987 |
| u / inner wall | .00961636 | 1.007 | .00961636 | 1.007 |
| u / outer wall | .0216536 | .957 | .00972430 | 1.023 |
| v / azimuthal | 5.63867e-4 | 2.011 | 1.28940e-6 | float noise |
| v / coarse-fine radial | .0136750 | .995 | 9.89811e-5 | 1.986 |
| v / inner wall | .0111045 | 1.007 | .0111045 | 1.007 |
| v / outer wall | .0376235 | .974 | .0112285 | 1.023 |
| u,v / axial cap | 0 | — | 0 | — |

The inner-wall numbers are exactly the same on hex and poly; hex outer-wall
errors are .00954260 for u and .0110183 for v, both order .996. CPU and Metal
RMS records match exactly for this diagnostic. This directly distinguishes
quadratic tangential reconstruction from a wall-normal-distance mistake.
All three levels, LSQ/exact variants and face types are in
[quadratic before](quadratic-baseline.json) and [after](quadratic-after.json),
with raw logs beside them. The probe is `bench_quadratic_viscous`.

The quick tier includes independent linear/quadratic pressure and original
primary/secondary diffusion oracles, plus new checks on Cartesian hex, tetrahedral
and all three poly annulus levels:

- Quadratic tangential projection reproduces the independent normal-foot finite
  difference and normal-central interpolation, including unequal momentum
  diagonals. It does not falsely demand a quadratic-exact one-sided wall normal
  derivative. Maximum float errors are below 5.8e-5 for the normal derivative,
  1.9e-6 for interpolation, and 2.2e-5 for the fitted gradient.
- For non-polynomial fields, the wall-gradient normal flux mismatch falls from
  .0948 (tet) / .365 (poly) to below 2.5e-7, while tangential derivatives change
  by less than 1.6e-7.
- Reversing cell numbering reverses physical face ownership while preserving the
  mesh invariant owner<neighbour. Corrections reverse sign and physical face
  velocities remain invariant within 2e-6. This guards against a rejected
  one-sided fit that passed the annulus but broke generated/imported cavity parity.
- An all-zero cell limiter restores the original face predictor exactly.

The general quadratic test includes axial curvature; axial cap error remains
constant under purely radial/azimuthal refinement because layer thickness is
fixed. Taylor–Couette itself is axially constant and has zero cap flux.

## Scope and unresolved limits

The correction is enabled for laminar momentum. SST, SA and realizable-k-epsilon
retain their existing diffusion and wall-function contracts. In particular this
work does not change V6's SST viscous operators or resolve its unmatched WFGC
setting. It does not establish Fluent same-mesh agreement, full FSAE-car
acceptance, arbitrary-polyhedral second-order pointwise truncation, or CUDA GPU
execution parity. Existing tetrahedral wall-load accuracy failures remain
visible; no original acceptance assertion or threshold was weakened.

Geometry fitting occurs only at setup. Every field-dependent reconstruction,
transport assembly and force operation runs through Device and generated
CPU/CUDA/Metal kernels. No app source, convergence control, pseudo time or
linear solver was changed by this lane. Main's independently developed
convergence changes are included by the final rebase and identified in metadata.
Historical diagnostics and rejected trials are preserved under [historical](historical/README.md);
their old `final` names do not identify this retained scheme.

## Full assertion comparison and runtime proof

The original baseline is `fb827e0d416746030eaf156f9f42fca87b712a31`. Numerical source is `d4fa853`, rebased onto main `587de30`. Machine: Apple M4 Pro / 24 GiB, macOS 26.6.2 arm64, AppleClang 21, C++20 Release. [Metadata](metadata.json) records source, generated-code, executable and fixture hashes. Every check set `KUTTA_FLUENT_EXAMPLES`; TMR grids were unzipped. The Linux CPU run additionally verified all 13 source-relative fixture hashes before executing.

All 24 slow cases and 10,884 assertions were compared by case, normalized source location, expression, INFO context and occurrence. No assertions are missing or added. Changed expanded values are retained even when outcomes match. The only context normalization is the known SA INFO ternary-pointer address, scoped to that case and exact message shape; a synthetic check confirmed distinct identities, changed-value detection, new-failure detection and missing-assertion detection.

| Slow case | Before failed | Final failed |
| --- | ---: | ---: |
| FV V3 Re1000 cavity refinement against Erturk | 0 | 0 |
| FV V3 Re100 and Re400 cavity reported deviations | 0 | 0 |
| FV V4 Fluent cavity matches its actual grid | 0 | 0 |
| FV V6 Fluent mixing elbow fixed point and cold start | 29 | 29 |
| FV V2 manufactured Navier Stokes hex refinement | 0 | 0 |
| FV V2 manufactured Navier Stokes tet refinement | 6 | 6 |
| FV V1 manufactured scalar diffusion and upwind refinement | 0 | 0 |
| FV SIMPLEC develops planar Couette flow with moving and slip walls | 0 | 0 |
| V5 TMR flat plate SST accuracy and three-grid trend | 6 | 6 |
| V5 TMR bump SST station accuracy and three-grid trend | 10 | 9 |
| FSAE V12 native Fluent polynomial fan import acceptance | 1 | 1 |
| V8 realizable KE scalable channel MKM DNS profile and momentum balance | 0 | 1 |
| V7 SA TMR flat plate and bump accuracy convergence and three-grid trend | 17 | 17 |
| FSAE V7.1 translating-wall Couette hex and poly | 0 | 0 |
| FSAE V7.2 rotating-wall Taylor Couette hex and poly | 2 | 0 |
| FSAE V8 Darcy Forchheimer pressure drop and oblique directions | 0 | 0 |
| FSAE V9 absolute-velocity MRF analytic solutions and frame invariance | 0 | 0 |
| FSAE V10 incompressible pressure and mass-flow inlets | 0 | 0 |
| FSAE V11 internal fan polynomial and porous jump pressure discontinuity | 0 | 0 |
| FSAE V7.3 SST moving-wall and translating-frame invariance | 0 | 0 |
| FSAE V9.2 conformal MRF interface preserves absolute Taylor Couette | 0 | 0 |
| FSAE V10.2 SST pressure and mass inlets match velocity-inlet transport | 0 | 0 |
| FSAE V9.3 rotation and translation preserve Taylor Couette pressure | 0 | 0 |
| Fluent research downloads retain published byte identities | 0 | 0 |

[Complete assertion comparison](comparison_final_source.json), [final XML](final-slow.xml), [final output](final-slow.log), and [original full baseline](baseline_all_slow.xml) preserve the evidence. Earlier full/run-attempt outputs are separate; they are not substituted for failures in this final run.

The final full tier is **not green**: Metal's unchanged k-epsilon channel stopped at 5300 steps / 182.5 seconds. Its new failed assertion is the time-limited convergence check. Original-channel trajectories agree, but Metal wall time varied enough to stop before the required 5,700 steps. A current-machine original-binary replay converged at 5,700 / 183.7 seconds; candidate attempts stopped at 4,300, 5,400 and 4,200 steps at the 180-second check. The legacy transport shader/signature has been restored verbatim, with projected diffusion on a separate entry point. This timing instability remains unresolved; no budget or convergence gate was changed. V1–V6 have no new failed assertions.

V2 finest tetrahedral errors improve from .0131268 to about .0118707 for velocity, .0764240 to .0558932 for pressure, and .0418577 to .0285646 for total wall force. Some order estimates decline: velocity 1.9778→1.9020, viscous wall load .8613→.8111, and total wall load .8718→.5468. The velocity gate remains passed and the same six wall-load assertions remain failed; lower finest-grid errors do not mean every measured slope improved. The comparator records all such values.

The apparent recovered SST-bump assertion is a budget endpoint effect (more CPU iterations), not evidence for this laminar reconstruction. SA/SST unconverged results remain diagnostics. Native polynomial fan import still lacks its V12 fixture/encoding acceptance.

## Cost and CUDA portability

One serial end-to-end sample uses the same 6,144-cell / 1,200-step scenario per backend. It includes setup, initialization, solve and final reductions; CPU also includes initial device discovery. It is not a kernel-throughput benchmark.

| Measurement | Before | After |
| --- | ---: | ---: |
| cpu scenario (seconds) | 45.85 | 57.14 |
| metal scenario (seconds) | 50.07 | 50.67 |
| Peak process RSS (MiB) | 54.39 | 121.41 |

[Timing and executable identities](performance.json) record the sample. Overlapping face clouds add storage and work; car-scale memory and throughput were not measured. No G4 speed claim is made.

CUDA 13.0.2 container / NVCC 13.0.88 compiles the complete C++/generated-kernel stack for `75-real;86-real;89-real;120`. A separate Linux CPU-only build passes 225 cases / 3,756,393 assertions with all 13 fixture hashes verified, including source-relative legacy readers. See [compile output](cuda-compile.log) and [accepted fixture-backed CPU output](cuda-fixture-quick.log). The first container CPU pass had source-relative fixture skips and is not the accepted fixture proof. No NVIDIA device was available, so CUDA GPU execution remains unverified. Native quick passes 225 cases / 7,413,872 assertions, compiled Metal presence and portability lint ([output](quick-ctest.log)).

## Reproduction

```sh
export KUTTA_FLUENT_EXAMPLES="${KUTTA_FLUENT_EXAMPLES:?set the external fixture root}"
python3 tools/prepare_tmr.py
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 6
ctest --test-dir build -LE slow -V
build/kutta-tests --test-suite=slow --success --reporters=xml --out=slow.xml
python3 tools/forensics/compare_slow_assertions.py baseline_all_slow.xml slow.xml
build/bench_taylor_couette_operators --baseline > operators-baseline.log
build/bench_taylor_couette_operators > operators-after.log
build/bench_quadratic_viscous --baseline > quadratic-baseline.log
build/bench_quadratic_viscous > quadratic-after.log
```
