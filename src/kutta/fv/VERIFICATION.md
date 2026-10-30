# Lane D verification

This branch preserves the Fluent-class discretization and checks its operator
consistency separately from integrated-load accuracy. The former tet pointwise
pressure ≥1.8 gate was mathematically incompatible with that operator. The new
wall-load ≥1.8 gates remain strict and fail on tets; V2 is therefore **not fully
passed**. The V3 acceptance criterion is now a refinement study, not agreement
of a finite grid with rounded grid-converged numbers.

The current build/test results below use this worktree, Apple M4 Pro (24 GB),
macOS 26.6.2, AppleClang 21, C++20 Release (`-O3 -DNDEBUG`), HDF5 and Metal enabled.
CPU denotes the serial reference backend. CUDA runtime is not verified here.
The V1, V4 and performance tables retained below originally used source
`5e527dc856280e6e69f436f0fc7438b5340a903d`; they are historical measurements,
not fresh performance claims. Current correctness runs are listed at the end.

## Update after review fixes (2026-10-30)

The automatic pseudo-time length scale now uses the domain volume, as Fluent's
conservative method does, instead of the mean cell volume (about N^(1/3) smaller).
Tests and `bench_fv` now use Fluent's default time scale factor 1; the former
factors 3 and 5 compensated for that error. Closed domains pin the cell containing
`reference_point` instead of internal cell 0. Tables below that report iteration
counts or wall times predate this change; their accuracy values are unchanged
to the printed digits except where listed here.

Fresh slow tier (same machine): 6/7 cases pass in 1,261 s (was 2,891 s), and in 1,012 s after rebasing onto the device-resident linear solvers (same iteration counts). The
only failures are the six retained tet wall-force ≥1.8 assertions; the measured
orders are the same as below.

| Case | Outer iterations before | after |
| --- | ---: | ---: |
| V3 Re=1000, 65² / 129² / 257² (CPU) | 1,400 / 1,900 / 3,300 | 700 / 600 / 1,000 |
| V3 Re=100 / Re=400, 129² | 400 / 800 | 1,000 / 500 |
| V4 Fluent mesh Re=100 / 400 / 1000 | 300 / 600 / 900 | 600 / 300 / 500 |
| V2 tet n=16 | 500 | 440 |

On the 3,550-cell Fluent cavity mesh (Metal, relaxation 1) the outer iterations
to the driver's steady criterion for time scale factors 0.1 / 0.3 / 1 / 5 were
700 / 300 / 600 / 2,200 at Re=100 and 2,700 / 1,000 / 500 / 700 at Re=1000: the
best factor depends on the flow, and Fluent's default is a reasonable middle.

## Mesh, gradients and flow checks

CPU and Metal pass analytic wall-distance tests on hex/tet boxes translated by
1e8, finite wall tiles with edge/corner nearest points, and no-wall domains.
The host BVH evaluates exact double point-to-triangle distances on each wall
face's triangle fan, then stores floats. Wall-normal distances are checked
against analytic projections; non-wall entries are zero.

The inverse-distance-squared LSQ reconstruction includes actual boundary face
values. Maximum affine-gradient component errors are 3.8147e-7 (60-cell hex),
8.58307e-7 (144-cell tet), and 5.20766e-5 (86,803-cell real Fluent poly box), equal
on CPU and Metal. On the poly mesh, single-precision field quantization is
amplified by short centroid spacing: the kernel agrees with an independent
host reconstruction of the rounded field within 2e-6, and every analytical error
lies inside its propagated input-rounding bound plus 2e-6. The limiter keeps
every cell-to-face reconstruction inside its stencil range. Fused vector
and separate scalar gradient kernels agree bitwise.

Rest, forward/reverse inlet–outlet flow, default density-aware normal backflow,
and moving-wall/symmetry checks pass on both backends. In planar Couette flow,
400 outer iterations give maximum velocity error 1.19209e-7 and scaled momentum
defect 4.04177e-8. Top-wall viscous force is -0.1 in x and coefficient is -2 for
the supplied reference values. Shifting reference pressure by 2 changes the
normal force by -0.2, as expected. Moving-wall normal motion is projected out.

`SteadySolver`, settings, patch-condition variants, `iterate(n)`, monitor snapshots,
and wall forces are exposed without a CLI; caller usage is in [README](README.md).
Solvers are selected in one place, `transport.hpp`: BiCGStab/Jacobi for transport
and AMG-PCG for pressure. Closed domains use symmetric reference-cell elimination;
pressure outlets remove the pin. `device_mesh.hpp` has no changes.

## V1 scalar manufactured solution

The field is sin(pi*x) sin(pi*y) sin(pi*z), diffusivity 1, with exact source and
Dirichlet face data. Convection speed is 4 in x. Each study runs 8/16/32 cells per
axis; tet boxes have six tets per brick. Reported orders are the finest 16→32
pair after three refinements, not a fit over all three levels.

| Mesh | Backend | Diffusion order | FOU order | SOU order |
| --- | --- | ---: | ---: | ---: |
| Hex | CPU | 2.00206 | 0.854671 | 1.91201 |
| Hex | Metal | 2.00192 | 0.854671 | 1.91201 |
| Tet | CPU | 1.93250 | 0.993216 | 1.94919 |
| Tet | Metal | 1.93243 | 0.993216 | 1.94917 |

The hard finest-pair gates are diffusion ≥1.9, SOU ≥1.8, FOU ≥0.8. All pass.
The coarser tet diffusion slope is 1.82407, so the whole sequence is not claimed
to have an order ≥1.9. Finest CPU L2 errors (diffusion/FOU/SOU) are
0.000284123/0.00630692/0.000249390 for hex, and
0.000288447/0.00412068/0.000286845 for tet. Finest outer counts are 3/3/6 and
19/21/18 respectively; Metal uses 17 rather than 21 for tet FOU, with the same
reported error. Scalar convergence on a poly-dual family was not measured.

## V2 Navier–Stokes manufactured solution

The smooth divergence-free polynomial velocity, cosine pressure, viscosity 0.1
and exact forcing are in `tests/fv/test_navier_stokes.cpp`. Pressure errors remove
a volume-weighted constant gauge offset. Hex studies are extruded 2D boxes
(8/16/32 by the same number in y, one z layer); tet studies are full 3D boxes
(4/8/16 per axis). Orders again use the finest pair.

| Mesh | Velocity order | Pressure order | Finest velocity L2 | Finest pressure L2 | Finest outer iterations |
| --- | ---: | ---: | ---: | ---: | ---: |
| Hex, CPU and Metal | 1.91567 | 1.82088 | 0.00418077 | 0.00300871 | 840 |
| Tet, CPU and Metal | 1.98247 | 0.942387 | 0.0122571 | 0.0689336 | 500 |

### Pressure-operator contract

For this exact MMS pressure, the cell-centred Gauss gradient
`sum_f p(face centroid) S_f / V` is first order on Freudenthal tets even when
face pressure is exact. Midpoint-face quadrature errors cancel across opposite
hex faces, giving second order there; tets do not have those opposite pairs.
Inverse-distance² LSQ and Gauss applied to reconstructed Second Order pressure
also remain first order on tets. The deliberately matched Fluent-class momentum
pressure force uses this same Gauss form. Requiring second-order pointwise tet
pressure from this family was an unattainable gate, not evidence of a Kutta bug.
No reconstruction, quadrature, limiter or solver discretization was changed.

`FV MMS pressure gradient operator convergence` independently computes the
exact-face Gauss operator in host double precision. It also runs the actual
LSQ -> `fv_pressure_faces` -> `fv_momentum_source` device kernels on exact
cell pressure, supplying exact boundary face pressure to isolate the operator
from boundary extrapolation error. The source kernel returns the negative
pressure gradient; that sign is removed before comparison. All errors are
volume-weighted vector L2 norms against the analytic cell-centre gradient.
The hex operator mesh is extruded 2D, matching V2; tets are full 3D.

The exact-face finest-pair gate is tet order in [0.9,1.1] and hex ≥1.9. The
actual tet momentum operator must be ≥0.9. The pre-existing LSQ order ≥0.9
check remains in place on both meshes. The V2 solver tet pressure order must
be ≥ its independently measured actual momentum-operator order minus 0.1,
using the same 8→16 tet refinement. Hex pressure remains ≥1.8. The existing
velocity gate ≥1.9 is retained on both meshes (stronger than ≥1.8).

| Mesh | Exact-face Gauss order | Actual momentum operator order | LSQ order | Solver pressure order | Solver consistency threshold |
| --- | ---: | ---: | ---: | ---: | ---: |
| Hex, 8→16 | 1.99791 | 1.36247 | 1.56286 | 1.82088 (16→32) | 1.8 |
| Tet, 8→16 | 0.999409 | 1.01406 | 0.994008 | 0.942387 | 0.91406 |

CPU and Metal agree to the printed precision. Only the exact-face hex Gauss
operator must be ≥1.9. The reconstructed operator has a lower measured slope;
the load checks below assess its effect on the converged solution.

### Integrated wall-load contract

The same converged V2 solutions call `SteadySolver::forces` on the lower half
of xmin (`x=0, 0<=y<=1/2`, full span). Boundary faces are split into named wall
patches without moving or reordering any face. This avoids zero net pressure
from whole-wall cancellation; the selected wall also has nonzero viscous force. Pressure is aligned with the same
volume-mean gauge offset used for its L2 error. Force errors are vector norms
per unit span; the hex span is 1/n and the tet span is 1.

With outward normal `n=(-1,0,0)`, wall-on-solid pressure traction is
`p*n=(-cos(pi*y),0,0)` and viscous traction is
`-mu*(grad u+grad u^T)*n=(0,-20*y²*(1-y)²,0)`. Exact antiderivatives give
pressure force/span `(-1/pi,0,0)` and viscous force/span `(0,-1/3,0)`.
No approximate reference quadrature is used. The total-force gate cannot hide
cancellation because pressure and viscous forces are also gated separately.

| Mesh | n | Pressure error/span | Viscous error/span | Total error/span |
| --- | ---: | ---: | ---: | ---: |
| Hex | 8 | 0.0397838 | 0.0713798 | 0.0351425 |
| Hex | 16 | 0.0117129 | 0.0197706 | 0.00914660 |
| Hex | 32 | 0.00310731 | 0.00506342 | 0.00225342 |
| Tet | 4 | 0.161730 | 0.222746 | 0.152052 |
| Tet | 8 | 0.0714325 | 0.133492 | 0.0778614 |
| Tet | 16 | 0.0358039 | 0.0727936 | 0.0407649 |

| Mesh | Pressure-force order | Viscous-force order | Total-force order | ≥1.8 gates |
| --- | ---: | ---: | ---: | --- |
| Hex | 1.91437 | 1.96517 | 2.02112 | pass |
| Tet | **0.996464** | **0.874869** | **0.933580** | **all fail** |

The three ≥1.8 tet force assertions fail on each backend. This MMS does not
establish second-order tet loads for FSAE use, so V2 remains failed.

## V3 cavity refinement

Re=1000 uses uniform 65², 129² and 257² cell grids against every point in the
supplied Erturk 601² u/v centreline tables. The error is the maximum absolute
deviation separately for u and v, with linear interpolation to reference
coordinates. Each finest-grid error must be smaller than the 129² error, and
its finest-pair observed order must be ≥1.6. Orders use the actual cell-width
ratio 257/129, not an assumed ratio of two. The 65² level and its slope are
also reported. A 129² second-order solution cannot reasonably be required to
fall inside the printed rounding intervals of a nearly grid-converged 601²
reference; those intervals have been removed as an acceptance criterion.

The solver settings are unchanged: unit lid speed, viscosity 1/Re, symmetry
front/back, momentum/pressure relaxation limits 1, pseudo-time scale 5, monitors
every 100 iterations. Convergence requires momentum defects <5e-7 and corrected
mass L1 <1e-5. The refinement study tightens the maximum velocity update over
100 iterations to 2e-6 (from the reporting driver's default 2e-5).

| Grid | CPU outer iterations | Max u error vs Erturk | Max v error vs Erturk | u order | v order |
| --- | ---: | ---: | ---: | ---: | ---: |
| 65² | 1,400 | 0.020743157625 | 0.030959975930 | — | — |
| 129² | 1,900 | 0.004344232038 | 0.007102702184 | 2.280873 | 2.147893 |
| 257² | 3,300 | 0.000607678559 | 0.001067561871 | 2.853709 | 2.749453 |

Metal matches the 65²/129² rows to every printed digit. At 257² its max errors
are 0.000610226658 for u and 0.001065764046 for v, giving finest-pair orders
2.847639 and 2.751898, after 3,700 outer iterations. Both backends pass.
Fine-grid single-precision results
are close but not bitwise identical.

All four finest-pair component orders exceed 1.6 and both errors decrease
on each backend. These slopes
measure error against the finite, rounded 601² reference; they do not establish
a scheme order above two.

CPU centre and primary-vortex values:

| Grid | u at (0.5,0.5) | v at (0.5,0.5) | Primary-vortex psi | Primary-vortex omega |
| --- | ---: | ---: | ---: | ---: |
| 65² | -0.063983514905 | 0.024281581864 | 0.116309238798 | 2.047185257478 |
| 129² | -0.062479555607 | 0.025318449363 | 0.118306555813 | 2.062453139863 |
| 257² | -0.062144596130 | 0.025664789602 | 0.118777723601 | 2.066107778952 |
| Richardson, order 2 | -0.062031779164 | 0.025781439698 | 0.118936416589 | 2.067338689919 |
| Erturk extrapolated | — | — | 0.118942 | 2.067213 |
| Botella/Peyret spectral | — | — | 0.1189366 | 2.067753 |

Metal Richardson values are u=-0.062034249242, v=0.025782229027,
psi=0.118937994671 and omega=2.067430924350, at vortex position
(0.530771922789, 0.565226268158). Its psi differs from Botella by +1.39e-6
and omega by -3.22e-4.

The CPU Richardson vortex position is (0.530771834413, 0.565226296680), compared
with Botella/Peyret (0.5308, 0.5652). Its psi differs from Botella by -1.83e-7
and its omega by -4.14e-4. Erturk's extrapolated row has no centre coordinates.
The raw 601² table centre velocities are -0.0620 and 0.0258.

Second-order Richardson estimates use the finest pair and actual width ratio.
The primary vortex is measured independently of centreline velocities: integrate
`psi_y=u` from the stationary bottom wall with trapezoids, fit the local minimum
with a quadratic including its mixed derivative, and interpolate `-curl(u)`
there. A cheap analytic-streamfunction test checks this diagnostic's position,
sign and second-order psi/vorticity errors. These comparisons report deviations
from the supplied primary-vortex CSV; they do not impose rounding tolerances.

Re=100/400 comparisons retain reported deviations and convergence assertions;
no unsupported scatter band is inferred from Ghia alone. Only the supplied
SOURCE.md's suspect Re=400 v point at x=0.90625 is omitted.

Fresh Re=100/400 results on both backends:

| Re | Outer iterations | Ghia max u deviation | Ghia max v deviation |
| --- | ---: | ---: | ---: |
| 100 | 400 | 0.005039 | 0.008967 |
| 400 | 800 | 0.002638 | 0.004046 |

## V4 real Fluent cavity mesh

The historical driver read the exact file
`/Users/ezzy/.Others/Kutta/validation/fluent_examples/lid_driven_cavity/lid_driven_cavity.msh`.
It has 3,550 cells (71×50); `top_moving_wall` selects the moving wall. The driver
also solves a uniform box at the same cell count and a generated box mapped to
the exact Fluent x/y grid lines. This distinguishes spacing effects from import
or connectivity errors. All imported and comparison cases converge on CPU and
Metal. The comparison below concerns centerline velocities, not full cell fields.

| Re | Imported outer iterations | Ghia max abs u deviation | Ghia max abs v deviation |
| --- | ---: | ---: | ---: |
| 100 | 300 | 0.00530634 | 0.00904564 |
| 400 | 600 | 0.00328528 | 0.00343363 |
| 1000 | 900 | 0.00894912 | 0.00409038 |

| Re | Matched-grid max u difference | Matched-grid max v difference | Uniform-grid max u difference | Uniform-grid max v difference |
| --- | ---: | ---: | ---: | ---: |
| 100 | 3.72529e-07 | 3.72529e-07 | 0.00252729 | 0.000836089 |
| 400 | 6.25849e-07 | 2.23517e-07 | 0.0088028 | 0.00545676 |
| 1000 | 4.02331e-07 | 3.42727e-07 | 0.0203209 | 0.012485 |

These values are the same on CPU and Metal. The largest matched-grid difference
is 6.25849e-7. Re=1000 imported-mesh deviations from Erturk are 0.0118548 (u) and
0.0212349 (v). Equal cell counts alone do not imply equal grid spacing.

The current slow test reads the worktree copy at
`validation/fluent_examples/lid_driven_cavity/lid_driven_cavity.msh`, whose
SHA-256 matches the historical external fixture: `9b09f1287db71978e10c67a528616c1c95118cfcf4f763abad047df565e6a6dd`.
It passes convergence and all centreline samples against the generated grid
with identical grid lines, with a maximum allowed deviation of 2e-6, at
Re=100/400/1000 on CPU and Metal.

## Performance

`bench_fv --size 100 --repeats 3` uses a 100³-cell moving-lid box, viscosity 0.01,
momentum relaxation 1 and pseudo-time scale 5. Mesh/solver setup and five warmup
steps are excluded. Outer step 6 is profiled with synchronization around every
operation; throughput measures asynchronous outer steps 7–9, with no periodic
monitor in that window. The 1M-cell case was not run to steady convergence.
Both backends' last pressure solve takes 11 iterations and reports convergence.

| Backend | Outer iterations/s | Wall time / outer |
| --- | ---: | ---: |
| CPU serial reference | 0.551709 | 1.81255 s |
| Metal, Apple M4 Pro | 3.09470 | 0.323133 s |

The profiler's 512-cell adapter check reproduces native Device velocity,
pressure and face flux bitwise on both backends. Timings are wall times including
dispatch/completion, not GPU event times. The following array-size model counts
the full size of each buffer argument once per call (copy counts read+write).
It is neither measured DRAM traffic nor a strict bound: boundary kernels touch
only part of some cell arrays, and gather reuse depends on cache. Consequently
modeled GB/s in the raw log can exceed hardware bandwidth.

FV kernels total **3,721.84 modeled bytes/cell/outer**, above the ~2 KB target.
All recorded operations, including Krylov/AMG, total 19,305.51 bytes/cell/outer.
The FV estimate at 240 GB/s implies 15.508 ms; the synchronized Metal FV kernels
sum to 21.588 ms. All synchronized operations sum to 1.941 s, illustrating the
intrusiveness of serial completion around the many inner-solver operations.
These measurements do not establish a speedup over Fluent.

| Kernel | Calls / outer | CPU µs / call | Metal µs / call | Logical bytes / cell / outer |
| --- | ---: | ---: | ---: | ---: |
| `fv_continuity` | 1 | 2055.8 | 396.8 | 68.12 |
| `fv_correct_flux` | 1 | 1768.2 | 785.8 | 52.24 |
| `fv_correct_velocity` | 1 | 723.4 | 734.3 | 44.00 |
| `fv_gradient` | 2 | 5162.9 | 1205.8 | 280.48 |
| `fv_limiter` | 3 | 7107.2 | 624.5 | 360.72 |
| `fv_momentum_faces` | 1 | 68861.1 | 889.5 | 185.32 |
| `fv_momentum_source` | 1 | 5786.4 | 760.1 | 128.48 |
| `fv_pressure_assemble` | 1 | 7533.7 | 901.0 | 160.00 |
| `fv_pressure_boundary` | 3 | 245.6 | 347.3 | 304.68 |
| `fv_pressure_faces` | 1 | 10055.1 | 690.6 | 136.48 |
| `fv_rhie_chow` | 1 | 15465.7 | 769.1 | 157.44 |
| `fv_transport_assemble` | 3 | 11218.8 | 1776.6 | 614.16 |
| `fv_transport_faces` | 3 | 18214.6 | 938.3 | 590.52 |
| `fv_vector_gradient` | 2 | 14128.4 | 898.2 | 345.44 |
| `fv_velocity_boundary` | 2 | 1032.2 | 196.4 | 293.76 |

A fresh `bench_bandwidth` run on the same source measures Metal triad 241.718
GB/s (best of seven, 256 MiB/array), add 241.658, copy 161.206, dot 238.897 and
sum 220.768 GB/s. Serial CPU triad is 126.982 GB/s. Metal single-dispatch
synchronization costs 90.96 µs; batched small dispatches cost 1.139 µs each.

## Test tiers and reproduction

The requested CMake configure and build succeed. Fresh numerical verification
uses source `f57c891`; the later `10b9a4f` only restores four existing LSQ
assertions and corrects a comment. A fresh rebuild, focused operator run and
quick suite pass with those assertions restored. No numerical code changed
after the full run began.

| Run | Cases passed / run | Failed assertions | Wall time |
| --- | ---: | ---: | ---: |
| Original baseline | 137 / 138 | 2, former tet pressure gate | 266.23 s |
| Final quick, excludes slow | 135 / 135 | 0 | 12.36 s |
| Full verification | 141 / 142 | 6, retained tet force gates | 2,890.75 s, 48m 10.75s |

The full run has exactly one failed case, V2 tet refinement: wall pressure,
viscous and total-force order each fail on CPU and Metal. Every other case
passes, including V1, both V3 studies, V4 and the Metal runtime/kernel checks.
The focused final operator run passes all 12 assertions; the final quick run
passes all 7,248,706 assertions.

| Slow case | Full-run wall time |
| --- | ---: |
| V1 scalar refinements | 107.476 s |
| V2 hex refinement | 44.318 s |
| V2 tet refinement | 112.646 s |
| Couette | 5.796 s |
| V3 Re=1000 refinement, both backends | 2,241.925 s |
| V3 Re=100/400 reports, both backends | 174.296 s |
| V4 imported/matched grids, both backends | 191.430 s |

The longest quick-tier case in the full run is external Fluent mesh acceptance
at 4.981 s. Every measured case over five seconds is in `slow`. The added
cavity studies account for the increase over the original ~270-second suite.

The `slow` doctest suite contains V1/V2 refinements, V3/V4 cavities and Couette
(the baseline Couette test exceeded five seconds). All other measured baseline
tests took <5 seconds. CI excludes `slow` on every job and retains its existing
GPU-presence exclusion. Full verification has no suite exclusion and continues
to return failure while the three tet force gates fail on each backend.

```sh
cmake -S . -B build -G Ninja && cmake --build build && ./build/kutta-tests
# Quick loop / CI (CI also excludes its GPU-presence test):
./build/kutta-tests --test-suite-exclude=slow
# Focused studies:
./build/kutta-tests --test-case='FV MMS pressure*'
./build/kutta-tests --test-case='FV V2*'
./build/kutta-tests --test-case='FV V3*,FV V4*'
# Optional standalone CSV checker, after generating all three grids:
./build/bench_fv --cavity 1000 --nx 65 --ny 65 --iterations 10000
./build/bench_fv --cavity 1000 --nx 129 --ny 129 --iterations 10000
./build/bench_fv --cavity 1000 --nx 257 --ny 257 --iterations 10000
python3 validation/cases/check_cavity_profiles.py
```

The generated Metal source is an array of raw-string parts, each ≤60,000 UTF-8
bytes, joined at runtime before appending the reduction kernels. Adjacent
literals would still produce one overlength literal. The current parts are 60,000 and 9,893 bytes. A fresh build emits no
`-Woverlength-strings` warning and the generated C++ compiles with
`-Werror=overlength-strings`. Generated kernel/transfer tests exercise CPU and Metal successfully.

Fresh evidence is in `build/baseline-tests.log`, `build/operator-tests.log`,
`build/v2-tests.log`, `build/quick-tests.log`, `build/full-tests.log`,
`build/quick-tests-final.log`, `build/operator-tests-final.log` and
`build/final-build.log`. These are local artifacts, not committed fixtures.
Only local commits are made; nothing is pushed.

## V6 mixing-elbow SST

**The requested same-settings accuracy gate is not achieved.** SST transport and
momentum are coupled, but the strict adapter rejects the original elbow's
`recon/cell-lsf?` (warped-face gradient correction, WFGC). Its hybrid initializer
is also rejected by `initial_fields()`. The slow test keeps both requirements
failing. The numbers below are an explicitly labelled **diagnostic**: WFGC is
disabled and cold start uses the case's stored standard initialization values.
They are not results for Fluent's complete original settings.

Measured 2026-10-29 in worktree `t3-265ec8e3`, branch
`t3/couple-sst-steady-solver`, on Apple M4 Pro / 24 GiB, serial CPU and Metal.
Numerical implementation is `bbc176c`, setup `c70561a`, stopping criterion
`764fa63`; subsequent changes add diagnostics, input guards and this report.
macOS HDF5 is 2.2.0. Float fields, double/compensated reductions, deterministic
gather kernels, and the existing test thresholds are retained.

External files (not committed), from `KUTTA_FLUENT_EXAMPLES/mixing_elbow`:

| File | SHA-256 |
| --- | --- |
| `mixing_elbow.cas.h5` | `8145813f50f989a521ed8d1b02c891c7ca22ccc5fafa35a4eeba011949f19cbc` |
| `mixing_elbow.dat.h5` | `863a813227cf4274a9ca24fdad43e6c1f9f132dd62fcabc39a3be8cf493cff42` |

The saved solver version is Fluent **22.2**. The mesh has 17,822 cells, 91,581
faces, two velocity inlets, a pressure outlet, symmetry, and 2,436 wall faces.
There are **zero report definitions**. The saved iteration counter is **absent**;
the last retained residual sample is iteration **294**, with 294 samples. Do not
substitute 294 for a separately verified save counter. Fluent wall time is absent.

Active material `air` has constant density **1.225 kg/m³** and dynamic viscosity
**1.7894e-5 Pa s**; gravity is off. Although energy is enabled, neither density
nor viscosity depends on temperature, and there are no other enabled couplings
that feed energy back into these incompressible flow equations. Omitting energy
does not change this case's converged flow.

### Settings, classifications and remaining limits

Every item originally flagged by `kutta inspect` is accounted for:

| Setting | Effect and evidence | Adapter action |
| --- | --- | --- |
| Missing `turbulence/wall-omega-treatment` | Changes wall equations if a different model is selected. The case explicitly has legacy `kw-wall-omega-treatment-r13? #t`; the saved-field forensics establish this correlation's wall-cell omega, shear and y+ on this exact case. | Use the established legacy correlation; no version-independent default. Reject tabulated treatment. |
| `kw-wall-omega-treatment-r13? #t` | Discrete wall closure; forensics §3/§7 and the saved-field SST test are the evidence. | Implemented by `SstModel`. |
| `kw-dissipation-newton? #t` | Convergence path: Newton versus Picard linearization of destruction has the same source zero. | Picard implicit sink; recorded in setup notes. |
| `sst-des-iddes-simplify_omit_fe? #t` | Inactive: case configuration explicitly has `rp-des? #f`. | Accept only the globally disabled DES problem. |
| `sst-linearize-neg-cross-diff? #t` | Linearization path; source sign still changes the equation. | Signed cross diffusion, negative part as an implicit sink, implemented in SST. |
| Informational k-epsilon flags | Inactive in the explicitly selected k-omega SST model. | Parser's model allowlist retains their provenance. |
| Energy | Decoupled constant-property flow as verified above. | Record omission. |

Other consequential case settings:

- LSQ cell gradients, standard cell-to-face limiter, SOU momentum/k/omega,
  second-order pressure, momentum-based Rhie-Chow: represented. WFGC is **on**
  (`recon/cell-lsf? #t`): this changes gradients and hence the converged equations.
  The [Fluent WFGC documentation](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_sec_solve_conv_stab.html)
  describes its purpose but gives no reproducible algorithm. The adapter errors
  rather than silently disabling it.
- Coupled is replaced by SIMPLEC; implicit relaxation uses the recorded numerical
  factors (momentum/pressure .5, k/omega .75). Conservative automatic global
  pseudo-time uses scale 1. These substitutions affect the convergence path.
  Rhie-Chow uses the trace/3 of physical component momentum diagonals, including
  wall drag and excluding relaxation/pseudo-time. The scalar trace choice is
  rotationally invariant; exact equivalence to Fluent's undocumented high-order
  Rhie-Chow correction remains unverified.
- `wall-distance/method 2`, quad correction on: finite-surface BVH wall distance
  and face-normal wall-function distances are used. At the saved fields, wall
  distance relative L2 error is **2.09324e-5** (0.00209324%), max absolute
  **1.16476e-5 m**. This is much smaller than the eddy-viscosity discrepancy.
- Backflow pressure code 0 is total pressure: on backflow faces,
  `p=p0-rho*|u_cell|²/2`, using the adjacent cell velocity (UG 8.4.9.1.3 and 8-88).
  Outflow remains static pressure. The [boundary-condition chapter](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_ug/flu_ug_bcs_sec_bound_cond.html)
  supplies the formulas; a CPU/Metal regression checks a nonzero gauge-pressure
  case. Backflow direction is normal and turbulence is intensity .05 / viscosity
  ratio 4. Prevent-reverse-flow is explicitly off.
- Inlets use `k=1.5*(U*I)²`, `ell=.07*Dh`,
  `omega=sqrt(k)/(.09^(1/4)*ell)`; viscosity-ratio specification instead uses
  `omega=rho*k/(mu*ratio)`. Cold inlet: speed 1, Dh .1016 m, k .00375,
  omega 15.72038792. Hot inlet: speed 1.2, Dh .0254 m, k .0054,
  omega 75.45786201. Both exactly match saved inlet fields. UG 8-73/75/79/80
  and the SVG exponent discrepancy are documented in the research report's
  authoritative corrections box.
- The case has `hyb-init? #t`. [Hybrid initialization](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_ug/flu_ug_sec_init_hybrid.html)
  solves auxiliary Laplace equations and patches other variables from recipes.
  It changes the starting fields and iteration count, not the converged equations.
  This initializer is not implemented. The diagnostic cold proxy uses stored
  uniform defaults: velocity (.3771192189, .06864234338, 2.999943839e-17) m/s,
  pressure 8.731149137e-11 Pa, k 5.657097174, omega 597.1623783. Those defaults
  are not described as the hybrid result.

### Fixed-point residual check

Initialization maps source zone IDs through the solver's retained RCM permutation.
It primes physical momentum coefficients/face fluxes without advancing the fields.
`monitor()` does not change the seeded u/p/k/omega. Momentum residuals use the
physical diagonal, turbulence residuals use their own physical diagonals and omit
fixed wall-cell omega rows. Fluent's coupled/pseudo-time normalization can differ;
these numbers compare observed residuals, not equivalent normalization operators.
Continuity uses each Kutta run's first-five-iteration normalization, learned after
its initial snapshot. Fluent's scale is .002614088497 kg/s.

| Equation | Fluent final scaled | Kutta CPU at saved fields | Kutta Metal at saved fields |
| --- | ---: | ---: | ---: |
| x-velocity | 1.39741825e-06 | 0.00116463854 | 0.00116463849 |
| y-velocity | 3.88742419e-06 | 0.00204521918 | 0.00204521912 |
| z-velocity | 1.00155813e-06 | 0.000516084586 | 0.000516084586 |
| k | 5.16072107e-05 | 0.0118328971 | 0.0118328977 |
| omega | 2.10806917e-05 | 0.0306121987 | 0.0306121985 |
| continuity | 0.000947130927 | 0.0393532417 | 0.0393531944 |

Initial Kutta continuity L1 is .002068248363 kg/s on CPU and .002068248380 kg/s
on Metal, versus Fluent's final 2.475884062e-6 kg/s. At those same saved fields,
Kutta's eddy-viscosity relative L2 discrepancy is **3.40282%**, max absolute
**5.09735e-5 Pa s**. The existing independent SST saved-field test already measures
scatter from inverse-distance-squared reconstruction in strain-limited cells.
Thus the saved Fluent state is observably not a fixed point of these reconstructed
coefficients, before taking any solver step.

### Iteration counts, timings and field differences

Stopping requires every scaled steady momentum/k/omega defect <1e-6 and
provisional continuity L1 below Fluent's recorded final unscaled continuity,
2.475884062e-6 kg/s. This uses the reference run's physical mass scale; an earlier
arbitrary 1e-7 kg/s target ran into float-level flux noise despite sub-1e-7 steady
transport defects. The accuracy thresholds below were never changed. Monitors
run every 50 iterations; time includes iteration/monitor work, excludes setup,
input and force post-processing. Docker's emulated CUDA build was running
concurrently; these are observed wall times, not an isolated speed benchmark.
There is no measured speedup against Fluent.

| Backend | Start | Outer iterations | Wall time (s) | Velocity L2 (%) | Pressure L2 (%) | k L2 (%) | omega L2 (%) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cpu | seeded | 300 | 35.7029 | 3.925067 | 3.147520 | 15.618063 | 8.596569 |
| cpu | cold-standard-proxy | 600 | 83.6633 | 3.926663 | 3.145914 | 15.620213 | 8.596902 |
| metal | seeded | 300 | 11.2184 | 3.925066 | 3.147465 | 15.618062 | 8.596579 |
| metal | cold-standard-proxy | 600 | 20.2913 | 3.926664 | 3.145997 | 15.620217 | 8.596901 |

Unweighted cell relative L2 is `sqrt(sum(|q-q_ref|²)/sum(|q_ref|²))`; velocity
uses vector norms. Maximum is an absolute cell error in the field's SI units;
relative maximum divides that maximum by the reference field's maximum magnitude,
not by a potentially zero local reference value.

| Backend / start | Field | Max absolute difference | Relative maximum (%) |
| --- | --- | ---: | ---: |
| cpu / seeded | velocity | 0.443521633 | 26.711847 |
| cpu / seeded | pressure | 0.121506399 | 10.601535 |
| cpu / seeded | k | 0.01788308 | 63.397187 |
| cpu / seeded | omega | 10826.0641 | 12.012038 |
| cpu / cold-standard-proxy | velocity | 0.443568948 | 26.714696 |
| cpu / cold-standard-proxy | pressure | 0.121516323 | 10.602401 |
| cpu / cold-standard-proxy | k | 0.0178871275 | 63.411536 |
| cpu / cold-standard-proxy | omega | 10826.2565 | 12.012251 |
| metal / seeded | velocity | 0.443521575 | 26.711843 |
| metal / seeded | pressure | 0.121507532 | 10.601634 |
| metal / seeded | k | 0.0178830818 | 63.397193 |
| metal / seeded | omega | 10826.0837 | 12.012059 |
| metal / cold-standard-proxy | velocity | 0.443569083 | 26.714705 |
| metal / cold-standard-proxy | pressure | 0.121512196 | 10.602041 |
| metal / cold-standard-proxy | k | 0.0178871312 | 63.411549 |
| metal / cold-standard-proxy | omega | 10826.2878 | 12.012286 |

### Pressure drop and wall forces

Pressure drop is the area-weighted static pressure over both inlets minus the
area-weighted outlet static pressure. Reference cold/hot/outlet pressures are
.5016180027 / 1.087188619 / 0 Pa. Reference combined drop is **.5351137855 Pa**.

| Backend / start | Cold inlet (Pa) | Hot inlet (Pa) | Outlet (Pa) | Drop (Pa) | Drop error (%) |
| --- | ---: | ---: | ---: | ---: | ---: |
| cpu / seeded | 0.510110193 | 1.11442965 | 0 | 0.544678447 | 1.787407 |
| cpu / cold-standard-proxy | 0.510144157 | 1.11448749 | 0 | 0.544713776 | 1.794009 |
| metal / seeded | 0.510109738 | 1.11442894 | 0 | 0.544677976 | 1.787319 |
| metal / cold-standard-proxy | 0.510140742 | 1.11448362 | 0 | 0.544710335 | 1.793366 |

Forces below are on the body, in N. Saved `SV_WALL_SHEAR` is force on the fluid;
its sign is reversed for body shear. Saved face `SV_P` is integrated with the
same mesh's double face areas, gauge/reference pressure zero. Kutta uses its
SST wall shear, not a molecular-viscosity gradient approximation.

| Reference zone | Part | Fx (N) | Fy (N) | Fz (N) |
| --- | --- | ---: | ---: | ---: |
| wall-elbow | pressure | 0.006362455524 | -0.005934901867 | 0.01552523764 |
| wall-elbow | shear | 0.0003481685928 | 0.0006880408909 | -6.734241096e-05 |
| wall-elbow | total | 0.006710624117 | -0.005246860976 | 0.01545789523 |
| wall-inlet | pressure | 7.988724203e-05 | -4.692876481e-06 | 0.001612299233 |
| wall-inlet | shear | 1.013114139e-06 | 7.091390076e-05 | 4.221385465e-07 |
| wall-inlet | total | 8.090035617e-05 | 6.622102428e-05 | 0.001612721372 |

| Backend / start | Zone | Part | Fx (N) | Fy (N) | Fz (N) | Vector error (%) |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| cpu / seeded | wall-elbow | pressure | 0.006429916884 | -0.005949214478 | 0.01528421649 | 1.408618 |
| cpu / seeded | wall-elbow | shear | 0.0003327747477 | 0.0006874346094 | -5.847051166e-05 | 2.296713 |
| cpu / seeded | wall-elbow | total | 0.006762691632 | -0.005261779868 | 0.01522574598 | 1.350647 |
| cpu / seeded | wall-inlet | pressure | 7.262223624e-05 | -5.181007285e-06 | 0.001601499837 | 0.806847 |
| cpu / seeded | wall-inlet | shear | 1.3234273e-06 | 7.261430428e-05 | 4.293509314e-07 | 2.437173 |
| cpu / seeded | wall-inlet | total | 7.394566354e-05 | 6.7433297e-05 | 0.001601929188 | 0.797971 |
| cpu / cold-standard-proxy | wall-elbow | pressure | 0.006430119715 | -0.005949366355 | 0.01528612308 | 1.398680 |
| cpu / cold-standard-proxy | wall-elbow | shear | 0.0003327702131 | 0.0006874739025 | -5.847263314e-05 | 2.296916 |
| cpu / cold-standard-proxy | wall-elbow | total | 0.006762889928 | -0.005261892453 | 0.01522765045 | 1.340429 |
| cpu / cold-standard-proxy | wall-inlet | pressure | 7.262826666e-05 | -5.181260033e-06 | 0.001601611789 | 0.800897 |
| cpu / cold-standard-proxy | wall-inlet | shear | 1.323552797e-06 | 7.261475904e-05 | 4.294074343e-07 | 2.437836 |
| cpu / cold-standard-proxy | wall-inlet | total | 7.395181945e-05 | 6.7433499e-05 | 0.001602041197 | 0.791975 |
| metal / seeded | wall-elbow | pressure | 0.006429920191 | -0.005949211152 | 0.01528422088 | 1.408598 |
| metal / seeded | wall-elbow | shear | 0.0003327746459 | 0.0006874345954 | -5.847054399e-05 | 2.296723 |
| metal / seeded | wall-elbow | total | 0.006762694837 | -0.005261776557 | 0.01522575033 | 1.350626 |
| metal / seeded | wall-inlet | pressure | 7.262227842e-05 | -5.181003054e-06 | 0.001601498779 | 0.806900 |
| metal / seeded | wall-inlet | shear | 1.323425071e-06 | 7.261429323e-05 | 4.293473738e-07 | 2.437158 |
| metal / seeded | wall-inlet | total | 7.394570349e-05 | 6.743329018e-05 | 0.001601928127 | 0.798025 |
| metal / cold-standard-proxy | wall-elbow | pressure | 0.006430107648 | -0.005949365789 | 0.01528597329 | 1.399470 |
| metal / cold-standard-proxy | wall-elbow | shear | 0.0003327700705 | 0.0006874736837 | -5.847260683e-05 | 2.296935 |
| metal / cold-standard-proxy | wall-elbow | total | 0.006762877719 | -0.005261892106 | 0.01522750068 | 1.341240 |
| metal / cold-standard-proxy | wall-inlet | pressure | 7.262727542e-05 | -5.181244506e-06 | 0.001601602589 | 0.801403 |
| metal / cold-standard-proxy | wall-inlet | shear | 1.323563213e-06 | 7.261472454e-05 | 4.294127328e-07 | 2.437791 |
| metal / cold-standard-proxy | wall-inlet | total | 7.395083863e-05 | 6.743348003e-05 | 0.001602032001 | 0.792482 |

### Gate outcomes and diagnosis

| Gate | Required | Observed | Outcome |
| --- | --- | --- | --- |
| Original case settings / initializer | Both represented | WFGC and hybrid initialization explicitly rejected | FAIL |
| Cell velocity relative L2 | <0.5% | 3.9251–3.9267% | FAIL, all four runs |
| Pressure drop | Within 0.5% | 1.7873–1.7940% | FAIL, all four runs |
| Integrated wall shear | Within 0.5% per zone | Elbow ~2.297%; inlet ~2.437–2.438% | FAIL, both zones/all runs |
| Integrated total wall force | Within 0.5% per zone | Elbow ~1.34–1.35%; inlet ~0.792–0.798% | FAIL, both zones/all runs |
| Fixed-point drift, max of u/p/k/omega relative L2 | <0.1% | 15.6181%, dominated by k | FAIL, both seeded runs |
| Diagnostic convergence | Steady defects <1e-6 and reference mass target | All four runs satisfy | PASS for the diagnostic equations |

The slow case exits nonzero with **28 retained failed assertions**. There are no
threshold exceptions, expected-failure decorators or skipped accuracy assertions.

A concrete discrepancy is reconstruction-dependent eddy viscosity: already at
the saved state mu_t differs by 3.40% L2 and k/omega residuals are .0118/.0306,
far above Fluent's final residuals. Wall distance has a much smaller global L2
error; its isolated nonlinear effect has not been measured.
Missing WFGC is a known spatial mismatch, and inverse-distance LSQ's effect on
strain-limited viscosity was independently established by the existing SST test.
The contribution of WFGC alone has **not** been isolated: no Fluent executable or
same-case WFGC-off reference was available. It would be unsupported to attribute
all drift to that one switch.

Other unresolved discrete details remain the wall-face omega reconstruction's
known ~0.42% systematic offset, unverified wall-cell production blending, and
Fluent's unpublished high-order Rhie-Chow correction. The treatment of the
isotropic `-2*rho*k/3` momentum stress is also unresolved in research A1; no
physical-pressure equivalence is claimed for it here. The data file contains
primary fields, not the gradients/source terms needed to settle these separately.
The constitutive options retain their documented/inferred provenance in
`turbulence/README.md`; the production limiter is inactive in the saved fields,
which does not verify its form during a cold solve. Hybrid initialization explains
only why the requested cold iteration-count comparison remains unverified.

Two implementation defects were resolved during verification: extrapolated
boundary feedback can preserve seeded gradients in deficient tet stencils (the
new constant-field CPU/Metal regression reproduces and removes it), and unscaled
identity wall-omega rows swamp the inexact Krylov residual (retaining the physical
row scale removes the artificial ~1e-3 interior-omega residual plateau). Neither
fix was made by relaxing an accuracy gate.

### Reproduction and repository checks

```sh
export KUTTA_FLUENT_EXAMPLES=/Users/ezzy/.Others/Kutta/validation/fluent_examples
cmake -S . -B build -G Ninja
cmake --build build
./build/kutta inspect "$KUTTA_FLUENT_EXAMPLES/mixing_elbow/mixing_elbow.cas.h5" \
  --data "$KUTTA_FLUENT_EXAMPLES/mixing_elbow/mixing_elbow.dat.h5"
./build/kutta-tests --test-suite-exclude=slow
ctest --test-dir build -LE slow --output-on-failure
# Expected nonzero exit while the accuracy/settings gates above remain unmet:
KUTTA_V6_ITERATIONS=1200 ./build/kutta-tests --test-case='FV V6*'
# Optional single backend; omitting it runs every compiled/present backend:
KUTTA_V6_BACKEND=metal ./build/kutta-tests --test-case='FV V6*'
```

Fixtures are located through `KUTTA_FLUENT_EXAMPLES` or `KUTTA_SOURCE_DIR`; absent
HDF5 or the external case/data pair produces a clean skip. The independent SST
saved-field test now honors the same environment variable. No Ansys data is added
to git. No push, deployment or remote-resource change is performed.

Final checks at code commit `94cdeca` (the final documentation-only commit
follows):

| Check | Result |
| --- | --- |
| Fresh macOS CMake/Ninja build | PASS |
| Default quick tier, CPU + Metal | 146/146 cases, 1,415,147 assertions; **8.82 s** |
| `ctest -LE slow`, external fixture environment enabled | 2/2 entries pass, including portability lint; doctest 17.74 s during the emulated CUDA build (earlier fixture-enabled run 13.11 s) |
| V6 slow case, CPU + Metal | Expected FAIL: 28 accuracy/settings assertions, 2,601 total assertions |
| CUDA 13.0.2 / Ubuntu 24.04 / HDF5 1.10.10 | Library, generated kernels and full test executable compile for `75-real;86-real;89-real;120` |
| CUDA build's CPU-oracle quick tier | 145/145 cases, 764,034 assertions pass; GPU-presence test excluded as in CI |
| Linux portability lint | PASS |
| CUDA GPU execution | UNVERIFIED: no NVIDIA device exposed |

CUDA configuration uses the exact CI command:

```sh
cmake -S /src -B /build/compile -G Ninja -DKUTTA_ENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES="75-real;86-real;89-real;120"
cmake --build /build/compile -j 6
/build/compile/kutta-tests --test-suite-exclude=slow \
  --test-case-exclude='compiled GPU backends*'
python3 tools/lint_portability.py
```

The container starts from `nvidia/cuda:13.0.2-devel-ubuntu24.04`, with CI's apt
packages `cmake ninja-build python3 libhdf5-dev zlib1g-dev git`. A temporary local
image cached those packages for the final incremental current-source check;
this does not change the CUDA/HDF5 toolchain. The original container's in-flight
build predated the last input guard and was superseded by that successful final
build, rather than being used as current-source proof.
Local evidence: `build/v6-evidence.log`, `build/v6-results.json`,
`build/ctest-final.log`, `build/quick-final-default.log`,
`build/quick-with-fixture.log`, `build/input-tests.log`, and
`build/cuda-evidence/{build,tests,lint}.log`. The complete V1–V4 slow ladder was not rerun: prior
retained failures above remain unresolved, and no claim is made that these
laminar long runs have been revalidated after the boundary-gradient change.

## Laminar regression repair after SST coupling

Code commit `bf4f2f5` restores the verified laminar reconstruction when SST is
disabled. `bbc176c` had enabled `boundary_least_squares()` for every flow solve,
omitting extrapolated wall/outlet/symmetry observations as appropriate to each
field. That changes the pressure-face reconstruction, Gauss momentum source,
Rhie–Chow pressure jumps, velocity reconstruction and pressure-correction
gradient. It is a discretization change, not just an SST initialization change.
The wall-momentum and effective-viscosity hooks already return without changing
laminar flow. This repair keeps the boundary-aware reconstruction for SST and
the full original stencil for laminar flow, including wall-force reconstruction.
It also copies the common laminar momentum diagonal/inverse directly instead of
rounding it through the SST trace/3 average. No Fluent-matching correction is
being introduced, and no existing threshold is changed.

The original `fbb5111` test binary reproduced n=4/8 nonconvergence after 2,000
steps (`u_L2=0.564466/1.30634`) and the n=16 AMG invalid-diagonal exception.
Controlled 384-cell tet experiments on CPU and Metal isolated the change:

| Reconstruction restored to main | Steps | Velocity L2 | Outcome |
| --- | ---: | ---: | --- |
| Pressure correction only | 200 (limit) | 0.566300 | Still fails convergence |
| Pressure field only | 140 | 0.226447 | Converges to a degraded solution |
| Full laminar gradient path | 160 | 0.177223 | Reference accuracy restored |

`FV laminar tet manufactured solution converges without SST` reuses the V2
analytical solution on 4³ bricks / 384 tetrahedra. It retains the existing V2
residual criterion, requires convergence within 200 steps, and checks velocity
L2 < .18 and gauge-adjusted pressure L2 < .28 on every available backend. The
pre-fix run fails both convergence and velocity accuracy on CPU and Metal in
2.88 s. The slow V2 studies retain their 2,000-step limit and all existing gates.

A fresh isolated archive/build of main `4e98429` supplies a direct V2 reference;
only the new test harness is copied into that archive. A separate field
comparison after 160 steps checks all 2,400 floats in u/v/w/p/face mass flux for
the 384-cell tet: fixed versus main is **bitwise equal on CPU and Metal**.
The reference's full V2 runs converge at hex n=8/16/32 in 140/320/840 steps and
tet n=4/8/16 in 160/220/420 steps, with only the six known tet wall-force-order
failures. The historical accuracy/iteration tables above predate some changes;
the direct main comparison is authoritative for this regression.

### Complete post-repair verification

The rebuilt Release binary at code commit `bf4f2f5` has SHA-256
`09d26f3ce46777a0a08b136416b3b1a02e7134a5fc0cb7f421e29da6088b54f3`.
Runs use this worktree on macOS 26.6.2 / arm64, CPU and Metal, with HDF5 enabled.
Both tiers use
`KUTTA_FLUENT_EXAMPLES=/Users/ezzy/.Others/Kutta/validation/fluent_examples`;
V6 and the independent saved-field SST quick check actually run.

| Check | Result | Wall time |
| --- | --- | ---: |
| `cmake --build build` | PASS, current source rebuilt | — |
| `python3 tools/lint_portability.py` | PASS | — |
| `./build/kutta-tests --test-suite-exclude=slow -d` | **147/147 cases; 1,449,526/1,449,526 assertions pass** | **8.62 s** |
| New 384-cell tet quick regression, both backends | PASS, 6 assertions | 3.09 s (included above) |
| `./build/kutta-tests -ts=slow -d` | **6/8 cases pass; 34 retained failures in 3,476 assertions** | **1,020.69 s** |

Every slow case ran; none skipped. The nonzero slow exit is exactly 28 V6
assertions and six tet force-order assertions. There are no convergence,
velocity/pressure-order, AMG or other new failures.

| Slow case | Main/reference outcome | Fixed outcome, CPU and Metal | Case time |
| --- | --- | --- | ---: |
| V1 scalar refinements | PASS; diffusion/FOU/SOU gates 1.9/.8/1.8 | PASS; finest-pair orders below | 68.08 s |
| V2 hex | PASS; fresh main 140/320/840 steps | PASS; same steps and printed solution/load values | 28.37 s |
| V2 tet | FAIL, six wall-force-order gates; fresh main 160/220/420 steps | Same six failures only; same steps and printed solution/load values | 71.23 s |
| Couette | PASS; documented max error 1.19209e-7, residual 4.04177e-8 | PASS; 400 steps, max error 1.78814e-7, residual 3.1049e-8, both backends | 4.60 s |
| V3 Re=1000 refinement | PASS; updated CPU counts 700/600/1000 | PASS; CPU 700/600/1000, Metal 700/600/900; finest errors/orders below | 581.09 s |
| V3 Re=100/400 | PASS; updated counts 1000/500 | PASS; same counts; max u/v deviations .005038/.008966 and .002639/.004047 | 78.22 s |
| V4 imported/matched Fluent cavity | PASS; updated Re=100/400/1000 counts 600/300/500 | PASS; same counts, every matched-grid sample within the existing 2e-6 gate | 73.26 s |
| V6 elbow (not present on main) | Prior SST branch: FAIL, 28 settings/accuracy gates | Same 28 failures; all recorded numerical diagnostics unchanged to printed precision | 115.76 s |

Except for the fresh V2 comparison, the reference numbers here are the
historical measurements retained above, not a second complete main slow run.
Small differences from those tables are reported explicitly rather than
claimed as bitwise agreement. The Couette max-error difference is 5.96e-8.

V1 finest-pair orders (diffusion / FOU / SOU):

| Mesh / backend | Documented reference | Fixed |
| --- | --- | --- |
| Hex CPU | 2.00206 / .854671 / 1.91201 | 2.002060 / .854671 / 1.912087 |
| Hex Metal | 2.00192 / .854671 / 1.91201 | 2.002005 / .854671 / 1.912087 |
| Tet CPU | 1.93250 / .993216 / 1.94919 | 1.931832 / .993209 / 1.949123 |
| Tet Metal | 1.93243 / .993216 / 1.94917 | 1.931099 / .993197 / 1.949137 |

V2 comparison checks all 12 solutions and 204 printed numeric values, including
velocity/pressure errors, forces and orders, against fresh main within the
displayed rounding. Fresh-main finest hex velocity/pressure errors are
.00418355/.00301253; the fixed slow log rounds them to .004184/.003013. Both
tet backends report .012256/.068933. Tet finest velocity/pressure orders are
1.982589/.942388; retained pressure/viscous/total wall-force orders are
.996553/.874871/.933505, all identical to fresh main at the printed precision.

V3 finest Re=1000 comparisons against Erturk:

| Backend | Documented u / v error | Fixed u / v error | Fixed finest u / v order |
| --- | --- | --- | --- |
| CPU | .000607678559 / .001067561871 | .000613706079 / .001062930888 | 2.839101 / 2.755230 |
| Metal | .000610226658 / .001065764046 | .000614044335 / .001062692470 | 2.838302 / 2.755555 |

All errors decrease on refinement and both orders exceed the unchanged 1.6
gate. V4 imported-mesh Ghia max u/v deviations at Re=100/400/1000 are
.005306/.009046, .003285/.003434 and .008948/.004081 on both backends, versus
the historical .00530634/.00904564, .00328528/.00343363 and .00894912/.00409038.

V6 comparison covers all four seeded/cold CPU/Metal iteration counts, all 16
field relative-L2 errors, all 24 wall-force vectors and relative errors, and all
four pressure drops. Every comparison matches the preceding V6 report within
its printed precision. Seeded/cold counts stay 300/600 on both backends;
velocity L2 remains 3.925066–3.926664%, and pressure-drop error remains
1.787319–1.794009%. Only measured wall times changed. No V6 discretization,
initializer, stopping criterion or accuracy threshold was edited.

Local evidence: `build/laminar-evidence-identity.txt`,
`build/laminar-tet-before.log`, `build/laminar-tet-quick-before.log`,
`build/laminar-tet-correction-only.log`, `build/laminar-tet-pressure-only.log`,
`build/laminar-tet-quick-after.log`, `build/laminar-field-compare.log`,
`build/laminar-main-v2.log`, `build/laminar-v2-compare.log`,
`build/laminar-v6-compare.log`, `build/laminar-couette-detail.log`,
`build/laminar-slow-final.log`, `build/laminar-quick-final.log`,
`build/laminar-build-final.log` and `build/laminar-portability-final.log`.
Temporary reference sources/builds and the field-comparison harness were
removed after retaining the results and field hashes. Commits are local only.
