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
