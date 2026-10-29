# Lane D verification

The implementation is committed locally, but the validation ladder is not fully
passed. V2 tetrahedral pressure remains first order, V3 fails a published-data
scatter check, and the array-size traffic estimate exceeds the ~2 KB target.
No acceptance threshold was relaxed.

Numerical source: `5e527dc856280e6e69f436f0fc7438b5340a903d`, branch
`t3/implement-fv-steady-solver`, based on `c517bb6`. Measurements were taken on
2026-10-29, Apple M4 Pro (24 GB), macOS 26.6.2, AppleClang 21, C++20, Release
(`-O3 -DNDEBUG`), HDF5 and Metal enabled. CPU means the serial reference backend.
The subsequent report commit changes documentation and the case-data checker.

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
| Hex, CPU and Metal | 1.91521 | 1.82043 | 0.00418213 | 0.00300986 | 860 |
| Tet, CPU and Metal | 1.98241 | **0.942386** | 0.0122576 | 0.0689337 | 520 |

The tet pressure assertion **≥1.8 remains failed** on both backends. Removing
boundary-adjacent cells gives finest pressure L2 0.0670728; it does not restore
second order. The maximum pressure error, 0.212136, is at a wall cell. Standalone
LSQ reconstruction of exact MMS pressure has tet order 0.994008. Momentum's
Gauss pressure force and Rhie–Chow's LSQ reconstruction differ, but replacing the
latter with the former worsened both orders. A wall zero-normal pressure
extrapolation experiment also worsened them. Both experiments were removed.
The converged operator retains the specified LSQ gradients, Second Order
pressure faces, momentum-based Majumdar flux and over-relaxed non-orthogonal
diffusion. The added pressure-update cap changes the path only. Root cause is
unresolved; the bounded investigation and rejected experiments are in
[V2-DIAGNOSIS](V2-DIAGNOSIS.md).

## V3 cavity at 129×129 cells

All three Reynolds numbers converge on CPU and Metal, with identical printed
results. The case uses unit lid speed, viscosity 1/Re, symmetry front/back,
momentum/pressure relaxation limits 1, pseudo-time scale 5, and monitors every
100 iterations. Convergence requires scaled momentum defects <5e-7, velocity
change over the last 100 iterations <2e-5 and corrected mass L1 <1e-5.
These are convergence tests; published-data accuracy is checked separately.

| Re | Outer iterations | Ghia max abs u deviation | Ghia max abs v deviation | Erturk max abs u deviation | Erturk max abs v deviation |
| --- | ---: | ---: | ---: | ---: | ---: |
| 100 | 400 | 0.00503866 | 0.00896686 | — | — |
| 400 | 800 | 0.00263811 | 0.00404582 | — | — |
| 1000 | 1300 | 0.00594781 | 0.0100679 | 0.00434731 | 0.00711588 |

Only the supplied SOURCE.md's suspect Ghia Re=400 v point at x=0.90625 is
excluded. All other points are retained. No arbitrary tolerance is described
as data scatter. `validation/cases/check_cavity_profiles.py` compares the
Re=1000 coordinates common to both papers, including their printed-value
rounding intervals. The center values are u=-0.0624667927623 and v=0.025315169245;
u lies outside the interval [-0.06205, -0.060795] on both backends, while v passes.
The checker reports 10 passing and 2 failing common-point checks.
The supplied Re=100/400 data do not contain an independent scatter band.
Thus V3's stated accuracy criterion is **not passed**, despite convergence.

## V4 real Fluent cavity mesh

The exact file read in place is
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

## Final checks, scope and reproduction

The requested CMake configure/build succeeds. The complete test binary reports
**76 cases: 74 passed, 2 failed; 7,168,742 assertions: 3 failed**. Two failed
assertions are the V2 tet pressure gate (one per backend). The other failure is
`tests/io/test_fluent_real.cpp:108`: the external inventory has an extra
`LICENSE-ansys-example-data`, with all 12 expected files present. The extra file
and I/O test were left untouched, outside Lane D ownership. Every other FV
case passes. The cavity common-point checker separately returns 1 for the
retained Re=1000 u scatter failure. Bench cavity exit 0 means convergence only.

Only the authorized FV, FV-kernel, FV-test, benchmark and case paths changed.
There are no declarations added to `device_mesh.hpp`, and no changes to runtime,
linear algebra, mesh, I/O or turbulence. Solver changes from the concurrent
linear-algebra branch are not integrated. CUDA compilation/runtime was not
verified in this lane. Commits are local; nothing was pushed.

```sh
cmake -S . -B build -G Ninja && cmake --build build && ./build/kutta-tests
./build/kutta-tests --test-case='FV V2*'
./build/kutta-tests --test-case='FV MMS pressure*'
./build/bench_fv --cavity 100 --nx 129 --ny 129 --iterations 6000
./build/bench_fv --cavity 400 --nx 129 --ny 129 --iterations 6000
./build/bench_fv --cavity 1000 --nx 129 --ny 129 --iterations 6000
python3 validation/cases/check_cavity_profiles.py
./build/bench_fv --cavity 1000 --fluent /Users/ezzy/.Others/Kutta/validation/fluent_examples/lid_driven_cavity/lid_driven_cavity.msh
./build/bench_fv --size 100 --repeats 3
./build/bench_bandwidth
```

Omit `--backend` to run CPU and Metal; `--backend cpu` or `--backend metal`
selects one. Repeat the Fluent command at Re=100/400 for the other V4 cases.
Generated evidence is under `build/lane-d/`: `tests-final.log`,
`performance-1M.log`, `bandwidth-final.log`, `profiler-oracle-final.log`,
`cavity129-final-Re{100,400,1000}.log`, `fluent-final-Re{100,400,1000}.log`,
`scatter-final.log`, `external-inventory.json`, and `results.html`/`results.svg`.
Profile CSVs are `build/fv-cavity-{cpu,metal}-Re*-*.csv`. These are local build
artifacts, not committed external fixtures.
