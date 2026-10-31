# Projected polyhedral faces for turbulent flow

**Result.** The turbulent momentum predictor and the omega, nu-tilde and
k/epsilon diffusion operators now use the projected normal-line face
reconstruction; SA also projects its mu+mu_t momentum diffusion. V6 (the Ansys
mixing elbow, SST, polyhedral, compared with Fluent's own saved fields)
improves in every field: omega error falls from 8.47% to 2.81%, velocity from
3.60% to 3.37%, pressure from 2.96% to 2.67%, k from 15.48% to 14.72%, and
the six wall-force errors fall by 28-80%. The inlet-outlet
pressure-drop error worsens from 1.19% to 1.97%. The laminar 3-ring fit could
not be reused: at FSAE-car scale it needs on the order of 100 GB. The turbulent
path uses a new compact fit (walk-count weights, per-cell moments) that is as
accurate on V6 and stores 40 floats per cell plus 30 per corrected face.

**Slow tier (this machine, same session conditions as far as possible):**
main `63287b9` 24 cases, 10,884 assertions, 68 failed; this branch `38c033e`
24 cases, 10,876 assertions, 64 failed. `compare_slow_assertions.py`
([comparison.json](comparison.json)) reports 5 recoveries (four V6
inlet-shear force gates, V7 SA bump Metal level-2 Cf at x=0.75), 8 missing
V6 REQUIREs (the seeded and hybrid solves converge 50 steps sooner, so the
progress loop runs one fewer chunk), and **one regression: V8's Metal 80-cell
channel missed its 180 s wall-clock budget** (5,600 steps in 181.4 s; main
5,700 in 168.5 s). V8 never constructs a reconstruction; it drives
`TransportEquation` with none and runs the same kernels as main, and its CPU
results are bitwise identical. In alternating V8-only runs main missed the
same budget once in three (157.8-182.7 s) and the branch twice in three
(169.1-180.2 s); I count it as the known shared-GPU timing flake, but by the
letter of the acceptance rule it is a new failed assertion in this tier. No
threshold, budget or gate was changed. The first version of this branch
(`87d62ce`, [its tier](first_candidate_comparison.json)) also failed the
level-2 x=0.632 Cf gates of V5 SST and V7 SA bump on CPU at the 600 s budget;
`38c033e` removes that cost and both pass (section 6).

Commits: `11b429c` (reconstruction), `a7b525c` (V6 operator probe), `87d62ce`
(SA wall forces), `38c033e` (cost and scope). Baseline: main `63287b9`.
Machine: Apple M4 Pro (14 cores), macOS 26.6.2, Apple clang 21, shared with
other lanes' slow tiers throughout. V6 fixture SHA-256: `mixing_elbow.cas.h5`
`8145813f…`, `.dat.h5` `863a8132…`.

## 1. What limits V6, measured

All V6 numbers are relative L2 errors against Fluent 22.2's saved `.dat.h5`
fields on the same 17,822-cell polyhedral mesh, with the documented `wfgc`
approximation (warped-face gradient correction off). Seeded and cold-hybrid
starts agree to about 1e-5 relative.

### Operator ablation (converged V6, Metal, seeded)

The first rows enable only the stated projected operators, using the laminar
lane's distance-weighted 3-ring fit, through a scratch environment switch
([patch](ablation/scratch_ablation_knobs.patch), logs in [ablation](ablation/)).
"Walls omitted" removes wall values from the velocity/omega fit samples.

| Projected operators | Velocity % | Pressure % | k % | Omega % | Drop error % |
| --- | ---: | ---: | ---: | ---: | ---: |
| none (main) | 3.602 | 2.962 | 15.48 | 8.475 | 1.187 |
| momentum predictor | 3.200 | 3.018 | 15.89 | 8.778 | 1.910 |
| momentum diffusion | 3.574 | 3.004 | 15.37 | 8.534 | 1.167 |
| both momentum | 3.169 | 3.063 | 15.78 | 8.805 | 1.886 |
| both momentum, walls omitted | 3.441 | 2.833 | 14.90 | 8.781 | 2.366 |
| k diffusion only | 3.600 | 2.960 | 15.48 | 8.474 | — |
| omega diffusion only | 3.607 | 2.963 | 15.14 | 4.061 | — |
| k + omega diffusion, walls omitted | 3.602 | 2.966 | 15.34 | 3.418 | — |
| all, walls omitted | 3.443 | 2.832 | 14.77 | 3.467 | — |

With the retained walk-weighted fit (section 3), walls omitted:

| Projected operators | Velocity % | Pressure % | k % | Omega % | Drop error % |
| --- | ---: | ---: | ---: | ---: | ---: |
| all | 3.363 | 2.688 | 14.71 | 2.821 | 1.987 |
| all but momentum predictor | 3.600 | 2.964 | 15.41 | 2.900 | 1.20 |
| all but k diffusion | 3.364 | 2.690 | 14.71 | 2.821 | 1.988 |
| all but momentum diffusion | 3.366 | 2.666 | 14.73 | 2.806 | 1.974 |
| **all but k and momentum diffusion (retained)** | **3.367** | **2.667** | **14.73** | **2.806** | **1.975** |

Findings:

- **Omega diffusion on skewed faces is the largest turbulent-operator error.**
  It alone halves the omega error.
- **The momentum predictor moves velocity, pressure, k and the wall forces,
  and the pressure drop the other way.** Without it the walk fit keeps the
  omega gain but loses the rest.
- **k diffusion and mu+mu_t momentum diffusion are marginal** (≤0.03 points
  in any field with either fit). They are left on the compact operator under
  wall functions: together they cost four of the eight fits per SST iteration
  (section 6). SA keeps projected momentum diffusion: its walls are resolved
  Dirichlet data and its wall traction is the projected wall flux.
- Including wall values in the fit helps velocity but hurts pressure and k.
  Kutta's SST and k-epsilon walls are wall-function closures: the wall model,
  not a polynomial, relates the wall to the first cell. Omitting wall samples
  under wall functions improves all four fields together and is the same rule
  SST's strain LSQ already applies.
- Per-side gradients or a gradient-of-gradient Hessian in the secondary
  diffusion term do **not** reproduce the omega gain (omega 8.86% and 9.23%,
  CPU, k/omega only). The gain needs a quadratic fit over a wide cloud.

### Where the error lives

`bench_v6_operators … cpu converge` splits squared field error by cell class:
wall-adjacent cells, and interior cells by their worst face skew (tangential
centroid offset over normal distance). Logs: [main](probe_main_cpu.log),
[this branch](probe_final_cpu.log).

| Field, class L2 % | wall (2,436) | skew<0.1 (2,858) | skew<0.3 (8,085) | skew≥0.3 (4,443) |
| --- | ---: | ---: | ---: | ---: |
| velocity, main | 9.73 | 3.62 | 3.29 | 2.94 |
| velocity, this branch | 8.93 | 3.51 | 3.15 | 2.56 |
| pressure, main | 3.01 | 3.30 | 3.02 | 2.75 |
| pressure, this branch | 2.81 | 3.06 | 2.66 | 2.48 |
| k, main | 13.6 | 11.2 | 17.2 | 20.0 |
| k, this branch | 12.3 | 11.1 | 16.4 | 18.9 |
| omega, main | 1.62 | 7.27 | 9.14 | **35.0** |
| omega, this branch | 1.52 | 7.01 | 8.77 | **8.74** |

On main, 94% of the squared omega error sat in the high-skew interior cells;
the projected omega diffusion removes most of it. Velocity, pressure and k
errors are nearly uniform across classes, including the least skewed cells.
**The remaining velocity/pressure/k error is therefore not a skew-truncation
error of these operators.** At Fluent's own saved state Kutta's k and omega
equations still have residuals about 150× and 1000× Fluent's final residuals
(raw k 3.81e-5 vs 2.50e-7; omega 1.92 vs 1.92e-3); the projected operators cut
the saved-state scaled residuals by 1-12% (u 1.3%, v 5.6%, w 12%, k 7.3%,
omega 7.8%). Kutta's eddy viscosity evaluated at Fluent's own k and omega
already differs by 1.03% L2 (10.5% max), unchanged by this work. The remaining
candidates, from earlier evidence in `src/kutta/fv/VERIFICATION.md`, are
inferences rather than measurements here: the unimplemented WFGC gradient
correction (which feeds strain, production and mu_t), Fluent's undocumented
momentum face interpolation (its saved flux fits an average of side
reconstructions; that formula gave 2.50% velocity but was rejected for V2 tet
order), and SST source evaluation.

Mass flux at the saved state, interior faces, relative to Fluent's `SV_FLUX`:
main 1.075% L2 (median .762%); this branch .965% L2 (median .837%).

## 2. Why the laminar fit could not be reused

The laminar `ProjectedFaceReconstruction` stores, for every corrected face, a
weight per functional for every sample in the 3-ring cloud around both cells.
On V6 that is 153 samples per face.

| V6 setup (Metal, converge probe) | Samples / face | Peak RSS | Velocity / pressure / k / omega % |
| --- | ---: | ---: | --- |
| main (no projection) | — | 180 MB | 3.602 / 2.962 / 15.48 / 8.475 |
| laminar fit, 1 ring | 19 | 411 MB | diverged |
| laminar fit, 2 rings | 64 | 1.02 GB | 3.716 / 3.042 / 15.69 / 3.982 |
| laminar fit, 3 rings | 153 | 1.65 GB | 3.443 / 2.832 / 14.77 / 3.467 |
| **walk-weighted moments (retained)** | — | **218 MB** | **3.367 / 2.667 / 14.73 / 2.806** |

The laminar rows used one reconstruction per equation. Even shared, the 3-ring
fit stores about 2.4 KB per corrected face for the three turbulent
functionals. An FSAE car is 20-40M poly-hexcore cells
(`docs/research/fluent-numerics-and-files.md`); its skewed poly and prism faces
alone would need tens to hundreds of GB, and stream that per field per
iteration against the ~2 KB/cell/iteration budget in `docs/ARCHITECTURE.md`.
Smaller clouds lose the accuracy or diverge. The laminar path is unchanged.

## 3. Compact walk-weighted reconstruction

`WalkFaceReconstruction` evaluates exactly the same functionals as the laminar
fit: the deferred diffusion correction `alpha[q(N*)-q(N)-q(P*)+q(P)]`, and the
normal-line predictor value and centroid link. Only the fit weights differ.
A sample's weight is its number of walks of length ≤3 (self-steps allowed)
from the owner cell plus from the neighbour, i.e. `(I+A)^3` from each side.
That makes every fit moment a sum of per-cell moments

```
M_c^(k+1) = M_c^(k) + sum_{d in faces(c)} T(x_d-x_c) [M_d^(k) + (q_d-q_c) G_d^(k)]
M_c^(k)   = sum_j walks_c^(k)(j) b(x_j-x_c) (q_j-q_c),   b = (x,y,z,x²/2,xy,xz,y²/2,yz,z²/2,1)
```

where `T` re-expresses a quadratic moment about another origin and `G^(k)`
are the same moments of the constant 1. Observed boundary values join their
owner cell at hop 0. Three gather kernels (`fv_moment_hop`) build `M^(3)`;
each face shifts its two cells' moments to the face centroid and applies its
stored 10-term functional rows. The normal matrix, its pseudo-inverse and the
functional rows are built once in double precision (multi-threaded; 1.2 s for
V6), with the laminar lane's rank and exactness checks. Faces on the normal
line at roundoff, rank-deficient faces and clouds with discontinuous wall
traces keep the original operators. Moment, geometry and coefficient rows are
cell- or face-major, so each neighbour gather touches one or two cache lines;
on the serial CPU backend that made the fits 38% cheaper than plane-major
storage, with identical results.

Device memory: 40 floats per cell (`G^0..G^3`), 30 floats per corrected face,
and 24 floats per cell of reusable moment scratch. One reconstruction serves
velocity and every turbulence scalar with the same observed boundary values.
On SST and k-epsilon those are the same (inlets and pressure inlets; walls are
omitted), so V6 builds one reconstruction, used by the momentum predictor and
omega.

Model rules:

- **SST, realizable k-epsilon:** wall values are omitted from velocity, omega
  and epsilon fits. Only the momentum predictor is projected for momentum;
  SST k diffusion stays compact. Wall faces have zero momentum diffusivity in
  the wall-function contract, so the wall-function rows, wall shear and force
  path are unchanged. Realizable k-epsilon projects k and epsilon diffusion;
  no case here measures that choice (V8 is a hex channel with no corrected
  faces).
- **SA:** walls are resolved Dirichlet data (u = wall velocity, nu-tilde = 0),
  sampled as on laminar walls, with the laminar discontinuous-trace guard.
  Momentum diffusion is projected, and wall traction keeps its normal
  derivative equal to the assembled projected wall flux (commit `87d62ce`), so
  SA forces and Cf stay flux-consistent.
- The projected transport kernel now applies a correction only on Dirichlet
  boundary faces, as the legacy kernel does; previously an outflow face of a
  pressure inlet with a fit could have received one.

Verification in the quick tier (`tests/fv/test_walk_projection.cpp`, CPU and
Metal, 1,758 assertions):

- Quadratic exactness on hex, tet and three poly annulus levels, with and
  without boundary samples: normal-flux error ≤5.8e-5 (bound 2e-4, the laminar
  test's), normal-line value ≤1.9e-6 (bound 2e-5). Hex faces are not corrected.
- An independent double-precision oracle enumerates each face's walk-weighted
  cloud and solves its normal equations by Gaussian elimination. The device
  moment path matches it to 3.9e-8 absolute (correction scale 0.10) and
  predictor values to 3.9e-7, on tet and poly meshes.
- Reversing every face's ownership reverses the corrections and preserves face
  values within 2e-6; an all-zero limiter restores the original predictor.

## 4. V6 before and after

Main `63287b9` and this branch (`38c033e`), both from their slow tiers. Errors
in %. Standalone runs of the branch: [CPU](v6_final_cpu.log), [Metal](v6_final_metal.log).

Seeded from Fluent's fields:

| Quantity | CPU main | CPU branch | Metal main | Metal branch |
| --- | ---: | ---: | ---: | ---: |
| steps to converge | 350 | 300 | 350 | 300 |
| velocity L2 | 3.602 | 3.367 | 3.602 | 3.367 |
| pressure L2 | 2.962 | 2.667 | 2.962 | 2.667 |
| k L2 | 15.48 | 14.73 | 15.48 | 14.73 |
| omega L2 | 8.475 | 2.806 | 8.475 | 2.806 |
| pressure drop (worse) | 1.187 | 1.974 | 1.187 | 1.975 |
| wall-elbow pressure force | 1.556 | 0.656 | 1.557 | 0.655 |
| wall-elbow shear force | 2.259 | 1.636 | 2.259 | 1.636 |
| wall-elbow total force | 1.506 | 0.600 | 1.506 | 0.599 |
| wall-inlet pressure force | 2.073 | 0.527 | 2.073 | 0.527 |
| wall-inlet shear force | 0.610 | 0.119 | 0.610 | 0.119 |
| wall-inlet total force | 2.070 | 0.527 | 2.070 | 0.527 |

Cold hybrid start:

| Quantity | CPU main | CPU branch | Metal main | Metal branch |
| --- | ---: | ---: | ---: | ---: |
| steps to converge | 550 | 500 | 550 | 500 |
| velocity L2 | 3.601 | 3.367 | 3.601 | 3.367 |
| pressure L2 | 2.962 | 2.666 | 2.962 | 2.665 |
| k L2 | 15.48 | 14.72 | 15.48 | 14.72 |
| omega L2 | 8.475 | 2.805 | 8.475 | 2.805 |
| pressure drop (worse) | 1.185 | 1.973 | 1.185 | 1.973 |
| wall-elbow pressure force | 1.558 | 0.655 | 1.557 | 0.655 |
| wall-elbow shear force | 2.259 | 1.635 | 2.259 | 1.635 |
| wall-elbow total force | 1.508 | 0.599 | 1.507 | 0.599 |
| wall-inlet pressure force | 2.074 | 0.527 | 2.074 | 0.527 |
| wall-inlet shear force | 0.610 | 0.119 | 0.610 | 0.119 |
| wall-inlet total force | 2.071 | 0.527 | 2.070 | 0.527 |

CPU and Metal agree to 1e-5 relative in every field and force. Failed V6
assertions: 29 on main, 25 on this branch, both backends together. The
four recoveries are the inlet shear-force gates (0.61% → 0.12%, gate 0.5%).
V6 records 8 fewer assertions because seeded starts converge in 300 steps,
not 350, so the 50-step progress loop runs one fewer chunk of its two
finite-velocity REQUIREs per backend and start mode.

## 5. Other verification cases

From the final slow tier against main's (both on this machine; logs and XML
in this directory). Max-error is each test's normalized frozen error over
stations, drag and Cp (lower is better); "Cf x=0.632" is the level-2 gate
that the first version of this branch failed on CPU.

| Case | Backend | Main | Branch |
| --- | --- | --- | --- |
| V5 SST bump max-error, levels 0/1/2 | CPU | 2.108 / 1.174 / 0.528 | 2.025 / 0.889 / 0.468 |
| V5 SST bump max-error, levels 0/1/2 | Metal | 2.117 / 1.186 / 0.635 | 2.033 / 0.896 / 0.349 |
| V5 SST bump level-2 Cf x=0.632 error (gate 0.1) | CPU | 0.0992 (2,300 steps) | 0.0899 (2,100 steps) |
| V5 SST bump level-2 Cf x=0.632 error | Metal | 0.0636 | 0.0414 |
| V7 SA bump max-error, levels 0/1/2 | CPU | 2.133 / 1.149 / 0.619 | 2.046 / 0.926 / 0.360 |
| V7 SA bump max-error, levels 0/1/2 | Metal | 2.131 / 1.149 / 0.669 | 2.047 / 0.932 / 0.278 |
| V7 SA bump level-2 Cf x=0.632 error (gate 0.1) | CPU | 0.0765 (2,600 steps) | 0.0966 (2,100 steps) |
| V7 SA bump level-2 Cf x=0.632 error | Metal | 0.0581 | 0.0467 |
| V7 SA bump level-2 Cf x=0.75 error (gate 0.1) | Metal | 0.105 (fail) | 0.069 (pass) |
| V5 SST / V7 SA flat plate level-2 Cf and Cd | both | 0.0085 / 0.0082; 0.0078 / 0.0072 | same to 3e-7 |

The V5 and V7 bump grids are curvilinear hexes, so most of their faces are
corrected. Their CPU level-2 runs stop for the 600 s budget; with that budget
the V7 CPU gate passed by 0.003 in the slow tier (2,100 steps) and by 0.011 in
the standalone runs of section 6 (2,200 steps). Cases with orthogonal hex or
tet meshes are unchanged: V1, V2 (hex and tet), V3, V4 and the FSAE V7-V11
cases have identical assertion outcomes. V2 tet keeps its 6 existing
failures on the laminar path.

V8 (realizable k-epsilon channel, 1D hex, no corrected faces, so no walk
reconstruction is built) gives bitwise identical CPU results. Its Metal
80-cell run has a 180 s wall-clock budget, which it missed in the final tier
(5,600 steps in 181.4 s against 5,700 steps in 168.5 s on main); every V8 Metal
grid ran about 18% slower in that tier than in main's, with the same numbers.
Alternating main/branch runs of V8 alone
([log](timing/v8_metal_alternation.log)), 80-cell Metal:

| Round | Main | Branch |
| --- | --- | --- |
| 1 | 5,700 steps, 178.4 s, pass | 5,600 steps, 180.2 s, fail |
| 2 | 5,700 steps, 157.8 s, pass | 5,700 steps, 169.1 s, pass |
| 3 | 5,400 steps, 182.7 s, fail | 4,400 steps, 180.2 s, fail |

## 6. Cost

### Level-2 TMR budget runs (serial CPU backend)

V5 (SST) and V7 (SA) bump level 2 have a 600 s wall-clock budget, and on the
serial CPU backend they stop for the budget, not for convergence, so their
level-2 Cf gates depend on how many steps fit in 600 s. The first version of
this branch (`87d62ce`) reached 1,700 and 1,900 steps there in the slow tier,
against main's 2,300 and 2,600, and failed the x=0.632 Cf gate on both cases
although it is more accurate per step. Commit `38c033e` cut that cost.
Runs below alternate main and branch, two rounds, with
`bench_tmr … bump cpu 2 3000 … default 600` and the scratch SA driver
[`timing/sa_bench.cpp`](timing/sa_bench.cpp), which calls the slow test's
`sa::run` ([summary](timing/summary.log)). Other lanes' slow tiers shared the
machine.

| Case | Main steps (s/step) | Branch steps (s/step) | Cf x=0.632 error, main → branch (gate ≤0.1) |
| --- | --- | --- | --- |
| V5 SST, round 1 | 2,400 (0.258) | 2,100 (0.297) | 0.0924 → 0.0898 |
| V5 SST, round 2 | 2,400 (0.251) | 2,100 (0.297) | 0.0924 → 0.0898 |
| V7 SA, round 1 | 2,600 (0.233) | 2,200 (0.285) | 0.0765 → 0.0890 |
| V7 SA, round 2 | 2,800 (0.219) | 2,200 (0.273) | 0.0665 → 0.0890 |

At equal step counts the branch is closer to the reference on both cases
(final slow-tier Metal runs, 3,000 steps, Cf error at x=0.632: SST 0.0636 → 0.0414, SA 0.0581 → 0.0467), so the per-step cost decides these gates. The SA margin is
narrower because SA keeps projected momentum diffusion. Both pass in every
run above and in the final slow tier (section 5), V7 by only 0.003 there.

### Where the time goes

V5 bump level 2, CPU, 100 steps, fit time measured with temporary timers:

| Version | Fits per step | Fit time per 100 steps | Solve time per 100 steps (main 19.5-21.8 s) |
| --- | ---: | ---: | ---: |
| `87d62ce`, plane-major rows | 8.65 | 9.9 s | 31.7 s |
| cell-major rows, min/max packed | 8.65 | 6.1 s | 27.7 s |
| plus compact k and momentum diffusion (`38c033e`) | 4.25 | 3.6 s | 24.5-26.1 s |

The rest of the gap is the pressure solve: over 400 steps the branch needs
14,294 AMG-PCG iterations against main's 11,502 (+24%), because the projected
predictor changes the provisional mass flux that the pressure equation
corrects. Per V6 solve (single runs from the two slow tiers): CPU
119 → 121 ms/step seeded, 101 → 119 hybrid; Metal 24.4 → 26.9
seeded, 24.0 → 25.9 hybrid. V6 converges in 50 fewer steps, so the
solve is faster overall except CPU hybrid (55.6 s → 59.4 s).

### Memory and setup

- Device: 40 floats per cell of geometry moments, 24 per cell of moment
  scratch, 30 per corrected face. V6 peak RSS 180 MB → 218 MB (Metal converge
  probe); the CPU probe, where device buffers are host memory, peaks at
  198 MB → 258 MB.
- Setup: the fits are built once in double precision on all host threads.
  V6 setup 1.2 s; V5 bump level 2 (56,320 hex cells) 0.70 s → 0.98 s.
  Host transient memory is about 124 bytes per face plus 320 bytes per cell.

## 7. What remains

- **Pressure drop.** V6's inlet-outlet pressure-drop error worsens from 1.19%
  to 1.97%. The ablation pins it on the momentum predictor (without it the
  drop error is 1.20%), which is also what improves velocity, pressure, k and
  every wall force. Not resolved.
- **The rest of V6's error is not in these operators.** Velocity, pressure
  and k errors are uniform across skew classes, and Kutta's k/omega residuals
  at Fluent's saved state remain 150× and 1000× Fluent's. Next candidates
  (inferred, not measured here): WFGC, Fluent's momentum face interpolation,
  SST source evaluation. V6 still fails its existing accuracy gates
  (velocity 0.5%, k 0.1%, drop 0.5%, forces 0.5%).
- **CPU per-step cost** is about 1.17× main on SST and 1.23× on SA (V5/V7
  bump, alternating runs). The level-2 Cf gates now pass with margins of 0.010
  (V5) and 0.003 (V7) on the serial CPU backend in the final tier; a slower or more
  heavily shared machine could flip them, as it can on main (main passed V5
  by 0.0008 in its tier).
- **V8 Metal timing.** The 180 s V8 budget failed in this branch's final
  tier and in one of three main runs; see section 5.
- **Car scale is extrapolated, not measured.** At 30M cells the device
  storage is about 7.7 GB plus 120 bytes per corrected face; host setup
  transients reach tens of GB, and setup time grows with the number of
  non-orthogonal faces (orthogonal hexcore faces are rejected before any
  fit). Measure on a real FSAE mesh before relying on it there.
- **CUDA runtime is unverified.** The CUDA build (CI flags, sm 75/86/89/120,
  `--fmad=false`) compiles and links in `nvidia/cuda:13.0.2-devel-ubuntu24.04`
  (arm64 on this Mac), and the quick tier passes there on the CPU backend
  (227 cases, 814,247 assertions, GPU-presence test excluded as in CI). No
  NVIDIA device was available, so CPU==CUDA parity was not run.
- **Unmeasured model choices.** Realizable k-epsilon projects k and epsilon
  diffusion and SA projects momentum diffusion; no polyhedral case with a
  reference covers either model.

## Reproduction

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build
export KUTTA_FLUENT_EXAMPLES=<fixtures>/validation/fluent_examples
# V6 per backend
KUTTA_V6_BACKEND=cpu build/kutta-tests -tc="FV V6 Fluent mixing elbow fixed point and cold start"
# Operator probe (h5py needed for the flux dump)
python3 -I tools/forensics/dump_interior_flux.py $KUTTA_FLUENT_EXAMPLES/mixing_elbow mixing_elbow flux.txt
build/bench_v6_operators $KUTTA_FLUENT_EXAMPLES/mixing_elbow flux.txt cpu converge
# Ablations: apply ablation/scratch_rule_knobs.patch to 87d62ce, then
# SCRATCH_NO_KDIFF=1 and/or SCRATCH_NO_MOMDIFF=1 (the distance-weight rows
# use ablation/scratch_ablation_knobs.patch on main 63287b9)
# Budget timing (alternate with a main build)
build/bench_tmr validation/data/nasa_tmr bump cpu 2 3000 out default 600
# Slow tier and assertion comparison
build/kutta-tests --test-suite=slow --success --reporters=xml --out=final_slow.xml
python3 tools/forensics/compare_slow_assertions.py baseline_slow.xml.gz final_slow.xml.gz > comparison.json
```
