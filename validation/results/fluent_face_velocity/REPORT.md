# Fluent's face velocity, the limiter and the V6 mixing elbow

**The default is unchanged.** One formulation closes most of the V6 gap.
It combines a Péclet-weighted side-reconstruction face velocity with Fluent's
differentiable limiter on momentum. On V6 it gets:

- velocity 3.37% → 1.39%;
- pressure 2.67% → 0.94%;
- k 14.7% → 8.4%;
- pressure drop 1.97% → 0.105%.

It also passes three V6 gates that fail on main. It keeps every V2 order
gate. It does not ship, because the rule is "at least as close to Fluent
on every converged reference gate", and on the k-epsilon elbow it moves
velocity, pressure, k and epsilon away from Fluent (all still inside
their gates).

The two leads resolve as follows.

- **Lead 2, face velocity against V2 order: resolved.** The side
  reconstruction loses V2 order only on the hex mesh. On tets it is
  better than main. On hex the loss comes from the no-slip wall layers
  in the viscous limit:
  - a quadratic wall-normal gradient does not fix it;
  - switching to the point form near walls makes it worse.

  Weighting the side correction by the cells' advective share of the
  momentum diagonal fixes it. The share is O(cell Reynolds number), so the
  correction is O(h³) as h → 0 and every V2 hex gate passes. Hex pressure
  order rises from main's 1.821 to 1.881, and tet stays within 0.004 of
  main. On V6, where cells are convection dominated, the face velocity
  becomes Fluent's side form.
- **Lead 1, limiter: partly resolved.** Main passes the V6 wall-inlet shear
  gate (0.12%) by cancellation. Face by face its shear is 4.5% off, and
  along the inlet pipe it is −2% to +4.5% off. Most variants that move the
  flux or the limiter toward Fluent fail that gate (0.9-3.3%), while the
  ones measured face by face improve it.
  - With the Péclet-weighted face velocity, Barth-Jespersen on momentum
    stalls Coupled on V6 and fails the inlet shear.
  - An ε-floor on Barth-Jespersen (Fluent's `recon/limit/epsilon`) is worse.
  - The differentiable limiter converges in 150 iterations. It halves
    every face-by-face wall error and removes k's +8.7% bias.
  - The smooth limiter has to act in momentum convection, not only in the
    flux.
  - The limiter rpvars are identical in all 55 corpus cases. So the data
    cannot tie a smoother limiter to V6 alone, and a global switch fails a
    k-epsilon gate (previous lane).
- **A Fluent setting that differs:** `pressure/ho-rc-scheme` is 2 in V6
  and 1 in the k-epsilon elbow and the SST nozzle. It is undocumented. It is
  the only rpvar found that could select the face-velocity form per case
  (section 5).

This commit adds no solver change. It adds the face-by-face
wall-shear tooling (`tools/forensics/wall_shear.py`, `KUTTA_V6_WALL_DIR`
in the V6 test) and a `cff.load` fix that lets the saved-flux forensics
read meshes whose zones interleave (the SST nozzle). The experiment
switches are stored as `evidence/scratch_face_velocity_knobs.patch`.

## Builds and data

- Main is `1bb54a9`, built from `git archive` in a separate tree.
- The candidate is this commit's tree.
- The variants are the candidate tree plus the scratch patch, selected
  by environment variables (section 8). The default path of that build
  is the candidate.
- All are Release (`-O3 -DNDEBUG`), on an Apple M4 Pro with macOS 26.6.2
  and Metal. The machine was shared.
- The Fluent pairs are the converged corpus references, read through
  `KUTTA_FLUENT_EXAMPLES`:
  - V6 is `mixing_elbow`, Fluent 22.2, SST, poly, Coupled.
  - The k-epsilon elbow is `elbow_param`, tet, first-order momentum,
    standard pressure.
  - The SST nozzle is `nozzle_para_named_1`.
  - The laminar cylinder is `cylinder_tutorial`, Fluent 18.
- All V6 numbers are CPU, seeded start, against Fluent's saved fields.
  For the formulation that matters (pd below), Metal and CPU agree to
  every printed digit, and the cold-hybrid start converges in 200
  iterations to the same fields within 3e-5.

## 1. What Fluent documents

Ansys Fluent Theory Guide 2024 R2:

- **23.4.2, eqs. 23-46 to 23-49.** Momentum-weighted face flux
  `Ĵ_f = ρ_f (a_p,c0 v_n,c0 + a_p,c1 v_n,c1)/(a_p,c0 + a_p,c1) + J_ho,ap`.
  The "high-order correction" `J_ho,ap` has no published form.
- **23.3.4 and 23.3.4.1.** The Standard limiter is Barth-Jespersen,
  cell-to-face by default.
- **23.3.4.3.** Standard limiters "tend to stall"; the Differentiable
  limiter is a modified Venkatakrishnan limiter (Wang). No coefficient is
  given.
- `pressure/ho-rc-scheme`, `recon/lim/method 6`, `recon/lim/sigp 0.05` and
  `recon/limit/epsilon 0.01` appear in the case files but are not
  documented. Their meanings below are inferences from data.

## 2. Lead 2: the side reconstruction and V2 order

Notation (VERIFICATION.md, "Why affine exactness did not protect V2"):

```text
w = aP/(aP+aN), rP = xf-xP, rN = xf-xN
point = w uP + (1-w) uN + (GP+GN)/2 . (w rP + (1-w) rN)
side  = w (uP+GP.rP) + (1-w)(uN+GN.rN)
```

Kutta's turbulent and laminar paths replace `point` by the projected
normal-line value on non-orthogonal faces. That is what main runs on V6
and on V2 tets.

### 2.1 Measured V2 orders (CPU, finest pair)

Columns: velocity, pressure, then wall pressure, viscous and total force.
The gates are u ≥ 1.9 everywhere; on hex p ≥ 1.8 and forces ≥ 1.8; on
tet p ≥ operator order − 0.1 = 0.914. Main fails tet u and the three tet
forces (known failures).

| face velocity | hex u / p / forces | tet u / p / forces |
| --- | --- | --- |
| main (point, projected on tets) | 1.916 / 1.821 / 1.915, 1.959, 2.012 | **1.815** / 0.967 / **1.052, 0.843, 0.623** |
| side everywhere | **1.825 / 1.571 / 1.204, 1.061, 0.811** | 1.916 / 1.028 / 1.078, 0.852, 0.540 |
| side, quadratic wall-normal gradient in wall cells | **1.838 / 1.564 / 1.171, 1.023, 0.770** | 1.918 / 1.027 / 1.081, 0.855, 0.544 |
| side, point on faces of wall cells | **1.785 / 1.368 / 3.081, 2.873, 0.684** | not run |
| side, point on faces within two layers | **1.782 / 1.308 / 1.667, 1.868, 2.355** | not run |
| **Péclet-weighted: point + β (side − point), β = (outP+outN)/(aP+aN)** | 1.917 / **1.881** / 1.980, 1.997, 2.010 | **1.813** / 0.963 / **1.054, 0.843, 0.621** |
| product weight β = (outP/aP)(outN/aN) | 1.917 / 1.830 / 1.921, 1.960, 2.008 | **1.805** / 0.941 / **1.074, 0.852, 0.612** |

`out` is a cell's convective outflow Σ max(F_out, 0); `a` is the physical
momentum diagonal. V2 has `limit=false`, so the limiter does not enter
these numbers. The slow tier repeats them on CPU and Metal (section 7).

What this establishes:

- **The order loss is the hex wall layer, as the 2026-10-29 analysis said.**
  On tets the side form passes the velocity gate that main fails.
- **It is not the wall cell's biased LSQ gradient.** Replacing that
  gradient's wall-normal component by the quadratic through the wall
  value, the cell and its opposite neighbour leaves hex at 1.56. (This
  holds on these orthogonal hexes; on general meshes it is not a
  candidate anyway.)
- **Switching formulas by wall layer makes it worse.** The switch puts an
  O(h²) flux-error jump on an interior face, which is an O(h) divergence
  there.
- **The Péclet weight preserves order by construction.** β is the
  advective share of the cells' momentum diagonal. At fixed Reynolds number
  it is O(ρ U h/μ) as h → 0. The correction β·(side − point) is then O(h³),
  below the scheme's truncation. On V6, β is near 1 in the core, so the
  face velocity is the side form that Fluent's saved flux carries. In the
  saved-flux test (`face_flux.py` section 3, mixing elbow) this weighting
  ("Momentum-Peclet") has:
  - median face error 1.02e-3 |F|, against 6.2e-4 for the pure side form
    and 7.7e-3 for point;
  - sum|div| 7.7e-4 kg/s, the lowest of every formula tested.

  The weighting is Kutta's construction for order, not a documented Fluent
  formula.

### 2.2 V6 with each face velocity (CPU, seeded)

"Converged" is the V6 test's own criterion: all scaled residuals below
1e-6 and continuity below Fluent's final raw value. Force columns are the
relative error of each wall zone's pressure, shear and total force. Bold
marks a gate that fails on main and passes here, or the reverse.

| run | limiter | iterations | velocity | pressure | k | omega | drop | elbow p / shear / total | inlet p / shear / total |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| main | BJ | 200 | 3.365 | 2.663 | 14.71 | 2.806 | 1.967 | 0.66 / 1.63 / 0.61 | 0.52 / 0.12 / 0.52 |
| side | BJ | **4000, stalls at 1e-4** | 2.723 | 1.689 | 11.50 | 2.619 | 1.041 | 1.02 / 1.50 / 1.05 | 1.11 / **1.57** / 1.11 |
| side, unlimited gradients | none in flux | 200 | 3.179 | 2.090 | 11.99 | 2.695 | 1.290 | 1.00 / 1.41 / 1.03 | 1.15 / **0.93** / 1.15 |
| Péclet | BJ | **1500, stalls at 3e-5** | 2.552 | 1.403 | 11.52 | 2.593 | 0.621 | 0.24 / 1.51 / **0.22** | 0.21 / **1.61** / **0.23** |
| Péclet, SIMPLEC override | BJ | 300 | 2.554 | 1.399 | 11.49 | 2.593 | 0.605 | 0.23 / 1.49 / **0.20** | 0.20 / **1.61** / **0.22** |
| Péclet | BJ, ε-floor 0.01 | **1500, stalls at 1e-5** | 3.511 | 1.630 | 14.93 | 2.700 | 0.680 | 0.31 / 2.43 / **0.27** | 0.50 / **1.95** / 0.51 |
| Péclet | BJ, ε-floor 0.05 | 150 | 4.311 | 1.926 | 17.02 | 2.946 | **0.486** | 0.28 / 2.81 / **0.21** | 0.45 / **3.28** / **0.48** |
| main face velocity | differentiable (VW K=0.05) | 150 | 2.280 | 1.926 | 9.88 | 2.583 | 1.699 | 0.48 / 0.89 / **0.44** | 0.20 / **1.61** / **0.21** |
| main face velocity | differentiable in the flux only | 200 | 3.416 | 2.560 | 14.53 | 2.798 | 1.889 | 0.70 / 1.67 / 0.64 | 0.27 / 0.09 / **0.27** |
| Péclet | differentiable in the flux only | 200 | 2.407 | 1.304 | 11.29 | 2.545 | 0.507 | 0.18 / 1.47 / **0.15** | 0.21 / **1.54** / **0.23** |
| side | differentiable | 150 | 1.299 | 0.995 | 8.34 | 2.335 | **0.107** | 0.49 / 0.51 / 0.51 | 0.67 / 0.19 / 0.67 |
| **Péclet (pd)** | **differentiable** | **150** | **1.390** | **0.936** | **8.40** | **2.337** | **0.105** | 0.06 / 0.56 / **0.08** | 0.22 / 0.15 / **0.22** |
| Péclet, point base instead of projected | differentiable | 150 | 1.386 | 0.988 | 8.19 | 2.334 | **0.229** | 0.09 / 0.61 / **0.07** | 0.52 / 0.33 / 0.52 |
| product weight | differentiable | 150 | 1.505 | 1.098 | 8.47 | 2.391 | **0.060** | 0.26 / 0.67 / **0.24** | 0.75 / 0.20 / 0.75 |

- Under V6's own Coupled scheme, any side-type face velocity with
  Barth-Jespersen stalls. The SIMPLEC override converges the same fixed
  point (Péclet: 300 iterations, identical fields).
- The ε-floor variant (Fluent's `recon/limit/epsilon 0.01` read as Wang's
  threshold on Barth-Jespersen) does not converge at 0.01 and is far from
  Fluent at 0.05.
- The smooth limiter has to act on momentum convection: limiting smoothly
  only in the flux keeps the inlet shear at 1.54%.
- pd is the best balance: every field, the drop and four of the six wall
  forces are the best or near the best of all runs. The wall-elbow shear
  (0.56%, gate 0.5%) still fails.

## 3. Lead 1: the inlet-wall cells

`tools/forensics/wall_shear.py` matches each Kutta wall face to Fluent's
saved `SV_WALL_SHEAR`. It reports:

- the integrated force error (what the gate checks);
- the face-by-face error;
- the wall-cell tangential velocity and k;
- five bands along the coordinate the zone spans most.

Full output: `evidence/wall_shear.txt`.

| run | inlet: integrated | inlet: per face L2 | inlet: wall-cell U_t L2 | inlet: wall-cell k L2 | elbow: integrated | elbow: per face L2 | k bias (volume-weighted) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| main | 0.12% | 4.49% | 4.18% | 8.0% | 1.63% | 8.80% | +8.7% |
| differentiable only | 1.61% | 3.58% | 2.88% | 5.2% | 0.89% | 6.15% | +3.0% |
| Péclet + BJ (SIMPLEC) | 1.61% | 4.19% | 4.25% | 8.6% | 1.49% | 7.14% | +6.3% |
| **pd** | 0.15% | **2.01%** | **1.91%** | **4.7%** | 0.56% | **4.18%** | **+0.6%** |

Inlet-pipe bands, Kutta/Fluent shear ratio from the hot inlet downstream:

| run | band 1 | 2 | 3 | 4 | 5 |
| --- | ---: | ---: | ---: | ---: | ---: |
| main | 0.985 | 0.979 | 1.030 | 1.045 | 0.986 |
| differentiable only | 0.999 | 1.007 | 1.038 | 1.058 | 1.009 |
| Péclet + BJ | 0.974 | 0.966 | 1.010 | 1.008 | 0.983 |
| pd | 0.990 | 0.996 | 1.019 | 1.023 | 0.997 |

What this shows:

- **Observed:** main's 0.12% is a cancellation of errors from −2.1% to
  +4.5%.
  - In every run the shear ratio tracks the wall-cell tangential velocity
    ratio to within 0.01. The hot pipe (Re about 2,000) sits in the
    viscous sublayer, where the shear is μU/y.
  - The wall-cell k is 17-51% high in bands 3-5 in every run, but it
    barely enters the shear there.
- **Observed:** the Fluent-ward variants fail the gate because they remove
  the upstream deficit at different rates from the downstream excess.
  pd improves every band.
- **Inferred:** Kutta's Barth-Jespersen limiter, combined with Fluent's
  face velocity, is farther from Fluent's converged state than
  Venkatakrishnan-Wang at K = 0.05. In that sense Lead 1 holds: Kutta's
  BJ over-limits.
  - At Fluent's saved state the two limiters fit about equally well
    (momentum residual 3.10e-4 vs 3.40e-4 on deep cells,
    `docs/research/fluent-momentum-forensics.md` 2.2). Both sit about 100
    times above Fluent's own residual floor, so the saved state does not
    separate them.
- **Not established:** that Fluent's Standard limiter is a smooth function.
  - The case says `recon/lim/type 0`.
  - Its limiter rpvars are identical across the whole corpus (section 5).
  - On the k-epsilon elbow with second-order pressure, the differentiable
    limiter alone fails the 1.7% drop gate (1.85%, previous lane).
  - Which stencil, face or cell form, or direction Fluent uses beyond
    TG 23.3.4 is undocumented. `recon/lim/method 6` and `sigp 0.05` are
    unexplained.

## 4. The other converged Fluent references

All slow-tier reference cases, CPU. Metal agrees to the digits shown.
Percentages.

| case, metric (gate) | main | Péclet + BJ | pd | differentiable only |
| --- | ---: | ---: | ---: | ---: |
| k-epsilon elbow velocity L2 (1.5) | 0.881 | 1.036 | 1.159 | 0.964 |
| pressure L2 (10) | 5.82 | 8.31 | 8.95 | 6.32 |
| k L2 (14) | 9.10 | 11.05 | 10.84 | 9.04 |
| epsilon L2 (27) | 17.6 | 20.6 | 20.4 | 17.5 |
| pressure drop (1.5) | 1.01 | **2.24** | 1.12 | 0.43 |
| wall force (0.6) | 0.316 | 0.250 | 0.310 | 0.249 |
| k-epsilon elbow, second-order pressure: velocity (1.6) | 1.05 | 1.22 | 1.33 | |
| pressure (8.7) | 5.99 | 7.91 | 8.31 | |
| k (13.1) / epsilon (25.4) | 8.71 / 16.9 | 10.7 / 20.2 | 10.5 / 20.0 | |
| pressure drop (1.7) | 1.11 | 0.63 | 0.22 | **1.85** (previous lane) |
| wall force (0.61) | 0.408 | 0.033 | 0.118 | |
| SST nozzle velocity (2.5) | 1.60 | 1.50 | 1.43 | 1.55 |
| pressure (0.9) | 0.527 | 0.593 | 0.165 | 0.114 |
| k (8) / omega (8.5) | 4.92 / 5.52 | 4.54 / 5.55 | 4.47 / 5.49 | 4.87 / 5.50 |
| pressure drop (0.8) | 0.462 | 0.555 | 0.070 | 0.020 |
| wall force (1.0) | 0.566 | 0.563 | 0.155 | 0.083 |
| cylinder velocity (0.097) | 0.0645 | 0.0595 | 0.0597 | |
| pressure (0.78) | 0.517 | 0.435 | 0.441 | |
| pressure drop (0.76) | 0.506 | 0.686 | 0.750 | |
| wall force (0.11) | 0.0706 | 0.0789 | 0.0750 | |
| Cd error (0.105) | 0.0695 | 0.0786 | 0.0750 | |
| Cl error / Cd (0.0185) | 0.0121 | 0.0071 | 0.0003 | |
| V4 cavity (763 assertions) | pass | not run | pass | |

- pd passes every gate in these cases.
- pd is farther from Fluent than main on:
  - every k-epsilon elbow field (both pressure schemes);
  - the cylinder's drop (0.750% against a 0.76% bound), wall force and Cd.
- So it does not meet the shipping rule.
- The cylinder's drop with pd sits 1.3% inside its bound. That is too
  close to call robust.

## 5. A Fluent setting that separates V6 from the others

Scanned with every `.cas.h5` in the corpus (55 readable):

- `pressure/ho-rc-scheme` is 2 in 18 cases, including V6, the axial fan,
  the rotating disk and the periodic cases.
  It is 1 in 36 cases, including the k-epsilon elbow and the SST nozzle.
  The Fluent 18 cylinder has no such rpvar.
- `recon/lim/type 0`, `recon/lim/method 6`, `recon/lim/sigp 0.05` and
  `recon/limit/epsilon 0.01` are the same in all 55.
- Saved-flux fits (`face_flux.py`; median face error / quiet-face median):

  | case (scheme) | side | point | Péclet |
  | --- | --- | --- | --- |
  | V6 (2) | 6.2e-4 / 3.3e-4 | 7.7e-3 / 3.5e-3 | 1.0e-3 / 4.4e-4 |
  | k-epsilon elbow (1) | 2.3e-3 / 1.4e-3 | 2.1e-3 / 1.3e-3 | 2.3e-3 / 1.4e-3 |
  | SST nozzle (1) | 3.4e-3 / 2.7e-3 | 4.9e-3 / 3.3e-3 | 3.4e-3 / 2.7e-3 |

  The k-epsilon elbow uses standard pressure and first-order momentum.
  The script's Rhie-Chow model assumes second-order pressure, so that row
  only says that nothing fits there. The nozzle's floor is eight times
  V6's and nothing tested explains it (R² ≤ 0.04).

**Inferred, not established:** `ho-rc-scheme 2` is the side-reconstruction
high-order Rhie-Chow term, and scheme 1 is something else. Only V6 is a
converged scheme-2 reference in the slow tier, so one case supports this.
A scheme-2 mapping alone would also not ship: with Barth-Jespersen, V6
stalls under Coupled and fails the inlet shear (section 2.2).

## 6. Decision and options

The default stays (Barth-Jespersen, point/projected face velocity). The
evidence supports one of these owner decisions:

1. **Recommended: keep the default and get one more converged
   `ho-rc-scheme 2` Fluent reference** (for example the rotating disk or a
   periodic case, saved converged) before mapping anything to that rpvar.
   With two scheme-2 references agreeing, map scheme 2 to pd and keep
   scheme 1 unchanged. Every scheme-1 reference is then bit-identical.
2. Ship pd for `ho-rc-scheme 2` now. V6 improves as above and nothing
   else in the slow tier changes. This ties a limiter change to a
   face-velocity rpvar with one case of evidence.
3. Ship pd globally. Section 7 lists every moved slow-tier value. This
   trades the k-epsilon elbow and cylinder agreement for V6 and the
   nozzle.

## 7. Verification

- **Quick tier** (`kutta-tests --test-suite-exclude=slow`, fixtures
  present, `KUTTA_FLUENT_EXAMPLES` set): 262 cases and 7,457,289
  assertions, all passed, 91.8 s, peak RSS 2.76 GB.

- **App** (`ctest --test-dir build/app -L app` with `KUTTA_REAL_CLI`,
  `KUTTA_ELBOW_CASE` and `KUTTA_FLUENT_EXAMPLES` set): 10 of 10 passed.
  Real-CLI ran 8 tests, none skipped. No solver result changed, so no
  fixture was re-recorded.
- **Portability lint** (`tools/lint_portability.py`): clean.

### Slow tier

Main, candidate and pd ran concurrently, each in thirds (doctest
`--first/--last` 1-11, 12-22, 23-33), with `-r=xml --success`. pd is the
patched tree with `KUTTA_EXP_FACE=side KUTTA_EXP_PECLET=1 KUTTA_EXP_DIFF=1`
forcing pd into every case. `tools/forensics/compare_slow_assertions.py`
compared each third with main (`slow-{main,cand,pd}-*.xml.gz`,
`slow-comparison-{cand,pd}-*.json.gz`). The XML writes paths as
`<checkout>`, `<main>`, `<variant>` and `<fixtures>`.

| build | assertions passed | failed | against main |
| --- | ---: | ---: | --- |
| main `1bb54a9` | 39,502 | 65 | the known 65 (V6 25, V2 tet 8, V5 plate 6, V5 bump 9, V7 SA 16, V12 fan 1) |
| candidate (this commit) | 39,502 | 65 | 0 moved values, 0 regressions, 0 recoveries |
| pd forced everywhere | 39,507 | 52 | 1,004 moved values, **0 regressions, 13 recoveries** |

**Candidate.** No asserted value moves in any third.
- The comparison reports one missing and one new assertion: V6's
  `CHECK(solved>0)`, which moved from line 226 to 246 because the
  `KUTTA_V6_WALL_DIR` block sits above it. It is the same assertion and it
  passes on both builds.
- V8 passed on all three builds.
- V5 bump's CPU level 2 stopped at the same iteration on main and the
  candidate this time.

**pd.** Its 13 recoveries:
- 12 on V6 (pressure drop, wall-elbow total and wall-inlet total, for both
  starts on CPU and Metal);
- 1 on V7 SA (bump CPU level 2, Cf 0.120 → 0.087, bound 0.1).

The V7 SA recovery comes from a run that stops on its 600 s wall-clock
budget, so it is not attributable to pd. V5 bump and V7 SA CPU level 2
move for the same reason.

Missing assertions are 8 per-iteration `REQUIRE`s, because V6 converges
in fewer chunks, plus the moved `solved>0`.

Every moved value, by case. The largest relative move is shown; the full
lists are in `slow-comparison-pd-*.json.gz`.

| case | moved | largest move |
| --- | ---: | --- |
| V3 Re1000 cavity vs Erturk | 8 | `eu[2]` 6.00e-4 → 6.71e-4 (+12%), all bounds pass |
| V4 Fluent cavity | 698 | interpolation residuals 0 → at most 7.5e-8, bound 2e-6 |
| k-epsilon elbow and SST nozzle | 24 | as section 4 (nozzle drop 0.462% → 0.070%) |
| k-epsilon elbow, second-order pressure | 16 | as section 4 (drop 1.11% → 0.22%) |
| laminar cylinder | 12 | as section 4 (Cl error 1.2e-4 → 3.0e-6 of Cd; drop 0.506% → 0.750%, bound 0.76%) |
| V6 mixing elbow | 42 | drop 1.97% → 0.105% (section 2.2) |
| V2 hex | 10 | as section 2.1 (pressure order 1.821 → 1.881) |
| V2 tet | 10 | pressure order 0.967 → 0.941 (bound 0.914) |
| SIMPLEC Couette | 6 | shift-invariance 3.3e-8 → 1.5e-8 |
| V5 TMR flat plate | 4 | Cd error 0.826% → 0.821% |
| V5 TMR bump | 13 | CPU level 2 (budget-limited) cp error 0.0023 → 0.0056, bound 0.03 |
| V7 SA flat plate and bump | 19 | CPU level 2 (budget-limited) cp error 0.0035 → 0.0071, bound 0.03 |
| FSAE V7.1-V11, V7.3, V9.2, V9.3, V10.2 | 142 | residual-level values (for example 1.8e-12 → 2.2e-12); none near a bound |

Not measured:
- The four FSAE aero cases (Tyrrell wing, Ahmed L1 and y+ ~ 1, S1223) and
  the two transition cases skip on every build here. Their generated meshes
  (`runs/fsae_aero/mesh`, `runs/transition`) are absent from this checkout,
  from the git-archive trees and from the fixture master.
- So pd on external aero is unmeasured in this lane. The previous lane
  measured the differentiable limiter alone there: wing lift +0.23%, Ahmed
  lift +11% against main's SIMPLEC.

Wall time: all nine processes ran at once and finished in about 60
minutes. Third 2 was the longest, as before.

## 8. Reproduce

```sh
# Variants: the candidate tree plus the scratch switches
git apply validation/results/fluent_face_velocity/evidence/scratch_face_velocity_knobs.patch
#   KUTTA_EXP_FACE=side           side blend over the face velocity (KUTTA_EXP_PROJ=0 skips the projected predictor)
#   KUTTA_EXP_PECLET=1|2          weight by (outP+outN)/(aP+aN) or (outP/aP)(outN/aN)
#   KUTTA_EXP_DIFF=1              differentiable limiter on momentum; KUTTA_EXP_FLUXDIFF=1 in the flux only
#   KUTTA_EXP_BJEPS=K             Barth-Jespersen with an epsilon floor K * component range
#   KUTTA_EXP_UNLIM=1, KUTTA_EXP_QUAD=1, KUTTA_EXP_SIMPLEC=1
# V6 with wall export, then the face-by-face split
KUTTA_V6_BACKEND=cpu KUTTA_V6_WALL_DIR=<dir>/pd KUTTA_EXP_FACE=side KUTTA_EXP_PECLET=1 KUTTA_EXP_DIFF=1 \
  kutta-tests -tc='FV V6 Fluent mixing elbow*'
python3 -I tools/forensics/wall_shear.py "$KUTTA_FLUENT_EXAMPLES/mixing_elbow" mixing_elbow <cache> <dir>/pd cpu-seeded
# V2 and the references with the same switches
KUTTA_V2_BACKEND=cpu KUTTA_EXP_FACE=side KUTTA_EXP_PECLET=1 kutta-tests -tc='FV V2 manufactured Navier Stokes hex*'
KUTTA_EXP_FACE=side KUTTA_EXP_PECLET=1 KUTTA_EXP_DIFF=1 kutta-tests -tc='Fluent standard k-epsilon elbow and SST nozzle*'
# Saved-flux fits on another pair (the nozzle needs the cff boundary-first fix)
python3 -I tools/forensics/face_flux.py "$KUTTA_FLUENT_EXAMPLES/flprj_3891/nozzle_para_named_1.cffdb/nozzle_para_named_1-Solve/Base DP" nozzle_para_named_1 <cache>
```

`evidence/runs.tar.gz` holds every V6, V2 and reference log named above.

## Open items

1. Whether `pressure/ho-rc-scheme 2` selects the side-reconstruction
   high-order Rhie-Chow term. This needs a second converged scheme-2
   reference.
2. What Fluent's Standard limiter is, in detail. Kutta's literal
   Barth-Jespersen with a Fluent-like face velocity is farther from Fluent's
   converged V6 than Venkatakrishnan-Wang. But the limiter rpvars do not
   differ between cases, and the differentiable limiter alone fails a
   k-epsilon gate.
3. Why every side-type face velocity moves the k-epsilon elbow (tet,
   first-order momentum, standard pressure) away from Fluent.
4. The k excess in the inlet-pipe wall cells (17-51% in the downstream
   bands, in every run).
5. The SST nozzle's saved-flux floor: 2.7e-3 |F| on quiet faces, which no
   tested formula explains.
6. The V6 wall-elbow shear: pd leaves it at 0.56% against a 0.5% gate.
