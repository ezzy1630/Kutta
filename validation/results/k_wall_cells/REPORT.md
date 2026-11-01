# k in V6's inlet-pipe wall cells

Open item 4 of `fluent_face_velocity/REPORT.md` reports Kutta's wall-cell k
17-51% above Fluent's in the downstream half of V6's hot inlet pipe
(`wall-inlet`). This lane asked why, and whether Kutta deviates there from
Fluent's documented SST and wall treatment.

**Kutta's SST and its ω-based wall treatment match Fluent's saved state
wherever the data can test them. One documented boundary condition was
missing, and it is now the default:** a symmetry plane has zero normal
gradient of every variable. Kutta's k and ω gradients omitted the
symmetry face. Where the hot pipe's wall meets `symmetry-xyplane`, ω came out
10-23% low and the wall-cell k up to twice Fluent's.

- **What the fix changes.**
  - V6 improves on every field, both starts and both backends (CPU and Metal agree):
    - velocity L2 3.3655% → 3.3616%;
    - pressure 2.663% → 2.656%;
    - k 14.71% → 14.48%;
    - ω 2.806% → 2.776%;
    - pressure drop 1.9667% → 1.9653%.
  - The hot-pipe wall-cell k ratio in bands 3-5 goes from 1.313/1.357/1.511
    to 1.275/1.299/1.366. The `wall-inlet` per-cell wall k L2 goes from 8.0% to 6.7%.
  - Every other converged Fluent reference (the k-epsilon elbow and SST nozzle,
    k-epsilon with second-order pressure, the laminar cylinder, the V4 cavity),
    V7 SA and FSAE V7-V11 are bit-identical.
  - The TMR SST, external-aero and transition cases move inside their bounds,
    most at round-off size. Ahmed moves the most: its Coupled/SIMPLEC drag
    agreement goes from 0.043% to 0.175%, bound 1% (section 7).
- **What is left.**
  - The excess away from the symmetry plane is flow-driven. With Fluent's
    velocity held fixed, Kutta's k and ω reproduce Fluent's wall-cell k there to
    0.96-1.06 (section 4).
  - At the corner, about half the excess remains (wall-cell k 1.29 with
    Fluent's velocity). No documented setting explains it (section 5).

Observed means measured here. Inferred means it follows from a measurement
plus an assumption, which is named. Undocumented means no Ansys document gives
the formula.

## Builds and data

- Main is `2a16e1f`, built from `git archive` in a separate tree. The
  candidate is this branch: `9f23db5` (tooling) and `60b2b32` (the fix).
- All builds are Release. They ran on an Apple M4 Pro, macOS 26.6.2, with
  Metal and CPU. The machine was shared.
- The Fluent pair is the corpus `mixing_elbow` (Fluent 22.2, SST, Coupled, poly).
  The slow test runs it with `--approximate wfgc`, as before.
- Every V6 number below is Metal, seeded start, unless marked otherwise. The
  slow tier repeats them on CPU to the digits shown.
- Isolation runs used environment switches in a scratch build:
  - freeze the flow at the seeded state;
  - also freeze ω;
  - scale k production;
  - variant gradient rows.

  They are in `evidence/scratch_k_isolation_knobs.patch`, which applies on
  `60b2b32`. They are not in the commit.

## 1. Settings and what Kutta does with each

`tools/forensics/settings.py` on V6. The table lists every rpvar that enters
the k or ω equation or the wall treatment.

| rpvar | V6 | Kutta | Evidence on V6 |
| --- | --- | --- | --- |
| `kw-sst-on?` | #t | SST | |
| `kw-low-re-mod?` | #f | α* = 1, β* = 0.09 (setup reads it) | μt = ρk/ω exactly in 62% of cells (`fluent-sst-forensics.md` 2.1). Wall-cell k balances with β* = 0.09 to 0.16% (section 3). |
| `kw-not-trans-sst-production-limiter?`, `sst-pk-factor` | #t, 10 | min(G_k, 10ρβ*kω) | Inactive in the pipe: max G_k/(10ρβ*kω) = 0.388 |
| `kw-not-trans-sst-kato-launder?`, `kw-vorticity-pk?` | #f, #f | G_k = μt S² | In limited cells, Fluent's balance needs S to 0.1% (section 3) |
| `kw-sst-strain?` | #f | limiter uses S | established before |
| `kw-shear-mod?` | #t | required #t; read as f_β* = 1 for SST (TG 4.4.3) | Standard k-ω's f_β* is refuted on V6 (section 3) |
| `kw-wall-omega-treatment-r13?` | #t | correlation wall ω, fixed in wall cells | Wall-cell k balance closes (section 3) |
| `sa-kappa`, `kw-beta-i` | 0.4187, 0.072 | wall κ, wall β_i | established before |
| `sst-limit-cross-diff-to-n-times-w-destruct`, `sst-linearize-neg-cross-diff?` | 1e6, #t | clip, implicit negative part | ω only; inactive at 1e6 |
| `curvature-correction?` | #f | rejected if #t | |
| `wall-distance/method`, `quad-correction?` | 2, #t | exact faceted distance | relative L2 2e-5 |
| `sst-pk-lam-factor`, `sst-po-lam-factor`, `kw-sst-min-visc-ratio`, `kw-sst-min-visc/apply-below-f2-value`, `kw-set-wall-w?`, `sst-v61?`, `turb-wf-ustar-limiter?` (coef 0.1) | 0, 0, 0, 0, #f, #f, #t | **not read** | Undocumented. All 55 corpus cases hold these same values, so the data cannot say what they do. u* matches Fluent without a limiter (`fluent-sst-forensics.md` 3.2). |

Ansys Fluent Theory Guide 2024 R2:

- SST, section 4.4.3: Γ_k = μ + μt/σ_k; G̃_k = min(G_k, 10ρβ*kω);
  Y_k = ρβ*kω with f_β* = 1 for SST.
- The k-ω wall boundary conditions use the ω-model near-wall treatment,
  equations 4-413 to 4-426:
  - zero-flux k;
  - blended wall-cell production (4-419 to 4-422);
  - the ω correlation (4-423 to 4-426).
- Fluent User's Guide, symmetry boundary conditions: zero normal velocity and
  zero normal gradients of all variables at a symmetry plane.

## 2. Isolating the k equation

The three runs below start from Fluent's saved state. Bands are the five
equal-face bands of `wall-inlet` along y, from the hot inlet downstream.
"Ratio" is `wall_shear.py`'s shear-weighted Kutta/Fluent wall-cell k.

| run (main, Metal) | band 1 | 2 | 3 | 4 | 5 | V6 k L2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| full solve | 0.992 | 0.898 | 1.313 | 1.357 | 1.511 | 14.71% |
| velocity, pressure and flux frozen at Fluent's state; k, ω solved | 1.036 | 1.014 | 1.122 | 1.135 | 1.212 | 6.88% |
| velocity and ω frozen; only k solved | 1.040 | 1.077 | 1.395 | 1.402 | 1.419 | 11.35% |
| same, k production × 0.99 | | 1.014 | 1.299 | 1.277 | 1.264 | 9.91% |
| same, k production × 0.98 | | 0.989 | 1.209 | 1.163 | 1.126 | 12.14% |

Observed:

- With Fluent's velocity and ω frozen, Kutta's k equation alone still runs
  about 40% high downstream. So the k equation disagrees with Fluent there.
- The hot pipe is near neutral: in layers 2-3, production ≈ dissipation
  (Re ≈ 2,000). A 1% change in production moves downstream k by 10-15
  points. Source-level differences of about 1% are enough to explain the excess.
- Solving ω as well removes most of it: ω compensates.

## 3. The k equation at Fluent's state, term by term

`tools/forensics/k_budget.py` evaluates Kutta's k-equation terms at Fluent's
saved fields with Kutta's operators. It uses Kutta's own strain from the V6
test (`<backend>-saved-sst.csv`). Fluent's residual there is about 0, so any
nonzero R = C + D − P + Y is a disagreement. Full output:
`evidence/k_budget.txt`.

- **Wall cells close.** Wall-cell R/Y is −0.0016 over every wall zone. In the
  pipe's wall cells it is within ±0.002 in bands 1-4 and +0.032 in band 5.
  - This is the full documented treatment:
    - zero-flux k;
    - fixed correlation ω;
    - harmonic production (4-422), using the saved wall shear and κ = 0.4187.
  - A wall cell's μt S² instead of 4-422 gives +0.11 to +0.81.
  - No wall production gives +0.47 to +1.06.
  - **Established on this case:** Fluent uses 4-419 to 4-422 in wall cells.
    That was open in `fluent-sst-forensics.md` 3.7.
- **Production is μt S² with the limiter's S.**
  - In μt-limited cells Fluent's S is exact: S_F = ρ a1 k/(μt F2).
  - Kutta's S/S_F has these medians:
    - whole domain: 1.007, 0.988, 0.996 and 1.004 in layers 2-5, 1.000 deeper;
    - pipe: 1.002.
  - Solving the balance for the S it needs, S_need = √((C+D+Y)/(μt V)),
    returns S_F to medians of 1.000, 1.001, 0.996 and 1.018 in the pipe's
    layers 2-5. So Kutta's convection and diffusion agree with Fluent's cell
    by cell.
- **The pipe's k residual is in unlimited near-wall cells.**
  - Over pipe bands 2-5, R/Y = −0.0197, a net source excess. With S_F in
    the limited cells it is −0.0145.
  - In the unlimited cells of layers 2-3 where production matters
    (P/Y ≥ 0.3), S_need/S_Kutta is 0.889 in the pipe (211 cells) and 0.985
    elsewhere (1,314 cells).
  - **Inferred:** Kutta's strain reconstruction is about 10% high in those
    pipe cells. Fluent does not save S there, and its gradient (WFGC, LSQ
    weights) is undocumented.
  - No convection or diffusion variant tried closes R. The variants were:
    unlimited SOU, first-order upwind, harmonic Γ, μt from face k and ω,
    σ_k fixed, and k-gradient wall rows. They move pipe R/Y between −0.010
    and −0.031.
- **f_β* = 1 holds for SST.** Low-F1 cells (∇k·∇ω > 0) show S_need/S_F
  of 0.90-0.95, which suggested standard k-ω's dissipation factor
  f_β* = (1+680χ²)/(1+400χ²).
  - Binned by its predicted value, Fluent's balance gives
    S_need/S_F = 0.987 / 1.022 (layers 2-4) and 0.997 / 1.017 (deeper)
    where f_β* would be 1.15-1.8.
  - Adding f_β* overshoots to 1.06-1.17.
  - **Established on V6:** Fluent's SST uses f_β* = 1, as TG 4.4.3 says,
    even with `kw-shear-mod? #t`.

## 4. Where the solved k–ω pair departs: the symmetry corner

The cells are split by hops from `symmetry-xyplane` (hop 1 has a symmetry
face) in pipe bands 3-5. The volume-weighted Kutta/Fluent ratio, k / ω:

| run | wall layer | hop 1 | hop 2 | hop 3 | hop > 3 |
| --- | --- | --- | --- | --- | --- |
| main, full solve | L1 | 2.047 / 1.019 | 1.472 / 1.008 | 1.046 / 1.000 | 1.097 / 1.004 |
| | L2 | 1.769 / **0.904** | 1.385 / 1.032 | 1.032 / 1.051 | 1.080 / 1.031 |
| | L3 | 1.376 / **0.774** | 1.266 / 0.971 | 1.110 / 1.091 | 1.080 / 1.017 |
| main, velocity frozen | L1 | 1.530 / 1.003 | 1.249 / 1.001 | 0.956 / 1.000 | 0.964 / 1.000 |
| | L2 | 1.739 / **0.911** | 1.323 / 1.024 | 0.956 / 1.036 | 0.976 / 1.017 |
| | L3 | 1.424 / **0.783** | 1.274 / 0.965 | 1.085 / 1.075 | 1.051 / 1.016 |
| velocity and ω frozen | L1 | 1.471 | 1.415 | 1.410 | 1.370 |
| **fix, full solve** | L1 | 1.659 / 1.015 | 1.457 / 1.007 | 1.094 / 1.000 | 1.116 / 1.004 |
| | L2 | 1.533 / 0.979 | 1.337 / 1.012 | 1.081 / 1.042 | 1.097 / 1.030 |
| | L3 | 1.282 / 0.869 | 1.231 / 0.952 | 1.123 / 1.073 | 1.091 / 1.015 |
| **fix, velocity frozen** | L1 | 1.286 / 1.001 | 1.231 / 1.001 | 0.998 / 1.000 | 0.981 / 1.000 |
| | L2 | 1.499 / 0.981 | 1.273 / 1.004 | 1.002 / 1.028 | 0.993 / 1.017 |
| | L3 | 1.320 / 0.874 | 1.236 / 0.944 | 1.098 / 1.057 | 1.061 / 1.015 |

Observed:

- With Fluent's velocity, Kutta's k–ω pair reproduces Fluent's wall-cell k
  away from the symmetry plane (0.96-1.06).
- The k excess sits in the two cells next to the plane, where ω is 9-22%
  low. Low ω raises μt = ρk/ω and P/Y, so k grows.
- With ω frozen, the excess is uniform around the pipe. So the corner is
  an ω-equation effect.

Kutta's k and ω gradients had no symmetry row. The LSQ omitted every boundary
row whose value is the cell's own, as extrapolated data. At a symmetry plane,
that row is not an extrapolation. It is the documented zero-normal-gradient
condition, and Kutta already uses it for velocity: the cell value at the
normal foot (1f052c3). The fix gives k and ω the same row.

- It also reaches the transition model's γ gradient, which shares the k LSQ.
- It also reaches ω's projected diffusion reconstruction, which copies the ω
  rows.
- Kutta's saved symmetry faces hold the cell value. So does Fluent's: saved
  symmetry-face k and ω equal the cell value, median ratio 1.000000.

At Fluent's saved state, ω's R/Y (R > 0 is an excess sink, so ω too low)
over every layer 2-4 cell with a symmetry face:

- before the fix: +0.174 (Σ|R|/Y 0.259);
- with it: +0.088 (Σ|R|/Y 0.199).

By hop, the first cell improves (L2 +0.332 → +0.156, L3 +0.966 → +0.574).
The second gets worse (L2 +0.041 → +0.129, L3 +0.200 → +0.382). The ω
residual at Fluent's state is large in near-wall layers everywhere (L4,
away from symmetry: −0.46), as in `fluent-momentum-forensics.md` section 6.
So this residual is a noisy check. The case for the change is these three:
the documented condition, the localisation above, and the solve results
(section 6).

## 5. Variants tried and rejected

| variant | V6 velocity / pressure / k L2 | drop | verdict |
| --- | --- | ---: | --- |
| main | 3.3655 / 2.663 / 14.71% | 1.9667% | |
| **symmetry rows for k and ω (shipped)** | **3.3616 / 2.656 / 14.48%** | **1.9653%** | every field closer |
| + k wall rows at the normal foot (∂k/∂n = 0) | 3.3841 / 2.676 / 14.46% | 1.9967% | velocity, pressure, drop worse |
| + strain LSQ symmetry rows at the normal foot | 3.3618 / 2.656 / 14.47% | 1.9665% | no gain |
| standard k-ω f_β* | | | refuted at the saved state (section 3) |
| face Γ_k harmonic, or μt from face k, ω | | | saved-state R unchanged (section 3) |

The remaining corner excess (wall-cell k 1.29 with Fluent's velocity) is not
explained by any documented setting tested here. **Inferred:** the cause is
gradient reconstruction in those cells (the strain LSQ, the omitted WFGC), but
the saved data cannot isolate it.

## 6. V6 after the fix

CPU and Metal agree to the digits shown (slow tier). Seeded start; the hybrid
cold start moves the same way.

| quantity (gate) | main | fix |
| --- | ---: | ---: |
| velocity L2 (0.5%, known failure) | 3.3655% | 3.3616% |
| largest field error, k (0.1%, known failure) | 14.713% | 14.478% |
| pressure drop (0.5%, known failure) | 1.9667% | 1.9653% |
| wall-elbow shear (0.5%, known failure) | 1.635% | 1.641% |
| wall-elbow total (0.5%, known failure) | 0.606% | 0.602% |
| wall-inlet shear (0.5%) | 0.118% | 0.152% |
| wall-inlet total (0.5%, known failure) | 0.522% | 0.534% |
| `wall-inlet` per-face shear L2 / wall-cell U_t L2 / wall-cell k L2 | 4.49 / 4.18 / 8.0% | 4.44 / 4.10 / 6.7% |
| iterations, seeded / cold | 200 / 250 | 200 / 250 |

- The wall-inlet integrated shear moves from 0.118% to 0.152%. Face by face it
  improves (per-face L2 4.49% → 4.44%). Main's 0.12% was already a
  cancellation (`fluent_face_velocity` section 3).
- Wall-inlet total moves the wrong way by 0.012 points, and wall-elbow shear
  by 0.006. Both already failed on main. No gate changes outcome.

## 7. Verification

**Quick tier.** `kutta-tests --test-suite-exclude=slow`, with
`KUTTA_FLUENT_EXAMPLES` and `KUTTA_RUNS` set: 262 cases and 7,457,397
assertions passed. No skip messages.

**App.** Built with `-DKUTTA_BUILD_APP=ON`.

- `mixing_elbow` and `naca4412_ground` (which has symmetry zones) were
  re-recorded with `app/scripts/capture_fixtures.py`. The wing case was
  regenerated with `tools/fsae_aero/naca_ground_case.py`.
- `elbow_param` (k-epsilon) and `cylinder_tutorial` (laminar) do not use this
  path and were not re-recorded.
- `ctest --test-dir build/app -L app` with `KUTTA_REAL_CLI`,
  `KUTTA_ELBOW_CASE` and `KUTTA_FLUENT_EXAMPLES` set: 10 of 10 passed.
- `-R realcli -V`: 8 passed, 0 skipped.

**Portability lint.** Clean.

**Slow tier.**

- Main `2a16e1f` and the fix ran concurrently with `KUTTA_RUNS` set, in
  doctest parts 1-7, 8-11, 12-20 and 21-33, with `-r=xml --success`.
- `tools/forensics/compare_slow_assertions.py` compared each part
  (`evidence/slow-comparison-*.json.gz`).
- The tooling commit adds ten lines to `test_fluent_elbow.cpp`. Main's V6 line
  numbers after the insertion were shifted by +10 before comparing, so each
  assertion keeps its identity.
- The first part 12-20 run found the decompressed TMR level-2 bump grid
  missing in both trees. It is gitignored and normally made by
  `tools/prepare_tmr.py`. It was symlinked from the fixture master
  (identical to the tracked `.gz`) and part 12-20 was rerun on both builds.
- Four cases first printed "skipping: missing ... (restricted reference
  data)" on both builds: the Tyrrell wing, S1223 (SST and transition) and
  T3A. Their gitignored reference data was symlinked from the fixture master
  into both trees. Those four cases were rerun on both builds and compared
  the same way. Their comparison JSONs are not committed, because the
  expanded values carry restricted reference data.

| build | cases run | skip messages | assertions passed | failed |
| --- | ---: | ---: | ---: | ---: |
| main `2a16e1f` | 33 | 0 | 201,025 | 62 |
| fix | 33 | 0 | 201,025 | 62 |

The totals combine the four parts and the four reruns. Against main: 0
regressions, 0 recoveries, 0 missing and 0 new assertions. The 62 failures are
main's known ones. 126 asserted values move:

| case | moved | values (main → fix; bound) |
| --- | ---: | --- |
| V6 mixing elbow, CPU and Metal, both starts | 44 | section 6. Plus each monitor's max velocity, 1.6666 → 1.6674 (bound 100) |
| Ahmed 35 deg, level 1 | 1 | drag guard 0.0198% → 0.0226% of 0.3531 (bound 2%) |
| Ahmed 35 deg at y+ ~ 1, SIMPLEC and Coupled | 4 | SIMPLEC drag guard 0.0034% → 0.0082% (2%); Coupled/SIMPLEC drag 0.043% → 0.175% (1%), lift 1.50% → 1.32% (3%); y+ mean 0.87960 → 0.87958 (1.5) |
| V5 TMR flat plate SST, level 2, CPU and Metal | 4 | cf and Cd errors move by at most 8e-6 (CPU cf 0.008535 → 0.008533; bounds 0.05, 0.07) |
| V5 TMR bump SST, level 2, CPU and Metal | 12 | cf errors move by at most 3e-5 and cp by at most 4e-5. CPU cf at 0.868 is 0.14533 → 0.14533 (fails on both; bound 0.1) |
| Tyrrell 026, SST | 3 | Cl errors move by at most 1e-6 (freestream 1.2357% → 1.2356%) |
| S1223, SST | 2 | CL error 0.0747027 → 0.0747029 (0.10); CD error 22.728% → 22.727% (25%) |
| T3A transition | 53 | laminar 12.0954% → 12.0958% (13%), turbulent 6.9437% → 6.9431% (7.5%), against the independent solution 1.6957% → 1.6960% and 1.3480% → 1.3480% (3%); settle drag 2.3e-7 → 6.1e-7 (5e-4); 48 data-range checks |
| S1223 transition | 3 | CL error 0.0988868 → 0.0988863 (0.10); CD error 16.2315% → 16.2312% (25%); settle lift 7.7e-5 → 7.4e-5 (5e-4) |

- Unchanged:
  - the k-epsilon elbow and SST nozzle, and k-epsilon with second-order pressure;
  - the laminar cylinder, V1-V4, V8 and V12;
  - V7 SA (its transport has its own gradients);
  - FSAE V7.1-V11.
- The Ahmed moves are the largest outside V6. That half-body has a symmetry
  plane. It has no Fluent reference; the bounds are Kutta's own regression
  guards. The 2D wing, S1223 and T3A move at round-off size.
- The S1223 transition CL error sits 0.0011 inside its 0.10 bound on both
  builds; the fix moves it 5e-7 away from the bound.
- Wall times are not compared. The eight parts ran concurrently on a shared
  machine. Every budget is now an iteration count, so the comparison does not
  depend on load.

## Reproduce

```sh
# Kutta's SST terms at Fluent's saved state, and solved fields
KUTTA_V6_BACKEND=metal KUTTA_V6_WALL_DIR=<dir>/fix kutta-tests -tc='FV V6 Fluent mixing elbow*'
# isolation runs (scratch build)
git apply validation/results/k_wall_cells/evidence/scratch_k_isolation_knobs.patch
#   KUTTA_EXP_FREEZE_FLOW=1 [KUTTA_EXP_FREEZE_OMEGA=1] KUTTA_V6_ITERATIONS=1500  frozen velocity [and omega]
#   KUTTA_EXP_PK_SCALE=0.99   KUTTA_EXP_NO_SYM_SCALAR=1 (main's rows)   KUTTA_EXP_WALL_K=1   KUTTA_EXP_SYM_STRAIN=1
python3 -I tools/forensics/k_budget.py "$KUTTA_FLUENT_EXAMPLES/mixing_elbow" mixing_elbow <cache> <dir>/fix metal wall-inlet <dir>/<run>/metal-seeded ...
python3 -I tools/forensics/wall_shear.py "$KUTTA_FLUENT_EXAMPLES/mixing_elbow" mixing_elbow <cache> <dir>/<run> metal-seeded
```

## Open items

1. The remaining corner excess: wall-cell k 1.29 and layer-3 ω 0.87 with
   Fluent's velocity. The likely cause is gradient reconstruction near the
   wall-symmetry corner (inferred, unisolated).
2. The flow-driven excess away from the symmetry plane, about +10% in the full
   solve against 0.96-1.06 with Fluent's velocity. It follows the momentum
   limiter and face-velocity items of `fluent_face_velocity/REPORT.md`.
3. Kutta's strain about 10% above Fluent's implied S in the pipe's unlimited
   near-wall cells (section 3). This is undocumented WFGC and LSQ weighting.
   ω compensates for it in the solved pair.
4. Seven undocumented rpvars (section 1) are not read. Every corpus case holds
   the same values, so nothing can test them.
