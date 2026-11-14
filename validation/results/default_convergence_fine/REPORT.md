# Default `kutta solve` on finer Ahmed meshes

**A default solve now stops by itself on every Ahmed mesh tested, with an
honest verdict.** Seven of eight meshes converge, including all three hex
levels. Native L1 poly stops `statistically_steady`, not converged, on a
±0.09% lift cycle whose mean is within 0.25% of its converged answer. Before
this change, the default reached its 1000-iteration limit on hex L2, hex L3
and the hybrid mesh, with hex L2 and L3 cycling.

Two changes do it:

1. **The starter template (`kutta new-case`) writes momentum relaxation 0.9**
   instead of Fluent's 0.7. Fluent's other defaults stay, including SIMPLEC
   and pressure 0.3. Cases that state their own relaxation, Fluent or Kutta,
   keep it.
2. **A statistically-steady stop.** It catches a bounded force cycle that
   never converges. The run stops unconverged and reports each force's
   Hann-weighted mean, min and max over the last two windows. A cycle that
   drifts, grows or decays never stops this way; a decaying one is still
   converging. Cases with a force report now default to a 5000-iteration cap,
   since they stop by themselves.

Also: a diverging solve now says so. It names the iteration, the last finite
residuals and the likely causes, instead of printing `solve: nonfinite input
vector`.

## Before and after (Metal, M4 Pro, default `kutta solve`)

Setup: the Ahmed 35° case of the car-mesh and tet/poly lanes. That is unit
speed and density, Re_L 2.78e6, Tu 1%, viscosity ratio 1 and A_ref 0.112032,
plus the methods the starter writes. Every mesh has the 15 µm first cell.
- **Before:** main `0ceba3c`. SIMPLEC, momentum 0.7, pressure 0.3, a
  1000-iteration default cap.
- **After:** this commit. Momentum 0.9, a 5000-iteration cap, the averaged
  stop.

Cd and Cl are the final-iteration coefficients. For a statistically steady
run they are the mean, with the range in brackets.

| mesh | cells | before: iterations, stop | before: Cd / Cl | after: iterations, stop | after: Cd / Cl |
| --- | ---: | --- | --- | --- | --- |
| hex L1 | 446,364 | 921, converged | 0.3545 / 0.2199 | 621, converged | 0.3546 / 0.2206 |
| hex L2 | 757,988 | 1000, iteration limit (cycling) | 0.3307 / 0.1472 | 561, converged | 0.3287 / 0.1400 |
| hex L3 | 1,443,692 | 1000, iteration limit (cycling) | 0.3157 / 0.1040 | 1621, converged | 0.3091 / 0.0854 |
| native L1 → poly | 86,215 | 861, converged (residuals waived) | 0.4692 / 0.5251 | 1600, statistically steady (not converged) | 0.4696 [0.4695, 0.4696] / 0.5267 [0.5263, 0.5272] |
| native L2 → poly | 133,976 | 661, converged | 0.4388 / 0.4293 | 881, converged (residuals waived) | 0.4392 / 0.4307 |
| native L3 → poly | 232,438 | 961, converged (residuals waived) | 0.3943 / 0.3115 | 621, converged | 0.3945 / 0.3122 |
| native, level-5 core → poly | 267,235 | 981, converged (residuals waived) | 0.3661 / 0.2529 | 921, converged (residuals waived) | 0.3664 / 0.2540 |
| hybrid: hex wall layer + level-5 poly core | 417,155 | 1000, iteration limit | 0.3622 / 0.2122 | 801, converged (residuals waived) | 0.3623 / 0.2126 |

Notes:
- **Hex L1:** stays within 0.5% (Cd +0.03%, Cl +0.32%). The rows that
  converged before move at most 0.43% (C5 lift).
- **How the after rows were run:** they are the before cases with the setup
  override `controls/mom-relaxation` 0.9. A case file that states momentum
  0.9 gives bit-identical forces and residuals (checked on hex L1 and native
  L3).
- **The averaged stop never fires on a converging run.** Force sampling
  leaves the solve bit-identical: hex L1 with the averaged stop at momentum
  0.7 gives exactly the old 921 iterations, Cd 0.354459, Cl 0.219861.
- **Native L1 at 0.9:** a stationary lift cycle of ±0.09% (0.5263–0.5272)
  that never settles to forces-settled's 0.1% range. Its mean is +0.06% Cd and
  +0.23% Cl from the converged momentum 0.8 result (0.4693 / 0.5255).
- **Before, hex L2 and L3:** the values are one point of a cycle.
  - L2 cycles over Cd 0.3267–0.3349 and Cl 0.1197–0.1647. It is still
    cycling at 4000 iterations (Cd 0.3349, Cl 0.1203).
  - L3 cycles over Cd 0.3129–0.3226 and Cl 0.0720–0.1148.
  - The tet/poly report's 0.3274 / 0.1590 (L2) and 0.3183 / 0.0938 (L3) are
    other points of the same cycles.
- **Before, hybrid:** it settles by forces at 1241 iterations, with
  continuity plateaued near 3e-3. That is past the old 1000 cap. The earlier
  report's "drifting at 2500" row was a `bench_fsae_aero` run, not
  `kutta solve`.
- **With only the averaged stop** (momentum 0.7):
  - hex L2 stops statistically steady at 1820: Cd 0.3309 (0.3267–0.3349), Cl
    0.1422 (0.1197–0.1647). Replaying the final rule on its per-iteration
    history gives the same stop;
  - hex L3 stopped at 2160: Cd 0.3177 (0.3129–0.3226), Cl 0.0931
    (0.0720–0.1148). This was measured with the first version of the rule,
    which also accepted a shrinking range, and was not rerun;
  - every other mesh is unchanged.

## Diagnosis: an unstable steady state, not a plateau

On hex L2 the default's "stall" is a self-sustained oscillation, not a
residual plateau.
- From about iteration 420, lift and drag swing nearly sinusoidally with a
  period of about 44 iterations. Two close lines, 43.5 and 46.2, make it beat
  over about 750 iterations.
- The swing is Cl ±15% and Cd ±1.2%. It is stationary from about iteration
  900 to 4000.
- Fluent-scaled continuity rises with the cycle to about 0.11; it does not
  plateau.
- Forces were recorded every iteration for this. A force evaluation does not
  change the solve: monitor interval 1 and 20 give bit-identical forces.

One switch at a time on hex L2:

| L2 run | behaviour | Cd | Cl |
| --- | --- | --- | --- |
| default: SIMPLEC, momentum 0.7, pressure 0.3 | cycle, period ≈ 44, stationary | 0.3267–0.3349 | 0.1197–0.1647 |
| pressure 1.0 (Fluent's SIMPLEC advice) | same period, growing faster | 0.329–0.338 at 840 | 0.118–0.173 at 840 |
| momentum 0.5, pressure 0.2 | period ≈ 80, still growing | ≈ 0.329 | 0.137–0.144 at 1300, growing |
| differentiable (smooth) limiter | same period-40 mode, smaller | ≈ 0.341 (different fixed point) | 0.1645–0.1700 |
| Coupled with pseudo time, CFL 200 (Fluent's steady default) | small cycle, not settled at 1500 | 0.3286–0.3289 | 0.1397–0.1419 |
| momentum 0.8 | cycle decays; converged at 1261 | 0.3287 | 0.1398 |
| momentum 0.9 | converged at 561 | 0.3287 | 0.1400 |

What this shows:
- **The steady RANS solution exists, but small pseudo-time steps can't reach
  it.**
  - Three different paths reach the same answer: SIMPLEC at 0.8, SIMPLEC at
    0.9, and Coupled's cycle centre.
  - Each SIMPLEC iteration acts like a local pseudo-time step proportional to
    α/(1−α): 1.0 at momentum 0.5, 2.33 at 0.7, 4 at 0.8, 9 at 0.9.
  - The cycle's period in iterations stretches as that step shrinks
    (80 iterations at 0.5, 44 at 0.7).
  - More relaxation slows the cycle but doesn't damp it, the signature of an
    unstable equilibrium rather than an iterative instability.
  - A larger step damps it. That is consistent with implicit pseudo-time
    marching suppressing an unstable mode once the step is large enough.
  - Coupled, with implicit steps at CFL 200, only shows its known small
    Barth-Jespersen limiter cycle
    ([coupled plateau](../coupled_plateau/REPORT.md)).
- **It is not the limiter.** The smooth limiter shows the same mode, smaller.
  It also moves the fixed point, so it is not a candidate.
- **It is not the stop criterion.** A ±15% lift swing is not a bounded
  near-converged state that a looser tolerance should accept.
- **Linear-solver tolerances were not varied.**
  - They are already tighter than Fluent's default AMG termination of 0.1:
    momentum 1e-2, pressure 1e-3.
  - The relaxation scaling above points elsewhere.

### Why 0.9 and not 0.8

Hex L3 is more unstable than L2, and of the two only 0.9 converges it:

| mesh | momentum 0.8 | momentum 0.9 |
| --- | --- | --- |
| hex L1 | 821, converged: 0.3545 / 0.2201 | 621, converged: 0.3546 / 0.2206 |
| hex L2 | 1261, converged: 0.3287 / 0.1398 | 561, converged: 0.3287 / 0.1400 |
| hex L3 | 1600, statistically steady: 0.3151 [0.3090, 0.3211] / 0.0912 [0.0612, 0.1212] | 1621, converged: 0.3091 / 0.0854 |
| native L1 → poly | 841, converged (waived): 0.4693 / 0.5255 | 1600, statistically steady: 0.4696 [0.4695, 0.4696] / 0.5267 [0.5263, 0.5272] |
| native L2 → poly | 801, converged (waived): 0.4392 / 0.4305 | 881, converged (waived): 0.4392 / 0.4307 |
| native L3 → poly | 601, converged: 0.3943 / 0.3116 | 621, converged: 0.3945 / 0.3122 |
| native, level-5 core → poly | 861, converged (waived): 0.3662 / 0.2533 | 921, converged (waived): 0.3664 / 0.2540 |
| hybrid | 941, converged (waived): 0.3622 / 0.2124 | 801, converged (waived): 0.3623 / 0.2126 |

- **At 0.8, hex L3's cycle grows to Cl 0.063–0.118 by iteration 700 and then
  holds.** Its statistically steady mean is +1.9% Cd and +6.8% Cl from the
  steady solution 0.9 converges to, with a ±33% lift range. That is a
  materially wrong answer, labelled honestly but still wrong.
- **At 0.9 the one unconverged mesh is native L1.** Its cycle is ±0.09% of
  lift and its mean within 0.25% of the converged 0.8 result.
- So 0.9 has the smaller worst case, converges every hex level, and is faster
  on the hex meshes (621 vs 821 iterations on L1, 561 vs 1261 on L2).
- **On hex L3 at 0.9, convergence is slow but steady.** A decaying lift
  oscillation narrows from a 1.3% range over 300 iterations at 1100 to 0.2% at
  1500, and the run converges at 1621.
  - The first version of the averaged stop accepted any range that did not
    grow. It stopped this run "statistically steady" at 1600 while it was
    still converging. The final rule requires the range to stay within ±10%,
    so a decaying cycle runs on.
- Values between 0.8 and 0.9 were not tried.

## What Fluent does by default

From the Fluent 2025 R2 User's Guide, §37.3, "Pressure-Based Solver
Settings":
- Coupled with a pseudo-time method "is the default for cases with" steady,
  single-phase flow, and "the coupled solver is usually recommended for
  steady-state simulations".
- On SIMPLEC: "the pressure-correction under-relaxation factor is generally
  set to 1.0". Also: "it is usually not necessary to reduce the pressure
  under-relaxation for SIMPLEC."
- The default under-relaxation factors are "near optimal for the largest
  possible number of cases".

Fluent 2024 R2 User's Guide §41.2: report definitions offer "Average Over" N
iterations, a running mean whose default 1 is the instantaneous value.
That is the Fluent user's tool for an oscillating steady monitor; the
statistically-steady stop automates it.

How Kutta relates:
- **Coupled:** Kutta's Coupled reproduces Fluent's default here. It lands at
  the right answer with a small cycle, but it limit-cycles on hex L1 (lift
  range 2.3% at 3000, per the car-mesh report). So it stays out of the
  template.
- **Pressure 1.0:** it makes the L2 cycle grow faster.
- **Momentum 0.9:** a departure from Fluent's SIMPLE/SIMPLEC dialog default
  of 0.7. The case states it, so a user sees it in `inspect` and in the app's
  Controls.

## The statistically-steady stop

`src/kutta/setup/convergence.*`, `app/solve.cpp`; contract in
[docs/app-protocol.md](../../../docs/app-protocol.md#convergence-criteria-v1-additive).
- **Sampling:** from iteration 2 × forces-settled window + 1 (601), every
  force report is evaluated at every iteration. Stride samples would alias a
  period-44 cycle.
- **The test:** at each monitor sample n, compare iterations n−999…n−500 with
  n−499…n. The run stops when, for every report:
  - the Hann-weighted means differ by ≤ 0.1% of the forces-settled scale;
  - the second range is within ±10% of the first, plus 0.1% of the scale.
- **Defaults:** `monitors/forces-averaged/{active,tolerance,range,window}`
  = true, 1e-3, 0.1, 500.
- **Outcome:** `stop_reason: "statistically_steady"`, `converged: false`.
  `convergence.averaged` holds the window, span, and per-report
  mean/min/max/drift/range_change plus coefficients.
- A run that hits its limit while cycling still reports that state, with
  `steady: false`. The app shows "Statistically steady; not converged".
- **Why Hann weighting:**
  - A plain window mean of a cycle is off by up to amplitude × period /
    (π × window) from the cut at its ends.
  - On L2's per-iteration history the plain means of adjacent 500-iteration
    windows never come within 2e-3 for lift.
  - The Hann-weighted ones reach 2e-5 to 5e-4.
- **Why 500:** a 1e-3 tolerance between adjacent windows admits at most 0.2%
  drift per 1000 iterations, tighter than forces-settled's 0.33%. At
  W = 300 lift needed until about 1500 just to get under the tolerance.
- **Why the range must not shrink:** a decaying cycle is converging. Hex L3
  at 0.9 is the case (above). A synthetic decaying cycle in
  `tests/setup/test_convergence.cpp` stops at 600 under a grow-only test and
  converges at 1941 under this one.
- **Accuracy on L2 (momentum 0.7):** the live run stops at 1820, exactly
  where the offline simulation of the rule on the per-iteration history does.
  Its means (Cd 0.33090, Cl 0.14224) are within 1e-4 and 2e-4 of the means
  over iterations 1321–2796.
- **The caveat this verdict carries:** a limit cycle's mean is not the
  steady answer, and its range is not a physical unsteadiness band.
  - On L2 at 0.7 the cycle mean is 0.7% (Cd) and 1.6% (Cl) above the steady
    solution that 0.8 and 0.9 converge to.
  - On L3 at 0.8 it is 1.9% and 6.8% above.
  - Treat a statistically steady result as an estimate. The cases'
    `forces.csv`, `loads.json` and `result.vtu` are the last iteration, one
    point of the cycle.

## Divergence message (`--unlimited`, `--poor-cell-quality 0`)

Both are real divergences, not a backend bug.
- **CPU, native L1 poly:** `--unlimited` is finite after iteration 1
  (residual 1.10) and nonfinite by iteration 30.
- **Metal, the same mesh:** `--poor-cell-quality 0` reaches Cl ≈ −6.5e9 by
  iteration 50.
- **`kutta solve` with `gradient_limiter: none` on hex L1:** x-momentum
  residual goes 1.9e-2 → 4.1 → 18 over iterations 2–4, then nonfinite at 9.
  The hex L1 mesh also has 2,127 cells at orthogonal quality ≤ 0.05. Unlimited
  second-order reconstruction on 15 µm × 40 mm wall-layer cells is unstable;
  Fluent's defaults keep the limiter and poor-mesh numerics on.
- **Why Metal fails one iteration earlier:** its reductions accumulate in
  compensated float pairs, since Metal has no FP64. So `dot(b,b)` overflows
  once |b| passes about 1.8e19; CPU and CUDA accumulate in double.
- **The fix:** a typed `NonfiniteError` from the linear solver's entry check.
  `kutta solve` turns it, and non-finite residuals or forces, into "the
  solution diverged at iteration N (…). Last finite residuals, iteration M:
  …". It then gives the limiter hint when limiting is off and the poor-cell
  count. Exit code 1, with a JSONL `error` event.
- Coupled's own finite checks under `src/kutta/fv` are unchanged and keep
  their message.

## Gates

All on the final tree, macOS (SDK 26.5), M4 Pro.
- **Quick tier** (`ctest -LE slow`): kutta-tests 328/328 cases, 7,996,959
  assertions, including CPU/Metal parity and the CLI protocol check
  (`tests/app/check_protocol.py`). `tools/lint_portability.py` passes.
- **App tests:** 13/13.
- **Slow tier:** nothing under `src/kutta/fv` or the kernels changed, so the
  full slow tier was not run.
  - The four slow CLI cases that convert, write and solve Kutta cases were
    run against main `0ceba3c` with the same fixtures. They are "Kutta cases
    converted from Fluent cases…", "The NACA ground-effect wing…", "Native
    FSAE Fluent physics…" and "Corpus FSAE physics…".
  - Result: 4/4 pass with 546 assertions, and
    `tools/forensics/compare_slow_assertions.py` finds 0 regressions and 0
    changed values.
  - Its 21 missing/new pairs are the same assertions moved two lines by the
    edit to `tests/app/test_kutta_case_cli.cpp`. "Native FSAE Fluent physics"
    now clears the starter's relaxation along with its coupling, because a
    Coupled pseudo-time case rejects classic under-relaxation.
- **Not run:** CUDA (no NVIDIA machine).

## Reproduce

Meshes are regenerated, not committed.
- Native L1 is sha256 `35cee20e…`.
- The hex L1 mesh is `runs/fsae_aero/mesh/ahmed35_L1_first15um.msh`.
- The hybrid replaced each prism of the level-5-core `kutta polyhedra`
  output with `subdivide.py`'s three hexes, using a scratch script that is
  not committed. It reproduces the tet/poly lane's 417,155 cells, with volume
  and body area unchanged.

```sh
uv run --with gmsh==4.15.2 --with numpy python tools/fsae_aero/ahmed.py --slant 35 --first 1.5e-05 --level 3 --out hex_L3.msh
uv run --with gmsh==4.15.2 --with numpy python tools/fsae_aero/ahmed.py --slant 35 --first 1.5e-05 --level 2 --cell-type tet-prism --out nat_L2.msh
kutta polyhedra nat_L2.msh --out nat_L2_poly.msh
kutta new-case hex_L3.msh --speed 1 --out ahmed_L3.kutta.json  # then the Ahmed setup above
kutta solve ahmed_L3.kutta.json --backend metal --out run       # converges at 1621
```
