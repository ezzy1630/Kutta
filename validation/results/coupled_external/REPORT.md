# Coupled on external aero: wing in ground effect, S1223, Ahmed body

Kutta's native Coupled solver diverged within 5-57 iterations on every FSAE
2D mesh, where SIMPLEC converged. It now converges on all seven wing and
S1223 points tried, to SIMPLEC's answer. It takes 700-1300 iterations against
SIMPLEC's 1800-2850, and 1.6-9x less wall time on Metal. The Ahmed level-1
hex mesh still diverges with Coupled (open item 1). The fixed point is
unchanged, SIMPLE and SIMPLEC are bitwise identical to main, and V6 takes the
same or fewer Coupled iterations.

- Branch `t3/coupled-external-convergence`, base main `2747dd2`. The commit
  carrying this report holds the change.
- Machine: M4 Pro, 24 GB, Metal and the serial CPU reference backend. Other
  lanes shared the CPU and GPU, and the FSAE runs went two at a time, so all
  times are contended. Iteration counts are deterministic.
- Evidence is in [`evidence/`](evidence), with home paths replaced by
  `<checkout>`. Meshes, cell fields and wall coordinates are not committed:
  they carry the restricted geometry.
- The diagnosis and ablation runs used a diagnostic build that is not
  committed. That build added environment switches for each part of the
  change, plus a dump of the block RHS, the update and the matrix diagonals.
  Every number for the final state comes from the committed code.

## Reproduction (main `2747dd2`)

Meshes are the FSAE lane's level-2 `wing2d.py` meshes (and `ahmed.py` level
1), run through `bench_fsae_aero --coupled` from a uniform freestream on Metal.
The table gives the last iteration completed before a non-finite residual or
linear-solver input (`evidence/reproduction.tar.gz`; Ahmed in
`evidence/startup.tar.gz`, `ahmed/main_*`):

| mesh | relaxation 0.7/1.0 (bench default on main) | 0.5/0.5 (Fluent's Coupled default) |
|---|---|---|
| Tyrrell 026 freestream | 7 | 17 |
| h/c 0.448 | 7 | 57 |
| h/c 0.224 | 8 | 12 |
| h/c 0.134 | 5 | 9 |
| h/c 0.090 | 7 | 8 |
| S1223 4.76 deg | 8 | 10 |
| S1223 10 deg | 10 | 10 |
| Ahmed 35 deg L1 | 2 | 3 |

SIMPLEC converges on all of them.

## Diagnosis

The candidate causes in the brief, in the order they were ruled in or out, on
the h/c 0.224 wing with main's bench defaults unless stated
(`evidence/startup.tar.gz`, `ruled_out/`, rerun for this report):

1. **Not the linear solver.** Solving every block system to 1e-4 instead of
   Fluent's 0.1 (about 90 BiCGStab iterations each) still diverges at
   iteration 8. The outer iteration is unstable even with near-exact solves.
2. **Not the pseudo-time step.** Kutta's global step follows Fluent (Theory
   Guide 23.6.2; checked by the previous lane): Δτ = 0.3 L / max(0.1 U_bc,
   U_mean, U_Δp). L = min(V^(1/3), largest extent), which is 12.8 c for the
   70c x 30c domain extruded at unit depth, so Δτ = 3.85 chord transit times.
   Smaller steps only delay the blow-up by a few iterations: time scale 0.1
   diverges at iteration 10 and Δτ = 0.001 at 15, against 8. Fluent's
   equations on the 23.6.2 page are images, so whether its 2D length scale
   uses the cube root of the unit-depth volume could not be checked from the
   text. Since a step thousands of times smaller still diverges, it does not
   change the conclusion.
3. **Not turbulence, the boundary rows or initialization.** Laminar flow at Re
   4.62e5 diverges at iteration 8, while laminar Re 1e2 and 1e4 stay bounded
   for 300 iterations. The freestream mesh fails the same way without a moving
   ground (table above). Hybrid initialization diverges at 10 and explicit
   relaxation 0.2/0.2 at 17. First-order convection also diverges (about
   iteration 40), so the second-order deferred correction is not the cause.
4. **The fixed point itself is stable.** Seeded from a converged SIMPLEC
   solution (freestream), Coupled stays bounded for 400 iterations with lift
   unchanged to 1e-5. The failure is on the start-up path.
5. **Sink cells make the momentum Jacobian singular.** The diagnostic build
   took one Coupled step of α = 0.05 from the iteration-6 state and dumped
   the block RHS, the update and the residual after the step. That gives the
   measured change of the true residual, which can be compared with the
   Jacobian's prediction cell by cell (`evidence/startup.tar.gz`,
   `jacobian_check*`). The assembled diagonal is a_P = Σ(D + max(F, 0)), so
   it carries the net outflow. In cells with inflow through nearly every face,
   a_P was 2e-4 to 3e-2 of Σ|a_nb|. Those were the five largest updates:

   | cell (x/c, y/c) | update \|Δu\| | a_P (with Δτ) | Σ\|a_nb\| | net outflow | continuity residual |
   |---|---|---|---|---|---|
   | (0.652, -0.051) | 66.9 U | 2.2e-5 | 6.1e-2 | -6.1e-2 | 1.16 |
   | (1.136, 0.044) | 43.5 U | 2.3e-5 | 8.0e-3 | -8.0e-3 | 0.15 |
   | (0.598, -0.061) | 31.7 U | 2.4e-5 | 3.3e-3 | -3.3e-3 | 0.063 |
   | (0.671, -0.046) | 26.9 U | 3.4e-5 | 1.5e-1 | -1.5e-1 | 2.95 |
   | (1.001, 0.058) | 23.7 U | 1.5e-5 | 4.9e-4 | -4.7e-4 | 0.009 |

   The next three (17-19 U, at the trailing edge) have normal diagonals. The
   worst continuity mismatch is another sink cell at (0.779, -0.023), with
   a_P = 2.6e-5 and Σ|a_nb| = 3.2e-3. Its continuity residual was 0.061 and
   the Jacobian predicted a change of +0.003, but the true change was -1.11.
   At that state, 5755 of the 195,581 cells had a_P < Σ|a_nb| and 34 had less
   than 0.2 of it. Overall, the true change was 1.9x the prediction for u and
   v and 5.1x for continuity. The residual keeps this term, and it vanishes at
   convergence. The Jacobian inherits it. SIMPLEC is shielded from it by its
   implicit relaxation (a_P/0.7). Coupled's relaxation is explicit and never
   touches the Jacobian.
6. **Recirculation needs Fluent's Courant-number relaxation.** With sink
   cells floored, the h/c 0.224 wing converged at Fluent's 0.5/0.5. Closer to
   the ground the runs then wandered without diverging. At h/c 0.134 lift swung
   between 0.06 and 0.98, and raw continuity stayed at about 0.1-0.24, against
   SIMPLEC's converged 1.07e-3. Where streamlines close (the separated lower
   surface, the blunt trailing edge's wake, the gap above the moving ground),
   upwind rows have zero row sums. Wherever the global step is small relative
   to the local convective time scale, the momentum block is singular there.
   Fluent's pressure-based Coupled algorithm relaxes its equations implicitly
   with a flow Courant number: a_P φ/α = Σ a_nb φ_nb + b + (1-α)/α a_P φ_old
   with (1-α)/α = 1/CFL. The guide calls this "equivalent to the
   location-specific time step" ([Fluent 12.0 Theory Guide 18.4.4, eqs.
   18.4-23 and 18.4-24][relax]). Kutta's Coupled had no such term. With
   Fluent's default CFL of 200, h/c 0.134 converges to lift 1.46736 (SIMPLEC
   1.46746) by iteration 800; CFL 50 gets there by iteration 900
   (`evidence/ablation.tar.gz`).
7. **The flux block did not follow the residual's face velocity.** The
   residual interpolates face velocity with Fluent's momentum weights,
   u_f = (a_P u_P + a_N u_N)/(a_P + a_N) ([Theory Guide 18.4.2, eq.
   18.4-6][cont]). That is `face_velocity_weights()` in `fv_momentum_faces`
   and the momentum face reconstruction. The Jacobian's continuity-velocity
   block used equal halves, so its own-cell part summed to zero in every
   interior cell. With sink cells floored and no Courant term, equal halves
   still diverge at 0.5/0.5, while the weighted block converges. With the
   Courant term, both converge (table below).
8. **Under-relaxation.** The FSAE bench passed SIMPLEC's factors (momentum
   0.7, pressure 1.0) to Coupled as explicit relaxation. Fluent's
   pseudo-transient Coupled factors are 0.5 and 0.5, which is what its
   Coupled cases save: the V6 elbow stores
   `pressure-coupled/{mom,pressure}/pseudo-explicit-relax 0.5`. Kutta reads
   those from a Fluent case, so only the library bench was affected.

[relax]: https://www.afs.enea.it/project/neptunius/docs/fluent/html/th/node374.htm
[cont]: https://www.afs.enea.it/project/neptunius/docs/fluent/html/th/node372.htm

**Not the cause, and unchanged:** the smoother and multigrid on stretched
cells. With the fix, each block solve reaches Fluent's 0.1 in 2-6 BiCGStab
iterations after the first hundred iterations (ablation traces). The
Schur-bounded smoother from the previous lane was not touched.

## The change

All three parts are Jacobian-only. The block RHS is still the steady SIMPLEC
residual (delta form), so the fixed point is unchanged.

1. **Momentum diagonal** (`coupled_copy_diagonal`):
   - the physical a_P is floored at Σ|a_nb| (Patankar's M-matrix rule);
   - the pseudo-time term is kept;
   - Fluent's Courant relaxation adds max(a_P, Σ|a_nb|)/CFL, with
     `CoupledOptions::courant = 200`, Fluent's default.

   At convergence, the amount the floor adds equals the continuity residual,
   which vanishes.
2. **Continuity-velocity block** (`coupled_velocity_block`,
   `coupled_fine_row`, `coupled_galerkin_fine`, `coupled_fine_schur`):
   - H_ij = wc s S_f a_j/(a_i + a_j), with the weights the residual uses
     (`face_velocity_weights()`; distance weighting reproduces the old
     halves).
   - The own-cell part is refreshed every solve into the slots the setup used
     to fill.
   - `coupled_galerkin_fine` hit Metal's 31-argument limit, so it now finds
     its diagonal entry by the existing row search instead of taking a `diag`
     array.
3. **FSAE bench** (`validation/cases/fsae_aero.hpp`):
   - Relaxation 0 now means the scheme's Fluent default: 0.5/0.5 for Coupled,
     and still 0.7/1.0 for SIMPLEC, so SIMPLEC runs are unchanged.
   - `history.csv` gains `continuity_raw`.

No new device buffers: `CoupledSolver::bytes()` is unchanged (elbow 393.8
B/cell, nozzle 390.5).

**Ablation** on the h/c 0.224 wing, from a uniform start, Metal, 800
iterations; the floor + CFL 20 run at 0.7/1.0 went to 4000
(`evidence/ablation.tar.gz`). The diagnostic build switched each part with
environment variables. In the file names, `dom` is the floor, `h05` the equal
halves, `a05` relaxation 0.5/0.5 and `cflN` the Courant number. Files without
`dom` or `nodom` that date from after the floor was added (`cfl10`, `cfl30`,
`cfl50`, `cfl200_a05`, `cfl50_a05`) have the floor on. The earlier `cfl20`
has it off.

| H block | sink floor | Courant term | explicit α (u/p) | outcome |
|---|---|---|---|---|
| halves (main) | - | - | 0.7/1.0 | diverges, iteration 8 |
| halves | - | - | 0.5/0.5 | diverges, 12 |
| weighted | - | - | 0.5/0.5 | diverges |
| halves | yes | - | 0.5/0.5 | diverges |
| weighted | yes | - | 0.5/0.5 | converges here; wanders closer to the ground (item 6) |
| weighted | - | CFL 200 | 0.5/0.5 | diverges |
| halves | yes | CFL 200 | 0.5/0.5 | converges |
| **weighted** | **yes** | **CFL 200** | **0.5/0.5** | **converges** (the change) |
| weighted | yes | - | 0.7/1.0 | diverges, 28 |
| weighted | yes | CFL 10, 20, 30 | 0.7/1.0 | converges |
| weighted | yes | CFL 50 | 0.7/1.0 | diverges |
| weighted | - | CFL 20 | 0.7/1.0 | diverges |

## Results

### Coupled vs SIMPLEC on the FSAE meshes

These runs used `run_fsae.sh`, then `summarize.py` (`evidence/fsae.tar.gz`,
`fsae_table.md`). They ran on Metal at level 2 (Ahmed level 1), from a uniform
freestream start, under the FSAE lane's conditions.

Both schemes stop on the bench criterion: Fluent-scaled residuals <= 1e-3, and
lift and drag settled to 5e-4 over 500 iterations (1000 for Ahmed). Coupled
stopped on that criterion at 900, 850, 750, 1300, 800 and 1150 iterations. The
h/c 0.448 run reached its 2000-iteration budget instead (next subsection).

The Coupled column is the first sample that meets all of these:
- the same settle rule;
- u, v, k and ω at or below 1e-3;
- raw continuity at or below SIMPLEC's final raw value.

That is the like-for-like comparison. Times are solve seconds at that
iteration, and the forces are each run's final values.

| case | cells | SIMPLEC its / s | Coupled its / s | speedup | SIMPLEC lift / drag | Coupled lift / drag | Δ lift | Δ drag |
|---|---|---|---|---|---|---|---|---|
| Tyrrell free | 185,964 | 1800 / 389 | 700 / 82 | 4.7x | 0.861225 / 0.016140 | 0.861297 / 0.016143 | +8.4e-5 | +2.2e-4 |
| h/c 0.448 | 227,294 | 2850 / 858 | 750 / 95 | 9.0x | 1.140300 / 0.021483 | 1.140385 / 0.021482 | +7.5e-5 | -2.1e-5 |
| h/c 0.224 | 195,581 | 2400 / 333 | 750 / 64 | 5.2x | 1.401467 / 0.030513 | 1.401503 / 0.030511 | +2.6e-5 | -6.6e-5 |
| h/c 0.134 | 179,126 | 2000 / 308 | 750 / 90 | 3.4x | 1.467457 / 0.045194 | 1.467365 / 0.045195 | -6.3e-5 | +6.5e-6 |
| h/c 0.090 | 177,312 | 2600 / 581 | 1300 / 162 | 3.6x | 1.410692 / 0.059040 | 1.410205 / 0.059060 | -3.5e-4 | +3.3e-4 |
| S1223 4.76 deg | 184,470 | 2150 / 332 | 800 / 109 | 3.1x | 1.584707 / 0.030069 | 1.584693 / 0.030069 | -8.4e-6 | +2.1e-7 |
| S1223 10 deg | 183,892 | 2100 / 198 | 1150 / 123 | 1.6x | 1.933725 / 0.047736 | 1.933738 / 0.047745 | +6.6e-6 | +1.9e-4 |
| Ahmed 35 deg L1 | 403,604 | 5000 / 683 † | fails at iteration 9 | | 0.209024 / 0.353173 | | | |

† Residuals met; lift still drifting 0.15 %, as in the FSAE lane's run.

Every force difference is inside the 5e-4 settle tolerance that stopped both
runs. On Metal, Coupled's cost per iteration (88-139 ms) is at or below
SIMPLEC's (94-301 ms) on these meshes. The wall-time gain therefore comes from
fewer iterations that are also cheaper.

**Fixed point, longer runs** (`evidence/fields.tar.gz`). SIMPLEC ran 4000
iterations and Coupled 2000, with `--settle 0` so neither stopped early. The
field columns are volume-weighted relative L2 differences.

| h/c | SIMPLEC lift / drag | Coupled lift / drag | Δ lift | Δ drag | velocity | pressure | k | ω | max wall Δcp |
|---|---|---|---|---|---|---|---|---|---|
| 0.224 | 1.4015132 / 0.0305099 | 1.4015109 / 0.0305115 | -1.6e-6 | +5.2e-5 | 1.5e-5 | 6.0e-4 | 5.5e-4 | 4.0e-5 | 3.3e-4 |
| 0.090 | 1.4105407 / 0.0590417 | 1.4104751 / 0.0590676 | -4.6e-5 | +4.4e-4 | 8.3e-5 | 4.4e-3 | 3.7e-3 | 2.0e-4 | 2.6e-3 |

At h/c 0.090, SIMPLEC's forces are still changing at the 5e-4 level at 4000
iterations, so that row bounds the remaining drift rather than a gap between
fixed points. The pressure norm is relative to a field whose far-field values
are near zero, which inflates it against velocity.

**CPU.** h/c 0.224 Coupled ran 300 iterations on the serial CPU backend. Lift
was 1.401545, against Metal's 1.401531 at the same iteration, in 223 s against
22 s.

### Fluent-scaled continuity asks more of Coupled

Fluent scales continuity by the largest raw value of the first five
iterations. On these meshes that scale is 25-112 for SIMPLEC, whose first
iterations violate continuity badly. For Coupled it is 0.48-1.07, because its
start-up already nearly satisfies continuity. The same 1e-3 scaled criterion
therefore asks Coupled for raw continuity 25-230 times lower than it asks of
SIMPLEC:

| case | SIMPLEC scale / final raw | Coupled scale / final raw |
|---|---|---|
| free | 42.3 / 8.9e-4 | 0.49 / 4.8e-4 |
| h/c 0.448 | 68.6 / 1.35e-3 | 0.55 / 6.2e-4 |
| h/c 0.224 | 82.4 / 1.13e-3 | 0.48 / 4.8e-4 |
| h/c 0.134 | 95.8 / 1.07e-3 | 0.52 / 3.5e-4 |
| h/c 0.090 | 112 / 8.3e-4 | 0.49 / 4.4e-4 |
| S1223 4.76 | 25.5 / 7.8e-4 | 0.99 / 5.2e-4 |
| S1223 10 | 36.9 / 1.16e-3 | 1.07 / 5.3e-4 |

Coupled ends every run with less raw imbalance than SIMPLEC. At h/c 0.448 it
ran out of budget at 2000 iterations, with scaled continuity 1.13e-3 and raw
continuity 6.2e-4, half of SIMPLEC's final value. The previous lane found the
same scale gap in Fluent's own Coupled histories. This is why `summarize.py`
compares raw continuity.

### V6 and the other Coupled cases, before and after

`validation/results/coupled_convergence/measure.sh`, two repetitions with the
order alternating (`evidence/v6_measure.tar.gz`, `v6_measure_summary.md`):

| criterion | case | backend | Coupled before | Coupled after | SIMPLEC |
|---|---|---|---|---|---|
| Fluent residuals | elbow | CPU | 160 / 15.1-15.3 s | 160 / 15.5-15.7 s | 240 / 26.2-26.4 s |
| Fluent residuals | elbow | Metal | 160 / 2.2-2.3 s | 160 / 2.2-2.3 s | 260 / 3.8-3.9 s |
| Fluent residuals | nozzle | CPU | 60 / 7.4 s | 60 / 7.7 s | 80 / 10.2-10.3 s |
| Fluent residuals | nozzle | Metal | 60 / 1.2-1.3 s | 60 / 1.2-1.3 s | 80 / 1.2-1.3 s |
| V6 (strict) | elbow | CPU | 230 / 24.0-24.1 s | **220** / 23.3-23.5 s | 500 / 62.2-62.4 s |
| V6 (strict) | elbow | Metal | 230 / 3.9-4.1 s | **220** / 3.7-3.8 s | 500 / 9.7-10.0 s |
| V6 (strict) | nozzle | CPU | 70 / 9.2-9.3 s | 70 / 9.4-9.5 s | 150 / 21.2-21.3 s |
| V6 (strict) | nozzle | Metal | 70 / 1.6 s | 70 / 1.6 s | 150 / 2.9-3.1 s |

The V6 slow test, which samples every 50 iterations, gives the same counts
before and after: 200 iterations seeded and 250 from a cold hybrid start, on
CPU and Metal. Times are contended. Main ran its Metal V6 while the FSAE batch
held the GPU; the candidate ran on a quieter GPU.

| V6 slow test | CPU seeded / cold | Metal seeded / cold |
|---|---|---|
| main | 200 / 250 its, 21.0 / 26.7 s | 200 / 250 its, 9.2 / 15.5 s |
| this branch | 200 / 250 its, 21.0 / 26.3 s | 200 / 250 its, 3.5 / 4.3 s |

**TMR flat plate**, level 2, Coupled, 600 iterations, Metal
(`evidence/plate.tar.gz`):
- Drag matches main to 2e-7 for SST (0.0028768362 vs 0.0028768368) and 3e-6
  for SA (0.0028803329 vs 0.0028803253).
- Raw continuity is 8.4e-8 to 9.0e-8 in both builds, and time per iteration is
  the same.

## SIMPLEC and SIMPLE are bitwise unchanged

`parity.sh <main kutta> <this kutta> <out>` runs `kutta solve --iterations 40`
and compares `fluxes.csv`, `forces.csv`, `loads.json`, `result.vtu` and
`residuals.csv` without its timing field (`evidence/parity.txt`):

```
elbow_simplec cpu identical Pressure-velocity coupling: SIMPLEC
elbow_param_simple cpu identical Pressure-velocity coupling: SIMPLE
elbow_simplec metal identical Pressure-velocity coupling: SIMPLEC
elbow_param_simple metal identical Pressure-velocity coupling: SIMPLE
```

The elbow runs with `--approximate wfgc` and the SIMPLEC override, and
elbow_param with `--approximate compressibility,ke-options`.
`steady_solver.cpp` changes only inside `coupled_step`.

## Tests

- **Quick tier** (`kutta-tests --test-suite-exclude=slow`, fixtures linked;
  `evidence/quick_tier.tar.gz`):
  - Main and this branch both pass all 255 cases / 7,453,826 assertions.
  - One earlier full run on this branch failed one assertion, in "Mixing
    elbow solve stops on its convergence condition and kutta aero reads its
    loads". The CLI emitted no iteration event, with exit code 0.
  - It passed in two doctest reruns, five direct runs of `check_aero.py`
    under load, and the final full run, so it is not reproducible. That log is
    kept as `final_first_run_with_flake.txt`.
- **`tests/fv/test_coupled.cpp`** passes unchanged: 2 cases, 491 assertions,
  the same as main.
- **Slow tests:** `FV V6*`, `Fluent standard k-epsilon*` and `Fluent laminar
  2D cylinder*`, main against this branch, compared with
  `tools/forensics/compare_slow_assertions.py` after the address and continuity
  normalization the previous report describes (`evidence/slow.tar.gz`):

  | | count |
  |---|---|
  | assertions passed / failed (both builds) | 19,791 / 25 |
  | regressions / recoveries / pass-fail changes | 0 / 0 / 0 |
  | missing / new assertions | 0 / 0 |
  | changed values | 52 |

  The 25 failures are the known ones on main.

  The largest changes are in V6's cold-start field differences against
  Fluent (`test_fluent_elbow.cpp:202`, gate 0.005): 0.0011679 → 0.0011884
  passes in both builds, and 0.0051758 → 0.0052429 fails in both. Every
  other change is under 0.5 % relative. Examples:
  - V6 velocity L2: 3.36608 % → 3.36552 %;
  - pressure-drop error: 1.9694 % → 1.9666 %;
  - nozzle velocity: 1.60340 % → 1.60343 %.
- `python3 tools/lint_portability.py` passes.
- **CUDA:** not compiled here; nvcc and CPU/CUDA parity run in CI.
  `clang++ -x cuda --cuda-device-only -nocudainc -nocudalib -fsyntax-only`
  parses `prelude.h` and `fv_coupled.kern.h`. The same command rejects an
  injected undefined symbol, so the check is real.

## Open items

1. **Ahmed L1 still diverges with Coupled.** It fails at iteration 9 (main:
   2). The trace shows the blow-up starting on the 35-degree slant, with a
   pressure spike there that alternates in sign and grows.
   - Nothing that changes the step helps (`evidence/startup.tar.gz`,
     `ahmed/`):
     - time scale 0.2 fails at 11;
     - relaxation 0.3/0.3 at 15;
     - hybrid start at 13;
     - CFL 20 at 11, and CFL 5 at 7;
     - Rhie-Chow block scaled by 4 at 18, and by 0.25 at 10.
   - That is the signature of the tetrahedral channel in
     [`../coupled/REPORT.md`](../coupled/REPORT.md#tetrahedral-channel-bench---case-channel).
     When momentum is held back, pressure is set by the Rhie-Chow term alone,
     whose reconstructed-gradient self-coupling has the wrong sign on
     tet-derived cells.
   - These hexes come from splitting tetrahedra and prisms (`subdivide.py`).
     That face reconstruction and pressure path belongs to another lane, so
     it was not changed here.
   - SIMPLEC converges on this mesh.
2. **No margin above Fluent's pseudo-transient relaxation.** At explicit
   relaxation 0.75/0.75 or 0.7/1.0, Coupled still diverges on the wing (by
   iteration 50; `evidence/reproduction.tar.gz`, `cand_urf07*`).
   - After the floor, the true residual change was 1.6-1.9x the Jacobian's
     prediction (u, v, continuity; `jacobian_check_floor.txt`). That matches
     the Picard linearization's missing ∂F/∂u term.
   - At h/c 0.224, CFL 20 made 0.7/1.0 converge, while 50 and 200 did not. A
     lower default Courant number would trade some of the 0.5/0.5 speed for
     that margin, and it would need V6 re-measured.
   - The Courant number is not yet a `SolverSettings` or case setting.
3. **Bench stopping rule for Coupled.** `bench_fsae_aero` stops on
   Fluent-scaled residuals, which ask Coupled for about half of SIMPLEC's raw
   continuity (see above). A raw-continuity criterion, as V6 uses, would
   compare the schemes on equal terms.
4. **CUDA build and parity** remain for CI.

## Reproduce

Run from the checkout root. You need:
- `KUTTA_FLUENT_EXAMPLES` set;
- the FSAE lane's meshes in `<mesh dir>` (`tools/fsae_aero/run.py` generates
  them);
- a Release build of main in `<main build>`.

```sh
validation/results/coupled_external/run_fsae.sh build/bench_fsae_aero <mesh dir> runs/fsae
python3 validation/results/coupled_external/summarize.py runs/fsae
build/bench_fsae_aero --mesh <mesh dir>/zerihan_h0.224_L2.msh --re 4.62e5 --downforce --moving-ground \
  --tu 0.001 --ratio 1 --backend metal --iterations 4000 --check 100 --settle 0 --cells \
  --out runs/fields/z0.224_simplec   # and --coupled --iterations 2000 --out runs/fields/z0.224_coupled
python3 validation/results/coupled_external/compare_fields.py runs/fields/z0.224_simplec runs/fields/z0.224_coupled
validation/results/coupled_external/parity.sh <main build>/kutta build/kutta runs/parity
validation/results/coupled_convergence/measure.sh <main build> build runs/measure 2
python3 validation/results/coupled_convergence/summarize.py runs/measure
```
