# Coupled convergence: Fluent's pseudo-time step and a stable block smoother

Coupled now reaches Fluent's convergence criteria in fewer iterations than
SIMPLEC on the V6 elbow (seeded and cold) and the SST nozzle, on CPU and Metal.
It takes less wall time everywhere except one tie: the nozzle on Metal under
Fluent's residual criteria. The fixed point is unchanged, and SIMPLE and
SIMPLEC are bitwise identical to main.

- Branch `t3/dca2ac99`, base `ccd9ed4` (main). Commits:
  - `2a7b5c7`: Fluent's pseudo-time step and a Schur-bounded smoother.
  - `57ca0c9`: parallel coarsest-level sweeps.
  - The report commit follows these.
- Machine: M4 Pro, 24 GB, CPU (serial reference backend) and Metal. Another lane
  shared the CPU and GPU throughout. Baseline and candidate runs alternate, and
  times are ranges over two repetitions. Iteration counts are deterministic.
- Evidence is in `evidence/`, produced by the scripts here, with home paths
  replaced by `<checkout>`.

## Diagnosis

### 1. The nozzle's "stall" is the continuity normalization, not convergence

Fluent scales continuity by the largest raw value of the first five
iterations. Coupled already satisfies continuity at iteration 1, so its scale
is about 127 times smaller than SIMPLEC's (`bench_convergence --case nozzle`,
main, CPU; `evidence/diagnosis.tar.gz: nozzle_normalization/`):

| nozzle, main | first-five scale | raw continuity at 100 its | iteration reaching the V6 criterion |
|---|---|---|---|
| Coupled | 5.07e-3 | 1.39e-6 | 100 (8.2 s) |
| SIMPLEC | 0.643 | 1.01e-4 | 150 (17.0 s) |

Coupled's raw continuity floor is lower than SIMPLEC's: about 1.3e-6 against
1.7e-6, both the float32 floor of the face-flux sums. Its Fluent-scaled
continuity still levels off at 2.5e-4 because of the smaller denominator. The
same holds in Fluent: its saved Coupled histories use scales of 8.83e-3 (nozzle)
and 2.61e-3 (elbow), comparable to Kutta Coupled's and not to SIMPLEC's. So
Kutta's scaled residual behaves like Fluent's own. Accuracy against Fluent was
already identical for both schemes.

### 2. The elbow is limited by the pseudo-time step, and Kutta's was too small

On main, the elbow's outer convergence responds to the pseudo-time step and
the explicit relaxation, not to the linear solve. Cold elbow, CPU, iterations
to the V6 criterion (`diagnosis/main_sweep/`):

| main setting | iterations |
|---|---|
| as saved (time scale 1, relaxation 0.5/0.5) | 325 |
| block-solve tolerance 0.1 → 0.01 | 325 |
| relaxation 0.75 / 1.0 | 225 / 175 |
| time scale 2 / 5 / 20 | 225 / 175 / 150 |
| SIMPLEC | 525 |

The time scale and relaxation are the case's saved settings, so Kutta must
honour them. The automatic step itself was wrong, though. Fluent's global time
step (Theory Guide 23.6.2, eqs. 23-102 to 23-105) is

`Δτ = 0.3 L / max(0.1 U_bc, U_domain, U_Δp)`

where U_domain is the arithmetic **mean** cell speed. Kutta used the
**maximum** cell speed. With the mean, Kutta's step at the elbow's saved state
is 0.0426537 s. The case stores Fluent's last automatic step,
`pseudo-auto-time-step 0.0426536236`, so the two agree to 2e-6 relative. Before
the fix the step was 0.0245 s (V6 seeded `V6_INITIAL pseudo_time_step`; also
`src/kutta/fv/VERIFICATION.md`). The documentation page also says that the
global time step method applies only to the Coupled (and density-based)
solvers. So SIMPLE and SIMPLEC keep the maximum speed, and their paths do not
change.

### 3. The correct step exposed a divergent smoother

With Fluent's step, the TMR flat plate (stretched boundary-layer grid,
`--level 2`) slowed from 400 to 700 iterations (SST). BiCGStab hit its
20-iteration cap on every outer iteration and reduced the residual only 0.5–1×
against the 0.1 target. Main shows the same at start-up (20 iterations, 0.29).
Sensitivity runs on the plate (`diagnosis/smoother_main/`):

- damping 0.5: healthy (2–10 iterations);
- damping 0.9, 2+2 sweeps or 3+3 sweeps: diverge.

A smoother that gets worse with more sweeps has the wrong diagonal. In each
cell's arrow block, the pressure row used only the compact Rhie-Chow diagonal
c. The u–p coupling to neighbours is skew (G_ji = −H_ij on the face parts), and
for that structure damped Jacobi converges only while
ρ(c⁻¹ H D⁻¹ G) < (2 − ω)/ω. On stretched cells at the larger step, that radius
is about 2, and ω = 0.7 sits on the limit.

### 4. On Metal, a single-thread coarsest solve dominated the block solve

Aggregation stalls near 100 rows: 134 on the elbow and 79 on the nozzle. That
is below the 192-row switch to `coupled_coarsest_solve`, a one-thread
Gauss-Seidel kernel that took about 1.4 ms per call on Metal, twice per
BiCGStab iteration. Timing around `CoupledSolver::solve` (temporary
instrumentation, `diagnosis/metal_profile/`) gave 4.95 s of the elbow's 7.06 s
in the block solve, about 4.8 ms per BiCGStab iteration. With parallel sweeps
it took 2.0 s (`diagnosis/metal_coarsest/`).

## The change

1. **Pseudo-time step** (`steady_solver.cpp`, `update_time`). For Coupled, the
   step uses the mean cell speed (one extra sum in the same host reduction as
   the maximum). The block scaling wc still uses the maximum. SIMPLE and
   SIMPLEC are untouched.
2. **Schur-bounded block smoother** (`fv_coupled.kern.h`). Every smoother adds
   a Gershgorin bound on the neighbour part of the Schur complement to each
   cell's pressure diagonal:
   `cs_i = Σ_{k≠i} Σ_c |H_ik,c| R_k,c / D_k,c`, where `R_k,c = Σ_j |G_kj,c|` is
   the full row of the pressure-force block, including the own-cell and
   least-squares Q terms. It applies on the fine level and on every Galerkin
   level. Only the smoothers use cs; J x, the residual and the fixed point are
   unchanged.
3. **Coarsest level**: 8 parallel block-Jacobi sweeps replace the single-thread
   Gauss-Seidel kernel, which is removed. This needs a stable smoother, which
   item 2 provides.
4. `bench_convergence --case nozzle`. Elbow and nozzle read their continuity
   target from the saved Fluent data.

Candidates the evidence ruled out:

- **Tighter block-solve tolerance**: no gain (325 → 325).
- **ILU or multicolour Gauss-Seidel**: in the 2×2 skew model, Gauss-Seidel
  converges only for σ² < 1, which is no better than damped Jacobi. The
  smoother needed a correct pressure diagonal, not a different sweep order.
- **A more consistent Jacobian**: after items 1–3 the outer contraction
  matches what the time step and relaxation allow (time scale 5 still gives
  150 iterations), so the Jacobian is not the limit.

Effect of the bound on the plate (`diagnosis/smoother_bound/`, CPU, 400 its):

- ω 0.7 (default): BiCGStab needs 2–9 iterations.
- ω 0.9: stable (was divergent).
- ω 1.0: still diverges. Jacobi on diffusion-dominated momentum rows has
  checkerboard eigenvalues near −1, so this is expected.
- 3+3 sweeps: no longer divergent, but up to 20 BiCGStab iterations, so some
  margin is left.

`diagnosis/bound_scale/` scales cs by 0, ¼, ½ and 1 with the Fluent step:

- **Outer iterations:** identical for every non-zero scale (elbow 240, nozzle 70).
- **Plate:** cs = 0 caps BiCGStab at 20 on every outer iteration; any non-zero
  scale brings it down to 3–9.
- **Elbow:** cs = 1 costs about one extra BiCGStab iteration per outer step
  (4–5 instead of 3).

The bound is kept at 1, its stable value, without tuning: FSAE prism layers
are more stretched than the TMR plate.

## Results

### Iterations and solve time to Fluent's convergence criteria

These are `kutta solve` runs with `--approximate wfgc` and the case's own
Fluent residual criteria (absolute 1e-3 on every Fluent-scaled residual). The
`--setup` override `monitors/check-convergence: true` switches checking on,
because both cases save it off. Convergence is checked every 20 iterations.
Times are the `done` event's elapsed solve time. "Before" is main `ccd9ed4`;
"after" is `57ca0c9`. Cold runs use the case's hybrid initialization; seeded
runs use `--data` with Fluent's saved solution.

| case | backend | Coupled before | Coupled after | SIMPLEC |
|---|---|---|---|---|
| elbow cold | CPU | 200 / 12.9–13.7 s | **160 / 12.1–12.6 s** | 260 / 22.6 s |
| elbow cold | Metal | 200 / 3.9–4.5 s | **160 / 2.1–2.4 s** | 260 / 3.5–4.1 s |
| elbow seeded | CPU | 140 / 8.9–9.1 s | **100 / 7.2 s** | 120 / 10.2–12.7 s |
| elbow seeded | Metal | 140 / 2.6–2.8 s | **100 / 1.33 s** | 120 / 1.60–1.63 s |
| nozzle cold | CPU | 80 / 6.0–6.1 s | **60 / 6.1–6.2 s** | 80 / 8.2–8.3 s |
| nozzle cold | Metal | 60 / 1.1–1.2 s | 60 / 1.1–1.3 s | 80 / 1.1–1.2 s |

On main, Coupled took more iterations than SIMPLEC from Fluent's solution
(140 vs 120) and was slower than SIMPLEC on Metal for the cold elbow. For
context, Fluent's own saved elbow run stopped at 294 iterations when its
continuity reached 9.5e-4.

The nozzle's Fluent-criteria Metal run is a tie at about 1.1 s. The run is
short, and on this 19 000-cell mesh Metal is bound by kernel launches. Coupled
costs 18–22 ms per outer iteration there against SIMPLEC's 14–15 ms. Its 60
iterations are mostly start-up, where BiCGStab needs 9–11 iterations. On the
elbow, Coupled's Metal cost per outer iteration (17–19 ms) is now below
SIMPLEC's (19–22 ms).

### Iterations and time to the V6 criterion

The V6 criterion is: every steady transport defect below 1e-6, and raw
continuity below the saved Fluent run's final value (2.476e-6 for the elbow,
7.45e-4 for the nozzle). Cold runs use `bench_convergence` (sampled every 10
iterations) and also appear in the V6 slow test (every 50).

| case | backend | Coupled before | Coupled after | SIMPLEC |
|---|---|---|---|---|
| elbow cold (bench) | CPU | 320 / 24.3–26.1 s | **240 / 20.5–20.6 s** | 510 / 53.5–53.7 s |
| elbow cold (bench) | Metal | 320 / 7.2–8.5 s | **240 / 4.0–4.6 s** | 510 / 9.5–11.1 s |
| nozzle cold (bench) | CPU | 100 / 8.3–8.6 s | **70 / 7.4–7.5 s** | 150 / 17.5 s |
| nozzle cold (bench) | Metal | 100 / 2.0–2.2 s | **70 / 1.5–1.6 s** | 150 / 2.7–3.0 s |
| V6 seeded (slow test) | CPU | 300 / 22.0 s | **200 / 16.1 s** | 350 / 38.6 s † |
| V6 seeded (slow test) | Metal | 300 / 8.2 s | **200 / 3.4 s** | 350 / 8.4 s † |
| V6 cold (slow test) | CPU | 350 / 24.9 s | **250 / 19.7 s** | 550 / 55.5 s † |
| V6 cold (slow test) | Metal | 350 / 9.4 s | **250 / 4.1 s** | 550 / 12.6 s † |

† SIMPLEC under the V6 slow-test protocol is quoted from
`validation/results/coupled/REPORT.md` (commit `63287b9`, contended). SIMPLEC
is bitwise unchanged since then, so its iteration counts still hold. Its times
are from that earlier run. The V6 test always runs the case's scheme.

### Other Coupled cases (Metal; `evidence/tmr.tar.gz`)

These use `bench_convergence --scheme coupled`, level 2. The criterion is
summarize.py's "transport+load" from the previous Coupled report: transport
defect ≤ 1e-5 and load change ≤ 0.1 % over 100 iterations. "Mean BiCGStab" is
the per-solve iteration count at monitor samples.

| case | before | after | mean BiCGStab before → after |
|---|---|---|---|
| plate SST | 400 / 12.2 s | **400 / 4.6 s** | 7.9 → 2.5 |
| plate SA | 300 / 7.2 s | **300 / 1.9 s** | 8.9 → 2.5 |
| bump SST, SA | unconverged at 1500 | unconverged at 1500 | 2.9 → 3.0, 2.0 → 1.9 |
| tet channel, 294 912 cells | unconverged at 400 | unconverged at 400 | 1.3 → 1.2 |

- **Plate:** 3000-iteration drag agrees with main to 3e-6 (SST). SA differs by
  1.8e-4 (0.00288033 vs 0.00288084). Both are stationary, but SA's settled
  drag depends weakly on the pseudo-time step (below). SIMPLEC's SA drag is
  itself still moving at 3000 iterations (0.00288186 at 2000, 0.00288047 at
  3000). Final Coupled is 4.9e-5 from that value; main's Coupled is 1.3e-4.
- **Bump and tet channel:** the bench gives both a fixed step (bump 0.01,
  channel 0.01), so only the smoother differs. Neither converges, before or
  after. At 1500 iterations the bump's drag differs between builds by 4e-4
  (SST) and 1.4e-4 (SA) on unsettled transients. The channel's raw continuity
  oscillates between 0.016 and 0.044 in both builds, as in the earlier Coupled
  report.

Robustness at larger steps, SA plate, Metal, user time scale
(`diagnosis/sa_plate_time_scale/`):

| time scale | main | final |
|---|---|---|
| 0.5 | not run | settles, drag 0.00288039 |
| 1 | settles, drag 0.00288084 | settles, drag 0.00288033 |
| 2 | does not settle (drag 0.00286, cf0 0.00209 at 3000) | not run |
| 4 | **diverges** at iteration 74 (non-finite) | settles by 300, drag 0.00288028, 4 BiCGStab its per solve |

## Fixed point and accuracy

The gates are unchanged and pass as they did on main:

- `tests/fv/test_coupled.cpp` passes (Coupled reaches SIMPLEC's fixed point;
  tet channel). The closed cavity now converges in 75 Coupled iterations
  instead of 175, and the open channel in 75 (was 75).
- **Slow tests**, V6 and "Fluent standard k-epsilon elbow and SST nozzle
  match converged Fluent fields", baseline `ccd9ed4` against final `57ca0c9`,
  with `tools/forensics/compare_slow_assertions.py` (`evidence/slow/`). Run:
  `kutta-tests -tc="FV V6*,Fluent standard k-epsilon elbow and SST nozzle*" -r=xml --success`.

  | | count |
  |---|---|
  | regressions / recoveries | 0 / 0 |
  | pass-fail changes among the 51 changed values | 0 |
  | coverage case: assertions, failed (both) | 17 044, 0 |
  | V6 failed, both (pre-existing: strict-setup rejection, Fluent-accuracy thresholds, initial continuity) | 29 |
  | missing assertions | 16 |

  The 16 missing assertions are V6's per-50-iteration loop REQUIREs
  (`isfinite(maximum)`, `maximum_velocity<100`). They run 2 fewer times per
  mode and backend because Coupled converges 100 iterations sooner.

  The coverage test's INFO prints a `const char*` as a pointer and embeds the
  run's continuity value, so assertion identities differ between any two
  binaries. Both XMLs were compared after
  `sed -E 's/0x[0-9a-f]+/0xADDR/g; s/continuity=[0-9.eE+-]+/continuity=X/g'`.
  The stored XMLs are raw; `comparison.json.gz` is the normalized comparison.
- **Field errors against Fluent**, final vs main:
  - V6 velocity L2: 3.60102–3.60106 % vs 3.60099–3.60102 %.
  - V6 pressure L2: 2.96225–2.96232 % vs 2.96219–2.96224 %.
  - V6 pressure drop: 1.1840–1.1843 % vs 1.1837–1.1839 %.
  - Nozzle (300 iterations): velocity 1.6154 %, pressure 0.5833 %, k 5.054 %,
    ω 5.546 %, drop 0.517 %, wall force 0.622 %. Identical to main to 4
    significant figures.
  - elbow_param (k-epsilon, SIMPLE path): unchanged.

## SIMPLEC and SIMPLE are bitwise unchanged

`parity.sh` runs `kutta solve --iterations 40` with both binaries. The elbow
and nozzle use SIMPLEC and SIMPLE forced by `--setup`, on CPU and Metal. It
compares `fluxes.csv`, `forces.csv`, `loads.json` and `residuals.csv` without
`seconds_per_iteration` (`evidence/parity.tar.gz`):

```
elbow simplec cpu identical       elbow simple cpu identical
elbow simplec metal identical     elbow simple metal identical
nozzle simplec cpu identical      nozzle simple cpu identical
nozzle simplec metal identical    nozzle simple metal identical
```

## Device memory per cell (`SteadySolver::coupled_solver_bytes`)

| mesh | before | after |
|---|---|---|
| elbow (polyhedral) | 389.0 | 393.8 |
| nozzle | 385.8 | 390.5 |
| TMR plate/bump (quads) | 339.9 | 344.5 |
| tet channel, 294 912 cells | 348.2 | 352.8 |

The +4.6–4.8 B/cell is the Schur bound: one float per fine cell (a 20th slot
in the packed cell array) and one per coarse row (a 4th slot in each level's
diagonal array). Its setup reuses the smoother scratch buffers, which are free
before the first sweep. The removed coarsest kernel owned no memory.

## Quick tier and checks

`kutta-tests --test-suite-exclude=slow` with every fixture linked
(`evidence/quick_tier.txt`):

| build | cases | assertions | result |
|---|---|---|---|
| main | 239 | 7 414 934 | all pass |
| final | 239 | 7 414 926 | all pass |

The −8 assertions are the cavity's Coupled loop running 4 fewer times per
backend.

- `tools/lint_portability.py` passes.
- **CUDA:** not compiled here (no NVIDIA toolchain). `clang++ -x cuda
  --cuda-device-only -nocudainc` with attribute stubs parses `prelude.h`,
  `fv_flow.kern.h` and `fv_coupled.kern.h` without diagnostics. That checks the
  dialect only; nvcc and CPU/CUDA parity run in CI.

## Reproduce

From the checkout root, with `KUTTA_FLUENT_EXAMPLES` set. The baseline build
is main's tree with only `bench/bench_convergence.cpp` from this branch, so
that its bench knows `--case nozzle`.

```sh
validation/results/coupled_convergence/parity.sh <main build>/kutta build/kutta <out>
validation/results/coupled_convergence/measure.sh <main build> build <out> 2
python3 validation/results/coupled_convergence/summarize.py <out>
```

## Open

- **Metal cost per iteration on small meshes.** On the nozzle start-up,
  Coupled costs 18–22 ms per outer iteration on Metal against SIMPLEC's
  14–15 ms. It is bound by kernel launches, roughly 90 per BiCGStab iteration
  counted from the V-cycle structure. This makes the nozzle's Fluent-criteria
  Metal run a tie. Batching dispatches in the Metal backend, or fewer launches
  per V-cycle, is the next lever. It was not attempted here.
- **SA's settled state depends weakly on the pseudo-time step.** On the SA
  plate the settled drag moves about 1e-5 per doubling of the step (time scale
  0.5 / 1 / 4: 0.00288039 / 0.00288033 / 0.00288028). The momentum and
  continuity fixed point cannot depend on the step (delta form), so the
  dependence is in the segregated SA update. SST shows 3e-6 between the two
  steps. That belongs to the turbulence lane. It does not explain the whole
  1.8e-4 shift between main and this branch. Main's SA drag was still drifting
  (1.3e-4 between iterations 1000 and 2000), and SIMPLEC is still moving at
  that level at 3000. In float32, the SA plate's settled drag is defined only
  to about 1e-4.
- **Fluent-scaled continuity floor.** The float32 floor of raw continuity
  (about 1e-6 here) sits at 2.5e-4 of Coupled's first-five normalization on the
  nozzle. A user criterion for continuity below about 3e-4 cannot be met
  there, by either scheme in Coupled's normalization.
- **Remaining smoother margin.** ω = 1 and 3+3 sweeps still lose stability on
  the plate. The defaults (ω 0.7, 1+2 sweeps) have margin up to ω 0.9.
- **k/ω relaxation.** k and ω converge last on the elbow. Kutta reads Fluent's
  `k/pseudo-relax` 0.75; the case also stores `k/pseudo-implicit-relax` 0.8,
  whose role in Fluent's pseudo-transient turbulence update is not verified.
  This is outside this lane (turbulence is segregated).
- **SIMPLEC's pseudo-time step.** SIMPLEC still uses the maximum cell speed.
  Fluent documents the global step for Coupled only, and SIMPLEC had to stay
  bitwise identical, so this was not changed.
