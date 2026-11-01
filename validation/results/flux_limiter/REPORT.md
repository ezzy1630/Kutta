# The limiter in the Rhie-Chow face velocity, and Fluent's differentiable limiter

**The default is unchanged. Fluent's saved fluxes say that its mass flux
carries the limited reconstruction, as Kutta's does. Routing an unlimited
reconstruction into the flux, or switching every case to a differentiable
limiter, fails Fluent-reference gates. Kutta now implements Fluent's
differentiable limiter as an opt-in. It removes the Coupled external-aero
plateau.**

- **What Fluent documents.** The momentum-weighted face velocity carries an
  undocumented "high-order correction" `J_ho`. Fluent does not say whether
  the correction is limited. The limiter is documented only for second-order
  face reconstruction. Fluent says the standard (non-differentiable)
  limiters "tend to stall" residual convergence, and offers a differentiable
  limiter for that ([section 1](#1-what-fluent-documents)).
- **What Fluent's data says.** On the mixing elbow's saved `SV_FLUX`, the
  limited minus unlimited reconstruction enters the flux at coefficient
  1.03. A limited face velocity fits the saved fluxes better than an
  unlimited one: 1.5 times better in the side-reconstruction form, 1.06
  times in Kutta's point form ([section 2](#2-fluents-saved-flux)).
- **Unlimited flux, measured anyway.** It converges the wing. Against the
  Fluent references it gives 6 regressions and 8 recoveries: the k-epsilon
  elbow pressure drop fails (1.01% to 1.94%, gate 1.5%) and the V6
  wall-inlet shear fails (0.12% to 0.78%, gate 0.5%). Rejected
  ([section 3](#3-arbitration-against-fluents-converged-fields)).
- **New opt-in: Fluent's differentiable limiter** (`methods/limiter`,
  Fluent `recon/lim/type 2`). It uses Venkatakrishnan's function with Wang's
  threshold, 0.05 times each velocity component's global range, on momentum
  only:
  - On the wing it brings Coupled from a 1.2e-3 limit cycle to raw
    continuity 7.7e-5, with Coupled and SIMPLEC lift within 1.1e-5
    relative.
  - On Ahmed 35 L1, Coupled meets the 1e-3 residual criterion (scaled
    continuity 6.0e-4).
    Coupled and SIMPLEC agree on Cd to 1.4e-4 and on lift to 0.08%.
  - All 7 FSAE 2D sweep cases converge under Coupled, including h/c 0.448,
    which main's Coupled never converges.
  - It changes the fixed point. Ahmed lift rises 11% against main's
    SIMPLEC, and the Fluent references move (section 3). So it stays opt-in.
- **Default (Standard) is bit-identical to main.** That holds on the wing
  (CPU and Metal, SIMPLEC and Coupled, histories `cmp`-identical). Across
  the slow tier, the only moved values come from one run that stops on a
  wall-clock budget (section 6). The quick tier, the app tests (including
  real-CLI) and portability lint all pass.

## Builds

- Main is `f5487d6`, built in a separate source tree whose tracked
  sources match `f5487d6` byte for byte.
- The candidate is this commit's tree. It differs from main only in the
  files this commit lists.
- Both builds are Release, the project's default (`-O3 -DNDEBUG`).
- Platform: macOS 26.6.2, Apple Metal, shared machine. Wall times say which
  runs overlapped.

## 1. What Fluent documents

Ansys Fluent Theory Guide 2024 R2:

- **23.4.2, eqs. 23-46 to 23-49.** The face mass flux is
  `J_f = Ĵ_f + d_f (p_c0 - p_c1)`. Its momentum part is
  `Ĵ_f = ρ_f (a_p,c0 v_n,c0 + a_p,c1 v_n,c1)/(a_p,c0 + a_p,c1) + J_ho,ap`.
  `J_ho,ap` is described only as "a high-order correction for non-uniform
  and/or non-orthogonal meshes". Its form, its gradient and whether it is
  limited are not stated.
- **23.3.4.** Gradient limiters prevent spurious oscillations in
  second-order face reconstruction. The default is Standard
  (Barth-Jespersen).
- **23.3.4.3.** Non-differentiable limiters "tend to stall the apparent
  residual's convergence". The Differentiable limiter is "a modified form"
  [692] (Z.J. Wang, IJNMF 33, 2000) of Venkatakrishnan's limiter [669]
  (AIAA-93-0880). Fluent does not give its formula or its threshold
  coefficient.

The Fluent 12 User's Guide (node798) numbers the TUI slope limiter as
0 standard, 1 multidimensional, 2 differentiable. Reading `recon/lim/type 2`
as differentiable is an inference from that numbering: no corpus case saves
2.

The documentation does not settle whether `J_ho` uses the limited
reconstruction, so the saved data decides it.

## 2. Fluent's saved flux

`tools/forensics/face_flux.py` rebuilds the mixing elbow's saved `SV_FLUX`
from the saved cell fields (`evidence/face_flux.txt`):

- The flux regression on all deep faces fits the coefficient of
  limited-minus-unlimited reconstruction at **+1.03** (section 2c). Fluent's
  flux carries the limiter's effect at full strength.
- Assembled formulas, median `|error|/|F|` on all interior faces
  (section 3):

  | face velocity | median | sum\|div\| (kg/s) |
  | --- | ---: | ---: |
  | side reconstruction, Barth-Jespersen | 6.20e-4 | 8.98e-4 |
  | side, differentiable K=0.01 | 6.54e-4 | 8.02e-4 |
  | side, differentiable K=0.05 | 6.98e-4 | 8.84e-4 |
  | side, unlimited | 9.48e-4 | 1.30e-3 |
  | point (Kutta's), Barth-Jespersen | 7.67e-3 | 2.01e-3 |
  | point, unlimited | 8.15e-3 | 1.83e-3 |
  | point, differentiable K=0.05 | 7.40e-3 | 1.87e-3 |

The limited forms fit better in both families. The data does not support
the premise that Fluent's flux avoids the limiter.

The side form fits about 12 times better than Kutta's point form. Kutta
already tried it and rejected it on 2026-10-29 because it loses V2 solution
order (`docs/research/fluent-flux-gradient-forensics.md`). I did not reopen
that here.

A probe of Kutta's own formula at Fluent's converged fields (CPU) reads
continuity about 800 times Fluent's, limited or unlimited:

| case | limited | unlimited |
| --- | ---: | ---: |
| elbow | 2.069e-3 | 2.058e-3 |
| nozzle | 3.28e-3 | 2.97e-3 |

The limiter is not what separates Kutta's flux from Fluent's.

## 3. Arbitration against Fluent's converged fields

The slow cases that compare against converged Fluent fields were run on CPU
and Metal:
- the k-epsilon elbow and SST nozzle;
- the k-epsilon elbow with second-order pressure;
- the laminar 2D cylinder;
- V6 mixing elbow;
- V4 cavity.

Each variant was a scratch build: an environment hook forced the variant
into every case. `compare_slow_assertions.py` compared each with main
(`evidence/fluent-refs-*.json.gz`). CPU and Metal agree to the digits
shown. The hooks are not in this commit.

Values are CPU. **Bold** marks a gate that fails where main passes.

| metric (gate) | main | unlimited flux | differentiable K=0.05 | differentiable K=0.01 | BJ with ε floor, K=0.01 |
| --- | --- | --- | --- | --- | --- |
| V6 velocity L2 (< 0.5%, known failure) | 3.37% | 3.64% | 2.28% | 2.33% | 4.30% |
| V6 largest field error, k (< 0.1%, known failure) | 14.7% | 14.7% | 9.88% | 8.75% | 17.8% |
| V6 pressure drop (< 0.5%, known failure) | 1.97% | 2.14% | 1.70% | 1.57% | 2.08% |
| V6 wall-elbow shear (< 0.5%, known failure) | 1.63% | 1.73% | 0.89% | 1.08% | 2.46% |
| V6 wall-inlet shear (< 0.5%) | 0.12% | **0.78%** | **1.61%** | **1.65%** | 0.36% |
| V6 wall-elbow total (< 0.5%, known failure) | 0.61% | 0.42% pass | 0.44% pass | 0.46% pass | 0.44% pass |
| V6 wall-inlet total (< 0.5%, known failure) | 0.52% | 0.35% pass | 0.21% pass | 0.25% pass | 0.86% |
| k-epsilon elbow velocity / pressure L2 | 0.88 / 5.82% | 0.93 / 6.18% | 0.96 / 6.32% | 0.88 / 5.82% | 1.00 / 6.70% |
| k-epsilon elbow pressure drop (< 1.5%) | 1.01% | **1.94%** | 0.43% | 0.86% | 0.39% |
| k-epsilon, second-order pressure, pressure drop (< 1.7%) | 1.11% | 0.12% | **1.85%** | 1.29% | **1.92%** |
| SST nozzle pressure L2 / pressure drop / wall force | 0.527 / 0.462 / 0.566% | 0.100 / 0.003 / 0.062% | 0.114 / 0.020 / 0.083% | 0.434 / 0.364 / 0.475% | 0.197 / 0.147 / 0.207% |
| cylinder pressure drop (< 0.76%) | 0.51% | 0.56% | 0.56% | 0.37% | **1.12%** |
| cylinder Cd error (< 0.105%) | 0.070% | 0.076% | 0.066% | 0.013% | 0.085% |
| assertions, CPU and Metal: regressions / recoveries | | 6 / 8 | 6 / 8 | 4 / 8 | 4 / 4 |

V6 runs two starts. The table shows the first; the second agrees to 3e-5
absolute. The recoveries are the V6 wall totals.

- Under the differentiable variants, V6 converges in fewer iterations. Its
  loop therefore runs 8 fewer per-iteration `REQUIRE`s (counted as
  "missing").
- The "ε floor" variant is Barth-Jespersen with Wang's threshold added to
  the allowed range, `ψ = (|Δ| + ε)/|d|`. It converges the wing (Coupled raw
  continuity 4.5e-5, lift 1.36621) but is worse on V6 and fails two gates.

**Decision.** No variant matches Fluent at least as well as main on every
gate, so the default stays Standard. The differentiable limiter is
Fluent's documented remedy for the stall and fits the saved flux about as well as
the others. It is added as an option. A case selects it with
`recon/lim/type 2`; the app selects it with `methods/limiter`.

**K = 0.05.** At 0.01, Ahmed Coupled did not converge (section 4). SU2's
default Venkatakrishnan-Wang coefficient is 0.05.

**Momentum only.** k and ω under the differentiable limiter stall. The
wing at K=0.05 shows it: whichever turbulence equation takes the
differentiable limiter stops at 2e-4 to 1e-3 (scaled residuals, final
iteration):

| limiter on k, ω | iterations | Coupled k / ω / continuity | SIMPLEC k / ω / continuity |
| --- | ---: | --- | --- |
| differentiable, differentiable | 4000 | **1.2e-3** / **2.3e-4** / 3.7e-4 | **1.2e-3** / **2.5e-4** / 2.2e-6 |
| Barth-Jespersen, differentiable | 3000 | 2.9e-6 / **2.3e-4** / 9.7e-5 | 2.2e-6 / **2.3e-4** / 2.9e-7 |
| differentiable, Barth-Jespersen | 3000 | **9.1e-4** / 1.6e-6 / 3.7e-4 | **8.8e-4** / 2.9e-6 / 1.9e-6 |
| Barth-Jespersen, Barth-Jespersen (shipped) | 4000 | 6.5e-6 / 8.3e-7 / 1.8e-4 | 6.0e-6 / 2.2e-7 / 6.9e-7 |

The setup note records that this is not Fluent's exact limiter: "Fluent
documents neither its formula nor its coefficient".

## 4. External aero

All runs use `bench_fsae_aero` on Metal. The sweep and settings are as in
the [coupled-plateau report](../coupled_plateau/REPORT.md).

**Zerihan wing, h/c 0.224, L1** (`--re 4.62e5 --downforce --moving-ground
--tu 0.001 --ratio 1 --settle -1`):

| run | iterations | scaled continuity | raw continuity | lift | drag |
| --- | ---: | ---: | ---: | ---: | ---: |
| main Coupled | 3000 | 2.7e-3, cycling | 1.15e-3 | 1.38335 | 0.032433 |
| main SIMPLEC | 4000 | 8.9e-7 | 1.07e-4 | 1.38446 | 0.032357 |
| unlimited flux, Coupled (scratch) | 3000 | 2.0e-4 | 8.7e-5 | 1.36633 | 0.032248 |
| differentiable K=0.05, Coupled | 4000 | 1.8e-4 | 7.7e-5 | 1.38760 | 0.032454 |
| differentiable K=0.05, SIMPLEC | 4000 | 6.9e-7 | 8.3e-5 | 1.38759 | 0.032453 |
| differentiable K=0.01, Coupled | 4000 | 1.1e-4 | 4.5e-5 | 1.40401 | 0.032775 |
| differentiable K=0.01, SIMPLEC | 4000 | 5.3e-7 | 6.4e-5 | 1.40403 | 0.032775 |

The differentiable limiter moves the wing's lift by +0.23% against main's
SIMPLEC at K=0.05, and by +1.4% at K=0.01.

**Ahmed 35, L1** (`--re 2.78e6 --length 1.044 --area 0.112 --walls body
--lift-z --tu 0.01 --window 1000 --settle -1`, 2000 iterations; the six
runs overlapped):

| run | continuity at 2000 | Cd | lift | peak memory |
| --- | --- | ---: | ---: | ---: |
| main Coupled | scaled 4.0e-3, plateau | 0.35392 | 0.21146 | 1.985 GB |
| main SIMPLEC | scaled 5.8e-6 | 0.35323 | 0.20931 | 1.703 GB |
| differentiable K=0.05, Coupled | scaled 6.0e-4 (criteria met), raw 2.1e-4 | 0.36315 | 0.23310 | 1.991 GB |
| differentiable K=0.05, SIMPLEC | scaled 1.1e-6 | 0.36310 | 0.23292 | 1.710 GB |
| differentiable K=0.01, Coupled | scaled 6.2e-3, not converged | 0.37761 | 0.26300 | 1.991 GB |
| differentiable K=0.01, SIMPLEC | scaled 1.6e-6 | 0.37763 | 0.25819 | 1.714 GB |

- With the differentiable limiter, Coupled and SIMPLEC agree on Cd to
  1.4e-4 and on lift to 0.08%. Main's pair disagrees by 0.2% and 1.0%.
- Against main's SIMPLEC, Cd moves +2.8% and lift +11%. The limiter is a
  discretization change, not only a convergence aid. Nothing here says
  which fixed point is closer to Fluent on this body.

**FSAE 2D sweep, L2** (`exp/sweep.sh` in `evidence/runs.tar.gz`; 7 cases in parallel; stops when
the criteria are met; 3000 iterations max). Main's Coupled column is copied
from the coupled-plateau report, whose code is main's head `f5487d6`:

| case | main Coupled: iterations, lift | main SIMPLEC: iterations, lift | differentiable Coupled: iterations, lift | differentiable SIMPLEC: iterations, lift |
| --- | --- | --- | --- | --- |
| wing, free air | 700, 0.86123 | 1800, 0.86123 | 700, 0.85906 | 1800, 0.85906 |
| wing h/c 0.448 | **not met at 3000**, 1.14036 | 2850, 1.14029 | 750, 1.13714 | 2850, 1.13712 |
| wing h/c 0.224 | 750, 1.40149 | 2400, 1.40143 | 700, 1.39821 | 2400, 1.39815 |
| wing h/c 0.134 | 800, 1.46714 | 2000, 1.46729 | 750, 1.46253 | 1850, 1.46263 |
| wing h/c 0.090 | 1800, 1.41056 | 2300, 1.41058 | 1050, 1.40514 | 2300, 1.40569 |
| S1223 α 4.76 | 850, 1.58469 | 2150, 1.58470 | 800, 1.58133 | 2150, 1.58132 |
| S1223 α 10 | 1200, 1.93387 | 2100, 1.93373 | 1300, 1.92854 | 2000, 1.92856 |

- With the differentiable limiter, Coupled meets the criteria on all 7
  cases.
- Its Coupled and SIMPLEC lifts agree to 3.9e-4 relative.
- Against main's SIMPLEC, lift moves -0.21% to -0.39% and drag by -0.26%
  to +0.11%.

## 5. Cost

**Wing h/c 0.224 L1, 2000 iterations, Metal, two runs each, nothing else on the GPU:**

| build, scheme | wall time | peak memory |
| --- | ---: | ---: |
| main, Coupled | 36.3-36.5 s | 0.420 GB |
| candidate Standard, Coupled | 36.3-36.7 s | 0.420 GB |
| candidate differentiable, Coupled | 58.0 s | 0.419 GB |
| main, SIMPLEC | 64.1-64.3 s | 0.393 GB |
| candidate Standard, SIMPLEC | 64.3-64.4 s | 0.393 GB |
| candidate differentiable, SIMPLEC | 63.2-63.4 s | 0.393 GB |

- Standard costs the same as main.
- Differentiable Coupled costs 60% more per iteration on this mesh.
  Standard Coupled takes 17-22 ms per iteration throughout. Differentiable
  Coupled takes 24 ms while converging and 29-30 ms once raw continuity is
  near 2e-4. The cost is the turbulence solves at a converged state, not
  the limiter. A scratch
  bench build logged each solve's BiCGStab iterations every 10th step
  (mean over iterations 1000-2000):

  | run | block (Coupled) or pressure (SIMPLEC) | k | ω |
  | --- | ---: | ---: | ---: |
  | Standard Coupled (cycling) | 6.0-6.5 | 7.2 | 2.9 |
  | differentiable Coupled (converged) | 4.7-4.9 | 12.3-12.4 | 16.8-17.1 |
  | Standard SIMPLEC (converged) | 11.0-11.9 | 24.6-25.6 | 19.6 |
  | differentiable SIMPLEC (converged) | 18.3-18.6 | 12.1-12.2 | 17.5-17.6 |

  - The block solve never reaches its 20-iteration cap, and needs fewer
    iterations under the differentiable limiter.
  - The k and ω solves, which run to a relative tolerance, cost what they
    cost in any converged run.
  - Standard Coupled is cheap only because it never converges.
  - On free-air L2, where Standard Coupled also converges, the
    differentiable limiter is not slower: 72-79 ms per iteration against
    70-98, and 115.6 s against 126.7 s for 1500 iterations.
  - The limiter's own work is one batched reduction folded into the
    time-step refresh, plus an exact refresh only when a component leaves
    its bound. Before that batching, the same run took 61.5 s instead of
    58.0 s.
- Time to the criteria is what users see, and there the differentiable
  limiter is far ahead. Over the sweep, differentiable Coupled took
  282-360 s; main's SIMPLEC took 629-786 s (both sets ran 7 at a time).
- Memory: the differentiable limiter allocates 6 floats (24 B) per cell of
  scratch, and only when selected. Ahmed peak RSS rose 0.3-0.4%.

## 6. Verification

- **Parity.** Standard is bit-identical to main on the wing L1 (CPU and
  Metal, SIMPLEC and Coupled; `history.csv` `cmp`-identical apart from
  wall-clock columns) after the final refactor.
- **Quick tier** (`kutta-tests --test-suite-exclude=slow`, fixtures present,
  `KUTTA_FLUENT_EXAMPLES` set, GPU idle, back to back):
  - candidate: 262 cases and 7,457,289 assertions, all passed, 88.5 s, peak
    RSS 2.60 GB;
  - main: 261 cases and 7,454,665 assertions, all passed, 86.2 s, 2.55 GB.

  The added case is `differentiable limiter bounds at zero threshold and
  passes small variations`.
- **New and changed tests:**
  - the differentiable-limiter test: field ranges, bound at zero threshold,
    fused equals staged, paired-with-zero-low equals unpaired, ψ > 0.97 at a
    large threshold;
  - the monitor-cadence test now covers SST with the differentiable limiter
    under SIMPLEC and Coupled;
  - the mixing-elbow settings test checks `methods/limiter` (editable,
    override, note).
- **App.** The inspect fixtures were re-recorded (`methods/limiter` added,
  `case/recon/lim/type` now consumed). `tst_protocol`'s editable count goes
  from 59 to 60, and it checks the new control. On the two SIMPLEC fixtures
  the `controls/courant-number` reason was stale since `f5487d6` and is now
  refreshed. `ctest --test-dir build/app -L app`: 10 of 10 passed, with
  real-CLI running (8 cases, none skipped).
- **Portability lint** (`tools/lint_portability.py`): clean.

### Slow tier

The tier was run in thirds (doctest `--first/--last` 1-11, 12-22, 23-33),
main and candidate concurrently, with `-r=xml --success`.
`tools/forensics/compare_slow_assertions.py` compared each third
(`slow-{main,cand}-*.xml.gz`, `slow-comparison-*.json.gz`; paths in the XML
are written `<checkout>` and `<main>`):
- 33 cases, **0 regressions, 0 recoveries, 0 missing assertions, 0 new
  assertions**.
- Main and candidate both pass 39,495 assertions and fail 65: the known 65
  (V6 25, V2 tet 8, V5 plate 6, V5 bump 9, V7 SA 16, V12 polynomial fan 1).
  V8 passed on both builds.
- 7 asserted values moved, all in `V5 TMR bump`, CPU level 2, none across
  its bound. That run stops on its 600 s wall-clock budget, checked every
  100 iterations: main stopped at 2100 iterations (601 s), the candidate at
  2200 (630 s). The three cf errors fell by 6.5-8.7%; two are known
  failures that stay failing (0.142 to 0.132, and 0.324 to 0.297). The
  three cp errors moved by at most 7.4e-4 absolute (bound 0.03), and the
  iteration count went from 2100 to 2200 (bound 3000). CPU levels 0 and 1
  (3000 iterations each) and all three Metal levels print identical
  residual changes and errors to every digit on both builds. This is a
  load-dependent stopping point, like V8's budget, not a solver change.

The first attempt at third 2 ran without the level-2 bump grid. Both
checkouts lacked the decompressed `bump_2levelsdown_353x161.p2dfmt`, so V5
bump and V7 SA raised "grid missing" on both builds. The missing grids are
gitignored. They were symlinked from the fixture master and the third was
rerun on both builds; the numbers above are from the rerun.

The candidate's slow tier ran on the binary from before the last edits.
Those edits touched comments only (`steady_solver.{hpp,cpp}`,
`fv_gradient.kern.h`, `fluent_case.hpp`). The quick tier and app tests ran
after all but the last of them, a comment in `steady_solver.hpp`. After that
edit the tree was rebuilt, and the limiter, monitor-cadence and settings
tests were rerun (5 cases, 9,088 assertions, passed).

| third | main wall / peak RSS | candidate wall / peak RSS |
| --- | --- | --- |
| 1-11 | 1282 s / 2.15 GB | 1279 s / 2.18 GB |
| 12-22 | 3008 s / 1.30 GB | 3050 s / 1.21 GB |
| 23-33 | 176 s / 0.159 GB | 176 s / 0.159 GB |

## Open items

1. **Kutta's Barth-Jespersen may over-limit relative to Fluent's Standard
   limiter.** Both differentiable variants improve V6 toward Fluent:
   - velocity 3.37% to 2.28%;
   - k 14.7% to 8.8-9.9%;
   - pressure drop 1.97% to 1.6-1.7%;
   - wall-elbow shear 1.63% to 0.9-1.1%;
   - both wall totals.

   The V6 wall-inlet shear (0.12% to 1.6%) and, at K=0.05, one k-epsilon
   pressure drop regress. A follow-up lane should start from the V6
   inlet-wall cells. If Fluent's Standard limiter differs from Kutta's in
   where it acts (wall cells, stencil range, boundary values), that is
   where it should show.
2. **The differentiable limiter is momentum only.** On k or ω it stalls
   that equation's residual at 2e-4 to 1e-3. The inferred cause is the
   turbulence floors. Fluent applies its limiter to every equation and
   documents no floor interaction.
3. **`recon/lim/type 2` = differentiable is inferred** from Fluent 12's TUI
   numbering. No corpus case saves it, and Fluent's coefficient is unknown.
4. **The side-reconstruction face velocity** fits Fluent's saved flux 12
   times better than Kutta's point form. It was rejected for V2 order on
   2026-10-29, and that conflict is still open.

## Reproduce

```sh
# saved-flux forensics (writes evidence/face_flux.txt)
python3 -I tools/forensics/face_flux.py "$KUTTA_FLUENT_EXAMPLES/mixing_elbow" mixing_elbow <cache>
# wing, Ahmed, sweep (scripts in evidence/runs.tar.gz: exp/*.sh)
bench_fsae_aero --mesh zerihan_h0.224_L1.msh --re 4.62e5 --downforce --moving-ground --tu 0.001 --ratio 1 --settle -1 --coupled --iterations 4000 --limiter differentiable --out w
bench_fsae_aero --mesh ahmed35_L1.msh --re 2.78e6 --length 1.044 --area 0.112 --walls body --lift-z --tu 0.01 --window 1000 --settle -1 --coupled --iterations 2000 --limiter differentiable --out ah
# slow tier per build in thirds, compared assertion by assertion
kutta-tests --test-suite=slow --first=1 --last=11 -r=xml --success > slow-cand-1-11.xml
tools/forensics/compare_slow_assertions.py slow-main-1-11.xml.gz slow-cand-1-11.xml.gz > slow-comparison-1-11.json
```

`evidence/runs.tar.gz` holds every history and log named here:
- `exp/diff`: wing;
- `exp/ahmed`;
- `exp/sweep`;
- `exp/timing3`: final-code timing;
- `exp/iters`: the linear-iteration logs (`LIN` lines);
- `exp/parity`;
- `exp/wc_*`: scratch flux variants.

It holds no meshes, no wall-force CSVs and no fixture data.
