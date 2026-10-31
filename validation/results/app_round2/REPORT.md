# Desktop app, round 2

Branch `t3/0dc6068c`, based on `868f478`. All numbers below are from this
checkout's build (`build-merge`, Release, Qt 6.12.0, M4 Pro, Metal) unless
stated otherwise. "main" means a build of `868f478` kept beside it.

## Result

- Last round's five open items are fixed. Two of them needed CLI changes,
  both v1-additive and covered by docs and CLI tests: a short label per
  approximation, and no z-velocity residual for 2D cases. Residuals now start
  at iteration 1: the CLI emits iterations 1-5 as residual-only events.
- The app now exposes the solver features from the brief: the intermittency
  transition model, the editable pressure scheme, Coupled availability with
  the Courant-number note, and the 2D ground-effect wing on the Aerodynamics
  page. Each was checked with the real CLI on real cases.
- CI now replays recorded runs of a standard k-epsilon case, a laminar 2D
  case and a transition case. The screenshot smoke test renders all three.
- `docs/app-screens` is regenerated (22 screens, light and dark) at 7.0 MB,
  down from 20.1 MB.
- The constraint holds: for the mixing elbow (`--approximate wfgc`, 40
  iterations), `kutta solve` writes fluxes.csv, forces.csv, loads.json and
  result.vtu byte-identical to main, on CPU and Metal, with and without
  `--progress jsonl`. No file under `src/kutta/fv` or `src/kutta/kernels` changed.

## What changed

### Open items from round 1

| Item | Before | Now |
| --- | --- | --- |
| Run History said "Completed" for runs that hit the iteration limit | `Completed` | `Iteration limit` (also in the status bar: "reached the iteration limit"); converged, stopped and failed runs keep their words |
| APPROXIMATE badge did not say which approximation | `APPROXIMATE` | `APPROXIMATE · WFGC off`, `APPROXIMATE · No compressibility correction, Documented k-epsilon only` in the viewport, charts, Forces and Aerodynamics. The CLI sends a `short` label per approximation in issues and `summary.json`, and runs keep the labels |
| 2D cases reported a z-velocity residual | CLI reported it, app hid it | The CLI omits z-velocity for `general/dimension` 2: iteration events, residuals.csv, the convergence inventory (`monitors/residual/*`) and `solved_residuals` |
| Residual charts started at iteration 20 | first point at the first monitor sample | The CLI emits iterations 1-5 as residual-only events. Details below |
| The nozzle's first contour view showed one face, obliquely | camera carried over from the mesh view | A results view whose items all lie in one plane opens face on to that plane. `Camera::fit` now uses the axes the view actually renders with, which fixes the zoom in Top views |

Iterations 1-5. The solver always monitors these (Fluent normalises
continuity by the largest of the first five), so the residuals already
existed and only the reporting was missing. Early events carry `residuals`
and an empty `reports`, and no `convergence`. They are not convergence
samples, and forces are not evaluated for them, because `solver.forces()`
writes scratch state. The cost was measured on the elbow, 40 iterations,
wall time main → branch: CPU 4.01 → 4.04 s, Metal 0.58 → 0.55 s, within
run-to-run noise. Outputs are byte-identical, as listed above. The protocol
doc states all of this (`docs/app-protocol.md`, "Iteration events"), and
`docs/app-protocol-requests.md` records the items as done.

### Solver features in the app

- **Transition model.** Models shows the transition model under the viscous
  model. The Intermittency Transition Model gets a description, its constants
  read-only (CTU1-3, CPG1-3 and limits, CSEP, σ), and pointers to where its
  scheme, relaxation and inflow intermittency live. Transition SST and the
  algebraic model show the CLI's rejection reason as an error note. The
  `intermit` residual, the Intermittency field and its scheme, relaxation and
  inlet value appear with Fluent's labels.
- **Pressure scheme.** `methods/pressure-scheme` (standard or second order)
  is editable through the CLI's inventory entry from main.
- **Coupled and Courant number.** Methods offers Coupled only with the
  pseudo-transient formulation. Controls shows the case's flow Courant number
  read-only (`controls/courant-number`, a new v1-additive inventory entry,
  editable false). Its note says Kutta's Coupled solver always uses Fluent's
  default of 200 as implicit relaxation on top of the pseudo time step.
- **FSAE ground-effect wing.** A new generator,
  `tools/fsae_aero/naca_ground_case.py`, builds the case: an inverted NACA
  4412 at 3° nose down, h/c 0.2 above a moving ground, chord 0.25 m, 15 m/s,
  SST with the intermittency model. Settings come from the public mixing
  elbow example, edited as a user would in Fluent (drag and downforce report
  definitions, reference values, a moving-wall ground). The mesh has 53,881
  cells. Results with the real CLI on Metal:
  - The run stopped at 500 iterations (22.6 s) when downforce and drag
    changed less than 0.1 % over 5 values.
  - Downforce 41.8 N and drag 0.8 N per metre of depth. CL −1.213, CD 0.023,
    L/D 52.4.
  - Continuity levels off near 2e-3, which is why the scenario stops on the
    reports rather than on residuals.
  - These numbers are not validated against an experiment.

  The Aerodynamics page treats a 2D case as a section:
  - Up is +y. Axes come from the case's force reports. A downforce report is
    often a lift with a downward vector, so a report fixes the up axis but
    not its sign.
  - The front axle, wheelbase and balance are hidden, and the page says that
    forces are for the reference depth.
- **2D framing and mesh.** A 2D external-flow case is framed on its body
  with a body length of flow around it, instead of a far field 30 chords out.
  The 2D plane's mesh draws on the neutral surface colour, not a zone colour.

### Defects found and fixed along the way

- The legend kept the previous field's range when the shown surfaces lacked
  the field (velocity on a wall showed the pressure range). The range now
  resets, the legend hides, and the subtitle says the field is not defined on
  the surfaces shown.
- The k-epsilon elbow's first results view was walls only, because the
  recording had no slice. It now records the app's default slice.
- The "Suggest from zone names" button had no icon (`wand-sparkles` was
  missing from the resources).
- `docs/app-screens/README.md` said the seeded elbow converges at 120. Main
  and this branch both converge at 100.

### CI coverage

- `app/scripts/capture_fixtures.py` replaces `capture_elbow_fixtures.py`. It
  records four cases from the real CLI and scrubs local paths. 2D views keep
  only the extrusion cap the app shows (that saved 2 MB).

| Fixture | Case | Model | Size |
| --- | --- | --- | --- |
| `mixing_elbow` | Ansys mixing elbow (MIT) | SST, Coupled, 3D | 3.86 MB (re-recorded with this CLI) |
| `elbow_param` | Ansys elbow_param | standard k-epsilon, SIMPLE | 0.67 MB |
| `cylinder_tutorial` | Ansys adjoint cylinder tutorial | laminar, 2D | 0.82 MB |
| `naca4412_ground` | generated wing, above | SST + intermittency, 2D | 1.64 MB |

  The fixture total is 6.99 MB, up from 3.86 MB.
- The fixture CLI serves any recorded case by file name. It replays early
  iterations, evaluates convergence only on monitor samples, and accepts the
  four approximations.
- New samples and scenarios: `param-run`, `cylinder-run`, `wing-models`,
  `wing-aero`, `wing-monitors`. The screenshot smoke test's default set adds
  `param-run`, `cylinder-run` and `wing-aero`.
- `tst_session::recordedCasesKeepTheirOwnEquations` runs each new recording
  through the process boundary. It checks residual sets, approximation
  labels, the planar flag, transition, iterations 1-5 then 20, the "Iteration
  limit" status, and, for the wing, up = +y and positive downforce from
  `kutta aero`.

### Screens

`app/scripts/capture_screens.sh` regenerated all 22 documented screens in
both themes, with the real CLI for the elbow, cylinder and wing (`car-*` stay
synthetic). They are downscaled to 1500 px wide and palette-reduced with
libimagequant: 20.1 → 7.0 MB. Lossless optimisation alone gave 21.5 MB.
Pillow's median-cut palette shifted contour colours visibly, so it was not
used. libimagequant output was compared against the originals on the
diverging Kutta − Fluent map and on text-heavy pages, with no visible
difference. New screens: `wing-models`, `wing-aero`, `cylinder-run`.

## Before and after

Before shots use main's app and CLI; after shots use this branch. All use
the real CLI, at a 1500×940 window, downscaled to 1200 px
(`screens/`).

| Surface | Before | After |
| --- | --- | --- |
| Run History, laminar 2D cylinder (status, residuals) | [light](screens/before-history-cylinder-light.png) · [dark](screens/before-history-cylinder-dark.png) | [light](screens/after-history-cylinder-light.png) · [dark](screens/after-history-cylinder-dark.png) |
| Monitors, elbow (badge, residuals from iteration 1) | [light](screens/before-monitors-elbow-light.png) · [dark](screens/before-monitors-elbow-dark.png) | [light](screens/after-monitors-elbow-light.png) · [dark](screens/after-monitors-elbow-dark.png) |
| Nozzle, first results view | [light](screens/before-nozzle-run-light.png) · [dark](screens/before-nozzle-run-dark.png) | [light](screens/after-nozzle-run-light.png) · [dark](screens/after-nozzle-run-dark.png) |
| Methods, elbow (Coupled note) | [light](screens/before-methods-elbow-light.png) · [dark](screens/before-methods-elbow-dark.png) | [light](screens/after-methods-elbow-light.png) · [dark](screens/after-methods-elbow-dark.png) |
| Models, elbow (transition row) | [light](screens/before-models-elbow-light.png) · [dark](screens/before-models-elbow-dark.png) | [light](screens/after-models-elbow-light.png) · [dark](screens/after-models-elbow-dark.png) |
| Controls, elbow (Courant number) | [light](screens/before-controls-elbow-light.png) · [dark](screens/before-controls-elbow-dark.png) | [light](screens/after-controls-elbow-light.png) · [dark](screens/after-controls-elbow-dark.png) |
| Models, 2D wing (intermittency) | new | [light](screens/after-models-wing-light.png) · [dark](screens/after-models-wing-dark.png) |
| Aerodynamics, 2D wing | new | [light](screens/after-aero-wing-light.png) · [dark](screens/after-aero-wing-dark.png) |
| k-epsilon elbow run (epsilon residual, two approximations) | new | [light](screens/after-residuals-param-light.png) · [dark](screens/after-residuals-param-dark.png) |

## Tests

| Check | Result |
| --- | --- |
| `ctest --test-dir build-merge/app -L app` with `KUTTA_REAL_CLI`, `KUTTA_ELBOW_CASE`, `KUTTA_FLUENT_EXAMPLES` | 10/10 passed |
| `ctest -R realcli -V` | ran, not skipped: 8 passed (inventory, live run + post-processing, convergence + aero, planar laminar cylinder, k-epsilon elbow, stop) |
| `./build-merge/kutta-tests --test-suite-exclude=slow` | 256 test cases, 7,453,841 assertions, all passed (includes `tests/app/check_solve.py`, `check_protocol.py`, `check_aero.py` against the CLI) |
| `python3 tools/lint_portability.py` | clean |
| Elbow byte identity vs main (CPU, Metal; plain and jsonl; 40 iterations) | fluxes.csv, forces.csv, loads.json, result.vtu identical |
| Screenshot smoke (6 scenarios × 2 themes) | passed, 22 s |

Idle repainting: nothing added animates. Progress still updates only on new
events; the new early events are five extra events per run.

## Open UX issues

- **Aero groups.** A zone named `wing` is suggested into "Body". The
  suggestions only know `front-wing`/`rear-wing` names. A lone wing should
  probably be its own group.
- **Fixture aero iteration count.** With the fixture CLI, Aerodynamics says
  "600 iterations" for a replay that stopped at 500, because the loads are the
  recording's final ones. The provenance says so; the real CLI is consistent.
- **Chart titles.** When approximation badges are long, the residual chart
  title is elided ("Scaled Resid…").
- **Full-bleed 2D contours.** The title and legend sit on the field and rely
  on a text halo. A translucent backdrop would read better.
- **Read-only items.** The Courant number and the transition model and its
  constants are read-only: Kutta uses 200 and the case's constants. Making
  them editable needs solver work, which is outside this lane.
- **Cylinder coefficients.** The Ansys cylinder tutorial's reference values
  make q·A = 1, so its "coefficients" equal its forces (Cd 1280). That is the
  case's own setting, shown faithfully. A hint when q·A = 1 might help.
- **2D framing.** It applies only to 2D external flow (walls that enclose
  the domain or a ground). The cylinder (inlet, outlet, symmetry sides) still
  opens on the whole domain.
- **No real-CLI wing test.** `tst_realcli` has no wing test, because the
  case is generated (gmsh, h5py) and not in `validation/fluent_examples`. CI
  covers it through the recording.
