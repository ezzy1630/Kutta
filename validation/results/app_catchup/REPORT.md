# Desktop app catch-up: coupling, convergence, new case types

The Qt app now covers the solver and CLI features added since its last pass:

- the three pressure-velocity coupling schemes
- why a run stopped, and which criterion it met
- k-epsilon, laminar and 2D cases
- Green-Gauss gradients
- opt-in for named approximations, with their measured effect

Every screenshot below is from the real CLI (`build-merge/kutta`) on corpus
cases. The app opened, edited and solved each one through its scripted
scenarios (`--scenario case-*`, documented in `docs/app.md`). Each surface is
shown in the light and dark themes at 1500 px.

Cases used:

| Case | File | Model |
| --- | --- | --- |
| mixing_elbow | `mixing_elbow/mixing_elbow.cas.h5` (+ `.dat.h5` where seeded) | SST k-omega, Coupled + pseudo-transient, 3D |
| elbow_param | `mixing_elbow/elbow_param.cas.h5` | standard k-epsilon, 3D |
| nozzle | `flprj_3891/.../nozzle_para_named_1.cas.h5` | k-omega, 3D |
| cylinder_tutorial | `adjoint_cylinder_tutorial/cylinder_tutorial.cas` | laminar, 2D, Green-Gauss cell based |

The cylinder runs without its `.dat`, which is corrupted upstream (an earlier
text-mode transfer). Kutta initializes it.

## What changed

### CLI (app-facing protocol fixes, v1-additive)

Commit `444f1ed`. Solver numerics are untouched: nothing under `src/kutta/fv`
or `src/kutta/kernels` changed. Every change only adds fields or entries, or
fixes an entry that was wrong. No existing field changed meaning.

- **Methods inventory.** New read-only entries describe the rest of Fluent's
  Methods panel: `methods/gradient`, `methods/pressure-scheme` and
  `methods/pseudo-transient`. `methods/epsilon-scheme` is the new k-epsilon
  dissipation scheme. A scheme entry is active only when the case's viscous
  model solves that equation.
- **Coupling options follow pseudo time.**
  `methods/pressure-velocity-coupling` offers `coupled` only when the case
  uses pseudo time or was saved with Coupled. Before this change, elbow_param
  offered Coupled and `solve` then rejected it.
- **Laminar cases.** The saved turbulence inputs and the k/dissipation
  relaxation are inactive and read-only. Before this change, the laminar
  cylinder's inlet k was editable.
- **`measured_effect`.** Each approximable issue carries its measured effect
  as a separate field, so the app can show it under its own heading. The
  explanation text is unchanged.
- **Issue labels.** `flow/scheme` issues are labelled "Pressure-velocity
  coupling".

Updated alongside: `docs/app-protocol.md`, `docs/app-protocol-requests.md`
(Done items 7–10), `tests/setup/test_settings.cpp` (new Methods test across
elbow_param and the cylinder) and `tests/app/check_protocol.py`.

### App

- **Methods page** (Fluent's layout).
  - *Pressure-Velocity Coupling* shows the scheme, Pseudo-transient
    formulation, and what will actually run. For example: "Saved with Coupled
    but without the pseudo-transient formulation, which Kutta's Coupled solver
    needs: this case runs SIMPLEC."
  - *Spatial Discretization* shows gradient, pressure, momentum, k,
    omega/epsilon and gradient limiting. A scheme appears only when the
    case's model solves that equation.
  - The run's actual coupling (from the `setup` event's
    `pressure_velocity_coupling`) appears in the console and on the Run page,
    and is saved with the run in Run History.
- **Convergence status** (Monitors and Run Calculation pages).
  - The title states the outcome: "Run 1 converged at iteration 120",
    "reached the iteration limit (300) without converging", or "stopped by you
    at iteration 102".
  - A note gives the reason and names the last residual to reach its
    criterion.
  - A table lists each condition's value against its criterion as met or not
    met, live during the run and after it.
  - The residual chart draws the criterion as a dotted line, and the stop
    iteration as a labelled vertical line.
  - The scaled residual definition sits under the criteria.
  - The Run History details add "Stopped" and "Coupling".
- **2D cases.**
  - From `general/dimension`, the app shows only the front plane of the
    one-cell extrusion, labelled "2D domain".
  - The camera looks down z with x right and y up. Left-drag and the arrow
    keys pan rather than orbit.
  - The mesh subtitle reads "2D · N cells".
  - The z-velocity residual is hidden from the chart, console table and
    criteria, and a note explains why.
  - General shows "Dimension 2D (planar)" with how Kutta solves it.
- **Approvals.**
  - Each approximation card shows its explanation, a *Measured effect*
    section and the evidence link.
  - The checkbox reads "Run with this approximation". The CLI flag is in its
    tooltip.
  - Before this change, the long checkbox text clipped the card's label.
- **Fields.**
  - Epsilon and y* have Fluent's labels.
  - Epsilon appears in the residual order.
  - Epsilon relaxation appears under Controls.
- **Fixtures.**
  - `app/scripts/capture_elbow_fixtures.py` re-recorded the elbow fixtures
    against the current CLI. The stale "SIMPLEC replaces Coupled" text is gone
    and the recordings carry the coupling and convergence objects.
  - The script now also captures `loads.json`.
  - The fixture CLI reports the coupling like the real CLI. It replays the
    recorded loads with a `FIXTURE REPLAY` provenance prefix.
  - The synthetic car's Methods entries match the real CLI's.
- **Idle cost.** Nothing new animates. Charts, guides and status tables
  repaint only when a progress event, a resize or the cursor changes them.

## Screenshots

| Surface | Light | Dark |
| --- | --- | --- |
| Methods, mixing_elbow (Coupled + pseudo-transient, LSQ, k-omega) | ![](methods-elbow-light.png) | ![](methods-elbow-dark.png) |
| Methods, elbow_param (k-epsilon; options SIMPLE/SIMPLEC only) | ![](methods-elbow-param-light.png) | ![](methods-elbow-param-dark.png) |
| Methods, cylinder 2D (saved Coupled without pseudo time → runs SIMPLEC; Green-Gauss; laminar) | ![](methods-cylinder-2d-light.png) | ![](methods-cylinder-2d-dark.png) |
| Approvals, elbow_param (compressibility approved, ke-options pending) | ![](approvals-elbow-param-light.png) | ![](approvals-elbow-param-dark.png) |
| Approvals, cylinder 2D (reverse-flow; Dimension 2D under Solver) | ![](approvals-cylinder-2d-light.png) | ![](approvals-cylinder-2d-dark.png) |
| Mesh, cylinder 2D (plan view of the 2D domain) | ![](mesh-cylinder-2d-light.png) | ![](mesh-cylinder-2d-dark.png) |
| Converged, mixing_elbow seeded from Fluent data (120 iterations) | ![](converged-elbow-seeded-light.png) | ![](converged-elbow-seeded-dark.png) |
| Converged, elbow_param k-epsilon (120 iterations; last to converge: continuity) | ![](converged-elbow-param-kepsilon-light.png) | ![](converged-elbow-param-kepsilon-dark.png) |
| Converged, nozzle (60 iterations) | ![](converged-nozzle-light.png) | ![](converged-nozzle-dark.png) |
| Iteration limit, cylinder 2D (300; continuity and x-velocity not met) | ![](iteration-limit-cylinder-2d-light.png) | ![](iteration-limit-cylinder-2d-dark.png) |
| Run History, cylinder 2D (Stopped, Coupling SIMPLEC) | ![](history-cylinder-2d-light.png) | ![](history-cylinder-2d-dark.png) |
| Stopped by the user, mixing_elbow (Stop pressed at iteration 100; the run finishes the current iteration and ends at 102) | ![](stopped-elbow-light.png) | ![](stopped-elbow-dark.png) |

The seeded elbow scenario also runs against the fixture CLI. It converges at
the same iteration (120) with the same residuals as the real CLI.

## Tests

All runs used this checkout and its `build-merge` Release build on an Apple M4
Pro.

- **App unit tests and screenshot smoke** (`ctest --test-dir build-merge/app -L
  app`): 10/10 passed. The real-CLI test skips here because it needs the
  corpus. New tests cover:
  - `tst_protocol`: Methods entries, dimension, measured effect, setup and
    recorded coupling, residual visibility
  - `tst_runhistory`: stop reasons in words, convergence chart marks, coupling
    persisted with the run
  - `tst_session`: a converged run names the conditions it met
  - `tst_viewformat`: a planar case shows its plane
  - `tst_settings`: updated for the re-recorded inventory
- **App against the real CLI** (`tst_realcli` with `KUTTA_REAL_CLI`,
  `KUTTA_ELBOW_CASE` and `KUTTA_FLUENT_EXAMPLES`): 8/8 passed. Two of these
  are new:
  - `planarLaminarCylinder` opens and solves the 2D laminar cylinder.
  - `kEpsilonElbow` solves elbow_param with its approximations.
- **Solver quick tier** (`./build-merge/kutta-tests
  --test-suite-exclude=slow`): 250 test cases, 7,417,230 assertions, 0 failed.
- **CLI protocol boundary** (`tests/app/check_protocol.py`, run by the
  `Solve CLI runs generated CFF cases…` test case, 1,267 assertions) and
  `tests/setup/test_settings.cpp` (5 cases, 69,223 assertions): pass
  against the final CLI build.
- **`python3 tools/lint_portability.py`**: exit 0.

## Open UX issues

- **2D residuals come from the CLI.** The CLI still reports a z-velocity
  residual (~1e-11) and criterion for 2D cases. The app hides them. Dropping
  them in the CLI for planar cases would match Fluent's list exactly (filed
  under Observations in `docs/app-protocol-requests.md`).
- **Approximation badges are generic.** Results say APPROXIMATE but not which
  approximation. A short per-approximation summary from the CLI is still an
  open request.
- **Charts show one point per monitor interval.** Residuals and report
  conditions are sampled at the monitor interval (20 on these cases). The
  chart starts at iteration 20, and convergence is detected at that
  granularity.
- **"Completed" in Run History.** The Status column says "Completed" for a
  run that hit the iteration limit. The details and the Monitors page say it
  did not converge, but the column could say "Iteration limit".
- **Nozzle default view.** The first contour view shows one oblique face. A
  better default surface or camera for this case wasn't attempted.
- **Cylinder convergence.** Continuity stays near 3e-2 after 300 iterations
  against its 1e-5 criterion. That is solver behaviour (the solver lane's
  remit), not something the app changes. The app reports it as an iteration
  limit with the unmet conditions.
- **Docs screenshots.** `docs/app-screens/` was not regenerated, so the
  screenshots there show the previous Methods and Monitors pages.
- **CI coverage of the new case types.** CI's fixture CLI covers the elbow and
  the synthetic car only. k-epsilon, laminar and 2D behaviour is exercised by
  `tst_realcli` and the CLI settings tests, which need the local corpus.
