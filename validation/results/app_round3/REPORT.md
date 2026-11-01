# Desktop app, round 3

Branch `t3/fb1927a8`, based on main `e308968`. All numbers are from this
checkout's Release build (Qt 6.12.0, M4 Pro, Metal) unless stated otherwise.
"main" means a build of `e308968` kept beside it.

## Result

- Round 2's five open UX items are fixed. Each is shown before and after in
  both themes below, with the real CLI (the fixture CLI for the replay item).
- The app now surfaces the two solver features merged since. A run's poor
  mesh numerics count (`poor_quality_cells`) shows in plain words beside the
  run's verdict, in Run History and in the console, only when it is not zero.
  The residual chart and convergence readout handle SIMPLE's residuals
  down to 1e-10.
- One CLI change, v1 additive and documented in `docs/app-protocol.md`:
  `summary.json` gains `poor_cell_quality` (the threshold the run used), so
  the app does not hard-code 0.05. The human output now prints that
  threshold as `0.05` rather than `0.050000000000000003`.
- The constraint holds: for the mixing elbow (`--approximate wfgc`, 40
  iterations), `kutta solve` writes fluxes.csv, forces.csv, loads.json and
  result.vtu byte-identical to main, on CPU and Metal, with and without
  `--progress jsonl` (16 of 16 files). No file under `src/kutta/fv` or
  `src/kutta/kernels` changed.

## What changed

### Open items from round 2

| Item | Before | Now |
| --- | --- | --- |
| A zone named `wing` was suggested into Body | `wing → Body` | `wing → Wing`. A new "Wing" group catches wings, flaps, elements, mainplanes, foils, slats, Gurneys and endplates whose names carry no front or rear. `fw-*`/`rear_wing_*` still go to Front wing / Rear wing |
| Fixture Aerodynamics said 600 iterations for a replay that stopped earlier | `Run 1 · 600 iterations · Converged` | The replay's iteration count, the same as its summary and done event |
| Residual chart title cut short by long badges | `Scal… [APPROXIMATE · No compressibility correction, …]` | The full title. Badges stay beside it when they fit, else wrap onto a row below it |
| 2D contour title and legend sat on the field with a text halo | halo only | Each label block (title, colour legend, zone key) sits on a plate of the viewport colour at 84 % opacity. Over empty background the plate matches the background and does not show |
| The cylinder opened on its whole domain | 30 m × 30 m domain, cylinder a dot | Framed on the cylinder with three diameters of wake |

Details:

- **Iteration count (fixture).** The fixture CLI copied the recording's
  `loads.json` verbatim, and `kutta aero` reads `iterations` from it. The real
  CLI writes the run's own count there. The replay now writes its own stop
  iteration too. The loads stay the recording's final values, which the
  provenance already says ("recorded wall loads of a 600-iteration Kutta
  run"). `tst_session` now checks that aero reports 300 for a 300-iteration
  replay of the 600-iteration recording.
- **Chart header.** Title and badge placement moved into `ChartLayout`
  (pure geometry, like the legend), so it is tested:
  `tst_plotting::chartHeaderKeepsTitleAndBadges` covers 5 widths × 3 heights
  × 3 badge sets. It checks that the title is never elided when it fits the
  chart alone, and that no badge overlaps the title, legend, plot or another
  badge.
- **2D framing.** Round 2 framed a 2D body only when the case also had a
  ground or tunnel wall. A 2D case is now external too when a body wall
  stays clear of the domain's edges in the plane (flow on every side).
  Framing gives one body length all round. The exception is when one side
  has more than 1.5 times the room of the side opposite, and both sides have
  more than a body length: that side is the wake and gets three lengths. The
  cylinder (10 D upstream, 20 D downstream) shows its wake. The wing (3.75 m
  ahead, 4.75 m behind, the ground just below) keeps round 2's frame. A
  first version took the space above the wing as its wake and put the wing
  at the bottom of the view; the both-sides rule fixes that, and
  `tst_viewformat::planarExternalFlowFramesTheBody` covers it, along with the
  cylinder and a channel, which still opens on the whole domain.

### Solver features

- **Poor mesh numerics.** `summary.json`'s `poor_quality_cells` and the new
  `poor_cell_quality` are parsed into the run summary. When the count is not
  zero, the app says "4 poor-quality cells (orthogonal quality ≤ 0.05) used
  low-order face pressures, as Fluent's poor mesh numerics does". This shows
  in the run's convergence box on Run Calculation and Monitors, as a "Mesh
  numerics" row in Run History, and in the console when the run finishes. A
  zero count, or an older CLI without the field, shows nothing. Of the app's
  real cases, only the 2D wing has such cells (4 of 53,881). The elbow,
  elbow_param and cylinder have none.
- **Residuals to 1e-10.** With the paired state, elbow_param's momentum
  residuals reach 1.6e-9 to 3.2e-9 in 300 iterations, and the cylinder's
  3.3e-10 to 5.2e-10 when it converges at iteration 1700. Two things did not handle that cleanly:
  - The log axis labelled every few decades counting from 1e0, so a 1e0 to
    1e-10 chart was labelled 1e0, 1e-4, 1e-8. The floor where the residuals
    end had no label. The label step now prefers one that divides the span
    (1e0, 1e-5, 1e-10). Every decade keeps its grid line.
  - Exponents were written three ways: ticks `1e-4`, criteria `1e-05`,
    values `5.18e-10` or `7.71e-06`. Residual values, criteria and the chart's
    criterion line now use one form, as the axis does: `7.71e-6 ≤ 1e-5`,
    `criterion 1e-5`. The console's iteration table keeps its fixed-width
    `e-05` columns.

### Fixtures

`elbow_param`, `cylinder_tutorial` and `naca4412_ground` are re-recorded
with this CLI (`app/scripts/capture_fixtures.py`, local paths scrubbed), so
CI replays the paired-state residuals and the wing's `poor_quality_cells`.
`mixing_elbow` runs Coupled, which the paired state does not change, so it
keeps its round 2 recording. Fixture total: 6.99 → 6.96 MB.

The wing now converges on downforce and drag at iteration 380 instead of
round 2's 500. The thin-layer poor-cell treatment from `05bb8da` changed its
history. Its coefficients barely moved: CL −1.213 → −1.214, CD 0.023 and
downforce 41.8 N unchanged.

### Tooling

`app/scripts/shrink_screens.py` is the libimagequant step round 2 ran from
an untracked environment: downscale, then a 256-colour palette with
Floyd-Steinberg dithering. It needs `pillow` and `imagequant` in a virtual
environment.

## Before and after

Before shots use main's app and CLI (`e308968`), after shots this branch.
All use the real CLI except the fixture row, at a 1500×940 window,
downscaled to 1200 px (`screens/`).

| Surface | Before | After |
| --- | --- | --- |
| Aerodynamics, 2D wing (Wing group, plates, framing kept) | [light](screens/before-aero-wing-light.png) · [dark](screens/before-aero-wing-dark.png) | [light](screens/after-aero-wing-light.png) · [dark](screens/after-aero-wing-dark.png) |
| Aerodynamics, 2D wing, fixture CLI (iteration count) | [light](screens/before-aero-wing-fixture-light.png) · [dark](screens/before-aero-wing-fixture-dark.png) | [light](screens/after-aero-wing-fixture-light.png) · [dark](screens/after-aero-wing-fixture-dark.png) |
| Run Calculation, cylinder (framing, plates) | [light](screens/before-run-cylinder-light.png) · [dark](screens/before-run-cylinder-dark.png) | [light](screens/after-run-cylinder-light.png) · [dark](screens/after-run-cylinder-dark.png) |
| Residuals, elbow_param (long badge, 1e-9) | [light](screens/before-residuals-param-light.png) · [dark](screens/before-residuals-param-dark.png) | [light](screens/after-residuals-param-light.png) · [dark](screens/after-residuals-param-dark.png) |
| Monitors, cylinder converged at 1700 (1e-10 axis, residual format) | [light](screens/before-monitors-cylinder-light.png) · [dark](screens/before-monitors-cylinder-dark.png) | [light](screens/after-monitors-cylinder-light.png) · [dark](screens/after-monitors-cylinder-dark.png) |
| Run History, 2D wing (poor mesh numerics) | [light](screens/before-history-wing-light.png) · [dark](screens/before-history-wing-dark.png) | [light](screens/after-history-wing-light.png) · [dark](screens/after-history-wing-dark.png) |

`docs/app-screens` is regenerated: all 22 screens plus a new
`wing-monitors` (the wing's report convergence with its poor-mesh-numerics
line), in both themes. The elbow, cylinder and wing use the real CLI; the
`car-*` screens stay synthetic. 46 files, 7.8 MB, up from 6.8 MB, through
`app/scripts/shrink_screens.py`.

## Tests

Run on the committed tree.

| Check | Result |
| --- | --- |
| `ctest --test-dir build/app -L app` with `KUTTA_REAL_CLI`, `KUTTA_ELBOW_CASE`, `KUTTA_FLUENT_EXAMPLES` | 10/10 passed |
| `ctest --test-dir build/app -R realcli -V` | ran, not skipped: 8 passed (inventory, live run + post-processing, convergence + aero, planar laminar cylinder, k-epsilon elbow, stop) |
| `./build/kutta-tests --test-suite-exclude=slow` with `KUTTA_FLUENT_EXAMPLES` | 259 test cases, 7,454,253 assertions, all passed, no fixture absent (includes `tests/app/check_solve.py`, which now checks `poor_quality_cells` and `poor_cell_quality`) |
| `python3 tools/lint_portability.py` | clean |
| Elbow byte identity vs main (CPU, Metal; plain and jsonl; 40 iterations) | fluxes.csv, forces.csv, loads.json, result.vtu identical (16/16) |

New app tests: `tst_plotting::chartHeaderKeepsTitleAndBadges` and the
1e-10 log axis case in `logAxisTicksOnDecades`;
`tst_viewformat::planarExternalFlowFramesTheBody`;
`tst_runhistory::poorMeshNumericsOnlyWhenUsed`; the lone-wing names in
`tst_aero::suggestedGroupsFollowZoneNames`. `tst_session` checks the
replay's aero iteration count and the wing's suggested group.

A first quick-tier run reported 1.65 M assertions. `tests/io/test_fluent_real.cpp`
reads the repository's `validation/fluent_examples`, not
`KUTTA_FLUENT_EXAMPLES`, and skips quietly when it is empty. The corpus and
the restricted `validation/data` sets are linked into this worktree from the
main checkout (gitignored, not committed).

Idle repainting: nothing added animates. The plates are painted with the
overlay, which repaints only when the view, data or size changes.

## Open UX issues

- **Wake heuristic.** The 2D frame finds the wake from where the domain has
  room, not from the inlet. A far field centred on the body would get one
  body length all round. Reading the flow direction from the inlet or the
  reference velocity would be more direct.
- **Wing and the ground.** The wing's frame is clipped at the ground, which
  is the bottom of the domain, so in a wide viewport an empty band shows
  below the field. That is round 2's behaviour, unchanged here.
- **Criterion label on the curve.** In the cylinder's converged chart the
  "criterion 1e-5" label sits on the continuity line where it crosses the
  criterion. This predates this round.
- **Plates in light mode** tint the field behind them slightly (white at
  84 %). That is readable in every shot checked, but a user might prefer a
  stronger plate on saturated fields.
- **Fixture provenance.** The re-recorded summaries say `kutta_commit`
  `e308968`. They were recorded with this branch's CLI, whose only change
  from `e308968` is the additive `poor_cell_quality` and the human-output
  formatting.
- **Carried over from round 2:** the Courant number and the transition
  model are read-only, the cylinder's q·A = 1 makes its coefficients equal
  its forces, and `tst_realcli` has no wing test (the case is generated).
