# Desktop app, round 4

Branch `t3/e52613ef`, based on main `1e5d5ae`. All numbers are from this
checkout's Release build (Qt 6.12.0, M4 Pro, Metal) unless stated otherwise.
"Main" means a build of `1e5d5ae`'s app kept beside it.

## Result

- **Getting Started card.** Fixed at every inspector width from 320 to 720 px
  and with long case names. Every scripted scenario now checks its inspector
  page at those widths before the screenshot, and the smoke test fails on
  overlapping, clipped or cut-short text. With main's card put back, the check
  reports the symptoms from the mesh_to_run screenshot: step 2's title is
  "24 px tall, its text needs 48", and "Contours" is "38 px wide, needs 69".
- **Paths.** The app never shows the home folder. It writes it as `~`
  everywhere, and narrow fields show `…/folder/file` with the full path in the
  tooltip. Every committed screenshot that showed a home path was regenerated
  through `app/scripts/capture_screens.sh` (docs) or the same scenarios
  (rounds 2 and 3, mesh_to_run) and shrunk with `shrink_screens.py`. That is
  72 files: 38 in docs/app-screens, 18 in round 2, 10 in round 3 and 4 in
  mesh_to_run. Before saving a shot, the app now refuses (exit 4) if any
  visible text contains the home folder.
- **New case.**
  - Fixture CI covers Save Case with setting changes.
  - The fixture serves a 3D starter with wheels.
  - The New Case page states its template rather than showing a disabled
    picker.
  - kutta new-case's notes sit in a "Before you run" box.
  - A wheel kept as a wall is marked "Stationary wall" in its row.
- **Found by the new checks and fixed.** These pages were wider than the
  narrowest panel, which squeezed them so text overlapped or was cut off:
  - the General page's notes button (this is what broke the card);
  - the settings forms, so Run Calculation;
  - the button rows on Run History, Compare and Design Studies;
  - the run badges on Aerodynamics;
  - the check box on Surface Integrals;
  - the New Case page's three action buttons.
- **No CLI changes.** `docs/app-protocol.md` is unchanged. No file under `src/`
  changed (`git diff 1e5d5ae -- src` is empty), so solver outputs are the
  same build of the same code as main. The byte comparison was not rerun.

## What changed

### Getting Started and narrow panels

The card's layout: each step's action now sits under its text instead of in
a third column. The result links (Contours · Forces · Run History) are one
wrapping line of links. Case and file names longer than 24 characters get
break points after `_` and `-`. That card change alone did not fix it. The
check showed the real cause: on the General page, the button "10 notes on
settings that do not affect this case" was 319 px wide. In a 320 px panel
that forced the page below its minimum width, and every wrapped label on it
was then laid out for a width it did not get. The count is now wrapping
text, and the button says "Show Notes".

The check (`ui/ScreenChecks`, run by `Inspector::layoutProblemsAtWidths`)
lays the current page out again at each width from 320 to 720 px in 20 px
steps, sizing it as the scroll area would. It reports:

- a control narrower than its minimum;
- a wrapped label shorter than its text needs at its width;
- a control extending past any widget that holds it;
- any two controls overlapping.

It runs before every scenario's screenshot. The screenshot smoke test
(`ctest -L app`) adds three scenarios to its set: `new-case-car` (the New
Case page), `new-case-run` (the card after a run, under the long case name
`fsae_front_wing_ground_effect_study_rev12_final`) and `new-case-save`.

Other pages it caught are fixed the same way:

- Settings form labels wrap. A method picker's label now spans its row with
  no alignment, so its height follows its width.
- Button and badge rows use a new `FlowLayout`, which wraps instead of
  squeezing.
- The Surface Integrals check box reads "Also compute on Fluent's field",
  with the data-file requirement in its tooltip.

Scripted windows now open without taking focus. During this round's
captures, typing elsewhere on this shared machine opened the command palette
in one shot.

### Paths

The rule is in `core/PathDisplay`, and `tst_pathdisplay` covers it.

- `homeRelativePath` and `withHomeAsTilde` write the home folder as `~`. A
  sibling folder whose name starts with home's is not home, and Windows
  separators are handled.
- `compactPath` keeps the file name and its two nearest folders.

Where it applies:

- **Page labels** (`kit::label`): free text, such as an error that names a
  file, shows from `~`.
- **Path fields** (`kit::pathLabel`): the case's file, data and mesh, the New
  Case mesh and a run's folder are compact, with the full path in the
  tooltip.
- **Status bar:** "Solver: ~/…/t3-e52613ef/build/kutta" and status messages.
- **Console:** whole paths from `~`, which still work in a shell.
- **About box:** the solver's description.

**Guard against regressions.** I chose an in-app check at capture time over
a lint rule, for these reasons:

- A PNG cannot be checked for text without OCR, and `lint_portability.py`
  skips binary files.
- The app knows exactly what it puts on screen. In every scripted scenario,
  it refuses to save a shot if the home folder appears in visible text
  (labels, buttons, fields, the console, table cells, tabs or the window
  title).
- The capture script also makes run folders under `/tmp`. On macOS `$TMPDIR`
  is a long per-user path that filled the console.
- The smoke test runs with its data under the build folder, which is under
  home. So CI exercises the check, and would fail if a page showed a raw path.

The before shots here were taken with main's app on a CLI copy and case
links under `/tmp` (as rounds 2 and 3 did), so they show main's raw paths
without a home folder in them.

### New Case from Mesh

- **Template:** the field that looked like a disabled input is now the
  template's name in bold ("FSAE external aero") above its description.
- **Choose Another Mesh** moved to the Mesh section. The actions are now
  Create Case and Cancel, which fit a 320 px panel.
- **Notes:** kutta new-case's notes (reference area, polyhedral WFGC, no
  inlet or outlet, stationary wheels) are in a "Before you run" box above
  Create Case. Before, they were faint text below the zone list.
- **Wheels:** the physics is unchanged; wheels are still stationary walls.
  Each wheel kept as a wall has an amber line in its row: "Stationary wall:
  this wheel will not turn (see Before you run, below)". The box adds where
  wall motion is set, because the app does not edit it yet: "Wall motion is
  set in the case file (a boundary's motion, docs/kutta-case.md); the app
  shows it but does not edit it yet."

### Fixture CLI

- `convert` writes these settings into a Kutta case where the real CLI
  writes them:
  - reference values;
  - the operating pressure;
  - a velocity inlet's speed, turbulence intensity and viscosity ratio.
- `inspect` of a Kutta case reads those settings back over the recording's
  inventory, with Fluent's defaults where the case leaves them out. A
  starter now shows its own area of 1 m² and its own speed, not the
  recording's 0.25 m² and 15 m/s.
- Other changes are refused by name, for example "…to save
  boundary/outlet/gauge-pressure."
- `new-case` serves `fsae_car*.msh`: the synthetic car's ten zones, roles by
  the CLI's suggestion rules, groups by the app's aero grouping, and the
  CLI's role-dependent notes in its words. The wheels note appears because
  `tyre-front` and `tyre-rear` stay walls.
- Tests: `tst_session::newCaseFromMeshCreatesRunsAndSaves` now saves a
  changed area and inlet speed, reopens, and checks the settings and the
  file, then checks the refusal. `newCaseNotesStationaryWheels` covers the
  car.

## Before and after

All shots are 1500×940 windows downscaled to 1200 px (`screens/`), except
`card-long-small`, which is 1200×800. Real-CLI rows use the wing mesh from
`tools/fsae_aero/naca_ground_case.py --write msh` (53,881 cells). Main's
fixture cannot save a changed setting, so that scenario fails on main:
"save: The fixture CLI cannot write setting changes into a case; point the
app at the real kutta CLI to save them."

| Surface | Before | After |
| --- | --- | --- |
| Getting Started after a run, real CLI | [light](screens/before-card-wing-light.png) · [dark](screens/before-card-wing-dark.png) | [light](screens/after-card-wing-light.png) · [dark](screens/after-card-wing-dark.png) |
| Same, long case name | [light](screens/before-card-long-light.png) · [dark](screens/before-card-long-dark.png) | [light](screens/after-card-long-light.png) · [dark](screens/after-card-long-dark.png) |
| Same, long name in a 1200×800 window | [light](screens/before-card-long-small-light.png) · [dark](screens/before-card-long-small-dark.png) | [light](screens/after-card-long-small-light.png) · [dark](screens/after-card-long-small-dark.png) |
| Case just opened: paths on General, console, status bar (real CLI, elbow) | [light](screens/before-open-elbow-light.png) · [dark](screens/before-open-elbow-dark.png) | [light](screens/after-open-elbow-light.png) · [dark](screens/after-open-elbow-dark.png) |
| New Case, mesh with a wheel zone (real CLI) | [light](screens/before-new-case-wheel-light.png) · [dark](screens/before-new-case-wheel-dark.png) | [light](screens/after-new-case-wheel-light.png) · [dark](screens/after-new-case-wheel-dark.png), end of page [light](screens/after-new-case-wheel-bottom-light.png) · [dark](screens/after-new-case-wheel-bottom-dark.png) |
| New Case, NACA wing (fixture) | [light](screens/before-new-case-fixture-light.png) · [dark](screens/before-new-case-fixture-dark.png) | [light](screens/after-new-case-fixture-light.png) · [dark](screens/after-new-case-fixture-dark.png) |
| New Case, synthetic car with wheels (fixture) | new | [light](screens/after-new-case-car-light.png) · [dark](screens/after-new-case-car-dark.png), end of page [light](screens/after-new-case-car-bottom-light.png) · [dark](screens/after-new-case-car-bottom-dark.png) |
| Case saved with reference area 0.5 m² and reopened (fixture) | fails (above) | [light](screens/after-save-fixture-light.png) · [dark](screens/after-save-fixture-dark.png) |

The wheel mesh is `tests/data/meshes/fsae_box.msh` with its `front_wing`
zone renamed `wheel_front`, in a scratch folder. The real CLI groups the zone
into Wheels and adds its wheel note.

### Regenerated screens

| Set | Files | Size (MB, before → after) |
| --- | --- | --- |
| `docs/app-screens` (all 23 screens, both themes; README updated) | 46 | 8.1 → 9.8 |
| `app_round2/screens` after-* | 18 | 3.0 → 4.6 (whole folder) |
| `app_round3/screens` after-* with the real CLI | 10 | 3.1 → 3.4 (whole folder) |
| `mesh_to_run/screens` real-CLI | 4 | 1.3 → 1.4 (whole folder) |

The before-* shots of rounds 2 and 3 carry no home path (they used `/tmp`
solvers and run folders) and are unchanged, as are the fixture-CLI shots.
The regenerated "after" shots show today's app, not those rounds' commits.
The rounds' features are all there, plus later changes such as compact
paths and wrapped button rows. The round reports' text is unchanged.

## Tests

Run on the committed tree.

| Check | Result |
| --- | --- |
| `ctest --test-dir build/app -L app` (fixture) | 11/11 passed, including `app-tst_pathdisplay` and the screenshot smoke test (9 scenarios × 2 themes, with screen checks, 35 s) |
| `ctest --test-dir build/app -R realcli -V` with `KUTTA_REAL_CLI`, `KUTTA_ELBOW_CASE`, `KUTTA_FLUENT_EXAMPLES`, `KUTTA_RUNS` | ran, not skipped: 8 passed (inventory, live run + post-processing, convergence + aero, planar laminar cylinder, k-epsilon elbow, stop, new case from a bare mesh run and saved, elbow saved as a Kutta case) |
| `./build/kutta-tests --test-suite-exclude=slow` with `KUTTA_FLUENT_EXAMPLES`, `KUTTA_RUNS` | 268 test cases, 7,457,460 assertions, all passed; the 35 skipped are the slow suite; no fixture reported absent |
| `python3 tools/lint_portability.py` | clean |
| Every after-shot captured this round (100: 22 here, 46 docs, 28 rounds 2 and 3, 4 mesh_to_run) | passed the screen checks; main's card code fails them |

The gitignored corpus and restricted `validation/data` files are linked into
this worktree from the main checkout, for the quick tier's corpus tests. They
are not committed.

Idle repainting: nothing added animates. `FlowLayout` lays out only when its
row's geometry changes.

## Open issues

- **New Case notes after a role change.** The page's own lines follow role
  edits immediately: the wall list, the "Stationary wall" marker and the
  wall-motion line. The CLI's notes come from reading the mesh and refresh
  only when the case is created. For example, "No zone is a velocity inlet"
  stays after you pick an inlet. This predates this round. Re-reading after
  each edit would cost a mesh read on large cars; an app-side rule, or
  structured notes from the CLI, would fix it.
- **Wall motion is read-only in the app.** Rotating wheels need an edit to
  the case file. The page says so.
- **What the screen checks cover.** They check the inspector page that a
  scenario ends on, not the charts, viewport overlays or table cells (the
  home-folder check does read table cells). The 320 px problems on Run
  Calculation show only in the check: resizing the window does not shrink
  the panel to its minimum, so a 1200×800 shot looks the same before and
  after.
- **Build paths in CLI errors.** Some CLI errors carry the build's source
  path, for example ".../src/kutta/io/fluent_msh.cpp:358: cannot open Fluent
  mesh". The app shows these from `~`, but the CLI still prints the build
  machine's layout. That is outside this round's app scope.
- **Screen set size.** docs/app-screens grew from 8.1 to 9.8 MB and round 2's
  folder from 3.0 to 4.6 MB at the same width and palette settings.
