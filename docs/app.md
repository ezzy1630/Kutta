# Kutta desktop app

A native Qt 6 workbench for Kutta in Fluent's vocabulary: open a Fluent case,
change what Kutta can change, run on the GPU with live residuals and contours,
then compare against Fluent's own solution and against other runs. Source in
`app/`; product brief in `app/PRODUCT.md`; screenshots in `docs/app-screens/`.

The app never links the solver. It runs the `kutta` CLI as a child process over
[the app protocol](app-protocol.md); requests for that protocol are in
[app-protocol-requests.md](app-protocol-requests.md). Dependencies and licences
are in [app-licensing.md](app-licensing.md).

## Build

Qt 6.12 with Svg and ShaderTools, CMake 3.26, a C++20 compiler.

```sh
# with the solver (the app finds and bundles the in-tree kutta CLI)
cmake -S . -B build -G Ninja -DKUTTA_BUILD_APP=ON -DCMAKE_PREFIX_PATH=$HOME/Qt/6.12.0/macos
# or the app alone
cmake -S app -B build-app -G Ninja -DCMAKE_PREFIX_PATH=$HOME/Qt/6.12.0/macos
cmake --build build-app
```

## Which solver it runs

In order: `--cli PATH`, `KUTTA_CLI`, the bundled CLI (`Kutta.app/Contents/Helpers/kutta`,
`Kutta\solver\kutta.exe`), the in-tree build, `kutta` on `PATH`. Only then the
fixture CLI (`app/tools/kutta-fixture`, data in `app/fixtures/`). The fixture
replays recorded runs of four cases (mixing elbow, k-epsilon elbow, laminar 2D
cylinder, 2D ground-effect wing with transition) and generates a synthetic car. Its results are
labelled FIXTURE REPLAY or SYNTHETIC wherever they appear. The status bar
names the solver in use.

## Scenarios and screenshots

`Kutta --scenario NAME --screenshot out.png [--theme dark|light] [--size WxH]`
drives the app to a state and saves the window, GPU viewports included.
`--list-scenarios` names them. `elbow-*` scenarios are meant for the real CLI
and the Ansys mixing elbow:

```sh
KUTTA_REAL_CLI=/path/to/kutta ELBOW_CASE=/path/to/mixing_elbow.cas.h5 \
  app/scripts/capture_screens.sh build-app/Kutta.app/Contents/MacOS/Kutta docs/app-screens
```

`case-*` scenarios take any case as the positional argument and are meant for
the real CLI: `case-mesh`, `case-approvals` (the first approximation
approved), `case-methods`, `case-run` (convergence checking on,
`KUTTA_SCENARIO_ITERATIONS`, default 300) and `case-stop` (stopped by the user
at that iteration). `KUTTA_SCENARIO_PAGE` shows another page at the end:

```sh
KUTTA_SCENARIO_ITERATIONS=1000 Kutta --cli build/kutta --scenario case-run \
  --screenshot cylinder.png validation/fluent_examples/adjoint_cylinder_tutorial/cylinder_tutorial.cas
```

`car-*` and `bench` use the fixture's synthetic car. `param-run` (standard
k-epsilon), `cylinder-run` (laminar 2D), `wing-models`, `wing-aero` and
`wing-monitors` (2D wing, intermittency transition, stopping when downforce
and drag settle) open the fixture's recorded samples, or the case given as
the positional argument. Scenarios use the fixture
unless `--cli` or `KUTTA_CLI` names a solver. The fixture can fail a run on
request (`KUTTA_FIXTURE_FAIL_AT=N`, with `KUTTA_FIXTURE_FAIL_MESSAGE` for an
error event, or none for a crash); `car-failed` uses it.

`new-case`, `new-case-run` and `new-case-save` start from a bare mesh (the
positional argument, else the fixture's NACA wing mesh): the New Case page
with a role and the speed changed, the created case after a run, and that case
saved with a changed reference area and reopened. `new-case-car` is the New
Case page for the fixture's synthetic car mesh (3D, with wheels).
`KUTTA_SCENARIO_CASE_NAME` names the created case and
`KUTTA_SCENARIO_SCROLL=bottom` shows the end of the inspector page.

Before its screenshot, every scenario checks the screen and fails (exit 4)
with each problem on stderr: the inspector page is laid out again at each
panel width from 320 to 720 px, and no text may overlap other text, be clipped,
or be shorter or narrower than it needs; and no visible text (labels,
buttons, console, tables, tabs) may contain the home folder.

Paths on screen never show the home folder: the app writes it `~`
(`core/PathDisplay.h`). Narrow fields (the case's file and mesh, the New
Case mesh, a run's folder, the solver in the status bar) show a compact
`…/folder/file` with the full path in the tooltip; the console and messages
keep the whole path from `~`, which still works in a shell.

## New cases from a mesh, and saving

File ▸ New Case from Mesh reads a mesh (`.msh`, `.msh.h5`, or the mesh of a
`.cas`) with `kutta new-case --json` and shows its boundary zones, each with
the role the CLI suggests from its name (inlet, outlet, ground, symmetry
words, else the mesh's zone type, else a car wall with its aero group) and why.
Pick a role per zone and the freestream speed. The template is stated, not
picked (FSAE external aero is the only one). The CLI's notes are in a "Before
you run" box above Create Case; a wheel that stays a wall is marked
"Stationary wall" in its row, since the template's walls do not move and the
app does not edit wall motion yet (set `motion` in the case file). Create Case writes the Kutta
case (`kutta new-case --out`) and opens it, after which the usual Setup and
Solution pages edit it and Run solves it ([docs/kutta-case.md](kutta-case.md)).

File ▸ Save Case writes the open case with the current changes as a Kutta
case (`kutta convert --setup`), in place for a Kutta case and as a new
`.kutta.json` for a Fluent case (Save as Kutta Case chooses the file). The
case is reopened from the saved file, so the changes become its own settings.
A saved case with different settings is a different case: its runs start a
new history, as a changed `.cas` does.

## Case types and convergence

The app follows the case's model and dimension as the CLI reports them:

- **Methods** lays out Fluent's panel: the pressure-velocity coupling
  (SIMPLE, SIMPLEC, Coupled; Coupled is offered only with the
  pseudo-transient formulation, and a case saved with Coupled but without it
  says it runs SIMPLEC), the gradient (least squares or Green-Gauss cell
  based), pressure interpolation and the schemes of the equations the model
  solves. Run History and Run Calculation state the coupling each run used.
- **k-epsilon and laminar** cases show their own residuals (`epsilon`, or
  none), relaxation factors, inlet turbulence inputs and fields.
- **Transition**: Models shows the transition model next to the viscous
  model. The Intermittency Transition Model (one transport equation for
  intermittency on top of SST) lists its constants; Methods, Controls and the
  inlets carry its scheme, relaxation and inflow intermittency, and runs add
  the `intermit` residual and the Intermittency field. Transition SST and the
  algebraic model are named with the CLI's reason Kutta does not solve them.
- **Coupled** is offered only with the pseudo-transient formulation. Controls
  shows the case's flow Courant number read-only: Kutta's Coupled solver uses
  Fluent's default, 200, as implicit relaxation on top of the pseudo time
  step.
- **2D cases** are shown as their plane, face on with x right and y up; drag
  pans rather than orbits. An external-flow case is framed on its body with
  a body length of flow around it. The CLI does not report z-velocity for 2D
  cases. Aerodynamics treats a 2D case as a section: up is +y unless a force
  report gives another axis, there is no vehicle or balance, and forces are
  for the reference depth.
- **Approximate runs** carry badges naming each approximation in a few
  words ("APPROXIMATE · WFGC off") from the CLI's `short` labels.
- **Run status** says how a run ended: Converged, Iteration limit, Stopped,
  Failed.
- **Approximations** show the CLI's explanation, its measured effect and
  evidence separately, each with its own opt-in.
- **Convergence**: the residual chart draws the run's criteria and marks the
  iteration it stopped at with the reason (converged, stopped by you,
  iteration limit). Monitors, Run Calculation and Run History say why a run
  stopped, with each condition's last value against its criterion and the
  residual that was last to converge. The Monitors page states Fluent's
  scaled-residual definition.

## First run and failures

Opening a case shows its boundary mesh at once (`kutta export-geometry`, cached
per case) and a four-step "Getting started" card on the General page: open,
approvals, run, results. Each step states what is true and offers the next
action; Hide removes it for good (View ▸ Show Getting Started restores it).

Nothing fails silently or behind a modal dialog. A case that cannot be opened
says why on the page, and the open case stays open. A failed, crashed, stopped
or interrupted run keeps a callout on Run Calculation with kutta's message and
advice (disk full, device out of memory, permissions, crash) and stays in the
top bar until the next run. A run whose record or `setup.json` cannot be
written is refused before kutta starts.

## Tests

- `ctest -L app -LE gui`: unit tests (protocol, settings, run history, plots,
  view format, session through the fixture CLI). They include a real-CLI probe
  and report written by `export-view --post` on the elbow (`app/tests/data/`).
- `tst_realcli` runs the session against the real CLI when `KUTTA_REAL_CLI`
  and `KUTTA_ELBOW_CASE` are set: full inventory, live run, post-processing
  with Fluent fields, mass balance, stop, a new case from a bare mesh run and
  saved with a change, and the elbow saved as a Kutta case. It skips otherwise.
- `app-screenshot-smoke` (label `gui`): renders scenarios through Metal or
  D3D11 (WARP on GPU-less Windows) and fails when a viewport showing data
  draws nothing, or when a screen check above fails. Its default set covers
  the elbow setup, the synthetic car, the recorded k-epsilon, laminar 2D and
  transition cases, and the new-case flow through the fixture (the car's New
  Case page, the wing's case after a run under a long case name, and that case
  saved with a change).

## Performance (M4 Pro, macOS, Metal)

Single viewport 1820×1126 physical pixels, 4× MSAA, camera orbiting 1° per
frame (`--benchmark N` on the `bench` scenario; `KUTTA_FIXTURE_CAR_DETAIL`
scales the synthetic car):

| Triangles | Mean frame | p95 | Rate |
| --- | --- | --- | --- |
| 0.32M | 8.33 ms | | 120 fps (vsync) |
| 8.03M | 11.99 ms | 12.56 ms | 83.4 fps |
| 11.56M | 16.66 ms | 17.28 ms | 60.0 fps |

Idle: 0 frames rendered and 0.01–0.02% CPU over 5 s with no input. Nothing
animates continuously; views repaint only on change.

On the real elbow (9,856 settings, 57 editable): `inspect` 0.6 s; the All Case
Settings tree builds in about 120 ms.

## Packaging

`cmake --build <build> --target package-app` writes `Kutta-<version>-macos.zip`
(macdeployqt, ad hoc signed, CLI and its HDF5 libraries inside the bundle) or
`Kutta-<version>-windows.zip` (windeployqt, CLI in `solver\`). The CLI comes
from `KUTTA_APP_SOLVER_DIR` or the in-tree build; without one, the fixture is
bundled. CI jobs `app-macos` and `app-windows` build, test, run the smoke test
and upload both packages; Windows bundles the `windows-cuda` job's kutta.exe.
