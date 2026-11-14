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
`--list-scenarios` names them. Without `--screenshot` the window stays open
on the scenario's state, with its throwaway settings, for driving by hand. `elbow-*` scenarios are meant for the real CLI
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
`KUTTA_SCENARIO_SCROLL=bottom` (or a fraction such as `0.5`) scrolls the
inspector page that far down.

`geometry-*` go from STL files (the positional arguments, else the fixture's
synthetic car STLs) through New Case from Geometry: `-import`, `-tunnel`,
`-refinement`, `-sizing`, `-meshing` (held at the layers stage with
`KUTTA_FIXTURE_MESH_HOLD=layers`), `-report`, `-case` (the New Case page from
the generated mesh), `-run` (a run and the Aerodynamics page) and `-variants`
(a rear-wing variant meshed, run and compared side by side). With the real
CLI (`--cli`) they mesh with `kutta mesh` (the meshing shot is taken as the
mesher reaches the layers stage) and take, for the real path:
`KUTTA_SCENARIO_MESH_SETUP` (a mesh.json loaded after import),
`KUTTA_SCENARIO_SPEED`, `KUTTA_SCENARIO_ROLES` (`ground=wall,...`),
`KUTTA_SCENARIO_WHEELS` (generate_geometry.py's `fsae-wheels.json`),
`KUTTA_SCENARIO_REFERENCE` (`area,length`), `KUTTA_SCENARIO_VEHICLE`
(`front-axle-x,wheelbase`) and `KUTTA_SCENARIO_ITERATIONS`; runs stop on the
mesher qualification's criteria (residuals 1e-3; drag and downforce within
5e-4 over 50 values after 500 iterations) or the iteration limit. When the
CLI cannot export the run's view, `geometry-run` prints why and goes on to
Aerodynamics (other scenarios fail at once with the CLI's error):

```sh
uv run tools/mesher/generate_geometry.py --out geom
KUTTA_SCENARIO_MESH_SETUP=geom/fsae-coarse.json KUTTA_SCENARIO_SPEED=20 \
  KUTTA_SCENARIO_WHEELS=geom/fsae-wheels.json KUTTA_SCENARIO_VEHICLE=-0.78,1.73 \
  KUTTA_SCENARIO_ITERATIONS=1000 Kutta --cli build/kutta --scenario geometry-run \
  --screenshot fsae.png geom/fsae.stl
```

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

## New Case from Geometry

File ▸ New Case from Geometry (Ctrl+Shift+N, or drop STL files on the window)
goes from STL files to a case without the terminal. The outline's Meshing
group holds the steps in Fluent Meshing's words:

- **Import Geometry**: STL files read in the app (ASCII with any number of
  named solids, or binary as one solid named after the file) and shown at
  once. Units are guessed from the size (above 20 "m" means mm) and can be
  changed. Each solid becomes a wall zone; names are editable and solids
  named wheel/tyre/tire are flagged Rotating wheel. Add, Replace… and Remove
  edit the file list.
- **Wind Tunnel**, **Local Refinement Regions**, **Sizing and Layers**:
  mesh.json's fields as numbers, drawn in the viewport as you type (tunnel
  faces translucent, tunnel and regions as outlines). The tunnel is the one
  `kutta mesh` builds: flow +x and +z up (the only directions v1 meshes),
  distances in body x lengths, the ground at z = 0 (geometry below it is
  refused before meshing), upper limits rounded out to whole base cells.
  Sizing and Layers also holds the freestream speed (the mesher's y+
  estimate and the case's starting speed) and **Load mesh.json…**, which
  takes a kutta mesh setup file's settings (keys it omits get the mesher's
  defaults).
- **Generate Mesh**: a cell estimate and the device memory it needs against
  the device the case will run on (a planning figure: surface bands by size,
  layers, regions and a graded far field, at about 2.6 kB per cell), then
  `kutta mesh` with its stages in the console, the page and the status bar,
  cancellable. An error shows the mesher's message and advice in place.
- **Mesh Report**, from mesh-report.json, each figure "not reported" when the
  report lacks it:
  - size, cells by type and the hex share of the core; the solver memory
    **estimate** (2.6 kB per cell, measured on Metal for kutta mesh meshes)
    with a warning above about 4.3 million cells, which probably will not
    fit a 12 GB GPU;
  - **slots and gaps** (`castellation/gaps`): closed slots (no fluid cells
    across) in an amber box, the narrowest open ones, and parts that only
    touch (a tyre on its plinth, width 0) listed apart;
  - **boundary layers**: generated of requested, prism cells, achieved first
    height, terminations, and coverage of the layered walls (full, thinned,
    none) as a bar; "none" when no layers were requested;
  - **distance from the geometry** (`snap/point_distance_m`): the wall's
    remaining error, p99 and max;
  - **first-cell y+ (estimate)**: p50, p05–p95, min–max at the mesher's
    speed and the wall area below 5, 5–30, 30–300 and above 300;
  - quality: orthogonal quality (min, mean, cells below 0.15 and 0.05, the
    report's histogram, the worst cell's position), maximum skewness,
    minimum volume, negative and nonconvex cells, and the adjacent volume
    ratio for core and layer faces (p50, p99, p999, max, faces above 20);
  - repairs, missing solids, the mesher's notes and the zones; the generated
    mesh's surfaces in the viewport. Set Up Case runs the New Case from Mesh
    flow on it with the rotating wheels and the speed carried over.

A geometry project is a folder (`geometry.json`, `mesh.json`, `mesh/`, the
case). **Duplicate as Design Variant** copies it beside itself (`name-v2`)
with the same STLs, settings and aero settings; swap or add an STL, generate,
and Set Up Variant Case copies the base case's setup onto the new mesh (zones
it lacks dropped, new solids added as car walls). Runs ▸ Design Variants lists
each variant's latest run with CL, CD, L/D, front % and their deltas against
the base, and shows two runs side by side.

Everything `kutta mesh` is, the app reaches through one adapter,
`app/src/core/MeshCli.cpp` (arguments, event and report names), coded to the
mesher's section of [app-protocol.md](app-protocol.md) ("Automatic STL
meshing"):

- `kutta mesh A.stl [B.stl ...] --setup mesh.json --out DIR --speed V
  --progress jsonl --force`; DIR holds `mesh.msh` and `mesh-report.json`.
- mesh.json lengths are metres; `units` (m or mm) scales the STL coordinates
  only. Projects saved before this (geometry.json `v` 1) wrote lengths in
  `units` and are converted when opened.
- The report is parsed strictly (a key of the wrong type is an error naming
  its path); `cells`, `faces`, `nodes` and `zones` are required.
- Tunnel zones: inlet, outlet, side-low, side-high, ground (bottom without
  ground), top, symmetry for a half car. These and bottom/symmetry are
  reserved: a solid with one of those names must be renamed first.
- `--case` is not used. `kutta new-case` on the mesh is the same
  `fsae_starter` path and writes the same case (checked: only the mesh path
  differs), and the New Case page adds what the app needs on top: a role per
  zone (a stationary ground for a wind-tunnel model), the wheels' geometry
  estimated from the mesh and confirmed before the case is written, and the
  CLI's notes before Create Case.
- Progress fractions can step back within a stage (layer retreat rounds); the
  bar holds its furthest point.
- Cancel closes stdin, then sends SIGINT after 1.5 s, then kills.

A run on a mesher mesh says how it ended in the CLI's terms: Converged,
Iteration limit, Stopped, or **Forces settled; residuals not converged** when
the CLI stopped it on settled forces with stalled residuals waived
(`residuals_waived`); Monitors heads such a run "not converged: forces
settled at iteration N", lists the forces-settled condition (each report's
range over the window) and marks waived residuals "not met, waived". The
design-study table says "forces settled, not converged" for it.

Volume-mesh slices are not shown yet (read-only on the Mesh Report page): that
needs a cut-plane export for a bare mesh from the CLI. The fixture CLI's `mesh`
meshes nothing: it speaks the same contract (arguments, setup rules, events,
zone names and every report key; `tst_meshflow` checks its report has the
shape of a real one, `app/tests/data/mesh-report-ahmed-coarse.json`), with
invented numbers that `geometry_method` says are synthetic, and the app labels
its meshes SYNTHETIC. It accepts only its synthetic car's STLs (`kutta-fixture
sample-stl DIR` writes them, standing on z = 0, with `rear-wing-v2.stl` for a
variant), and its runs follow the parts' moves. `tst_realmesh` meshes a box
with the real `kutta mesh` when the app is built with the solver (or
`KUTTA_REAL_CLI` is set). Reproducible car geometry for the real mesher comes
from `uv run tools/mesher/generate_geometry.py --out DIR`.

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
  view format, session through the fixture CLI, STL reading, mesh.json and the
  mesh report, and geometry to a variant case through the fixture's `mesh`).
  They include a real-CLI probe and report written by `export-view --post` on
  the elbow (`app/tests/data/`).
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
  saved with a change), and New Case from Geometry through the fixture's
  `mesh` (import, meshing in progress, the report, and a design variant run
  beside its base).

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
