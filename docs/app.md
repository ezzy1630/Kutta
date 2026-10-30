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
replays a recorded elbow run and generates a synthetic car. Its results are
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

`car-*` and `bench` use the fixture's synthetic car. Scenarios use the fixture
unless `--cli` or `KUTTA_CLI` names a solver. The fixture can fail a run on
request (`KUTTA_FIXTURE_FAIL_AT=N`, with `KUTTA_FIXTURE_FAIL_MESSAGE` for an
error event, or none for a crash); `car-failed` uses it.

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
  with Fluent fields, mass balance, stop. It skips otherwise.
- `app-screenshot-smoke` (label `gui`): renders scenarios through Metal or
  D3D11 (WARP on GPU-less Windows) and fails when a viewport showing data
  draws nothing.

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

On the real elbow (9,844 settings, 41 editable): `inspect` 0.6 s; the All Case
Settings tree builds in about 120 ms.

## Packaging

`cmake --build <build> --target package-app` writes `Kutta-<version>-macos.zip`
(macdeployqt, ad hoc signed, CLI and its HDF5 libraries inside the bundle) or
`Kutta-<version>-windows.zip` (windeployqt, CLI in `solver\`). The CLI comes
from `KUTTA_APP_SOLVER_DIR` or the in-tree build; without one, the fixture is
bundled. CI jobs `app-macos` and `app-windows` build, test, run the smoke test
and upload both packages; Windows bundles the `windows-cuda` job's kutta.exe.
