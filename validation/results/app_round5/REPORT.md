# Desktop app, round 5

Based on `c0da65a`. Implementation: `9f5b857` (CLI diagnostics) and
`74d75dc` (app and additive inspect fields). Release build, Qt 6.12.0,
Apple M4 Pro, macOS, Metal. No changes under `src/kutta/fv` or the kernels.
All commits are local; nothing was pushed or published.

## Result

- A wall zone has a **Rotate as a wheel** control. Its radius, centre X/Y/Z
  and axle X/Y/Z are editable. Enabling it starts with radius zero, which
  blocks Create until the user supplies a positive radius. No geometry is
  inferred from a zone name, and no geometry estimate is applied.
- The page states the rule: angular speed is freestream speed divided by
  radius; the wheel's bottom moves along +x with the moving ground. The app
  normalizes the axle and chooses its sign (+z in 2D, toward -y in 3D).
  The actual rolling axis and angular speed are shown before creation.
  The axle must be across the flow; a 2D axle must be out of the plane.
  Zero/nonfinite axes, nonfinite coordinates and nonpositive/nonfinite radii
  are rejected inline, and Create is disabled. Session validation also
  refuses invalid geometry before starting the CLI.
- Creation passes one `--wheel ZONE=RADIUS,X,Y,Z,AX,AY,AZ` per enabled wheel.
  Changing a zone away from Wall removes its wheel motion. The reopened
  Boundary Conditions page shows a rotating wall, angular speed in rad/s,
  origin and direction. Unknown values from an older CLI remain unknown.
- Notes refresh from the real CLI after a 300 ms debounce. Edits immediately
  replace old advice with a checking message. Only one note query can run
  at a time; old replies cannot overwrite a newer draft or resurrect a
  discarded one. Numeric edits and note replies update their own widgets,
  preserving keyboard focus. There is no recurring idle refresh.
- `KUTTA_REQUIRE` now includes a source basename and line, rather than the
  build machine's directory. For example, a missing mesh reports
  `fluent_msh.cpp:358: cannot open Fluent mesh: /tmp/kutta-round5-missing.msh`.
  Both human diagnostics and JSON errors use this path-free source location.
- The wheel boundary screenshot exposed a retained-settings button that was
  four pixels too wide at 320 px. Its label is shorter, with the zone in its
  tooltip. Wall motion, roughness and rotation speed are not repeated in
  the fixture's Other section.

## Evidence

The real-CLI app test starts from committed
`tests/data/meshes/fsae_wheel.msh` (216 cells, planar). It enters radius
0.2 m, centre (0.4, 0.25, 0), axis (0, 0, -3), and speed 20 m/s. The saved
case and a fresh inspect both report 100 rad/s about (0, 0, 1), with absolute
motion. The created case runs three iterations successfully. This proves
case creation and motion configuration, not converged wheel aerodynamics.

The synthetic car test configures both tyre zones, checks the reversed axle
sign, changes roles, rejects invalid geometry, creates and reopens the case,
and verifies both wheels' angular speeds. A quiet interval produces no note
refresh signals. A scripted idle measurement of the New Case editor rendered
zero viewport frames over 4.75 seconds (1.26% process CPU); that viewport is
empty before a case is created, so this is not a heavy-scene benchmark.

Before the fix, `newCaseNotesStationaryWheels` failed after changing Inlet
to Wall because the inlet warning never arrived. The regression now checks
both adding and removing that warning, after the current query finishes.
The baseline missing-mesh command also reproduced the full source-path leak;
the rebuilt CLI retains the filename and message without that directory.

`inspect --json` adds rotating-wall `angular_speed`, `axis_origin` and
`axis_direction` condition fields. This only changes inspection output.
`docs/app-protocol.md` documents the addition. No recorded app case fixture
needed re-recording; the synthetic CLI now accepts and preserves wheel flags.

## Screenshots

Every image was inspected in light and dark for clipping, overlap, home paths
and wording. The scenario checks also test inspector widths 320–720 px and
reject visible home-directory paths. Captures use the actual field controls,
including keyboard input, rather than directly injecting wheel state.

These 14 PNGs were captured at 2× scale and reduced to 1500 px wide with
`uv run --no-project --with pillow --with imagequant python3 -I app/scripts/shrink_screens.py`.
Total size: 8.9 MB before, 2.1 MB after. The two baseline wheel shots use the
existing `new-case` scenario, which also marks Top as Wall; the new wheel
scenario keeps Top as Symmetry. The car boundary captures hide Ground in
the viewport to show the vehicle and wheels; the boundary itself is retained.
The default ground-visible dark preview was observed obscuring the synthetic
car; its rendering cause is outside these wheel-configuration changes.

| Surface | Before | After |
| --- | --- | --- |
| New Case, committed 2D wheel mesh, real CLI | [light](screens/before-wheel-light.png) · [dark](screens/before-wheel-dark.png) | [light](screens/after-wheel-light.png) · [dark](screens/after-wheel-dark.png) |
| New Case, synthetic car, end of page | [light](screens/before-car-light.png) · [dark](screens/before-car-dark.png) | [light](screens/after-car-light.png) · [dark](screens/after-car-dark.png) |
| Zero-axis error and disabled Create | unavailable | [light](screens/after-wheel-invalid-light.png) · [dark](screens/after-wheel-invalid-dark.png) |
| Created wheel boundary, real CLI | unavailable | [light](screens/after-wheel-boundary-light.png) · [dark](screens/after-wheel-boundary-dark.png) |
| Created front tyre boundary, synthetic car | unavailable | [light](screens/after-car-boundary-light.png) · [dark](screens/after-car-boundary-dark.png) |

## Verification

| Required check | Result |
| --- | --- |
| Configure Release, Ninja, app ON, Qt 6.12.0; build | passed |
| `ctest --test-dir build/app -L app` | 11/11 passed, including screenshot smoke (11 scenarios × 2 themes) |
| `ctest --test-dir build/app -R realcli -V` | 11 Qt test entries passed, 0 failed, 0 skipped; all real-CLI scenarios ran |
| `build/kutta-tests --test-suite-exclude=slow` | 277 test cases, 7,457,592 assertions passed; 37 slow tests excluded, no missing-fixture skip messages |
| `python3 tools/lint_portability.py` | clean |
| Solver parity vs `c0da65a`, elbow, 40 iterations, `--approximate wfgc`, CPU and Metal, plain and JSONL | fluxes.csv, forces.csv, loads.json, result.vtu and view-input.bin byte-identical in all four runs; residual histories identical excluding elapsed time |
| New Case idle check | zero rendered viewport frames over 4.75 seconds |

Test environment: `KUTTA_REAL_CLI=<checkout>/build/kutta`,
`KUTTA_ELBOW_CASE=<fluent-examples>/mixing_elbow/mixing_elbow.cas.h5`,
`KUTTA_FLUENT_EXAMPLES=<fluent-examples>`, `KUTTA_RUNS=<runs>`; the external
corpus and generated fixtures are the ones specified in the task. Ignored
fixtures were linked for reads only, then those links were removed. No
restricted source data is committed. [Check output](checks.txt) and
[physical-output hashes](parity.json) retain the verification evidence.

## Limits

- Geometry entry is manual; surface-derived estimates are not implemented.
- CUDA and Windows were not available for runtime verification. The CPU and
  Metal paths, macOS Qt rendering and the portability lint were exercised.
- The synthetic car is demonstration data. No real 3D car's wheel geometry
  or converged aerodynamic loads were validated here.
- The design-context tool reported stale platform/build-path metadata in its
  existing project context. That metadata was left unchanged in this scoped
  implementation.
