# Desktop app, round 6

Base: `bc9ad1c`. Implementation: `78a1d0963d6b05809a480a5baa58be7a13d2dc93`. Release build, Qt 6.12.0,
Apple M4 Pro, macOS 26.6.2, Metal. All commits are local.

## Ground and camera

The baseline reproduces the ground-dominated view in both themes. The light
capture hides most of the car; the dark capture leaves a small car on a dense,
bright road mesh. The failure starts in boundary classification and camera fit:
the synthetic road is 9 × 5 m, inside a 20 × 8 × 4 m tunnel. The old classifier
required a plane to span 98% of the domain in two axes, so this road was treated
as body geometry. Every visible vertex then contributed to fitting the camera.
The road also used the opaque, depth-writing body pipeline and its dense edges.
Filled surfaces used slope-scaled depth bias, which displaces steep triangles
according to their screen-space slope and can put low car surfaces behind the
road or neighboring surfaces. Both sides of each surface were already drawn;
back-face culling was disabled.

Mesh views now fit visible body walls. Broad partial boundary planes and combined
tunnel shells are recognized geometrically; inlet, outlet and symmetry zones do
not control the body fit. The road remains visible by default, rendered first as
12% neutral context, without depth writes or dense edges. Other domain surfaces
use the same treatment when shown. Opaque body surfaces retain normal depth
testing; their edge bias is one constant depth unit, without a slope term.
The key labels these domain surfaces as context. Every zone remains available in
the surface checklist; a boundary-only selection still fits what is selected.

Opening another mesh also resets its camera before choosing the initial view.
The new regression opens the elbow, waits for its mesh, then opens the car. Before
the reset it fails with `new car retained the previous case's camera`. Afterward
it opens the car isometrically, shows all zones through the actual checklist and
fits the body. The existing car-boundary scenario no longer hides Ground.

## Wheel estimates

The existing `export-geometry` command already reads a bare Fluent mesh. New Case
now calls it once in a temporary directory before publishing the draft, and reads
each zone's exported surface. No CLI or solver implementation change was needed.
The documented protocol now states this bare-mesh use and its 2D extrusion.

The estimate uses an area-weighted surface covariance to find the thin principal
axis of the largest connected surface. This matters for the synthetic car: each
tyre zone contains two disconnected, coaxial tyres, whose track width exceeds
their diameter. The whole zone is then measured in that principal frame: its
bounding midpoint gives the centre, and half the larger transverse extent gives
the radius. Zone names do not participate. Float-export noise is removed and
proposals retain six significant digits. A 2D mesh uses XY geometry, Z = 0 and an
out-of-plane axle; the display's artificial unit-depth extrusion is ignored.

Selecting **Rotate as a wheel** fills the editable fields and shows
“Estimated from the tyre surface — check it”. **I checked this wheel geometry**
must be checked before Create or a wheel-bearing CLI request can proceed.
Editing a confirmed field clears confirmation. Estimates never enable rotation.
If geometry cannot be exported, the manual editor remains available with an
explicit unavailable-estimate message and the existing geometry validation.

| Input | Surface estimate | Observed result |
| --- | --- | --- |
| Committed `tests/data/meshes/fsae_wheel.msh` | Radius 0.2 m; centre (0.4, 0.25, 0); axle (0, 0, 1) | Unconfirmed creation refused; confirmed case saved, reopened and ran three iterations at 20 m/s, with 100 rad/s rotating-wall motion |
| Synthetic front tyre zone | Radius 0.232996 m; centre (-0.78, 0, 0.23); axle (0, -1, 0) | Both coaxial tyres measured together; confirmation required |
| Synthetic rear tyre zone | Radius 0.233 m; centre (0.8, 0, 0.23); axle (0, -1, 0) | Both coaxial tyres measured together; confirmation required |

The synthetic surface reaches 0.185 + 0.048 = 0.233 m, slightly beyond its nominal
0.23 m tyre radius. The estimate preserves the surface measurement. An unnamed,
tilted pair of coaxial cylinders also passes a geometry-only regression test.

## Verification

| Required check | Result |
| --- | --- |
| Release configure, Ninja, app ON, Qt 6.12.0; complete build | Passed; final app and real-CLI test targets rebuilt after their last changes |
| `ctest --test-dir build/app -L app` | 11/11 passed, including screenshot smoke (15 scenarios × 2 themes) |
| `ctest --test-dir build/app -R realcli -V` | 11 Qt entries passed, 0 failed, 0 skipped; every real-CLI scenario ran |
| `build/kutta-tests --test-suite-exclude=slow` | 277 cases and 7,457,592 assertions passed; 37 slow tests excluded, no missing-fixture skips |
| `python3 tools/lint_portability.py` | Clean |
| Populated car mesh idle measurement, both themes | Zero viewport frames during each approximately 5-second interval; 359,256 visible triangles, 598 × 446 logical viewport at 2× scale |

The idle check uses the app's existing `--benchmark 24` frame-counter and process
CPU measurement. Dark: 4.751 s, 0.73% process CPU. Light: 5.249 s, 1.49% process
CPU. These are single process samples, with other work on the machine; the idle
frame count is the redraw gate. [Raw measurements](idle.json).

Tests used `KUTTA_REAL_CLI=<checkout>/build/kutta`, the specified external elbow
case and Fluent corpus, and the specified external runs folder. Ignored fixture
files were linked for reads only and the links were removed afterward. No
restricted source data is included. No files under root `src/` changed, including
`src/kutta/fv` and `src/kutta/kernels`; the conditional solver-output parity check
therefore does not apply. [Sanitized check output](checks.txt).

## Screenshots

All retained captures were inspected in light and dark for overlap, clipping,
home paths and wording. They use temporary input and data folders. The scripted checks also cover inspector widths 320–720 px. Long
forms scroll normally; the 2D estimate capture includes every field, the unchecked
confirmation, its explanatory message and disabled Create button. The car's
manual-entry baseline was captured after the scenario typed the wheel values;
the baseline source starts the editor at radius zero.

| Surface | Before | After |
| --- | --- | --- |
| Default car mesh with Ground visible | [light](screens/before-car-light.png) · [dark](screens/before-car-dark.png) | [light](screens/after-car-light.png) · [dark](screens/after-car-dark.png) |
| Synthetic car wheel editor | [light](screens/before-wheel-light.png) · [dark](screens/before-wheel-dark.png) | [light](screens/after-wheel-light.png) · [dark](screens/after-wheel-dark.png) |
| Real 2D wheel estimate, awaiting confirmation | — | [light](screens/after-wheel-2d-light.png) · [dark](screens/after-wheel-2d-dark.png) |
| Car with all domain zones shown after opening the elbow | — | [light](screens/after-domain-context-light.png) · [dark](screens/after-domain-context-dark.png) |
| Car viewed from below the visible ground | — | [light](screens/after-below-ground-light.png) · [dark](screens/after-below-ground-dark.png) |

The 14 retained PNGs were captured at 2× scale and reduced to 1500 px wide with
`uv run --no-project --with pillow --with imagequant python3 -I app/scripts/shrink_screens.py FILES`.
The screenshot smoke test covers 15 scenarios in both themes, including default
car mesh, all domain zones, below-ground view, unconfirmed estimates and created
wheel boundaries.

## Limits

Windows/D3D and CUDA runtime verification were unavailable. The 3D car is
synthetic; no team's real 3D tyre mesh or converged car aerodynamics was tested.
Surface-derived dimensions remain proposals that the user must review.
