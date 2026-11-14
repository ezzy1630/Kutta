# Desktop app: geometry to aero results

A team can go from car STL files to a mesh, a run and its aerodynamics inside
the app, without the terminal, on the real `kutta mesh`. This round aligns
the app with the mesher that merged in fb44554 (contract: the "Automatic STL
meshing" section of [app-protocol.md](../../../docs/app-protocol.md);
qualification: [mesher REPORT](../mesher/REPORT.md)) and runs the whole path
in the app on the mesher's synthetic Ahmed body and synthetic FSAE car.

## Real path in the app (Apple M4 Pro, Metal)

Geometry from `uv run tools/mesher/generate_geometry.py --out geom`, the
mesher qualification's coarse setups (`ahmed-coarse.json`, `fsae-coarse.json`)
loaded with Sizing and Layers ▸ Load mesh.json…, then Generate, Set Up Case
(New Case from Mesh: wheels from `fsae-wheels.json`, moving ground), Run and
Aerodynamics. Runs stop on the qualification's criteria (scaled residuals
≤ 1e-3; drag and downforce changing < 5e-4 over 50 samples after iteration
500) together with the CLI's forces-settled condition, or the iteration cap.
The app runs the case `kutta new-case` writes (SIMPLEC); the qualification
ran Coupled, so the coefficients are not directly comparable. GPU and CPU
were shared with other lanes' solves, so times are contended.

| | Ahmed coarse | FSAE coarse |
| --- | --- | --- |
| geometry | 1 solid, 6,672 triangles, 1.044 m | 19 solids, 2,652 triangles, 3.106 m |
| mesh | 144,877 cells (93,313 hex, 51,090 prism, 474 poly); 1.9–3.0 s | 1,312,313 cells (524,228 hex, 782,016 prism, 6,069 poly); 62–107 s |
| app's estimate before meshing | 141 k | 1.41 M |
| solver memory (estimate, 2.6 kB/cell) | 377 MB | 3.4 GB; fits 12 GB |
| slots | none closed; body–ground 3 cells across 64.8 mm | **rear-wing-element-1–2 closed** (12.8 mm, 0 cells); front-wing-element-2–3 1 cell across 20.1 mm; 4 tyre–ground contacts listed apart |
| layer coverage full / thinned / none | 94.9 / 1.5 / 3.5 % | 76.2 / 16.9 / 6.9 % |
| snap error p99 / max | 0.12 mm / 7.84 mm | 2.76 mm / 13.5 mm |
| y+ estimate p50; area 30–300 | 108 at 40 m/s; 93 % | 52.5 at 20 m/s; 84.5 % |
| case | `kutta new-case`: SIMPLEC, SST, moving ground, 40 m/s, A = 0.112032 m² | `kutta new-case`: SIMPLEC, SST, moving ground, four wheels rotating, 20 m/s, A = 1 m² |
| iterations; stop | 2,721; **forces settled, residuals waived** (CLI `converged: true`, `residuals_waived: true`) | 2,000; **iteration limit**, not converged |
| final residuals | continuity 1.36e-3 (criterion 1e-3), others met | continuity 3.3e-2, k 1.0, ω 5.9 |
| coefficients | Cd 0.438, CL +0.422 (downforce −46.3 N, drag 48.0 N) | CL −1.463 (downforce 358.5 N), CD 1.332, 23.7 % front |
| solve wall time | 94 s (0.035 s/iteration); 202 s and 422 s in earlier identical runs | 467 s (0.23 s/iteration); 692 s in an earlier identical run |
| app status | "Forces settled; residuals not converged"; Monitors "Run 1 not converged: forces settled at iteration 2721", continuity "not met, waived" | "Iteration limit"; Monitors "reached the iteration limit (2000) without converging" |

Both runs compare with the mesher qualification (Coupled, 5,000 iterations):
Ahmed Cd 0.454 / CL +0.451 there (stationary ground) against 0.438 / +0.422
here (moving ground, SIMPLEC); FSAE CL 1.460, CD 1.330, 23.3 % front there
against 1.463, 1.332, 23.7 % here. Neither run is converged: Ahmed's
continuity stalls just above 1e-3 and the CLI waives it once the forces have
settled; FSAE's k runs away in the stagnant cells at the closed rear-wing slot,
as the qualification found. The app says so on every page that shows the
run (top bar, Monitors, Aerodynamics, run history, residual chart mark).

`kutta export-view` fails on both runs ("slice intersection is not a closed
polygon", view.cpp:115), so the viewport shows that error instead of a field;
forces, coefficients and the breakdown come from `kutta aero` and are
unaffected.

## What changed in the app

- **One adapter, the real contract.** `app/src/core/MeshCli.cpp` passes
  `--setup --out --speed --progress jsonl --force` and parses
  `mesh-report.json` strictly: `v` 1, `cells`, `faces`, `nodes` and `zones`
  are required; a key of the wrong type is an error naming its path; a
  missing or null optional key shows "not reported". The model follows the
  real report: `cell_types`, `quality` (orthogonal histogram as bins,
  cells below 0.05/0.15, worst cell, skewness, volume ratio core/layers),
  `castellation.gaps`, `layers.coverage`, `snap.point_distance_m`,
  `estimated_first_cell_y_plus`, `repairs`, `warnings`,
  `missing_surface_zones`, `stage_seconds`.
- **Units.** mesh.json lengths are metres; `units` scales STL coordinates
  only. geometry.json is v2 (the mesh setup in metres plus the speed); v1
  projects are converted when opened.
- **Tunnel.** Sized as the mesher sizes it (body x lengths, ground at z = 0,
  upper limits rounded out to whole base cells); `tst_realmesh` checks the
  app's tunnel against the real report's bounds. Zones are inlet, outlet,
  side-low, side-high, ground (bottom without ground), top, symmetry; solids
  with those names must be renamed. Setup ranges match the mesher's parser.
- **Mesh Report.** Closed slots in an amber box (and on Generate Mesh),
  narrowest open gaps, touching contacts apart; layer coverage full / thinned
  / none as a bar ("none" when no layers were asked for); snap error p99 and
  max; y+ bands (below 5, 5–30, 30–300, above 300) at the mesher's speed;
  the solver memory **estimate** at about 2.6 kB per cell with a warning
  above about 4.3 M cells (12 GB GPU) and a check against the selected
  device. SYNTHETIC appears only when the fixture made the mesh.
- **Cell estimate.** Recalibrated on the real coarse meshes (octree levels
  quantized to base/2ⁿ): within 5 % of both given the layered wall area;
  from the whole STL area, which is what the app has before meshing, it reads
  141 k against 144,877 (Ahmed) and 1.41 M against 1,312,313 (FSAE, +8 %).
  `tst_geometry` holds both within 10 %.
- **Honest stops.** The CLI can now stop a run as `converged` on settled
  forces with stalled residuals waived (`residuals_waived`). The app shows
  that as "Forces settled; residuals not converged", never Converged; the
  run history mark, the stop reason, Monitors ("not met, waived") and the
  forces-settled condition's per-report range follow.
- **`--case` is not used.** `kutta new-case` writes the same case (only the
  mesh path differs), and the New Case page adds what the app needs: a role
  per zone, the wheel axes estimated and confirmed, and the CLI's notes.
- **Fixture.** `kutta-fixture mesh` emits exactly the real contract (strict
  setup keys, the same events and every report key, `geometry_method` saying
  SYNTHETIC); `tst_meshflow` compares its report's shape with a real one
  (`app/tests/data/mesh-report-ahmed-coarse.json`). It is for tests and
  screenshots only. `app/scripts/make_test_stl.py` is gone
  (`generate_geometry.py` replaces it).

## Screens (`docs/app-screens/`)

Real `kutta mesh` (dark, window 1500×940):

| Screen | Ahmed | FSAE |
| --- | --- | --- |
| Import Geometry | `mesher-ahmed-import` | `mesher-fsae-import` |
| Wind Tunnel | `mesher-ahmed-tunnel` | `mesher-fsae-tunnel` |
| Sizing and Layers (coarse setup loaded) | `mesher-ahmed-sizing` | `mesher-fsae-sizing` |
| Generate Mesh, held at the layers stage | `mesher-ahmed-meshing` | `mesher-fsae-meshing` |
| Mesh Report: size, memory, slots | `mesher-ahmed-report` | `mesher-fsae-report` |
| Mesh Report: layers, snap, y+, quality | | `mesher-fsae-report-mid` |
| Monitors after the run | `mesher-ahmed-monitors` | `mesher-fsae-monitors` |
| Aerodynamics | `mesher-ahmed-aero` | `mesher-fsae-aero` |

I inspected each one. The import to report shots come from the `geometry-*`
scenarios with `--cli build-app/kutta` (the env vars in `docs/app.md`); the
Monitors and Aerodynamics shots are window captures of the `geometry-run`
scenario's window after its run, reached through the command palette.

The fixture's `geometry-*` screens (dark and light) were regenerated on the
new contract and still show the SYNTHETIC path for the pages the real runs
do not reach (refinement regions, variants).

## Checks (this branch, macOS, Metal)

- `ctest` in `build-app/app`: 14 of 14 pass, including `app-screenshot-smoke`
  (label `gui`) and `tst_realmesh`, which meshes a box in mm with the real
  `kutta mesh` through `MeshWorkflow`, checks the app's tunnel against the
  report's bounds and sets up the case (it runs, not skipped, when the app is
  built with the solver). New or changed: `tst_geometry` (metres in mesh.json,
  setup ranges, tunnel sizing, the estimate against both real meshes, the
  mesh command, a real report read and strict-parse errors), `tst_meshflow`
  (the fixture's report has the shape of the real one), `tst_runhistory`
  (forces settled with waived residuals is not converged).
- `tst_realcli` with `KUTTA_REAL_CLI` and the Ansys mixing elbow: 9 pass, 2
  skipped (they need the other Fluent examples).
- Root quick tier, `./build/kutta-tests --test-suite-exclude=slow`: 333 test
  cases pass, 37 skipped. `src/` is unchanged.
- `tools/lint_portability.py` passes.
- The fixture's `geometry-*` scenarios pass in dark and light; the real-path
  shots above were taken through the app on this Mac and inspected.
- Both real runs were driven in the app end to end (import, load the coarse
  mesh.json, generate, New Case from Mesh with the wheels, run, Monitors,
  Aerodynamics); times are on a GPU shared with other work.

## Contract issues found (not changed here; owned by the mesher and CLI lanes)

1. **`export-view` fails on mesher meshes**: "slice intersection is not a
   closed polygon (view.cpp:115)" on both the Ahmed and FSAE runs, so no
   field can be shown in the app after a run on a `kutta mesh` mesh.
2. **`inspect` and `new-case` count cell types differently from the mesher**:
   the Ahmed mesh the mesher reports as 93,313 hex / 51,090 prism / 474 poly
   is "136,151 hex, 8,725 poly, 1 tet" to `inspect`, and `new-case` notes a
   tetrahedron.
3. **`proximity_min_size` 0 is rejected** though the contract says 0 means
   min_size / 4, and the error prints a range that does not say why.
4. **Progress fractions step back** within the layers stage (0.80 → 0.78),
   against "monotonically increasing" in the contract. The app holds the bar
   at its furthest point.
5. **Contacts are reported as gaps**: a tyre on the ground is a gap of width 0
   and 0 cells, the same shape as a closed slot apart from the width. The app
   separates them by width; a `kind` (or leaving contacts out) would be
   clearer.
6. **The volume-ratio histogram counts faces under the key `cells`.**
7. **`new-case` (and `mesh --case`) writes SIMPLEC**, while the qualification
   solves Coupled, so the app's numbers and the qualification's are not like
   for like.
8. **No stationary-ground role**: `new-case` can make the ground a moving wall
   or a plain wall, and a plain wall is added to the car's force reports
   (Ahmed Cd 1.19 that way). A wind-tunnel model like the qualification's
   Ahmed run needs a stationary ground that stays out of the reports.
9. **`converged: true` with `residuals_waived: true`** is a converged stop by
   the CLI's word on an unconverged solution. The app now handles it, but
   other consumers that read only `converged` or `stop_reason` will call such
   a run converged.
10. FSAE `skewness_max` is 1.36 (> 1), which a skewness should not exceed.

## Open problems

- No volume-mesh slice: it needs a cut-plane export for a bare mesh from the
  CLI; `export-geometry` draws boundary zones only.
- New Case from Mesh has no stationary-ground role that stays out of the car's
  force reports, so the Ahmed runs here use a moving ground, unlike the
  qualification's stationary one.
- Refinement boxes and the tunnel are edited as numbers; there is no direct
  manipulation in the viewport yet.
- Each variant case asks for the approximation approval again (approvals are
  per case).
