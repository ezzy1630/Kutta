# Desktop app: geometry to aero results

A team can now go from car STL files to a run, its aerodynamics and a
design-variant comparison inside the app, without the terminal. The meshing
step runs against the `kutta mesh` contract from the mesher lane. That mesher
has not landed yet, so every mesh here comes from the fixture CLI's stand-in,
and those screens are marked SYNTHETIC.

## The path

| Step | Where | Screens (`docs/app-screens/`) |
| --- | --- | --- |
| Import STLs (ASCII multi-solid or binary, mm/m guessed), solids as named wall zones, wheels flagged | File ▸ New Case from Geometry, drag and drop, Meshing ▸ Import Geometry | `geometry-import-{dark,light}.png` |
| Wind tunnel, ground, half car; refinement boxes; surface sizes, curvature, proximity, layers | Meshing ▸ Wind Tunnel / Local Refinement Regions / Sizing and Layers; drawn in the viewport as you type | `geometry-tunnel-*`, `geometry-refinement-*` |
| Cell and memory estimate against the selected device, Generate with stages, progress and Cancel in page, console and status bar | Meshing ▸ Generate Mesh, Mesh ▸ Generate Mesh (Ctrl+Shift+M) | `geometry-meshing-*` |
| Mesh report: cells by type, quality histograms, worst cells, repairs, zones, y+; generated surfaces in the viewport | Meshing ▸ Mesh Report | `geometry-report-*` |
| Case from the mesh (New Case from Mesh flow, wheels rotating, case file chosen for the project), run, Aerodynamics | Set Up Case, then the existing pages | `geometry-case-*` |
| Duplicate as design variant, swap the rear-wing STL, mesh, Set Up Variant Case (base case settings copied), run, compare with ΔCL, ΔCD and Δ front % side by side | Mesh ▸ Duplicate as Design Variant; Runs ▸ Design Variants | `geometry-variants-*` |

With the synthetic car, the larger rear wing (`rear-wing-v2.stl`) gives
ΔCL −0.251 (−9.5 %), ΔCD +0.036 (+10.8 %) and moves the front share from
56.5 % to 44.9 %. These are fixture numbers that show the comparison works;
they are not CFD results.

## Design

- One adapter: `app/src/core/MeshCli.cpp` spells every flag, event, file
  and report key of `kutta mesh`. `MeshWorkflow` owns the project folder
  (`geometry.json`, `mesh.json`, `mesh/`, the case). STLs are read and drawn
  in the app (`core/Stl.cpp`); the CLI runs only to mesh.
- The case is written with the existing `kutta new-case`, so wheel axes are
  estimated and confirmed as for any mesh. A variant copies its base case's
  setup onto the new mesh: zones the new mesh lacks are dropped, and new
  solids are added as car walls.
- No CLI or solver changes. The fixture gained `mesh` (it validates the
  solids, emits the progress events, writes a SYNTHETIC report and an
  estimate-sized mesh, and honours cancel) and `sample-stl`. Its runs follow
  each part's transform, so a swapped wing changes the forces.
- Assumptions about the contract are listed in `docs/app.md` under "New Case
  from Geometry". The main ones: mesh.json lengths are in `units`, tunnel
  distances are in body lengths, a binary STL is one solid named after the
  file, renamed solids are passed as ASCII copies, `--case` is unused, report
  keys are read tolerantly, tunnel zones are named inlet/outlet/ground/top/sides,
  and `new-case` accepts `mesh.msh`.

## Checks (this branch, macOS, Metal)

- `ctest -L app -LE gui`: 12 of 12 suites pass. New: `tst_geometry` (STL
  ASCII and binary, welding, crease split, truncation, rename round trip,
  mesh.json in mm, setup problems, tunnel bounds, estimate, event and report
  parsing) and `tst_meshflow` (import, generate, stale, reopen, cancel,
  mesher error with advice, case from mesh, variant with a swapped wing, mesh
  colours).
- `tst_realcli` with the real CLI built from the base commit (0747475) on the
  Ansys elbow: 11 of 11 pass.
- `app-screenshot-smoke` (label `gui`), now including geometry-import,
  -meshing, -report and -variants: passes (38 shots).
- All nine `geometry-*` scenarios pass in dark and light at 1500×940. I
  inspected each shot.
- End to end through the real UI with cua-driver (button presses and
  menus, not scenario calls), starting from a fresh data folder:
  1. Outline ▸ Meshing ▸ Import the Synthetic Car, then Mesh ▸ Generate Mesh.
  2. The report opened; Set Up Case, confirm both wheels, Create Case.
  3. Approve the approximation and Run.
  4. Duplicate as Design Variant, Generate, Set Up Variant Case, approve, Run.
  5. Compare Design Variants and Side by Side: identical geometry gave
     ΔCL = ΔCD = 0, as it should.
- `tools/lint_portability.py` passes.

## Stubbed pending `kutta mesh`

- Meshing is fixture-only. The fixture accepts only its synthetic car's six
  solids, so adding a new part (as opposed to swapping one) is covered by code
  and docs but not demonstrated end to end.
- The mesh report's numbers are invented; the cell counts follow the app's own
  estimate.

## Open problems

- No volume-mesh slice. The Mesh Report says so (read-only). It needs a
  cut-plane or section export for a bare mesh from the CLI; the
  `export-geometry` command draws boundary zones only.
- The cell and memory estimate is a planning figure (±50 % at best). Calibrate
  it against real `kutta mesh` reports once they exist.
- Refinement boxes and the tunnel are edited as numbers; there is no direct
  manipulation in the viewport yet.
- The native file dialog could not be driven from cua-driver while the app was
  in the background, so the file pick itself was checked by scenario (the same
  `importGeometry` call) rather than by clicking.
- Each variant case asks for the approximation approval again (approvals are
  per case). That is safe, but it is friction when iterating.
- `make_test_stl.py` writes a closed car for trying the real mesher; it has
  not been meshed, because the mesher does not exist yet.
