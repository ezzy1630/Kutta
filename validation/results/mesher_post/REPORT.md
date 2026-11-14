# Post-processing, ground and run verdict on `kutta mesh` meshes

The app runs STL → `kutta mesh` → `kutta new-case` → `kutta solve`, but on a
mesher mesh it could show no field afterwards, put a plain-wall ground into the
car's forces, and called a run "converged" whose residuals were waived. This
round fixes all three in the CLI. Everything here was run on synthetic
geometry from `tools/mesher/generate_geometry.py` meshed with its coarse
setups (Ahmed: 144,877 cells; FSAE: 1,312,313 cells), Apple M4 Pro, Metal,
on a Mac shared with other lanes.

## 1. `export-view` on mesher meshes

**Before.** `kutta export-view RUN --slice 0,0,0.1,0,0,1` after a 30-iteration
Ahmed solve: exit 1, `view.cpp:115: slice intersection is not a closed
polygon`. Every slice the app asks for failed, so no contour, slice or vector
view; forces were unaffected.

**Cause.** Not the hanging faces themselves. Snapping leaves the nodes of an
octree level scattered within about ±2e-9 m of it, while the slice tolerance
on this mesh is 1.2e-9 m. The slicer classified nodes three ways (above, on,
below), so one face mixed "on" and "below" nodes; its triangulated faces then
produced segments to both a node and a crossing point 2e-9 m away, which the
distance-based point merge did not join, and one point ended up with three
segments. A cell with any face inside the tolerance also skipped its cut
entirely, which left holes (0.031 m² at z = 0.4) where only some 2:1 subfaces
of a coarse cell were within it.

**Fix** (`src/kutta/io/view.cpp`).
- Each node is above the plane (distance > tolerance) or not; there is no
  "on" class. Crossings lie only on edges whose ends differ and are identified
  by that edge, not by distance, so on a closed cell every crossing has
  exactly two segments and the rings always close.
- Faces of more than three nodes are cut as the triangle fan about their node
  average, the solver's own face geometry. A first version that joined the
  crossings within a face by straight chords overlapped its neighbours where a
  2:1 face's mid-edge node is not collinear (a concave pentagon, mid-node
  1e-4 m off the line); the fan follows the face, and both cells sharing a
  face cut the same surface.
- A face lying in the plane keeps the documented rule (emitted once, owner's
  values) only when it is the cell's whole section; a face only partly within
  the tolerance is cut like any other. A ring that touches itself (a fan centre
  on the plane) is split into simple rings.
- A cell whose faces do not close, or whose section cannot be triangulated
  (a tangled cell), is skipped and counted in the slice descriptor's
  `skipped_cells` (v1 additive) instead of failing the export.

Iso-surfaces had the same failure mode (`src/kutta/post/iso.cpp`): FSAE's
`post operation "wake": iso decomposition has no positive cell kernel at cell
3683` failed the whole export. A cell without a positive kernel or with an
inverted decomposition tetrahedron is now skipped whole (nothing partial or
inverted is emitted) and counted in the iso descriptor's `skipped_cells`.

**After.**

| check | Ahmed coarse | FSAE coarse |
| --- | --- | --- |
| planes at every octree level (6.25 mm Ahmed, 31.25 mm FSAE) in x, y, z, each also 1 µm off | 491 planes, 0 cells skipped | 208 planes, 11 cells skipped in total (tangled cells, below) |
| slice area vs the plane 1 µm away | equal except at horizontal body facets; each jump equals those facets' area (z = 0.1: 2.1156e-4 m² both; z = 0.325: 1.3814e-3 vs 1.3816e-3 m²) | jumps of 4–8e-3 m² at car facets (not itemised) |
| slice area vs domain | z = 0.4 (above the body): 48.7600 m² = 10.6 × 4.6 m; z = 0 and the top: the same | |
| overlap between triangles (exact clipping) | about 1e-8 m² per slice (float32 roundoff) | |
| slices, iso-surfaces (Cp0 = 0, Q = 100, u = 0), streamlines, pathlines, surface-seeded traces, line and section probes, integrals, in one `--post` batch | all written; 0 cells skipped; 4.3 s | all written; slices 0 skipped, each iso-surface 104 cells skipped of 1.31 M; 19.7 s |

Checked on the 30-iteration runs and again on the converged stationary-ground
Ahmed run; a rendering of its z = 0.1 and y = 0 slices shows no hole.

Regression tests: `tests/io/test_view.cpp` meshes the mesher's box with prism
layers and a refinement box (2:1 hanging faces), slices exactly through an
octree level, jitters that level's nodes by ±3 tolerances and slices again
(and the plane ±1e-4 m and with the normal reversed), and slices a diagonal
plane through lattice vertices; each must close, skip nothing and give the
exact area (tunnel minus the 0.36 m² body). Against the old `view.cpp` it
throws `slice intersection is not a closed polygon`. A second test checks the
owner rule on a coincident interior face and that a cell whose faces do not
close is skipped and counted. `tests/post/test_iso_probe.cpp` checks that a
C-shaped cell with no kernel is skipped and counted.

## 2. Stationary ground

**Before.** `new-case` offered the ground `moving-ground` or a plain `wall`;
a plain wall is a car wall, so it joined the drag and downforce reports.
Rerun here: Ahmed with `--role ground=wall`, Cd 1.192, downforce coefficient
1.65 (the ground's shear and pressure counted as the body's).

**Fix.** New role `stationary-ground` (v1 additive, appended to `roles`):
`{"type":"wall","motion":"stationary"}`, no slip at zero velocity, never a
report zone. Ground/road/belt zones are still suggested `moving-ground`, the
FSAE default. A zone that looks like ground but is given the plain `wall` role
now gets a note saying it is in the car's reports and naming the two ground
roles. Documented in `docs/app-protocol.md` and `docs/kutta-case.md`.

**After** (Ahmed coarse, `new-case --speed 40`, reference area 0.112032 m²,
SIMPLEC, Metal; all three converged with residuals met):

| ground | report zones | iterations | Cd | CL (up) |
| --- | --- | --- | --- | --- |
| `--role ground=stationary-ground` | body | 621 | 0.4535 | +0.448 |
| `--role ground=moving-ground` (default) | body | 641 | 0.4376 | +0.422 |
| `--role ground=wall` (before) | ground, body | 561 | 1.1917 | −1.65 |

The stationary-ground Cd agrees with the mesher qualification's Ahmed run on a
stationary ground (Coupled, 5,000 iterations: Cd 0.454, CL +0.451). It is the
coarse mesh's answer (y+ about 108, 145 k cells), well above the experimental
Ahmed drag; that gap is mesh resolution, not the ground.

## 3. Run verdict

**Before.** A run stopped on settled forces with stalled residuals had
`converged: true`, `stop_reason: "converged"` and only
`convergence.residuals_waived: true` to tell it apart.

**Fix.** `done` and `summary.json` add `verdict` (v1 additive):
`converged`, `forces_settled_residuals_waived`, `statistically_steady`,
`iteration_limit` or `stopped`. A diverged iteration's `error` event adds
`verdict: "diverged"`; other errors have none. `converged` and `stop_reason`
keep their meaning. `tests/app/check_protocol.py` drives every value through
the real CLI: a waived run (unreachable continuity criterion, forces-settled
window 1), a run diverging under an absurd inlet speed, and the existing
converged, statistically steady, iteration-limit and stopped runs.

## Mesher findings (for the mesher lane; not changed here)

- **Ledges on the Ahmed body.** The snapped body has small horizontal facets
  at octree levels (2.1e-4 m² at z = 0.1, 1.4e-3 m² at z = 0.325, more at
  0.1625, 0.225, 0.2625, 0.2875). A wall-parallel slice shows them as jumps.
- **A non-manifold fin on the Ahmed body.** The z = 0.2 body section is three
  open curves meeting at (−0.026, −0.1945) and (−0.0135, −0.2025), 8 mm
  outside the body's side: the snap error maximum (7.84 mm, one cell size) in
  the mesh report.
- **Tangled cells on the FSAE car.** 104 cells cannot be decomposed for an
  iso-surface (no positive kernel or an inverted tetrahedron), and a few more
  give self-intersecting slice sections; the three slice cells examined have
  inverted fan tetrahedra (faces folded about the cell centre). The export
  skips them; the solver still uses them.

## App changes needed (not made here; `app/` belongs to the app lane)

- Read `verdict` from `done`/`summary.json` (and from the `error` event)
  where `Protocol.cpp` and `RunHistory.cpp` now combine `converged`,
  `stop_reason` and `residuals_waived`; keep the combination as the fallback
  for runs from older CLIs, which have no `verdict`.
- New Case from Mesh: label `stationary-ground` ("Stationary ground") in
  `roleLabel` (`app/src/ui/NewCasePage.cpp`). The role list comes from the
  CLI, so the option already appears, under its raw name.
- Geometry flow: the Wind Tunnel page has no ground choice;
  `tunnelFaces` (`MeshSetup.cpp`) always gives the ground `moving-ground`. Add
  a moving / stationary choice there, pass `--role ground=stationary-ground`
  to `new-case` for a stationary one, and let `MeshWorkflow.cpp` treat that
  role as the ground (it maps only `moving-ground` to `zone:ground`).
- Show `skipped_cells` when it is non-zero on a slice or iso-surface ("N
  tangled cells could not be drawn").

## Gates (this branch)

- Quick tier, `./build/kutta-tests --test-suite-exclude=slow` with the Fluent
  examples and generated runs from the main checkout (CPU/Metal parity
  included): 340 test cases pass, 37 skipped, 7,997,146 assertions.
- App tests, `ctest --test-dir build-app/app`: 14 of 14 pass.
- `tools/lint_portability.py` passes.
- Nothing in `src/kutta/fv` or the kernels changed, so no slow tier.
