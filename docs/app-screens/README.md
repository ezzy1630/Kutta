# Desktop app screens

Captured with `app/scripts/capture_screens.sh` from the real `kutta` CLI
(app round 4, based on `1e5d5ae`, Metal on an M4 Pro) and the Ansys mixing elbow
(`validation/fluent_examples/mixing_elbow`, with its `.dat.h5` as the Fluent
reference where noted). Each run uses the approved warped-face gradient
approximation, so results are marked APPROXIMATE. The `car-*` screens are the
exception: the fixture CLI's synthetic FSAE car, marked SYNTHETIC, because no
real car case is on hand. Their numbers show the pages, not CFD results. The
`geometry-*` screens are New Case from Geometry from the fixture's synthetic
car STLs through the fixture's `mesh`, which speaks the real `kutta mesh`
contract but meshes nothing (SYNTHETIC report and mesh). The `mesher-*`
screens are the same path on the real `kutta mesh`: the mesher's synthetic
Ahmed body and FSAE car from `tools/mesher/generate_geometry.py`, meshed with
its coarse setups, solved and post-processed in the app (dark only; numbers in
`validation/results/app_geometry_flow/REPORT.md`; the
`mesher-ahmed-stationary-*` and `mesher-ahmed-slice` screens are a later
round on a stationary ground, numbers in `docs/app.md`). The
`wing-*` screens are a real Kutta run of the 2D ground-effect wing from
`tools/fsae_aero/naca_ground_case.py`, `cylinder-run` the Ansys adjoint
cylinder tutorial. Window 1500×940, downscaled to 1500 px wide and reduced to
256-colour palettes with libimagequant (Floyd-Steinberg dithering) by
`app/scripts/shrink_screens.py` to keep the set small. Every screen but `mesher-*` exists as
`-dark.png` and `-light.png`. Paths show from `~` and compact, never the home
folder; the app refuses to save a shot whose visible text contains it.

| Screen | File | What it shows |
| --- | --- | --- |
| Welcome | `start` | The getting-started steps, with Open Case on the first |
| Case not opened | `open-failed` | Why the file could not be opened, in place of a dialog |
| Case just opened | `elbow-open` | The case's mesh before any run (zones coloured, polygon edges) and the getting-started steps |
| Boundary condition | `elbow-setup` | Methods in Fluent's words, only the fields the method uses, units, modified markers, per-field reset |
| Command palette | `palette` | Cmd/Ctrl-K over commands, pages and settings |
| All case settings | `settings` | 9,857 settings grouped Fluent-style; Kutta's controls first, raw settings folded |
| Live run | `elbow-live` | Contours updating during the run, residuals at the monitor interval, Stop |
| Finished run | `elbow-run` | Velocity on the symmetry plane, residual history |
| Kutta vs Fluent | `elbow-fluent` | Kutta's field beside Kutta − Fluent on a diverging map |
| Post-processing | `elbow-post` | Pathlines from both inlets, an iso-surface, Δ pressure with the section, wall pressure along the section against Fluent |
| XY plot | `elbow-xy` | Velocity across the outlet leg, Kutta solid and Fluent dashed |
| Surface integrals | `elbow-integrals` | Area-weighted pressure on the ports, Kutta and Fluent columns |
| Run history | `elbow-history` | Two runs with their settings and residuals overlaid |
| Compare runs | `elbow-compare` | Settings diff, force deltas and the Run 2 − Run 1 field |
| Failed run | `car-failed` | kutta's error, what to do about it, and the failed run in the top bar (synthetic car, fixture CLI) |
| Convergence | `elbow-monitors` | Fluent's residual criteria read from the case; the seeded elbow stopped at iteration 100 because every checked residual reached its criterion |
| Report convergence | `car-monitors` | Residuals plus downforce and drag report conditions (change over the last 5 values below 0.1 %), and which ones stopped the run (synthetic car) |
| Transition model | `wing-models` | SST with the Intermittency Transition Model and its constants; the 2D wing's mesh framed on the body |
| 2D wing aerodynamics | `wing-aero` | A section's CL, CD, L/D and downforce for the reference depth, up along +y, the run stopped when downforce and drag settled; the zone `wing` suggested into a Wing group; title and legend on backing plates over the full-bleed field |
| 2D wing convergence | `wing-monitors` | Downforce and drag conditions met at the same sample, and the 4 poor-quality cells the solver treated with low-order face pressures (Fluent's poor mesh numerics) named beside the verdict |
| Laminar 2D | `cylinder-run` | The cylinder framed with three diameters of its wake, three residuals from iteration 1, the APPROXIMATE badge naming the approximation |
| Aerodynamics | `car-aero` | CL, CD, L/D and balance; ClA/CdA, % viscous drag, pitching moment, CoP from the front axle; breakdown by front wing, rear wing, wheels and body; vehicle and zone groups (synthetic car) |
| Design study | `car-study` | An inlet speed sweep with the reference velocity following: CL, CD, L/D, % front and CoP per speed, MODIFIED variants marked, CD vs speed plotted (synthetic car) |
| Import Geometry | `geometry-import` | Two STL files read in the app and shown at once: units guessed (mm), six solids as named wall zones, tyres flagged Rotating wheel, Replace/Add for parts |
| Wind tunnel | `geometry-tunnel` | Tunnel extents in car lengths, moving ground, half-car option; the tunnel drawn round the car as you type |
| Refinement regions | `geometry-refinement` | Body and wake boxes as numbers in mm and as outlines in the viewport |
| Sizing and layers | `geometry-sizing` | Surface sizes, cells per gap and the smallest gap cell, layers, the freestream speed with a y+ hint, Load mesh.json… |
| Generate Mesh | `geometry-meshing` | Cell and device-memory estimate, then `mesh` held at Add Boundary Layers: stages, progress, console, status-bar Cancel |
| Mesh Report | `geometry-report` | Cells by type, the solver memory estimate with the 12 GB warning, slots and gaps (narrowest open, touching parts apart), layer coverage, snap error, y+ bands, quality; the generated mesh's surfaces in the geometry's colours (SYNTHETIC) |
| Case from the mesh | `geometry-case` | New Case from Mesh on the generated mesh, the case file chosen for the project, tyres already rotating with the estimated axle to check |
| Design variants | `geometry-variants` | A rear-wing variant meshed and run: CL, CD, L/D, front % and ΔCL, ΔCD, Δ front % against the base, the two runs side by side |
| Stationary ground | `mesher-ahmed-stationary-tunnel` | The Wind Tunnel page's ground motion set to Stationary for the Ahmed body (a wind-tunnel model; cars keep Moving), with the one-line reason (real `kutta mesh`) |
| Converged on its verdict | `mesher-ahmed-stationary-monitors` | That Ahmed run (stationary ground, 144,877 cells, SIMPLEC, Metal) converged at iteration 1521 by the CLI's `verdict`: residuals, both report conditions and settled forces met; Cd 0.4535, CL +0.448 |
| Contour slice | `mesher-ahmed-slice` | Velocity magnitude on the run's y = 0 plane with the body, framed on its wake, on the Contours page (0 cells skipped, so no note) |
