# Desktop app screens

Captured with `app/scripts/capture_screens.sh` from the real `kutta` CLI
(app round 2 on `868f478`, Metal on an M4 Pro) and the Ansys mixing elbow
(`validation/fluent_examples/mixing_elbow`, with its `.dat.h5` as the Fluent
reference where noted). Each run uses the approved warped-face gradient
approximation, so results are marked APPROXIMATE. The `car-*` screens are the
exception: the fixture CLI's synthetic FSAE car, marked SYNTHETIC, because no
real car case is on hand. Their numbers show the pages, not CFD results. The
`wing-*` screens are a real Kutta run of the 2D ground-effect wing from
`tools/fsae_aero/naca_ground_case.py`, `cylinder-run` the Ansys adjoint
cylinder tutorial. Window 1500×940, downscaled to 1500 px wide and reduced to
256-colour palettes with libimagequant (Floyd-Steinberg dithering) to keep the
set small. Every screen exists as `-dark.png` and `-light.png`.

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
| 2D wing aerodynamics | `wing-aero` | A section's CL, CD, L/D and downforce for the reference depth, up along +y, the run stopped when downforce and drag settled |
| Laminar 2D | `cylinder-run` | The cylinder's plane, three residuals from iteration 1, the APPROXIMATE badge naming the approximation |
| Aerodynamics | `car-aero` | CL, CD, L/D and balance; ClA/CdA, % viscous drag, pitching moment, CoP from the front axle; breakdown by front wing, rear wing, wheels and body; vehicle and zone groups (synthetic car) |
| Design study | `car-study` | An inlet speed sweep with the reference velocity following: CL, CD, L/D, % front and CoP per speed, MODIFIED variants marked, CD vs speed plotted (synthetic car) |
