# Desktop app screens

Captured with `app/scripts/capture_screens.sh` from the real `kutta` CLI
(`t3/fluent-workbench-ux`, Metal on an M4 Pro) and the Ansys mixing elbow
(`validation/fluent_examples/mixing_elbow`, with its `.dat.h5` as the Fluent
reference where noted). Each run uses the approved warped-face gradient
approximation, so results are marked APPROXIMATE. `car-failed` is the one
exception: the fixture CLI's synthetic car, made to fail part-way as a full
disk would. Window 1500×940, downscaled to 1500 px wide. Every screen exists
as `-dark.png` and `-light.png`.

| Screen | File | What it shows |
| --- | --- | --- |
| Welcome | `start` | The getting-started steps, with Open Case on the first |
| Case not opened | `open-failed` | Why the file could not be opened, in place of a dialog |
| Case just opened | `elbow-open` | The case's mesh before any run (zones coloured, polygon edges) and the getting-started steps |
| Boundary condition | `elbow-setup` | Methods in Fluent's words, only the fields the method uses, units, modified markers, per-field reset |
| Command palette | `palette` | Cmd/Ctrl-K over commands, pages and settings |
| All case settings | `settings` | 9,844 settings grouped Fluent-style; Kutta's controls first, raw settings folded |
| Live run | `elbow-live` | Contours updating during the run, residuals at the monitor interval, Stop |
| Finished run | `elbow-run` | Velocity on the symmetry plane, residual history |
| Kutta vs Fluent | `elbow-fluent` | Kutta's field beside Kutta − Fluent on a diverging map |
| Post-processing | `elbow-post` | Pathlines from both inlets, an iso-surface, Δ pressure with the section, wall pressure along the section against Fluent |
| XY plot | `elbow-xy` | Velocity across the outlet leg, Kutta solid and Fluent dashed |
| Surface integrals | `elbow-integrals` | Area-weighted pressure on the ports, Kutta and Fluent columns |
| Run history | `elbow-history` | Two runs with their settings and residuals overlaid |
| Compare runs | `elbow-compare` | Settings diff, force deltas and the Run 2 − Run 1 field |
| Failed run | `car-failed` | kutta's error, what to do about it, and the failed run in the top bar (synthetic car, fixture CLI) |
