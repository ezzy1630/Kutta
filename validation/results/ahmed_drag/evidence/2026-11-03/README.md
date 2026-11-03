# Evidence for the 2026-11-03 investigation

These are Kutta run receipts, solver monitors, integrated loads and aggregate
comparison errors. They contain no restricted reference curves, full fields,
generated meshes or slow-tier XML. `<checkout>`, `<reference-checkout>` and
`<baseline>` replace machine-specific roots.

The compact drag matrix uses the frozen `e014e36` solver and the bench adapter
at `1c36f4d`, with explicit controls recorded in each receipt. The historical
settling matrix uses a different mesh and force window; its files are under
`settling/`. A failed run's last forces are invalid. A finite run is accepted
only when both its residual and force-window gates pass.

`summary.json` compares common supported Cp taps across the finite compact
runs. `cut_band.json` integrates clipped wall polygons over the first 40 mm in
x after the roof edge, which is about 48.8 mm along the slant. Regional loads
and Cp are final-iterate diagnostics. The legacy `edge_pressure_drag_final` field in `summary.json` is the old centroid-selected band; use `cut_band.json` for the exact geometric loads. The band file also integrates the
reversed wall-shear area in each exact clipped band; this is a discrete
near-wall-flow sign diagnostic, not a measured separation-onset coordinate. A drag average from an unconverged steady
trajectory is not a physical time average or an accepted steady solution.

The settling monitor arrays have columns `[iteration, Cl, Cd, max_residual,
max_relative_force_change]`; `null` means that a full force window was not yet
available. A continuation suppresses early stopping while retaining the same
6,000-iteration cap and gate evaluations. Its source adaptation is recorded in
`history/`. The candidate build receipt pins the header patch and binaries
used for verification. Parity records separate the rejected default proposal, explicit restoration
controls and the final automatic defaults.

The rejected inner-layer geometry experiment is archived in
`mesh/plane_layer_experiment.py.gz`. To reproduce it, restore the script to the
path in `mesh/plane_layer_method.json` in the pinned checkout, generate the
compact baseline with its recorded mesher settings, and run the script with
the mesh-generation Python environment. It checks an exact baseline replay
and preserves the wall and exterior geometry. It is an archived experiment,
not an installed mesher option.

The three `plot_*_2026_11_03.py` scripts in the report directory recreate
the figures using only these archived aggregates and own-run monitors.
They do not need the private reference fixtures. The complete slow-tier XML
and comparator output remain local; the exported comparison removes INFO
contexts and measured query coordinates. Repeated T3A coverage checks retain
only their changing own-solution upper bound and count.

The proposed automatic pressure relaxation 0.5 passed the existing slow gates
but was not retained. Earlier stopping made the transition S1223 last-2,000
average include startup and worsened its reported drag error. The final
automatic numerical defaults remain unchanged; pressure 0.5 is an explicit
Ahmed control. `verification/pressure_proposal.patch.gz` records the exact
rejected header patch against `bc81ad8` used by the full slow-tier binary.
Final default parity and the final quick tier are recorded separately.

The complementary full-stack plane-offset experiment is archived in
`mesh/full_plane_layer_experiment.py.gz`, with its source hash and restore
path in `mesh/full_plane_layer_method.json`. Unlike the inner-layer morph,
it moves the cap and incident exterior points while keeping all boundary
points and connectivity fixed. Neither recipe is an installed mesher mode.

`analysis/span_pressure.py.gz` and `analysis/span_loads.py.gz` reproduce the
span aggregates using their recorded restore paths. Tap coordinates and
reference curves remain local; exported comparisons contain aggregate errors
and integrated loads only.
