# Desktop app: requests for the CLI protocol

`docs/app-protocol.md` is owned by the CLI. These are the app's proposed
changes and observations, against `t3/implement-cli-app-protocol` at `57fcdc3`,
measured on the Ansys mixing elbow (17,822 cells, M4 Pro, Metal). The app
already works with the protocol as it is; none of these block it.

## Done (2026-10-30)

All four requests are in the CLI and used by the app; `docs/app-protocol.md`
has the wire format.

1. **Labels for enum values.** Every inventory entry carries `option_labels`
   (parallel to `options`) and `value_label`, in Fluent's words: the cold
   inlet's turbulence specification reads "Intensity and Hydraulic Diameter",
   not 3. Entries also carry `active`, `depends_on` and `active_when`, so the
   app hides the fields a method does not use and updates them as the method
   changes. The turbulence specification method is editable on k-omega inlets
   and outlets, offering only methods whose values the case saved.
2. **Lenient derived fields.** `"fields_policy": "available"` skips fields a
   run cannot provide and lists them in the manifest's `unavailable_fields`.
   The app asks for it and says in the console which were skipped; against an
   older CLI it falls back once to the CLI's own selection.
3. **Incremental post-processing.** Unchanged base surfaces, slices and
   operations are reused from a per-run cache keyed by their inputs; the
   manifest reports `post_cache` `{reused, computed, seconds}`. On the elbow,
   adding a Q-criterion iso-surface to a 35-unit batch takes 0.23 s instead of
   3.0 s; re-exporting an unchanged batch takes 0.02 s. Warm output is
   byte-identical to cold.
4. **Run controls do not mark results MODIFIED.** The `setup` event and
   `summary.json` list them under `run_controls`; `modified` covers
   case-changing overrides only. The app treats the iteration count and the
   monitor interval the same way.

Also new: `kutta export-geometry CASE --out DIR [--json]` exports a case's
boundary without a run (0.16 s on the elbow), and zones carry a KUTL edges
file of their polygon edges, so the app shows the mesh as soon as a case
opens, with face edges rather than triangle diagonals.

## Done (aero and convergence)

The CLI now gives the numbers an FSAE aero engineer designs with;
`docs/app-protocol.md` has the wire format.

5. **Convergence criteria.** The case's residual criteria and report
   convergence conditions (Fluent's lag residual over the last N values) are
   `monitors/` run controls, off by default. Iteration events carry each
   condition's state; `done` and `summary.json` say which conditions stopped
   the run. The seeded elbow with checking on stops at iteration 120 on its
   own residual criteria (2.8 s on Metal).
6. **Wall loads and `kutta aero`.** `solve` writes `loads.json` (per-wall
   pressure and viscous force and moment). `kutta aero` turns it into
   coefficients, a group breakdown, balance and centre of pressure in
   milliseconds, so editing groups or the vehicle never re-solves. On the
   elbow its totals match `forces.csv` to 4e-9 relative.

## Resolved since

- **The monitor interval changed the solution.** Found while testing request
  4: 40-iteration elbow runs at interval 20 and interval 1 differed by up to
  0.0084 m/s in u (of 1.61). The automatic pseudo time step was refreshed only
  in `SteadySolver::monitor()`, so it followed the monitor cadence. It is now
  refreshed every outer iteration, and runs at intervals 1 and 20 write
  bit-identical results on CPU and Metal
  (`tests/fv/test_monitor_cadence.cpp` covers SST, SA and realizable k-ε).

## Observations

- **Progress cadence.** `iteration` events arrive at the monitor interval (20
  on the elbow). The residual chart and console table are built for that:
  points connect across gaps, and the console prints one row per event. The
  app doesn't need finer events by default. The case's own interval is
  editable (`run/monitor-interval`, Monitors page) for anyone who wants every
  iteration. Measured cost on the elbow, one run each: 400 iterations in 8.2 s
  at interval 20 and 9.6 s at interval 1 (+17%). Most of that is likely the
  per-event residual reduction on a small mesh. The cost is worth knowing
  before considering an independent `--progress-every`.
- **Inventory size.** `inspect --json` returns 9,844 settings (41 editable) in
  0.6 s. The app groups them Fluent-style, keeps retained raw settings in a
  collapsed folder per group, and builds the searchable tree in about 120 ms.
- **Approximation text.** The `APPROXIMATE RESULT` warning is about 600
  characters. The app shows it in full in the console (wrapped) and in tooltips,
  and as an APPROXIMATE badge on every result. A short machine-readable summary
  per approximation (for example `{"approximation":"wfgc","short":"WFGC off"}`)
  would make badges more specific.
- **Stop.** Closing stdin through QProcess's pipe stops the run after the
  current iteration with outputs written; verified by the app's real-CLI test.
- **Report sampling.** Report conditions use values at the monitor interval,
  so with an interval above 1 they are stricter than Fluent's per-iteration
  check for the same N. A per-iteration force reduction would match Fluent
  exactly, at the cost measured above for residuals.
- **Residual criteria in other layouts.** Only the 7-element residual rows
  parse; elbow1's 10-element rows leave checking off with a null criterion,
  and the user must enter one to check that equation.
- **Wording.** Fixed: the protocol doc describes the app as a native Qt
  child-process host, not a "sidecar".
