# Desktop app: requests for the CLI protocol

`docs/app-protocol.md` is owned by the CLI. These are the app's proposed
changes and observations, against `t3/implement-cli-app-protocol` at `57fcdc3`,
measured on the Ansys mixing elbow (17,822 cells, M4 Pro, Metal). The app
already works with the protocol as it is; none of these block it.

## Wanted

1. **Labels for enum values.** Raw settings carry codes such as
   `boundary/cold-inlet/turbulence-specification = 3`. An `options` entry of
   `{"value": 3, "label": "Intensity and Hydraulic Diameter"}` (or a parallel
   `option_labels`) would let the app show Fluent's words instead of numbers.
   Editable enums (`methods/*-scheme`) are fine as they are.

2. **Lenient derived fields in `export-view --post`.** One unavailable name in
   `fields` fails the whole export (`unknown or unavailable post-processing
   field`). The app wants a fixed useful set on every surface (Cp, total
   pressure coefficient, skin friction, vorticity magnitude, x velocity)
   without first exporting to read `available_fields`. Today it retries without
   `fields`. Proposal: `"fields_policy": "available"` to skip unavailable ones
   and list them in the manifest.

3. **Incremental post-processing.** Every operation republishes the whole
   `view/`, so adding one iso-surface recomputes every trace, probe and report
   in the batch. That's fine on the elbow (0.2 s plain; 5.2 s for 30 surface
   traces, an iso-surface, two probes and a report); it won't be on an
   18M-cell car. Proposal: reuse unchanged operations' outputs from the
   previous view, keyed by the operation's JSON.

4. **Run controls should not mark results MODIFIED.** `run/iterations` and
   `run/monitor-interval` in a setup file set `modified: true` and
   `MODIFIED RESULT`, but they don't change the case being solved. The app now
   passes the iteration count with `--iterations`, so it no longer triggers
   MODIFIED. A changed monitor interval still does.

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
- **Inventory size.** `inspect --json` returns 9,844 settings (26 editable) in
  0.6 s. The app groups them Fluent-style, keeps retained raw settings in a
  collapsed folder per group, and builds the searchable tree in about 120 ms.
- **Approximation text.** The `APPROXIMATE RESULT` warning is about 600
  characters. The app shows it in full in the console (wrapped) and in tooltips,
  and as an APPROXIMATE badge on every result. A short machine-readable summary
  per approximation (for example `{"approximation":"wfgc","short":"WFGC off"}`)
  would make badges more specific.
- **Stop.** Closing stdin through QProcess's pipe stops the run after the
  current iteration with outputs written; verified by the app's real-CLI test.
- **Wording.** The protocol doc still says "sidecar" in the stdin paragraph;
  the app is a native Qt child-process host.
