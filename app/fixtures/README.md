# Fixture CLI data

`kutta-fixture` (app/tools/kutta-fixture) speaks the same protocol as the
`kutta` CLI from this data, so the app, its tests and CI run without a solver
build or Fluent case files. The real CLI always takes precedence when the app
can find one; every fixture result is labelled FIXTURE REPLAY or SYNTHETIC.

- `cases/mixing_elbow/`: a real Kutta run of the Ansys mixing elbow (MIT),
  recorded on an M4 Pro with Metal. `inspect*.json.qz` are the CLI's inventory
  (9,844 settings), `progress-*.jsonl` its event streams, `run-cold` and
  `run-seeded` the finished runs' outputs and views (without and with the
  Fluent data file), `live/iter-*` three live surface snapshots. Replay changes
  nothing when settings are overridden; the run says so. The recording has no
  volume snapshot, so streamlines, iso-surfaces and line probes on the elbow
  need the real CLI; sections and surface integrals work from its surfaces.
  It predates convergence criteria and `loads.json`, so it also covers the
  app's paths for an older CLI (no Monitors controls; Aerodynamics asks for a
  re-solve).
- `fsae_car*`: synthetic, generated in code (`SyntheticCar.cpp`), to exercise
  external-aero screens and the 10M-triangle viewport benchmark. Nothing about
  it is a CFD result.
- `devices-*.json`: device inventories for macOS and an RTX Windows machine.

`.qz` files are `qCompress`ed to keep the repository small.
