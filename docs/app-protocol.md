# Desktop app ↔ solver protocol (v1)

The desktop app (`app/`, Tauri + TypeScript) never links the solver. It runs the
`kutta` CLI as a sidecar process and talks to it through files and line-delimited
JSON on stdout. The CLI stays the single source of truth for what a case means;
the app only presents it. Every message carries `"v": 1`.

## Commands

| Command | Output |
| --- | --- |
| `kutta devices --json` | one JSON object: `{"v":1,"devices":[{"backend":"cuda","name":"NVIDIA GeForce RTX 4070","memory_bytes":12884901888,"default":true}, ...]}` |
| `kutta inspect <case> [--data <dat>] --json` | one JSON object: case summary (below) |
| `kutta solve <case> [--data <dat>] [--approximate a,b] [--backend b] [--iterations n] --out DIR --progress jsonl` | JSON lines on stdout while running; files in DIR |
| `kutta export-view DIR [--slice origin,normal]...` | writes `DIR/view/` (below) from a finished run |

Human-readable output without `--json` / `--progress` is unchanged.

## Case summary (`inspect --json`)

```json
{"v":1,
 "file":{"path":"...","sha256":"...","fluent_version":"22.2"},
 "mesh":{"cells":17822,"faces":91581,"nodes":0,"cell_types":{"tet":0,"hex":0,"poly":17822},
         "bounds":[[xmin,ymin,zmin],[xmax,ymax,zmax]]},
 "zones":[{"id":34,"name":"wall-elbow","type":"wall","faces":2436,"condition":{...}}],
 "models":{"viscous":"k-omega-sst","energy":true,"energy_coupled":false},
 "materials":[{"name":"air","density":1.225,"viscosity":1.7894e-5}],
 "reports":[{"name":"drag","kind":"drag","zones":["car"],"direction":[1,0,0]}],
 "support":{"status":"approximate",            // "supported" | "approximate" | "unsupported"
            "issues":[{"setting":"recon/cell-lsf?","label":"Warped-face gradient correction",
                       "class":"approximable","approximation":"wfgc",
                       "explanation":"…plain language…","evidence":"src/kutta/fv/VERIFICATION.md#…"}]},
 "data":{"present":true,"fields":["SV_U","SV_P"],"iterations_hint":294}}
```

`class` is `unsupported` (blocks a run), `approximable` (runs only when the
user opts in with `--approximate <approximation>`) or `informational`.

## Progress events (`solve --progress jsonl`)

One JSON object per line, in order:

```json
{"v":1,"event":"setup","cells":17822,"backend":"metal","device":"Apple M4 Pro","approximations":["wfgc"],"init":"hybrid"}
{"v":1,"event":"iteration","iter":20,"residuals":{"continuity":1e-3,"x-velocity":1e-4,"y-velocity":1e-4,"z-velocity":1e-4,"k":1e-3,"omega":1e-3},
 "reports":{"drag":{"force":0.0123,"coefficient":0.31}},"sec_per_iter":0.042,"elapsed":0.9}
{"v":1,"event":"warning","message":"…"}
{"v":1,"event":"done","converged":false,"iterations":600,"elapsed":16.9,"out":"DIR"}
{"v":1,"event":"error","message":"…","issues":[...]}
```

`iteration` events are emitted at the monitor interval. Report values are
optional per event. stderr carries human-readable diagnostics only. Exit
codes: 0 success, 1 rejected case or solver failure, 2 usage error. The app
stops a run by closing stdin or sending SIGINT / CTRL_BREAK; the CLI then
finishes the current iteration, writes outputs and emits `done`.

## View data (`export-view`)

Small, GPU-friendly files so the viewer never loads the volume mesh:

- `view/manifest.json`: zones, fields (name, unit, min, max, location),
  slices, bounds and, if a Fluent `.dat.h5` was supplied, `"reference": true`.
- `view/<zone>.surface.bin`: triangulated boundary surface per zone
  (little-endian header + float32 positions, uint32 indices, float32 per-vertex
  fields, area-weighted from faces).
- `view/slice-<n>.bin`: the same format for each requested plane slice.
- Reference fields from the Fluent data file are exported alongside, as
  `<field>@fluent` and `<field>@delta`, so the app can show Kutta vs Fluent.

The binary layout is specified next to its writer in `src/kutta/io/` and is
versioned by the manifest.

## Version 1 details

The exact view binary layout, field names/units, interpolation and missing-data
semantics are in [the writer's format specification](../src/kutta/io/VIEW_FORMAT.md).
The app should use manifest descriptors rather than assuming every zone has
all fields. `--slice` takes **one** argument with six comma-separated finite
numbers, `ox,oy,oz,nx,ny,nz`. The normal must be nonzero and is normalized.

Solve persists `view-input.bin`; export reads that snapshot from the completed
run, without needing the source case/data files. `summary.json` is the run's
completion marker and is written after the snapshot and other outputs close.
Old runs without that snapshot must be solved again. `--data` both seeds the run
and supplies the Fluent comparison; v1 has no separate comparison-only flag.
A view retains any approximation notice as `provenance`.

`devices` includes only compiled backends with detected devices. Inventory does
not compile kernels. `ordinal` is an additive device identity within a backend;
`--backend` currently selects that backend's runtime default device. Exactly one
entry has `default: true` (CUDA, then Metal, then CPU, as in `--backend auto`).
`memory_bytes` is total physical RAM for CPU, recommended maximum working set
for Metal, and total VRAM for CUDA; zero means capacity could not be queried.
A detected device can still fail at solver setup, for example on kernel compile
or allocation. Such failure emits `error`.

Unknown optional inspect values are JSON null; missing values never imply a
Fluent default. `poly` counts every cell other than a topological tet or hex.
`energy_coupled: false` means energy is off or active constant material laws
prove decoupling; null means coupling cannot be excluded. The support status is
computed by the same setup, initialization and report validation as solve.
When `--data` is supplied, saved fields replace initialization validation.
Informational issues do not block solve; the only approximable setting is the
explicit setup allowlist (currently `wfgc`). An unparsed WFGC flag is unsupported,
not approximable. A downforce report is normally a lift report whose direction
points down; its sign and reference scaling are preserved.

On valid runs `setup` is the first event, followed by warnings and iterations;
`done` is last and means outputs are closed. Rejected/setup-failed runs may emit
just `error`. No `done` follows an error. Usage errors return 2 and also emit
`error` when `--progress jsonl` is present. The iteration's `reports` object is
keyed by exact Fluent report name; values include projected total `force` in N
and `coefficient` when finite positive density, speed and area references exist.
Directions are used as saved, without changing their sign or magnitude.

Cancellation is checked between **outer iterations**, including inside a monitor
interval. The final partial interval emits an iteration event when at least one
iteration ran. A request during setup completes initialization before writing
a zero-iteration run. `done` and `summary.json` add `stop_reason: "stopped"` or
`"iteration_limit"`. `converged` remains false until the adapter exposes case
convergence criteria. JSONL mode watches stdin EOF; keep the sidecar's stdin
pipe open while it runs. Any received bytes are ignored. Human mode retains
normal redirected-stdin behavior and responds to SIGINT / CTRL_BREAK. Neither
signal handler runs solver or file I/O. SIGKILL/forced process termination cannot
produce final outputs.

## Settings inventory and setup overrides (v1 additive)

`inspect --json` also returns `settings`, an array with these exact fields:

```json
{"path":"boundary/velocity-inlet-5/velocity-magnitude","label":"Velocity magnitude","group":"Boundary Conditions","value":0.4,"unit":"m/s","kind":"number","options":[],"min":0,"max":1000000000000000,"editable":true,"reason":""}
```

Groups use Fluent's vocabulary: General, Models, Materials, Cell Zone Conditions,
Boundary Conditions, Reference Values, Methods, Controls, Report Definitions,
Monitors, Initialization, Run Calculation. `kind` is `number`, `integer`, `enum`,
`bool`, `vec3` or `string`. `options` is empty when unrestricted; `min`/`max` are
inclusive and null when unbounded. For vec3 they constrain each component.
Unavailable values are null. Read-only entries have a plain `reason`; fields
which would need solver changes say `not yet configurable in Kutta`.

Paths are stable lowercase slash-separated identities. Boundary and cell-zone
paths use the exact Fluent zone name, materials and reports use their names.
Lowercase ASCII letters, digits, `-`, `_`, `.` remain literal in each segment;
all other UTF-8 bytes, including uppercase letters, `/` and `%`, are percent
encoded with lowercase hex. This avoids collisions between distinct names
instead of guessing a slug. For example `Inlet A` is `%49nlet%20%41`.
The inventory includes retained raw settings under their source keys when no
canonical adapter control replaces them. Compound unparsed Scheme values are
read-only strings, preserving their original structure and meaning. Absent
required adapter controls also appear, with null values.

Only entries with `editable: true` may be changed. Initial editable controls are:
active inlet speed/components/direction, active inlet/outlet turbulence values,
outlet gauge pressure, active constant material density/viscosity, used reference
area/density/speed/pressure, implemented momentum/k/omega upwind schemes, gradient
limiting, pseudo-time scale and active equation relaxation, standard numeric
initialization values when not seeded, force-report directions, and run iteration
and monitor intervals. Changing the initialization method, operating pressure,
unused reference length/viscosity, model selection, wall motion, turbulence
specification or convergence criteria is read-only. A displayed Fluent value is
not evidence that Kutta can edit it.

```sh
kutta solve case.cas.h5 --setup setup.json --out run --progress jsonl
```

```json
{"v":1,"overrides":{"boundary/velocity-inlet-5/velocity-magnitude":0.8,"run/iterations":400,"run/monitor-interval":10}}
```

Setup files are strict JSON objects, at most 1 MiB. Version must be 1, the
`overrides` object is required, and duplicate keys and unknown top-level
properties are rejected. Every override is checked against the same inventory
for identity, editability, kind, options and range. Directions must also be
nonzero. All entries are validated before any are applied; invalid entries yield
exit 1 and an `error` event, with each offending path in `issues[].path` (also
`setting` for the existing support-issue shape). No solver initialization starts.
The ordinary strict case/setup validation still runs on the resulting values;
no override can opt into an unsupported physics model or approximation.

Overrides apply to in-memory configuration, preserving the original case file
and its SHA256. They apply after CLI defaults, so `run/iterations` in a setup file
wins over `--iterations`. `setup` and `summary.json` echo the applied object as
`overrides` and add `modified: true` when it is nonempty. Modified VTU/CSV headers
and view provenance carry `MODIFIED RESULT`, separately from approximation state.
`fluxes.csv` gives each boundary zone's integrated outward mass flow in kg/s,
using the solver's final mass flux, and carries the same provenance notice.
