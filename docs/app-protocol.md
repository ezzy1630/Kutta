# Desktop app ↔ solver protocol (v1)

The native desktop app (`app/`, Qt 6 + C++) never links the solver. It runs the
`kutta` CLI as a child process and talks to it through files and line-delimited
JSON on stdout. The CLI stays the single source of truth for what a case means;
the app only presents it. Every message carries `"v": 1`.

## Commands

| Command | Output |
| --- | --- |
| `kutta devices --json` | one JSON object: `{"v":1,"devices":[{"backend":"cuda","name":"NVIDIA GeForce RTX 4070","memory_bytes":12884901888,"default":true}, ...]}` |
| `kutta inspect <case> [--data <dat>] --json` | one JSON object: case summary (below) |
| `kutta solve <case> [--data <dat>] [--approximate a,b] [--backend b] [--iterations n] [--live-view DIR] [--live-every N] --out DIR --progress jsonl` | JSON lines on stdout while running; files in DIR |
| `kutta export-view DIR [--slice origin,normal]... [--post spec.json]` | writes `DIR/view/` (below) from a finished run |
| `kutta export-geometry CASE --out DIR [--json]` | writes `DIR/geometry/` (case surfaces without a run); with `--json` one `geometry` event line |
| `kutta aero DIR [--aero spec.json] [--json]` | aero coefficients, groups, balance and centre of pressure of a finished run from its `loads.json`; writes `DIR/aero.json`, with `--json` one `aero` event line |

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
Approximable issues also carry `measured_effect` (v1 additive): the measured
effect alone, as `summary.json` records it per approximation. Their
`explanation` still ends with it, so older readers lose nothing. They also
carry `short` (v1 additive), a few words naming what the run does differently,
for result badges: `wfgc` "WFGC off", `compressibility` "No compressibility
correction", `ke-options` "Documented k-epsilon only", `reverse-flow` "No
reverse-flow suppression". Each entry of `summary.json`'s `approximations`
carries the same `short`.

## Progress events (`solve --progress jsonl`)

One JSON object per line, in order:

```json
{"v":1,"event":"setup","cells":17822,"backend":"metal","device":"Apple M4 Pro","approximations":["wfgc"],"init":"hybrid","pressure_velocity_coupling":"coupled"}
{"v":1,"event":"iteration","iter":20,"residuals":{"continuity":1e-3,"x-velocity":1e-4,"y-velocity":1e-4,"z-velocity":1e-4,"k":1e-3,"omega":1e-3},
 "reports":{"drag":{"force":0.0123,"coefficient":0.31}},"sec_per_iter":0.042,"elapsed":0.9}
{"v":1,"event":"warning","message":"…"}
{"v":1,"event":"snapshot","run_id":"run-123-0","sequence":1,"iter":20,"elapsed":0.9,"manifest":"/absolute/DIR/live/run-123-0/snapshot-1/manifest.json","snapshot_seconds":0.006}
{"v":1,"event":"done","converged":false,"iterations":600,"elapsed":16.9,"out":"DIR"}
{"v":1,"event":"error","message":"…","issues":[...]}
```

`iteration` events are emitted at the monitor interval (the monitor samples)
and, v1 additive, at iterations 1 to 5. The solver monitors those anyway,
because Fluent scales continuity by its largest value in the first five
iterations, so reporting them costs no solver work and changes no output; the
residual history starts at iteration 1, as in Fluent. An iteration 1-5 event
that is not also a monitor sample has `residuals`, `sec_per_iter` and
`elapsed`, an empty `reports` object and no `convergence`: it is not a
convergence sample. Residual names follow
the solved equations: `continuity`, `x-velocity`, `y-velocity` and, except for
planar 2D cases, `z-velocity`, then `k`
and `omega` for SST or `k` and `epsilon` for standard k-epsilon; laminar runs
have none of these turbulence entries. SST with Fluent's Intermittency
Transition Model adds `intermit` (v1 additive, Fluent's name) after `omega`.
A planar 2D case (`general/dimension` 2) has no `z-velocity` residual or
`monitors/residual/z-velocity/*` controls, matching Fluent's 2D residual list;
Kutta's one-cell extrusion keeps w at round-off. Readers use the names
present. `residuals.csv` uses the same columns and rows.
Report values are
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

Each exported zone also has `view/zone-N.edges.bin`, the zone's polygon edges
(no triangulation diagonals) as KUTL polylines over the zone's own vertex
array, described by the zone's `edges` object (`file`, `primitive:"polylines"`,
`vertices`, `indices`, `line_count`, `fields:[]`). Live snapshots omit it.

## Case geometry (`export-geometry`, v1 additive)

```sh
kutta export-geometry case.cas.h5 --out DIR --json
```

```json
{"v":1,"event":"geometry","manifest":"/absolute/DIR/geometry/manifest.json","zones":6,"triangles":16978,"seconds":0.15}
```

Reads only the case mesh (no data, no volume snapshot, no solver) and publishes
`DIR/geometry/` by staging and renaming, like `view/`. Its manifest is a normal
view manifest with `geometry_only:true`, `reference:false`, empty `provenance`,
`bounds`, every boundary zone (interior zones excluded) with its surface and
edge files and empty `fields`, `slices:[]`, `fields:[]` and
`case:{"sha256":"…"}`. `manifest` in the event is an absolute UTF-8 path.
Errors print one `error` event (with `--json`) and exit 1; usage errors exit 2
and leave any previous `geometry/` untouched. On the mixing elbow (M4 Pro) it
takes 0.15 s and 62 MB peak memory; the case hash runs concurrently with the
mesh read.

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
Informational issues do not block solve; the only approximable settings are the
explicit setup allowlist (`wfgc`, `compressibility` for the k-epsilon
dilatation term, `ke-options` for undocumented k-epsilon solver switches,
`reverse-flow` for pressure outlets of Fluent 18.0 and earlier cases, which do
not record `prevent-reverse-flow?`). Each
approximable issue carries its measured effect and evidence. An unparsed WFGC
flag is unsupported, not approximable. A downforce report is normally a lift report whose direction
points down; its sign and reference scaling are preserved.

`setup` and `summary.json` carry `pressure_velocity_coupling` (v1 additive):
`"coupled"`, `"simplec"` or `"simple"`, the algorithm that actually ran. Each
of Fluent's SIMPLE, SIMPLEC and Coupled schemes runs natively; PISO is rejected.
`methods/pressure-velocity-coupling` switches between the three; Coupled needs
the pseudo-transient formulation, and a saved Coupled case without it runs
SIMPLEC with a convergence-path note. Its `options` therefore omit `coupled`
for a case without pseudo time unless the case was saved with Coupled.

`summary.json` carries `poor_quality_cells` (v1 additive): the number of cells
at or below orthogonal quality 0.05 whose second-order face pressures use the
cell value, following Fluent's poor mesh numerics (User's Guide 36.22.5).
Zero on a well-conditioned mesh and with the standard pressure scheme.
`poor_cell_quality` (v1 additive) is that threshold as the run used it (0.05;
negative when the treatment is off). The app names the count beside the run's
convergence and in Run History only when it is not zero.

The Methods group (v1 additive) follows Fluent's Methods panel. Every case
lists `methods/momentum-scheme`, `k-scheme`, `omega-scheme` and
`epsilon-scheme`; a scheme is `active` only when the case's viscous model
solves that equation (k and omega for k-omega models, k and epsilon for
k-epsilon models, momentum only for laminar), and only active saved schemes
are editable. Read-only entries describe the rest of the panel:
`methods/gradient` (`least-squares-cell-based`, `green-gauss-cell-based` or
`green-gauss-node-based`, labelled "Least Squares Cell Based", "Green-Gauss
Cell Based", "Green-Gauss Node Based"), `methods/pressure-scheme` (`standard`
or `second-order`; other saved schemes keep their Fluent name and are
unsupported) and `methods/pseudo-transient` (bool, whether the steady case
uses the pseudo-transient formulation).

`controls/courant-number` (v1 additive, read-only) is Fluent's Flow Courant
Number as saved (`pressure-coupled/courant-number`, null when absent). Kutta's
Coupled solver always applies Fluent's default, 200, as implicit relaxation on
top of the pseudo time step (Theory Guide 18.4.4); its `reason` says so. With
pseudo time it depends on `methods/pressure-velocity-coupling`
(`active_when` `["coupled"]`); without pseudo time Kutta never runs Coupled,
so it is inactive.

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
a zero-iteration run. `done` and `summary.json` add `stop_reason: "stopped"`, `"iteration_limit"`
or `"converged"`; `converged` is true only for the last, when convergence
checking is on (see Convergence criteria below). JSONL mode watches stdin EOF only when stdin is a pipe
at startup; the app keeps the child process's stdin pipe open while it runs. Any received
bytes are ignored. A file, `/dev/null`, a terminal or absent stdin never stops a run. Human mode retains
normal redirected-stdin behavior and responds to SIGINT / CTRL_BREAK. Neither
signal handler runs solver or file I/O. SIGKILL/forced process termination cannot
produce final outputs.

## Settings inventory and setup overrides (v1 additive)

`inspect --json` also returns `settings`, an array with these exact fields:

```json
{"path":"boundary/velocity-inlet-5/velocity-magnitude","label":"Velocity magnitude","group":"Boundary Conditions","value":0.4,"unit":"m/s","kind":"number","options":[],"option_labels":[],"value_label":"","min":0,"max":1000000000000000,"editable":true,"reason":"","active":true,"depends_on":"boundary/velocity-inlet-5/velocity-specification","active_when":[0,2]}
```

Groups use Fluent's vocabulary: General, Models, Materials, Cell Zone Conditions,
Boundary Conditions, Reference Values, Methods, Controls, Report Definitions,
Monitors, Initialization, Run Calculation. `kind` is `number`, `integer`, `enum`,
`bool`, `vec3` or `string`. `options` is empty when unrestricted; `min`/`max` are
inclusive and null when unbounded. For vec3 they constrain each component.
Unavailable values are null. Read-only entries have a plain `reason`; fields
which would need solver changes say `not yet configurable in Kutta`.

`option_labels` is parallel to `options` (`[]` when there are no options).
`value_label` gives Fluent's words for the current `value` of an enum or code
when known, else `""`; zone, material and report `type` labels are the Fluent
type string. Integer codes stay integers: `turbulence-specification` 0 "K and
Omega" ("K and Epsilon" for k-epsilon models), 1 "Intensity and Length Scale",
2 "Intensity and Viscosity Ratio", 3 "Intensity and Hydraulic Diameter";
`velocity-specification` (velocity inlets) 0 "Magnitude and Direction",
1 "Components", 2 "Magnitude, Normal to Boundary". Retained raw settings have
no labels.

`active` says whether the case as saved uses the setting, as Fluent's panel
would show it. When `depends_on` names a controlling enum path, `active` is
exactly whether that setting's current `value` is in `active_when`; otherwise
`depends_on` is null, `active_when` is `[]` and `active` is static. Turbulence
values depend on the zone's `turbulence-specification`: intensity [1,2,3],
length scale [1], viscosity ratio [2], hydraulic diameter [3], k [0], omega
[0] for k-omega models, epsilon [0] for k-epsilon models. Inlet magnitude
depends on `velocity-specification` [0,2], components [1], direction [0].
Velocity inputs on non-inlets, gauge pressure on non-outlets and the turbulence
specification of zones without one (walls, symmetry) are inactive.

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
active inlet speed/components/direction, inlet/outlet turbulence values and
turbulence specification (k-omega models), outlet gauge pressure, active
constant material density/viscosity, reference area/density/speed/pressure and
length (length scales moment coefficients in `kutta aero`),
implemented momentum/k/omega/epsilon upwind schemes, gradient limiting, the
pressure interpolation (`methods/pressure-scheme`, options `standard` and
`second-order`), the pressure-velocity coupling (`methods/pressure-velocity-coupling`, below), pseudo-time
scale and active equation relaxation, standard numeric initialization values
when not seeded, force-report directions, and run iteration and monitor
intervals, and the convergence controls under `monitors/` (below). Changing the
initialization method, operating pressure, reference viscosity, model
selection, velocity specification or wall motion is read-only. A displayed Fluent value is not
evidence that Kutta can edit it.

The turbulence specification's `options` are restricted to the methods whose
required values were saved as constants or resolved expressions (plus the
current method), so switching never runs on a guessed value. Turbulence values
a method of the case's model can use stay editable while inactive, so the app
can switch method and set its values in one setup file; an override of an
inactive value is accepted and has no effect, as Fluent keeps such values.
Values no method of the model reads (epsilon for k-omega) are read-only. A
laminar case reads no turbulence inputs: its saved turbulence specification
and values, and its k and dissipation relaxation entries, are inactive and
read-only.

Transition models (v1 additive). `models/transition` is a read-only enum on
every case: `none`, `gamma-transport-eqn` (Fluent's Intermittency Transition
Model), `gamma-algebraic` or `transition-sst`. Only `none` and
`gamma-transport-eqn` with SST run; the others are unsupported issues at
`case-config/rp-transition-model` with a plain reason. With the
intermittency model, `inspect --json` adds `"transition":"gamma-transport-eqn"`
to `models` (absent otherwise), and the settings add
`methods/intermittency-scheme` (first or second order upwind),
`controls/intermit-relaxation`, and `boundary/<zone>/intermittency` (0 to 1,
the inflow value) on inlets and outlets. `result.vtu` and `export-view` add a
cell field `intermittency` (dimensionless gamma). The model's constants are
read from the case and checked; see `src/kutta/fv/turbulence/transition.md`.

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
`overrides`, and its run controls (`run/iterations`, `run/monitor-interval`)
again as `run_controls`. Switching `methods/pressure-velocity-coupling` changes
the algorithm only: the case's saved scheme still selects which pseudo-time
controls and relaxation keys apply (`controls/mom-relaxation` stays the saved
scheme's value), and a warning reports the convergence-path change. Both
algorithms solve the same discrete equations. `modified` is true only when a case-changing override
was applied; run controls alone never mark a result modified. Modified VTU/CSV
headers and view provenance carry `MODIFIED RESULT` with the case-changing
overrides only, separately from approximation state. The monitor interval is
still part of the run's recipe: the solver's monitor refresh currently changes
the unconverged iterate path, so runs compared at equal iteration counts should
use equal intervals.
`fluxes.csv` gives each boundary zone's integrated outward mass flow in kg/s,
using the solver's final mass flux, and carries the same provenance notice.

## Volume post-processing (v1 additive)

```sh
kutta export-view DIR --post post.json [--slice ox,oy,oz,nx,ny,nz]...
kutta export-view DIR --streamlines seeds.json
```

`--post` batches streamlines/pathlines, iso-surfaces, XY probes and surface
reports into the same atomic `view/` publication as ordinary surfaces/slices.
`--streamlines` accepts either a complete batch or one trace object with `v:1`.
Only one specification may be passed. The CLI reads the volume snapshot; the app
reads the manifest, compact geometry, probe CSVs and reports. No volume mesh is
sent to the app. Invalid request syntax is exit 2; unavailable fields, mismatched
frames or run-data failures are exit 1; an operation's failure message starts
with `post operation "<name>":`. A failed export preserves the old view.

```json
{"v":1,
 "fields":["pressure_coefficient","total_pressure_coefficient","vorticity","q_criterion","helicity","velocity_x","velocity_y","velocity_z","skin_friction"],
 "streamlines":[{
   "name":"inlet-rake","mode":"streamline",
   "seeds":{"type":"rake","start":[0,0,0],"end":[0,1,0],"count":20},
   "direction":"forward","max_steps":10000,"max_length":10,"max_time":100,
   "step_length":0.01,"relative_tolerance":0.00001,"absolute_tolerance":0.0000001,
   "fields":["velocity_magnitude","pressure","total_pressure_coefficient"]}],
 "isosurfaces":[{"name":"wake","field":"total_pressure_coefficient","value":0,
                  "fields":["pressure","velocity_magnitude"]},
                 {"name":"vortices","field":"q_criterion","value":100}],
 "probes":[{"name":"wake-line","type":"line","start":[0,0,0],"end":[1,0,0],"count":100,
            "fields":["pressure_coefficient","velocity_x"]},
           {"name":"wing-section","type":"section","zones":["wing"],
            "origin":[0,0,0.5],"normal":[0,0,1],"fields":["pressure_coefficient"]}],
 "integrals":[{"name":"aero","zones":["wing"],"fields":["pressure_coefficient","skin_friction"],
               "origin":[0,0,0],"drag":[1,0,0],"lift":[0,0,1],"side":[0,1,0],
               "pressure_plane":{"origin":[0,0,0],"normal":[1,0,0]}}]}
```

Names identify operations; filenames use ordinals, so arbitrary names cannot
escape `view/`. Names must be unique across a batch. Specifications are strict
JSON, at most 1 MiB and 256 operations. Unknown keys, duplicate keys/names,
nonfinite values, degenerate planes/directions, and invalid ranges are rejected.
Counts and maximum steps are positive integers no larger than 1,000,000.

The global `fields` list adds selected fields to the ordinary surfaces and
slices. When omitted/empty, all available derived fields below are added. The
manifest's `available_fields` describes every sampleable primitive, derived and
Fluent/delta field, with unit, components and `volume` or `boundary` location.
Boundary-only fields cannot be used on volume iso-surfaces or volume line probes.
Missing source/reference data is reported rather than invented. Raw saved fields
are also selectable by their exact exported name.

Top-level `"fields_policy"` is `"strict"` (default: an unavailable name fails
the export) or `"available"`: unavailable names in the global `fields` list and
in every operation's `fields` list are skipped, as are boundary-only names on
traces, iso-surfaces and line probes. An iso-surface's own `field` stays strict.
The manifest's `unavailable_fields` lists each skipped name with its operation
name, null for the global list, in batch order:
`[{"name":"skin_friction","operation":null},{"name":"x","operation":"wake"}]`.
It is present, possibly empty, whenever a post request is given.

### Reconstruction and tracing

Post-processing uses unlimited inverse-square least-squares cell reconstruction
with face-neighbour cell values and available boundary values. A rank-aware
pseudoinverse preserves resolved directions on one-layer meshes. This is a
post-processing reconstruction; it does not alter the solver's transport limiter.
Velocity derivatives use this same reconstruction. Initial point location uses a
cell BVH; subsequent integration stages walk across the actual owner/neighbour
faces. Face intersections use their triangulated 3D geometry. Boundaries terminate
traces; no periodic wrapping, wall reflection or particle-force model is implied.

Trace seed variants are:

```json
{"type":"rake","start":[0,0,0],"end":[0,1,0],"count":20}
{"type":"plane","origin":[0,0,0],"u":[0,1,0],"v":[0,0,1],"nu":10,"nv":10}
{"type":"surface","zones":["inlet"],"count":100,"distribution":"uniform","offset":0}
{"type":"surface","zones":["inlet"],"count":100,"distribution":"random","random_seed":42,"offset":0.001}
```

Rakes include endpoints; a count of one chooses the midpoint. Plane `u` and `v`
span the entire grid; an axis with count one chooses its midpoint. The total
plane-grid count cannot exceed 1,000,000. Surface sampling is weighted by triangle
area. `uniform` is deterministic stratified area selection with low-discrepancy
barycentric samples; `random` uses a reproducible seed. `offset` is a nonnegative
inward distance in metres. Seeds outside the volume are retained as one-point
traces with `stop_reason:"outside"`. A zero-speed seed terminates as `stagnation`.

Integration uses RK4 step doubling, with adaptive step control and cell-to-cell
walking at each stage. A zero/omitted `step_length` chooses 0.4 times the local
cube-root cell volume; otherwise it is the maximum proposed travel distance in
metres. Tolerances bound the estimated spatial integration error. Defaults are
10,000 steps, 1,000 m maximum length, 1,000,000 s maximum travel time, relative
error 1e-5 and an absolute tolerance scaled to mesh size. `direction` is `forward`,
`backward` or `both`; both emits two independently identified curves per seed.
Each vertex includes requested fields, `arc_length` (m) and physical `time` (s).
Line metadata gives seed ID, direction, accepted steps, length, time and stop
reason: `boundary`, `outside`, `stagnation`, `max_steps`, `max_length`, `max_time`,
`time_range`, `missing_velocity` or `step_underflow`. Boundary endpoints are snapped
within the integration tolerance; failed integration is never extended outside.

For a steady run, `mode:"pathline"` and `mode:"streamline"` follow the same
velocity field. Time-dependent pathlines additionally accept:

```json
{"name":"particles","mode":"pathline","start_time":0,
 "frames":[{"time":0,"run":"run-000"},{"time":0.1,"run":"run-001"}],
 "seeds":{"type":"rake","start":[0,0,0],"end":[0,1,0],"count":20}}
```

Frame paths are relative to the specification file. At least two strictly
increasing frame times must bracket `start_time`. Every frame must have exactly
the same mesh coordinates/connectivity and cell order; counts alone do not prove
identity. Reconstructed fields are interpolated linearly in time, and integration
stops at the saved time range. Only two frames are resident. These are massless
advected pathlines, not inertial discrete particles.

### Iso-surfaces and XY probes

Iso-surfaces use marching tetrahedra, decomposing each cell into its centre and
triangulated faces. Shared mesh-node values are inverse-distance-weighted means
of neighbouring cell reconstructions, making face edges consistent between cells.
Concave faces are ear-clipped in their dominant projection. For warped faces,
the exporter selects a shared diagonalization that favours positive tetrahedra
for both neighbouring cells. If necessary, it finds an interior point in the
triangulated cell kernel instead of using the Fluent face-pyramid centroid.
A cell without a positive kernel fails explicitly; inverted or overlapping
tetrahedra are never silently accepted. Vertices
are welded by topological edge identity, duplicate triangles are removed, and
winding points toward increasing iso value. A scalar `field` uses its value;
a vector uses magnitude by default, or `component:0|1|2` for x/y/z. Requested
vertex fields use the same tetrahedral interpolation. A constant-valued volume
at the iso value has no unique two-dimensional level surface and emits no faces.
No crossing also produces an empty geometry, with null field ranges.

Line probes sample exactly `count` positions along the line, including endpoints.
They reconstruct fields within each located volume cell. Outside points retain
their positions, `valid:0`, and NaN values; valid segments remain separate curves.
Section probes intersect the selected boundary zones with a plane and interpolate
the zone's face-area-weighted vertex values onto the intersection. Segments are
connected into ordered curves; closed curves include their closing point and
are oriented with the plane normal. Multiple zones retain separate curves and
zone identity. Coplanar triangles contribute their outer perimeter. Each curve
has an `arc_length` axis. CSV columns include sample, x/y/z, arc length, valid and
requested fields (vectors use `.x`, `.y`, `.z` suffixes, avoiding collisions
with scalar `velocity_x` fields); missing values are empty CSV cells. This
supports Cp against chord coordinate or distance around a wing section.

### Derived fields and references

Reference density, velocity, pressure, area and length, and actual material
density/viscosity are saved after `solve --setup` overrides are applied. No case
file needs to be reopened. Private snapshot v2 adds these values and the solver's
boundary mass fluxes; the reader still accepts v1 snapshots. For old snapshots,
`reference` in a post request may supply missing values. It may also change the
normalizing reference values for an existing run. It cannot change a known
solved density or viscosity; use `solve --setup` for those physical changes.
The view records the values used and post-reference overrides in provenance.

Let `q_ref = 0.5 rho_ref U_ref^2`, `G_ij = dU_i/dx_j`, `S=(G+G^T)/2` and
`Omega=(G-G^T)/2`. Additional fields are:

| Name | Unit | Definition |
| --- | --- | --- |
| `pressure_coefficient` | 1 | `(p - p_ref) / q_ref` |
| `total_pressure` | Pa | incompressible gauge `p + 0.5 rho |U|^2` |
| `total_pressure_coefficient` | 1 | `(total_pressure - p_ref) / q_ref`; static reference pressure |
| `velocity_x`, `velocity_y`, `velocity_z` | m/s | reconstructed velocity components |
| `vorticity` | 1/s | `curl(U)`, xyz vector |
| `vorticity_magnitude` | 1/s | `|curl(U)|` |
| `q_criterion` | 1/s^2 | `0.5 (||Omega||_F^2 - ||S||_F^2)` |
| `helicity` | m/s^2 | `U dot curl(U)` |
| `wall_shear`, `wall_shear_magnitude` | Pa | saved wall stress on fluid from the active closure, vector/magnitude |
| `y_plus` | 1 | saved SST u-star wall diagnostic |
| `y_star` | 1 | saved k-epsilon wall-function y* diagnostic |
| `skin_friction`, `skin_friction_vector` | 1 | wall shear magnitude/vector divided by `q_ref`, on fluid |

The existing `@fluent` and `@delta` conventions apply. Derivatives are reconstructed
from each solution independently; nonlinear delta fields are Kutta's derived
quantity minus Fluent's derived quantity. No derivative of a difference is
substituted for a difference of nonlinear diagnostics.

### Surface reports

Each integral request emits each exact selected zone and a combined total. Empty
`zones` selects all boundary zones. Outputs include area, signed outward mass
flow, absolute mass flow, field area integrals, area-weighted and absolute-mass-
flow-weighted averages, coverage, pressure/viscous/total force and moment, split
lift/drag/side coefficients, vector moment coefficients and centre of pressure.
Directions are normalized; defaults are drag x, lift y, side z. The actual axes
and moment origin are echoed. Undefined zero-flow averages and coefficients with
missing references are JSON null. NaN reference samples carry explicit coverage.

New snapshots use the solver's saved mass flux. Old snapshots may instead use
`rho U_boundary dot S_face`, and reports identify that source. Pressure force on
the body is `(p-p_ref) S_face`; viscous wall load is minus saved stress-on-fluid
times area. Non-wall viscous traction uses the reconstructed symmetric velocity
gradient and effective dynamic viscosity. Moments use the face centroid relative
to `origin`. Force coefficients divide by `q_ref A_ref`, moment coefficients by
`q_ref A_ref L_ref`. These are discrete face-field integrals, consistent with the
saved solution, rather than a re-solve of boundary layers.

Centre of pressure uses pressure force/moment only. It is the closest point of
the pressure-resultant line to the moment origin unless `pressure_plane` is
specified, in which case it is that line's intersection with the plane. The
report includes the residual pressure couple and `zero_moment`; a general 3D
load can retain a couple parallel to its force. Zero pressure force or a
nonintersecting parallel plane produces a null location and an explicit status.
These choices follow [Fluent's force/moment conventions](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/flu_th_sec_report_force_moment.html)
and [absolute mass-flow weighting](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/flu_th_sec_compute_surfint.html).

### App artifacts

The manifest adds `streamlines`, `isosurfaces`, `probes`, `reports`,
`available_fields` and `reference_values`. Every geometry descriptor has its
operation `name`, primitive, file, vertex/index counts, field offsets/units/ranges
and operation-specific metadata. Iso-surfaces use KUTV triangles. Traces and
probes use KUTL polylines in the same little-endian float32/uint32 family; its
layout is documented [next to the writer](../src/kutta/io/VIEW_FORMAT.md).
Probe descriptors link a CSV, and `reports_file` links `view/reports.json`, whose
root is `{"v":1,"reports":[...]}`. All outputs carry the run/view provenance.

### Incremental post-processing

`export-view --post` reuses unchanged work from `DIR/post-cache/`. Units are the
base zone surfaces (with the global derived fields), each slice and each
operation. A unit's key is a SHA-256 of a cache format version, the CLI version
and commit, the run's snapshot identity (absolute run path, `view-input.bin` size and mtime plus
the `summary.json` bytes; the snapshot itself is not hashed), the post
`reference` overrides, `fields_policy` and the unit itself: the operation's
canonical parsed JSON with its effective name (a pathline also includes each
frame run's identity), a slice's exact plane, or the global `fields` list. A hit
is hard-linked (copied if linking fails) into the new view under the ordinal
filename this batch assigns, so inserting or reordering operations still reuses
them. When every unit hits, the snapshot is not even read. Otherwise the
snapshot is read and one volume reconstruction serves every miss; derived
surface/slice fields are computed only when the base or a slice missed. A warm
publication is identical to a cold export of the same batch except
`post_cache`. After a successful publication, entries it did not use are
pruned; a failed export keeps the old view and the cache.

The manifest adds `"post_cache":{"reused":R,"computed":C,"seconds":S}`: units
reused and computed (base, slices and operations) and the export's wall time
up to publication. Mixing elbow (17,822 cells, M4 Pro), 30 surface-seeded trace
operations, a pressure iso-surface, line and section probes and an all-zone
report: cold 2.8–2.9 s (35 computed); adding a Q-criterion iso-surface 0.23 s
warm (1 computed, 35 reused) versus 3.0 s cold; the unchanged batch 0.02 s.


## Live surfaces during solve (v1 additive)

```sh
kutta solve case.cas.h5 --out RUN --progress jsonl --live-view RUN --live-every 20
```

`--live-view DIR` enables surface updates under `DIR/live/`. It works in human
mode too; JSONL mode additionally emits `snapshot` events. `--live-every N` is a
positive integer **minimum iteration gap**, default **20**, and requires
`--live-view`. Publication uses the solver's existing refreshed monitor states
(the first five iterations, scheduled monitor iterations, and the terminal
iteration). Between those states, boundary pressure is a pressure-correction
scratch buffer, so exporting it would show incorrect contours. No extra monitor
work is introduced for live view. The first snapshot follows at least N completed
outer iterations, at the next refreshed state; the 20-iteration default matches
the ordinary monitor interval. A shorter/stopped run may produce no live snapshots.
For more frequent refreshed states, use the supported `run/monitor-interval`
setup control. There is no forced final live write. Normal final run outputs
and the `done` completion contract stay the same.

Snapshots contain the ordinary export-view boundary surfaces: pressure, Cp when
references exist, wall stress/y+ and non-wall velocity, with available Fluent and
delta fields, units and ranges. They use the identical KUTV v1 layout and surface
interpolation. Each is a coherent completed-iteration state, carrying the run's
approximation/override provenance. Capture borrows the solver mesh and downloads
velocity/boundary fields, without copying mesh topology, writing a volume snapshot
or computing volume post-processing. Extra derived fields, slices, probes and
traces remain available through `export-view` on the completed run.

`DIR/live/manifest.json` is a small, atomically replaced pointer:

```json
{"v":1,"run_id":"run-123-0","sequence":2,"iter":40,"elapsed":1.7,
 "manifest":"run-123-0/snapshot-2/manifest.json"}
```

Resolve its `manifest` **relative to `DIR/live/`**. That file is a normal
export-view manifest; resolve its zone filenames relative to its own directory.
Each publication first closes all files in a temporary snapshot directory, then
renames the directory to `run_id/snapshot-SEQUENCE/`. Finally it closes a temporary
pointer file and renames it over `live/manifest.json`. A reader sees the previous
complete snapshot or the new complete snapshot, with no missing-pointer interval.
Only after the pointer swap does the solver emit the `snapshot` event, whose
`manifest` is an **absolute UTF-8 path** and whose `snapshot_seconds` measures
capture and publication. `sequence` starts at 1 and increases without gaps for
each `run_id`; a new solve gets a new run ID, even when it reuses DIR. The setup
event includes that ID in `live_view`, so the app can discard a previous run's
contours before the new run publishes its first snapshot.

Published snapshots are immutable and retained. Readers may finish reading an
older sequence while the solver publishes another. A publication failure leaves
the previous pointer/assets readable, removes unpublished staging data and emits a
`warning`; the solve continues, since live view is a preview. One solver owns each
DIR. The caller owns cleanup of retained run directories after readers close them.
On Windows the pointer swap uses replacement rename, with a bounded retry for
readers briefly holding the pointer file open; consumers should close it after
reading, and retry when it cannot be opened for a moment during a swap. The
pointer is never written in place, so a failed or empty read is always transient. Atomicity
means visibility to concurrent readers on the same filesystem, not a guarantee
against storage failure or power loss.

Snapshot work has a **5% solve-time budget**. The first write calibrates its cost.
Later writes wait until a refreshed state, the iteration gap and the budget
permit publication:
`total_snapshot_seconds + 1.25 * previous_snapshot_seconds <= 0.05 * solve_seconds`,
where solve seconds are loop wall time minus snapshot time. This reserves a 25%
cost margin; unexpectedly slow storage or the first calibration can exceed the
target temporarily. Delayed updates skip intermediate states and use the current
iteration, rather than queueing stale snapshots. Snapshot cost includes field
downloads, interpolation, binary/manifest I/O and publication. There is no export
work when live view is disabled.

The `setup` event adds optional `live_view` with `directory` (absolute live root),
`run_id`, `every`, `unit:"iterations"` and `budget_fraction:0.05`. `done` and
`summary.json` echo those fields and add `snapshots`, `snapshot_time_seconds`,
`solve_time_seconds` and `overhead_fraction` (snapshot/solve time; null for zero
solve time). Ordinary iteration `sec_per_iter` excludes measured snapshot time;
progress `elapsed` and run `wall_time_seconds` include it. Cancellation finishes
any active snapshot publication and the current iteration before normal output
closure. Live views are intermediate solutions, not convergence certificates.

## Convergence criteria (v1 additive)

Fluent stops a steady run when its residual criteria and report convergence
conditions are met. `kutta solve` does the same when convergence checking is
on. It is **off by default**, so an existing script still runs to its
iteration limit; the app's Monitors page turns it on. Every control below is a
settings-inventory entry in group `Monitors`, editable through `solve --setup`.
They are run controls like `run/iterations`: they say when to stop, not what
case is solved, so they never mark a result MODIFIED and are echoed under
`run_controls`. (`setup::run_control` covers the `run/` and `monitors/`
prefixes.)

| Path | Kind | Default | Meaning |
| --- | --- | --- | --- |
| `monitors/check-convergence` | bool | `false` | Stop when the conditions below are met |
| `monitors/condition` | enum `"all"`, `"any"` | `"all"` | "All conditions are met" / "Any condition is met" |
| `monitors/residual/<equation>/check` | bool | case value | Fluent's Check Convergence for that residual |
| `monitors/residual/<equation>/absolute-criterion` | number ≥ 0 | case value | Fluent's Absolute Criteria (scaled residual) |
| `monitors/report/<report>/active` | bool | `false` | Use this report definition as a convergence condition |
| `monitors/report/<report>/stop-criterion` | number > 0 | `0.001` | Fluent's Stop Criterion (relative change) |
| `monitors/report/<report>/previous-values` | integer ≥ 1 | `5` | Np: report values to include (Fluent's Use Iterations, counted in monitor samples) |
| `monitors/report/<report>/ignore-before` | integer ≥ 0 | `100` | Fluent's Ignore Iterations Before |

`<equation>` is each residual name the solver reports (`continuity`,
`x-velocity`, `y-velocity`, `z-velocity` except for planar 2D cases, then `k`
and `omega` or `epsilon` for turbulent cases, and `intermit` with the
intermittency transition model). Residual `check` and
`absolute-criterion` come from the case's `residuals/settings` when its rows use
the 7-element layout Kutta parses, otherwise `check` is `false` and the
criterion is null. A raw row that also lists equations Kutta does not solve
(energy, for example) stays read-only. A checked residual with a null criterion
is exit 1 with an issue at `monitors/residual/<equation>/absolute-criterion`,
whether or not checking is on. Report entries exist for each force-type
report definition (`force`, `drag`, `lift`, `downforce`); other report kinds are
listed read-only. Report defaults are Kutta's, not read from the case.

Evaluation happens at every monitor sample (each `iteration` event at the
monitor interval, including the final partial interval; not the extra
iteration 1-5 events):

- **Residuals** are one condition: met when at least one equation is checked
  and every checked equation's scaled residual is at or below its criterion.
- **A report** with value history m (the projected total `force` in N at each
  monitor sample, in order) is met at sample n when n's iteration is greater
  than `ignore-before`, at least Np+1 samples after `ignore-before` exist, and
  `max over k=1..Np of |m(n-k+1) - m(n-k)| / |m(n)|` is **below** the stop
  criterion (Fluent's lag residuals; [Fluent 24.2 User's Guide §36.16](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_convergence_conditions.html)).
  m(n) = 0 is never met, and its value is null. Coefficients differ from forces by a constant factor,
  so the relative change is the same. Kutta evaluates reports at the monitor
  interval, so consecutive values are `run/monitor-interval` iterations apart;
  with interval 1 Np is exactly Fluent's Use Iterations, with a larger interval
  the check is stricter than Fluent's for the same Np.
- **Combination**: the conditions are the residual condition (when any residual
  is checked) and each active report. `all` needs every condition met at the
  same sample, `any` needs one. Checking on with no condition emits the warning
  "Convergence checking is on but no residual is checked and no report
  condition is active; stopping at the requested iteration limit." and runs to
  the iteration limit.

When met, the run stops after that iteration and writes its outputs normally.
`done` and `summary.json` then carry `"converged": true` and
`"stop_reason": "converged"`; this wins over a user stop or the iteration
limit reached at the same sample. Otherwise `converged` stays false. Human
output (no `--json`/`--progress`) evaluates the same conditions at the same
samples.

The `setup` event adds `convergence` (always present):

```json
"convergence":{"enabled":true,"condition":"all",
  "residuals":{"continuity":0.001,"x-velocity":0.001},
  "reports":[{"report":"drag","stop_criterion":0.001,"previous_values":5,"ignore_before":100}]}
```

`residuals` lists checked equations only. With checking on, each `iteration`
event adds the state of every condition at that sample, residuals first and
then active reports in report order:

```json
"convergence":{"met":false,"conditions":[
  {"kind":"residuals","met":false,"equations":{"continuity":{"value":0.0021,"criterion":0.001,"met":false}}},
  {"kind":"report","report":"drag","value":0.0042,"criterion":0.001,"met":false}]}
```

A report's `value` is null until it has enough samples. `done` and
`summary.json` add `convergence`: the setup object plus `met_by` (the met
conditions at the stopping sample, as `"residuals"` or `"report:<name>"`, empty
unless converged) and `iteration` (the stopping iteration, null unless
converged). `summary.json` sets `convergence_criteria_available: true`.
Without checking, the warning reads "Convergence checking is off
(monitors/check-convergence); stopping at the requested iteration limit."

## Aero loads and coefficients (v1 additive)

Every successful `solve` also writes `DIR/loads.json`: the force and moment on
each wall zone, integrated from the same final boundary pressure and wall
shear stress that `export-view` reports use, about the origin of the case
coordinates:

```json
{"v":1,"origin":[0,0,0],"units":{"force":"N","moment":"N m","area":"m^2"},
 "reference":{"area":1.0,"density":1.225,"velocity":16.0,"length":1.53,"pressure":0.0},
 "report_zones":["monocoque","front-wing"],
 "zones":[{"name":"front-wing","zone_id":12,"area":0.91,
           "pressure":[fx,fy,fz],"viscous":[fx,fy,fz],
           "pressure_moment":[mx,my,mz],"viscous_moment":[mx,my,mz]}],
 "provenance":"","approximate":false,"modified":false,"iterations":600}
```

Pressure force is `(p - p_ref) S_face` on the body; viscous load is minus the
stress on the fluid times face area, as in surface reports. `reference` holds
the values the run solved with (case values after setup overrides; null when
absent). `report_zones` is the union of the zones of the case's force report
definitions, in report order (empty when there are none). The per-zone totals
equal the wall rows of `forces.csv` to floating-point summation order. Solve
deletes an earlier `loads.json` and `aero.json` in the output directory when it
starts, so neither can describe a different solution.

`kutta aero DIR [--aero spec.json] [--json]` reads only `loads.json`, so it
takes milliseconds and never re-solves; changing groups, axes or vehicle
geometry is post-processing. It writes `DIR/aero.json` (staged and renamed) and
with `--json` prints `{"v":1,"event":"aero", ...the same object}`. A run
without `loads.json` (solved before this addition) is exit 1 with an `error`
event; invalid specs are exit 2.

The spec is strict JSON (at most 1 MiB, unknown or duplicate keys rejected);
every key is optional, and `v`, when given, must be 1:

```json
{"v":1,
 "axes":{"drag":[1,0,0],"lift":[0,0,1]},
 "moment_centre":[0,0,0],
 "car_zones":["monocoque","sidepods","front-wing","rear-wing","tyre-front","tyre-rear"],
 "groups":[{"name":"Front wing","zones":["front-wing"]},{"name":"Wheels","zones":["tyre-front","tyre-rear"]}],
 "vehicle":{"front_axle":[0,0,0],"wheelbase":1.53}}
```

- `axes.drag` is the free-stream direction (default +x), `axes.lift` points up
  (default +z). Both are normalized and must be nonzero and orthogonal
  (|cos| ≤ 1e-6). The side axis is `lift × drag` (+y for the defaults).
- `moment_centre` (m, default origin) is where moments and moment
  coefficients are taken.
- `car_zones` defaults to `report_zones`, or every wall zone when the case has
  no force reports, so the total matches the case's own force reports. Names
  must be wall zones in `loads.json`.
- `groups` have unique nonempty names; their zones must be car zones and each
  zone may be in at most one group. Car zones in no group are listed in
  `ungrouped`.
- `car_zones` and each group's `zones` are nonempty and list a zone once.
- `vehicle.front_axle` is the front tyre contact point (m) and `wheelbase` (m,
  > 0) the distance to the rear contact point along the drag axis; a vehicle
  needs both. Without a vehicle there is no balance.

`aero.json`:

```json
{"v":1,"provenance":"…","approximate":false,"modified":false,"iterations":600,
 "reference":{"area":1.0,"density":1.225,"velocity":16.0,"length":1.53,"dynamic_pressure":156.8},
 "axes":{"drag":[1,0,0],"lift":[0,0,1],"side":[0,1,0]},"moment_centre":[0,0,0],
 "rows":[{"name":"front-wing","kind":"zone","zones":["front-wing"],
          "force":{"pressure":{"drag":12.1,"lift":-180.2,"side":0.3},
                   "viscous":{"drag":1.9,"lift":-0.1,"side":0.0},
                   "total":{"drag":14.0,"lift":-180.3,"side":0.3}},
          "downforce":180.3,
          "moment":{"pitch":-41.0,"roll":0.2,"yaw":-0.1},
          "coefficients":{"cd":0.089,"cl":-1.150,"cs":0.002,
                          "cd_pressure":0.077,"cd_viscous":0.012,"cl_pressure":-1.149,"cl_viscous":-0.001,
                          "cla":-1.150,"cda":0.089,
                          "cm_pitch":-0.171,"cm_roll":0.001,"cm_yaw":0.0},
          "downforce_share":0.31}],
 "ungrouped":["monocoque"],
 "lift_to_drag":3.1,
 "balance":{"status":"defined","front_axle":[0,0,0],"wheelbase":1.53,
            "front_load":250.0,"rear_load":330.0,"front_percent":43.1,
            "cop_from_front_axle":0.870}}
```

- Rows: each car zone (`kind:"zone"`, in `loads.json` order), then each group
  (`"group"`, in spec order), then the car total (`"total"`, name `"Total"`).
  Forces are projected on the axes in N; `downforce` is `-lift`.
- Moments (N m) are about `moment_centre` and projected: `pitch` on the side
  axis, `roll` on the drag axis, `yaw` on the lift axis. With the defaults a
  positive pitch raises the front (the nose, upstream) of the car.
- Coefficients divide forces by `q A` and moments by `q A L`, with `q = ½ ρ V²`
  from the run's reference values. Every key of `coefficients` is always
  present. Force coefficients are null when density, velocity or area is
  missing or not positive; moment coefficients are also null without a
  positive length. `cla` and `cda` are lift and drag over `q` in m² (`cl·A`,
  `cd·A`): they need density and velocity but not the area. `downforce_share` is the row's share of the total
  downforce (null when the total is zero).
- `lift_to_drag` is total downforce over total drag (null for zero drag).
- `balance`: with F the total force and M_C the total moment about the moment
  centre C, the moment about the front contact point P_f is
  `M_f = M_C + (C - P_f) × F`. The rear tyres carry
  `rear_load = (M_f · side) / wheelbase` and the front
  `front_load = downforce - rear_load` (N, positive pushing the tyres down),
  so drag acting above the ground is included. `front_percent` is
  `100 front_load / downforce` and `cop_from_front_axle` is
  `wheelbase · rear_load / downforce`, the distance behind the front axle at
  which the total downforce alone would give the same split. `status` is
  `"defined"`, `"needs_vehicle"` (no vehicle in the spec; every other balance
  field, including `front_axle` and `wheelbase`, null) or
  `"no_downforce"` (total downforce ≤ 0: loads given, percentage and CoP null).
