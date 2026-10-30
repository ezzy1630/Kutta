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
{"v":1,"event":"snapshot","run_id":"run-123-0","sequence":1,"iter":20,"elapsed":0.9,"manifest":"/absolute/DIR/live/run-123-0/snapshot-1/manifest.json","snapshot_seconds":0.006}
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
convergence criteria. JSONL mode watches stdin EOF only when stdin is a pipe
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
constant material density/viscosity, used reference area/density/speed/pressure,
implemented momentum/k/omega upwind schemes, gradient limiting, pseudo-time
scale and active equation relaxation, standard numeric initialization values
when not seeded, force-report directions, and run iteration and monitor
intervals. Changing the initialization method, operating pressure, unused
reference length/viscosity, model selection, velocity specification, wall
motion or convergence criteria is read-only. A displayed Fluent value is not
evidence that Kutta can edit it.

The turbulence specification's `options` are restricted to the methods whose
required values were saved as constants or resolved expressions (plus the
current method), so switching never runs on a guessed value. Turbulence values
a method of the case's model can use stay editable while inactive, so the app
can switch method and set its values in one setup file; an override of an
inactive value is accepted and has no effect, as Fluent keeps such values.
Values no method of the model reads (epsilon for k-omega) are read-only.

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
again as `run_controls`. `modified` is true only when a case-changing override
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
frames or run-data failures are exit 1. A failed export preserves the old view.

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
| `wall_shear`, `wall_shear_magnitude` | Pa | saved SST stress on fluid, vector/magnitude |
| `y_plus` | 1 | saved SST u-star wall diagnostic |
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
