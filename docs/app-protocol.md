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
| `kutta mesh GEOM.stl [MORE.stl ...] [--setup mesh.json] --out DIR [--case] [--progress jsonl] [--json] [--force]` | STL to Fluent mesh and quality report (v1 additive, below) |
| `kutta new-case MESH [--template fsae-external] [--speed M/S] [--role ZONE=ROLE]... [--wheel ZONE=R,X,Y,Z,AX,AY,AZ]... [--out CASE.kutta.json [--force]] [--json]` | a starter Kutta case for a bare mesh (v1 additive, below); without `--out` the case is printed |
| `kutta convert CASE [--setup overrides.json] --out CASE.kutta.json [--force] [--json]` | a case's settings, with overrides applied, written as a Kutta case (v1 additive, below) |

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

Rotating walls additionally carry `condition.angular_speed` (rad/s),
`condition.axis_origin` (metres) and `condition.axis_direction` (v1 additive).
Unknown values are null; older CLIs may omit them.

Every command that takes a case also takes a Kutta case, `*.kutta.json` (a
mesh plus a JSON setup, [docs/kutta-case.md](kutta-case.md)). For a Kutta case
`file` gains (v1 additive) `"format":"kutta"` (`"fluent"` otherwise),
`setup_sha256` and `"mesh_file":{"path":...,"sha256":...}`; `file.sha256` is
the case identity, the SHA-256 of the setup's and the mesh's digests, and
`fluent_version` is null. `--data` is rejected: a Kutta case has no Fluent
data file. `summary.json` likewise gains `case_format`, `setup_sha256`,
`mesh_file` and `mesh_sha256` (null for a Fluent case), and its `case_sha256`
is the same identity.

## Kutta cases from a mesh (v1 additive)

`kutta new-case MESH --json` reads the mesh and prints one object; with
`--out` it also writes the case there (refusing to replace a file without
`--force`) and stores the mesh path relative to it:

```json
{"v":1,"template":"fsae-external","speed":15,
 "mesh":{"path":"/abs/car.msh","dimension":3,"cells":70,"faces":269,"nodes":140,"polyhedral":false},
 "zones":[{"id":6,"name":"ground","mesh_type":"wall","kind":"boundary","faces":24,
           "role":"moving-ground","suggested_role":"moving-ground",
           "reason":"name contains \"ground\"; the road moves with the freestream","group":null},
          {"id":10,"name":"front_wing","mesh_type":"wall","kind":"boundary","faces":10,
           "role":"wall","suggested_role":"wall","reason":"a wall of the car","group":"Front wing"}],
 "roles":["velocity-inlet","pressure-outlet","moving-ground","wall","symmetry","interior","stationary-ground"],
 "notes":["Reference area is Fluent's default, 1 m^2; ..."],
 "case":{...the Kutta case...},
 "out":"/abs/car.kutta.json"}
```

`suggested_role` comes from the zone's name (inlet, outlet, ground/road/belt,
symmetry/top/side/farfield words), else from the mesh's zone type, else wall;
`role` is the suggestion unless `--role ZONE=ROLE` overrides it. `roles` lists every
role in this order; `stationary-ground` (v1 additive, appended) is a no-slip
wall at rest, `{"type":"wall","motion":"stationary"}`, for a fixed-floor wind
tunnel. A ground/road/belt zone is still suggested `moving-ground`, the FSAE
default; neither ground role is ever a force-report zone. A zone suggested as
ground but given the plain `wall` role joins the car's reports, and the notes
say so. Wall zones
carry an aero `group` (Front wing, Rear wing, Wing, Undertray / diffuser,
Wheels, Body), the app's grouping. The `fsae-external` template: flow along
+x at `--speed` (default 15 m/s), up +z (+y in 2D), k-ω SST, air, velocity
inlet and pressure outlet with 1 % intensity and viscosity ratio 10, a
moving ground at the freestream speed, symmetry planes, stationary car walls,
second-order k and ω, SIMPLEC pressure-velocity coupling (v1 additive: Kutta's
Coupled solver limit-cycles on 3D external aero, so its forces do not settle;
SIMPLEC converges the same discretization), reference velocity and density of
the freestream, and drag (+x) and downforce (lift along −z, −y in 2D)
coefficient reports on every car wall. On a 3D mesh with polyhedral cells the
template also sets `methods/warped_face_gradient_correction` false (v1
additive): Fluent's default turns it on there and Kutta does not implement it,
so the case states what is solved and needs no `--approximate wfgc`. A mesh
with tetrahedra adds a note recommending `kutta polyhedra`.
See [car mesh convergence](../validation/results/car_mesh_convergence/REPORT.md).

`--wheel ZONE=RADIUS,X,Y,Z,AX,AY,AZ` is repeatable (v1 additive).
It requires an exact wall-zone name, a positive wheel radius in metres, its
axis origin in metres, and a nonzero signed axis. The saved `case` contains a
rotating wall with `angular_speed = speed / radius` in rad/s, the supplied
axis/origin, and `relative:false`. For flow +x and up +z, use axis `0,-1,0`
so the wheel's bottom moves along +x with the ground. Origins are supplied,
not inferred from names. Unspecified wheels remain stationary and the notes
say so. Duplicate wheels, invalid geometry and non-wall targets are errors.
The translating ground keeps the existing speed and direction convention.

The version-1 case object also accepts `pressure-inlet`, `mass-flow-inlet`,
`porous-jump`, named `cell_zones` with constant `frame_motion` / `porous`, named constant
`materials` and per-zone material assignments (uniform properties required),
`pressure_reference_point`, `relative_frame_origins_absolute`, SST
`cross_diffusion_clip`, configurable `models.intermittency` constants, and
k-ε `source_modification` / `time_marching`. See the typed properties and
limits in [the case format](kutta-case.md). These additions use the shared
Fluent settings and solve protocol; no new progress events are required.

`kutta convert CASE --out OUT --json` (CASE a Fluent or Kutta case;
`--setup` applies a setup-overrides file first, as `solve --setup` does)
prints `{"v":1,"out":"/abs/OUT","replaced_input":false,"omitted":[{"path":...,"reason":...,"class":"informational"}],"run_controls_not_saved":[...]}`.
It exits non-zero, writing nothing, when an `unsupported` setting cannot be
expressed; `informational` omissions do not change the solution. Failures of
either command with `--json` print one `error` event, as `solve` does, whose
`issues` name each rejected setting. The app
saves an edited case with it.

## Tetrahedra to polyhedra (`polyhedra`, v1 additive)

`kutta polyhedra MESH --out OUT.msh [--force] [--json]` is Fluent's
tet-to-polyhedra conversion (`mesh/polyhedra/convert-domain`): each vertex of
the tetrahedral region becomes one polyhedron, the median dual (a quarter of
every tet around the vertex, cut at edge midpoints, triangle centroids and
tet centroids). Prisms, hexes and pyramids are kept; their faces on the
converted region gain the edge midpoints, so the mesh stays conforming.
Volume, boundary zones and their types are preserved. A vertex whose dual
cell wraps a sharp convex edge so that its centroid falls behind one of its
faces keeps one hex quarter per tet instead (`split_vertices`). Steady solves
on tetrahedral regions oscillate instead of converging; the converted mesh
converges (see the report linked above). The output is a Fluent ASCII `.msh`
(one cell zone per input cell zone); internal face zones are not supported.
With `--json` it prints one object:

```json
{"v":1,"input":"/abs/car.msh","out":"/abs/car-poly.msh","input_cells":130331,
 "tetrahedra":55371,"polyhedra":11149,"split_vertices":7,"kept":74960,"cells":86215}
```

Failures print one `error` event, as `new-case` does.

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

`iteration` events are emitted at the monitor interval (the monitor samples),
v1 additive at iterations 1 to 5, and v1 additive at the forces-settled
condition's confirming sample, one iteration after a monitor sample (see
Convergence criteria). The solver monitors those anyway,
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
kutta export-geometry mesh.msh --out DIR --json
```

```json
{"v":1,"event":"geometry","manifest":"/absolute/DIR/geometry/manifest.json","zones":6,"triangles":16978,"seconds":0.15}
```

Accepts a Fluent bare mesh (`.msh`, `.msh.h5`), a Fluent case or a Kutta case.
New Case uses this existing command before a case exists to estimate wheel
geometry from each zone's surface, independently of its name. The proposal is
not a boundary condition until the user enables rotation and confirms it.
For 2D inputs, exported surfaces have an artificial unit-depth extrusion;
estimates use the original XY extent, centre Z = 0 and axle along Z.

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
Each slice descriptor carries `skipped_cells` (v1 additive): cells left out of
that slice because their faces do not close or their section cannot be
triangulated (tangled cells), instead of failing the export. It is 0 on a
valid mesh; the app may say "N cells could not be sliced" when it is not.
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

`methods/limiter` (v1 additive) is Fluent's limiter type (`recon/lim/type`),
`standard` or `differentiable`, labelled "Standard" and "Differentiable". It is
editable and depends on `methods/gradient-limiting` (`active_when` `[true]`).
Standard is Barth-Jespersen, Fluent's default and the type every corpus case
saves. Differentiable is a Venkatakrishnan limiter with Wang's threshold, 0.05
times each velocity component's range. It acts on momentum only (turbulence
keeps Standard) and adds a convergence note. It changes the fixed point: choose
it to converge a run that stalls on the Standard limiter's switching (Fluent
Theory Guide 23.3.4.3), not to reproduce a Standard-limiter Fluent run. Fluent
documents neither its formula nor its coefficient, and a case that saves type 2
is read with Fluent 12's numbering; validation/results/flux_limiter/REPORT.md
records the evidence.

`controls/courant-number` (v1 additive; editable since the Coupled
convergence change) is Fluent's Flow Courant Number as saved
(`pressure-coupled/courant-number`, null when absent). Kutta's Coupled solver
applies it as implicit relaxation, a_P(1 + 1/CFL), on top of the pseudo time
step (Theory Guide 18.4.4); a case that saves none runs Fluent's default, 200.
It is a convergence control: the fixed point does not depend on it. With
pseudo time it is editable (`min` the smallest positive float, `max` 1e15) and
depends on `methods/pressure-velocity-coupling` (`active_when`
`["coupled"]`). Without pseudo time Kutta never runs Coupled, so it is
inactive and read-only, and its `reason` says the case runs SIMPLEC. Fluent
itself documents the Flow Courant Number only for Coupled with the pseudo time
method off; with Fluent's global pseudo time step it uses the pseudo-time
explicit relaxation factors instead (User's Guide 36.3), which Kutta also
honours (`controls/mom-relaxation`, `controls/pressure-relaxation`).

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
a zero-iteration run. `done` and `summary.json` add `stop_reason: "stopped"`, `"iteration_limit"`,
`"converged"` or (v1 additive) `"statistically_steady"`; `converged` is true
only for `"converged"`, when convergence checking is on (see Convergence
criteria below). A run whose iteration diverges (a non-finite linear system,
residual or force) ends with an `error` event whose message begins "the
solution diverged at iteration N" and names its last finite residuals. That
`error` event adds `"verdict":"diverged"` (v1 additive); other errors carry no
`verdict`.

`done` and `summary.json` add `verdict` (v1 additive), the run's outcome in one
field, which readers should prefer to `converged` and `stop_reason`:

| `verdict` | `converged` | `stop_reason` | meaning |
| --- | --- | --- | --- |
| `converged` | true | `converged` | every convergence condition met, residual criteria included |
| `forces_settled_residuals_waived` | true | `converged` | the forces settled over twice the window while residuals stalled above their criteria (`residuals_waived: true`); not converged by residuals |
| `statistically_steady` | false | `statistically_steady` | forces cycle about a steady mean |
| `iteration_limit` | false | `iteration_limit` | the iteration limit came first |
| `stopped` | false | `stopped` | stopped on request |
| `diverged` (error event only) | | | the iteration went non-finite; no `done` or `summary.json` |

`converged` and `stop_reason` keep their meaning: a forces-settled stop with
waived residuals is still `converged: true`, so a reader of `converged` alone
cannot tell it from residual convergence. JSONL mode watches stdin EOF only when stdin is a pipe
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
implemented momentum/k/omega/epsilon upwind schemes, gradient limiting and its
limiter type (`methods/limiter`), the pressure interpolation
(`methods/pressure-scheme`, options `standard` and `second-order`), the pressure-velocity coupling (`methods/pressure-velocity-coupling`, below), pseudo-time
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
A cell without a positive kernel, or with an inverted decomposition
tetrahedron, is left out whole (no overlapping or inverted tetrahedra are
emitted) and counted in the iso-surface descriptor's `skipped_cells` (v1
additive); earlier CLIs failed the export on such a cell. Vertices
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
on. Since v1 additive, checking is **on by default for a case with a force,
drag, lift or downforce report** (a car run), together with Kutta's
forces-settled condition below; a case without such a report keeps checking
off and runs to its iteration limit, as before. Every control below is a
settings-inventory entry in group `Monitors`, editable through `solve --setup`.
They are run controls like `run/iterations`: they say when to stop, not what
case is solved, so they never mark a result MODIFIED and are echoed under
`run_controls`. (`setup::run_control` covers the `run/` and `monitors/`
prefixes.) Without `--iterations`, `run/iterations` defaults to 5000 for a case
with a force report (a cap: such runs stop by themselves) and 1000 otherwise.

| Path | Kind | Default | Meaning |
| --- | --- | --- | --- |
| `monitors/check-convergence` | bool | `true` with a force report, else `false` | Stop when the conditions below are met |
| `monitors/condition` | enum `"all"`, `"any"` | `"all"` | "All conditions are met" / "Any condition is met" |
| `monitors/residual/<equation>/check` | bool | case value | Fluent's Check Convergence for that residual |
| `monitors/residual/<equation>/absolute-criterion` | number ≥ 0 | case value | Fluent's Absolute Criteria (scaled residual) |
| `monitors/report/<report>/active` | bool | `false` | Use this report definition as a convergence condition |
| `monitors/report/<report>/stop-criterion` | number > 0 | `0.001` | Fluent's Stop Criterion (relative change) |
| `monitors/report/<report>/previous-values` | integer ≥ 1 | `5` | Np: report values to include (Fluent's Use Iterations, counted in monitor samples) |
| `monitors/report/<report>/ignore-before` | integer ≥ 0 | `100` | Fluent's Ignore Iterations Before |
| `monitors/forces-settled/active` | bool | `true` with a force report | Kutta's forces-settled condition (not a Fluent setting) |
| `monitors/forces-settled/tolerance` | number > 0 | `0.001` | Largest relative range of each force report over the window |
| `monitors/forces-settled/window` | integer ≥ 1 | `300` | Window, in iterations |
| `monitors/forces-averaged/active` | bool | `true` with a force report | Kutta's statistically-steady stop (not a Fluent setting, not convergence) |
| `monitors/forces-averaged/tolerance` | number > 0 | `0.001` | Largest relative change of each force report's mean between two consecutive windows |
| `monitors/forces-averaged/range` | number ≥ 0 | `0.1` | Largest relative change (growth or decay) of each report's range from the first window to the second |
| `monitors/forces-averaged/window` | integer ≥ 1 | `500` | Averaging window, in iterations |

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
  is checked), each active report and forces-settled when active. `all` needs
  every condition met at the same sample, `any` needs one. Checking on with no
  condition emits the warning
  "Convergence checking is on but no residual is checked and no report
  condition is active; stopping at the requested iteration limit." and runs to
  the iteration limit.

- **Forces settled** (v1 additive) watches every force, drag, lift and
  downforce report at once (not only Fluent-active ones). At a sample n it
  holds when, for each report, the monitor samples from iteration n−window to
  n cover the whole window (one at or before n−window) and their
  max − min ≤ tolerance × scale, where scale is the window's mean |value|, or
  1 % of the largest report's mean |value| when that is larger (so a report
  near zero, such as the side force of a symmetric car, is judged against the
  main loads). The settled forces must then be confirmed: when they are what
  stands between a sample and convergence, the next iteration is sampled too
  (off the monitor stride, so a cycle whose period divides the monitor
  interval cannot look steady) and the condition is met only if the forces
  are still settled there. That sample is an ordinary `iteration` event.
  Fluent's lag-residual report conditions are unchanged; this condition is
  stronger against slow drift and oscillation, because it bounds the range
  over the window rather than consecutive differences.
- **Stalled residuals** (v1 additive): residuals that stall above their
  criteria (the external-aero continuity plateau, k or ω) must not hold a run
  whose loads have stopped changing. When forces-settled is active and the
  forces have also been settled over twice the window (confirmed the same
  way), the residual condition counts as met; its object then carries
  `"waived": true` (its `met` stays false), `met_by` omits `"residuals"`, and
  `done`/`summary.json` report `"residuals_waived": true`. Forces that still
  drift or oscillate beyond the tolerance never stop a run as converged.
- **Statistically steady** (v1 additive; not a convergence condition). Some
  steady RANS solutions are unstable: the iteration settles into a bounded
  cycle that no amount of iterating converges (the hex Ahmed body's wake on
  finer meshes under SIMPLEC, validation/results/default_convergence_fine/REPORT.md). From
  iteration 2 × `forces-settled/window` + 1 on, every force report is
  evaluated at every iteration (this does not change the solve). At a monitor
  sample n, the samples of iterations n−2·window+1 … n−window and
  n−window+1 … n are compared for each report: the change of its
  Hann-weighted mean (w = sin²(π(k+½)/window), which suppresses the error a
  partial cycle at the window ends puts in a plain mean) must be
  at most `tolerance` × scale (scale as for forces settled), and the second
  window's max − min within `range` × the first's of it, plus
  `tolerance` × scale. When every report passes and the convergence
  conditions are not met, the run stops with `"stop_reason":
  "statistically_steady"` and `"converged": false`. A cycle that drifts or
  grows never stops this way; the run then ends at its iteration limit. A
  decaying cycle is still converging and does not stop this way either. This is
  Kutta's automation of what Fluent users do with a report definition's
  "Average Over" running mean.

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
  "reports":[{"report":"drag","stop_criterion":0.001,"previous_values":5,"ignore_before":100}],
  "forces":{"tolerance":0.001,"window":300,"reports":["drag","downforce"]}}
```

`residuals` lists checked equations only; `forces` (v1 additive) is present
only while the forces-settled condition is active. With checking on, each `iteration`
event adds the state of every condition at that sample, residuals first and
then active reports in report order:

```json
"convergence":{"met":false,"conditions":[
  {"kind":"residuals","met":false,"equations":{"continuity":{"value":0.0021,"criterion":0.001,"met":false}}},
  {"kind":"report","report":"drag","value":0.0042,"criterion":0.001,"met":false},
  {"kind":"forces","met":false,"settled":true,"criterion":0.001,"window":300,
   "reports":{"drag":0.00031,"downforce":0.00072}}]}
```

A report's `value` is null until it has enough samples. The `forces`
condition's `reports` give each report's relative range over the window
(null until the window is covered); `settled` is the window test alone and
`met` adds the confirmation. `done` and
`summary.json` add `convergence`: the setup object plus `met_by` (the met
conditions at the stopping sample, as `"residuals"` or `"report:<name>"`, empty
unless converged), `iteration` (the stopping iteration, null unless
converged) and, v1 additive, `residuals_waived` and `averaged`. `averaged`
is null for a converged run and for a run that stopped before averaging had
two complete windows. Otherwise it is the averaging state at the stopping
sample (`"steady": true` when that is why the run stopped, `false` for a run
that reached its iteration limit or was stopped while its mean drifted or its
range grew or decayed):
`{"window","from","to","steady":true,"reports":{"<report>":{"mean","min","max",
"drift","range_change","met","coefficient":{"mean","min","max"}}}}`, forces in N over
iterations `from`…`to` (both windows), with `coefficient` when the reference
values define one; `mean` is Hann-weighted over both windows and
`range_change` is the second window's range over the first's, minus 1
(negative when the cycle decays). While averaging has two complete windows, each `iteration`
event's `convergence` object also carries that `averaged` state (without
`coefficient`). The setup object adds, v1 additive, `"averaging":{"tolerance",
"range","window","start"}` while the stop is active. `summary.json` sets
`convergence_criteria_available: true`.
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

## Automatic STL meshing (`mesh`, v1 additive)

```sh
kutta mesh car.stl front-wing.stl --setup mesh.json --out car-mesh --case --progress jsonl
kutta solve car-mesh/case.kutta.json --backend metal --iterations 2000 --out car-run
```

`kutta mesh GEOM.stl [MORE.stl ...] [--setup mesh.json] --out DIR [--case]
[--speed M/S] [--wheel ZONE=R,X,Y,Z,AX,AY,AZ]... [--progress jsonl] [--json]
[--force]` is additive to v1. Omitting setup uses SI FSAE defaults: +x flow,
+z up, 20 m/s for optional case/y+ estimates, 0.25 m base size, 0.008–0.03 m
surface sizes, 20° feature angle, three cells across gaps down to 2 mm cells,
four prism layers (2 mm first height, growth 1.3), body-relative refinement
boxes (below) and 10 million maximum cells. These are starting sizes, not an
accuracy or cell-count guarantee; on the synthetic FSAE car they give about
4.2 million cells and an estimated 8.9 GB of solver memory
(`estimated_solver_memory`), with every wing slot open. Budget exhaustion
fails rather than silently coarsening the mesh. (v1 additive change of
defaults: the minimum surface size was 0.006 m and no boxes were added; that
gave 4.4 million cells.)

Setup is strict JSON (1 MiB maximum): unknown keys, duplicate keys, wrong types,
invalid ranges and unsupported directions fail with qualified paths such as
`mesh/surface/min_size`. Missing properties use defaults; `v:1` is required.
All lengths in setup are metres; `units` scales **STL coordinates only**.

```json
{"v":1,"units":"m","flow":"+x","up":"+z",
 "tunnel":{"upstream":3,"downstream":6,"side":2,"top":2,"ground":true,"symmetry":null},
 "base_size":0.25,
 "surface":{"min_size":0.008,"max_size":0.03,"curvature_degrees":20,"proximity_cells":3,"proximity_min_size":0},
 "refinement_boxes":[{"name":"wake","min":[1,-1,0],"max":[5,1,1.5],"size":0.03}],
 "layers":{"count":4,"first_height":0.002,"growth":1.3,"zones":["*"]},
 "max_cells":10000000}
```

Tunnel distances are multiples of the body's x length, measured from its bounds.
The ground is z=0 (geometry below it is rejected). Without ground the lower
boundary extends `side` body lengths below the body. `symmetry:"y"` keeps y≥0.
Tunnel upper limits round outward to a whole base cell; lower limits are exact.
Octree sizes are powers of two below the base size, so the achieved minimum can
be smaller than `min_size`. Refinement boxes use world coordinates. **When the
setup omits `refinement_boxes`** (v1 additive), two body-relative boxes are
added, because drag depends on the volume cells near the body and in the wake,
not only on surface sizing: `body`, the body bounds extended 0.1 L upstream,
sideways and up (0.25 L downstream, down to the ground), at 2 × `max_size`; and
`wake`, from the body's rear to 1 L behind it, 0.05 L wider than the body and
0.1 L higher, at 4 × `max_size` (L is the body's x length). An explicit list,
including `[]`, is used as given. `castellation/refinement_boxes` reports the
boxes used and the octree cell each achieves. Face, edge
and corner neighbours are balanced to at most 2:1. Curvature sizes facets by
their radius times the feature angle. Gaps (to another facet across the fluid,
or to the ground) are sized to `proximity_cells` cells across, refining below
`min_size` down to `proximity_min_size` (0, the default, means `min_size/4`,
2 mm at the default sizes; an explicit 0 is accepted, a v1 additive fix).
Next to layered
walls, cells are at most four last-layer thicknesses, so the layer-to-core
size jump stays small. Tessellate curved STL surfaces finely enough to describe
their curvature.

Binary and ASCII STL, including multiple ASCII named solids, are accepted. Binary
files use the file stem as their wall-zone name; ASCII solids use their names
(or file stem when unnamed). Repeated names receive `-2`, `-3`, … suffixes.
Reserved tunnel names (`inlet`, `outlet`, `ground`, `top`, `bottom`, `symmetry`,
`side-low`, `side-high`) cannot be body-zone names. Normals are reconstructed.
Welding within 1e-7 of the body diagonal, degenerate/duplicate triangle removal,
orientation repairs, open/nonmanifold edges and conflicts are reported. Large
holes are not filled. Independent majority ray queries classify the union of
connected shells without a flood fill, so a small crack cannot flood the domain.
Intersecting shells are treated as a solid union; nested shells are not cavities.

The mesh follows snappyHexMesh. Octree cells inside the body or crossed by an
STL facet are removed (castellation), keeping the fluid connected to the
tunnel. Wall points then snap onto the STL: each moves to the least-squares
intersection of the tangent planes at its faces' nearest surface points, which
captures feature edges and corners without a feature-line input. The motion
decays over three point layers into the volume. Wherever a cell would fall
below orthogonal quality 0.15, the motion retreats locally (halved, then
zero), so the remaining geometric error is reported rather than hidden
(`snap/point_distance_m`). Slots narrower than about two cells close; see
`castellation/gaps` for each zone pair's width and minimum fluid cells across.
When `proximity_cells` is 3 or more (the default), no fluid cell is left with
the body on two opposite sides (v1 additive): such one-cell notches, where a
wing slot meets an endplate for example, are removed repeatedly and counted
in `castellation/thin_gap_cells_removed`. k-ω diverged in them at the default
FSAE size. With `proximity_cells` 1 or 2 the setup asks for thin slots, so
nothing is removed (the count is 0).

Prism layers are inserted by pushing the snapped mesh off the selected walls
(`layers/zones`, default all) along point normals and filling the gap with one
stack per wall face. A stack collapses where any cell would fall below quality
0.15 or where retreat thins it below a quarter of its target. Where layers end,
the termination points split so the unlayered wall stays on the surface; a
termination whose split fails quality stays lifted instead (a step), and
`layers/terminations` counts both. `count:0` disables layers. The y+ estimate
uses a flat-plate friction correlation and actual wall-adjacent centroid
distances, not solved shear.

Successful output contains `DIR/mesh.msh` (Fluent ASCII polyhedra) and
`DIR/mesh-report.json`. Report keys include `v`, `cells`, `faces`, `nodes`,
`cell_types` (`hex`, `hex_without_hanging_nodes`, `poly`, `prism`,
`hex_fraction_of_core`), `zones` (`name`, `type`, `faces`), `body_bounds`,
`bounds`, `quality`, `repairs`, `castellation`, `snap`, `layers`,
`missing_surface_zones`, `estimated_first_cell_y_plus`, `stage_seconds`,
`wall_seconds`, `peak_memory_bytes`, `estimated_solver_memory` (v1 additive),
`warnings`. `cell_types` counts cells by how the mesher built them: an octree
cube, even with 2:1 hanging faces, counts as `hex` and each layer cell as
`prism`, although a layer cell on a quadrilateral wall face has six
quadrilateral faces. `cell_types_by_shape` (v1 additive: `tet`, `hex`, `poly`,
`method`) counts them as `kutta inspect` and `new-case` do, by face shape alone:
four triangles a tet, six quadrilaterals a hex, anything else (a cube with
hanging nodes, a merged cell, a layer cell on a triangle) a polyhedron; its
counts are the ones those commands print. Fluent stores every cell as a
polyhedron.
`quality` includes minimum/mean orthogonal quality, `orthogonal_histogram`,
`cells_below_0_05`/`_0_15`, `worst_cell`, skewness, the adjacent volume
ratio's distribution split into core and layer faces (`max`, `p50`, `p99`,
`p999`, histogram), minimum volume and nonconvex cells (a count, not an error).
`skewness_max` is Fluent's equiangle skewness, 0 (ideal) to 1, from each
face's corner angles against a regular polygon's (`skewness_method`; nodes
where the face boundary turns by under 30°, hanging nodes, are not corners),
with `faces_skewness_above_0_95` (v1 additive). Before this fix it held an
unbounded centroid-offset ratio, which exceeded 1 in thin, slightly sheared
layer cells; that value is now `face_skewness_max` (v1 additive): the distance
from a face centroid to where the owner–neighbour link crosses the face plane,
over the link length (OpenFOAM's face skewness). Each volume-ratio histogram
bin counts internal faces in `faces` (v1 additive); its `cells` key carries the
same number and is a deprecated alias. `castellation/gaps` entries add
`contact` (v1 additive): true where the shells touch (width 0, a tyre on the
ground), so a contact is not read as a closed slot (a positive width with zero
cells across).
`snap` reports snapped/partial points, feature points, removed and merged cells
and the point-to-surface distance (absolute and over the local cell size).
`layers` reports `requested`, `generated`, `stacks`, `cells`, `coverage`
(`full`, `thinned`, `none` area fractions of the selected walls),
`first_height_achieved`, `terminations` and `collapse_rounds`. The y+ estimate
gives area-weighted percentiles and area fractions below 5, 5–30, 30–300 and
above 300. Closure, manifold edges, positive volumes and an orthogonal floor
of 1e-8 are enforced before writing; snapping and layers keep every cell at
0.15 or above. Peak memory is process footprint, not an estimate of solver GPU
memory. Unavailable memory readings are null. `estimated_solver_memory`
(v1 additive) estimates the solve's device memory as cells times bytes per
cell measured on Metal for this mesher's car meshes: `bytes` and
`bytes_per_cell` (2,100) for the default steady solve (SIMPLEC, SST),
`coupled_bytes` and `coupled_bytes_per_cell` (2,600) for Coupled, and
`warning_above_bytes` (10 GB). It is not measured on CUDA and excludes a CUDA
context's own memory. Above 10 GB `warnings` says a 12 GB GPU may not fit the mesh.

`--case` calls the same `fsae_starter` path as `new-case`: inlet, outlet, moving
ground, symmetry sides/top and named body walls, with the same wheel syntax.
The resulting `case.kutta.json` explicitly disables warped-face gradient
correction; the solver can read it directly.
A separately generated `new-case mesh.msh` retains that command's existing
polyhedron defaults and needs its documented WFGC setting/approximation choice.
Reference area remains the starter's 1 m² default; set the actual reference area
before comparing coefficients. Wheel rotation is never inferred from a name.

Existing named outputs require `--force`; other files in DIR are retained. Work
is staged before publication. On success, JSONL emits monotonically
nondecreasing `fraction` progress (a later message in a retreat or collapse
round repeats the fraction reached rather than reporting a lower one) with `stage` in `read|repair|octree|trim|layers|quality|write`:

```json
{"v":1,"event":"mesh-progress","stage":"octree","fraction":0.15,"message":"Refining surface, slots and volume boxes"}
{"v":1,"event":"mesh-done","cells":34178,"faces":305262,"nodes":100051}
```

The final `mesh-done` includes the full report. `--json` alone prints only that
final object. Errors exit nonzero and emit
`{"v":1,"event":"error","message":"...","advice":"..."}` with either JSON
option; stderr can contain the human diagnostic. `mesh-done` means mesh validity,
not aerodynamic convergence or force accuracy. Cancellation can terminate the
process; a staged scratch directory may need removal after a hard kill.

Generate reproducible inputs with
`uv run tools/mesher/generate_geometry.py --out scratch-geometry`.
The synthetic FSAE wheels stand on 10 mm contact plinths (a tangent or flattened
tyre leaves a wedge of cells walled on two or three sides, in which k-omega
diverged), and every other part clears the tyres by at least 50 mm;
`fsae-wheels.json` gives each wheel's rolling radius, centre and axis for
`--wheel`. Generated STLs/meshes are local artifacts, not repository
fixtures.
