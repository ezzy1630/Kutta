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
