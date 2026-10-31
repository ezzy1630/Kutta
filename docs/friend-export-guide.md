# What to send from a Fluent run

Thanks for helping. Nothing here changes how you run Fluent. We only need copies of files from a run you've already finished. Your car files stay private and are used only to compare our solver against yours.

## From one finished run

1. **The final case and data files**: `<name>.cas.h5` and `<name>.dat.h5` from the same save. Please send the final pair, not autosaves (those look like `<name>-1500.cas.h5`).
2. **Every `*.out` file** in the run folder. These are the lift/drag/force histories Fluent writes for each report definition.
3. **The console log, if one exists**: a `.trn` file, or the batch or cluster log. If there isn't one, just tell us:
   - the machine (CPU, GPU, RAM);
   - how many cores or which GPU Fluent ran on;
   - roughly how long the run took.

**If you use Workbench:** zip the Fluent folder of that system (`<project>_files/dp0/FFF/Fluent/`). That contains everything above.

**Not needed:** the separate `.msh.h5` (the case file already contains the mesh), screenshots, or report printouts.

## Optional, only if it costs you nothing

- **Turn on the transcript:** File → Write → Start Transcript at the start of a run. This records hardware and per-iteration timing, which gives us a fair speed comparison.
- **Note the case:** a one-line note of what it is (e.g. "full car, 15 m/s, rotating wheels, ride height baseline").

## How the files are used

- **Same-mesh comparison:** we load your mesh and settings and solve them with Kutta. Then we compare downforce, drag, balance and surface pressure with your Fluent results on the identical mesh.
- **Fair timing:** we time Kutta on your own machine against your Fluent run, and report both with the hardware named.
- **Privacy:** nothing is uploaded or shared.

## If Kutta rejects a setting

Keep the original Fluent files. By default, `kutta solve` lists unsupported
settings and exits without solving. It never silently disables warped-face
gradient correction (WFGC, `recon/cell-lsf?`).

For an explicitly approximate diagnostic, you can permit WFGC on to run with
WFGC off:

```sh
kutta solve car.cas.h5 --data car.dat.h5 --approximate wfgc --out kutta-approximate
```

`--approximate` takes comma-separated names; **`wfgc` is the only allowed name**.
All other unsupported settings still fail. Unknown, empty or duplicate names
fail too. Missing/unparsed WFGC settings still fail. An opt-in on a case with
WFGC already off records no approximation. The original case is never changed.

The [V6 initialization comparison](../src/kutta/fv/VERIFICATION.md#v6-initialization-comparison)
measures the WFGC-disabled mixing-elbow diagnostic against Fluent: CPU hybrid
velocity L2 **3.366748%**, pressure L2 **2.665514%**, k L2 **14.722368%**, omega
L2 **2.805466%**, and pressure-drop error **1.973008%**. These include other
solver discrepancies; they do not isolate WFGC's effect or bound errors on a car.
The original same-settings accuracy gate remains unmet.

Approximate runs print **APPROXIMATE RESULT** at start and end. `summary.json`
records `approximate=true` and each approximation's setting, Fluent value,
Kutta action, measured effect and documentation path. `forces.csv` and
`residuals.csv` each start with a `# APPROXIMATE RESULT` comment (skip `#` lines
in CSV readers); `result.vtu` has the same warning in its XML header comment.
Retain these markers with exported results so they remain labelled as approximate.
Without an applied approximation, the CSV/VTU headers are unchanged and the
summary records `approximate=false` with an empty `approximations` list.
