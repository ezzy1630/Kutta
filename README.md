# Kutta

A local, GPU-native CFD solver for car aerodynamics. Its first goal is to reproduce Ansys Fluent's steady incompressible k-ω SST results on the same mesh, within 1–2% on downforce and drag, much faster on one consumer GPU. The first targets are an NVIDIA RTX GPU on Windows and Apple Silicon on macOS.

**Status: early development.** Done so far:

- device runtime: CPU reference, Metal, CUDA;
- mesh I/O: Fluent `.msh` / `.cas` / `.cas.h5`, NASA Plot3D, VTU;
- linear solvers, including aggregation AMG;
- Fluent-style SST constitutive kernels and smooth-wall closure, verified against
  the saved mixing-elbow fields on CPU and Metal.

SST transport, effective momentum stress and wall shear are coupled into the
steady solver. The strict Fluent setup adapter rejects unsupported spatial
settings. Hybrid initialization solves the velocity-potential and pressure
Laplace problems using stored case controls. V6 retains failing accuracy gates;
its explicitly labelled elbow diagnostic disables Fluent's undocumented
warped-face correction.
See [the measured results](src/kutta/fv/VERIFICATION.md#v6-mixing-elbow-sst)
and the gates in [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

FSAE wall motion, porous resistance, absolute-velocity MRF and additional flow
boundaries now have [dedicated implementation and analytic gates](docs/research/fsae-physics.md).
Native rotating-wall and fan polynomial import remain held by failing V7.2/V12
gates. These tests do not establish Fluent same-mesh agreement.

## Build

```
cmake -S . -B build -G Ninja && cmake --build build && ./build/kutta-tests
```

- **CUDA:** add `-DKUTTA_ENABLE_CUDA=ON` (CUDA 13.0, RTX 20-series and newer, driver R580+).
- **HDF5:** Fluent HDF5 import needs HDF5 1.10.2 or newer.

For quick local checks (and CI), exclude the `slow` doctest suite:

```sh
./build/kutta-tests --test-suite-exclude=slow   # or: ctest --test-dir build -LE slow
```

The build command above runs **all** tests, including the V1/V2 refinement
studies and V3/V4 cavity validation. Slow gates can intentionally remain red
when a measured accuracy target is unmet; see the current
[verification results](src/kutta/fv/VERIFICATION.md).
The external-aero and transition slow tests skip unless their generated
meshes exist: `uv run tools/generate_fixtures.py` writes them under
`$KUTTA_RUNS` (default `<checkout>/runs`). See
[Generated fixtures](src/kutta/fv/VERIFICATION.md#generated-fixtures-and-iteration-budgets).

The SST defaults include a switchable signed
cross-diffusion cap at `1e6*rho*beta*omega^2`, inferred from Fluent's rpvar name,
and `k_floor=1e-14`. Wall output includes a finite implicit shear coefficient
that **replaces** the wall-face viscous flux, plus a separately documented
wall-face omega value for gradient/reconstruction. See the
[SST model contract and measured errors](src/kutta/fv/turbulence/README.md).

## Run a Fluent case

```sh
./build/kutta solve car.cas.h5 --data car.dat.h5 --backend cuda --iterations 2000 --out results
```

Omit `--data` for the case's selected standard or hybrid initialization. `--backend auto` (the
default) tries CUDA, Metal, then CPU. Unsupported settings are listed together
and cause a nonzero exit before solving. No case settings are silently replaced;
convergence-path changes already allowed by the adapter are printed.

`--coupled-edge-storage recompute` selects the memory-limited Coupled mode: it
recomputes five edge coefficients instead of caching them, saving device memory
at a throughput cost. It also shares strain scratch across velocity components,
including with SIMPLEC. The default, `cached`, keeps the ordinary fused vector fit.
Both modes preserve the same
float rounding and numerical outputs. This run option also applies to Kutta
case JSON inputs; it is recorded in `summary.json`.

The CLI prints Fluent-scaled residuals and average seconds per iteration at the
solver's monitor interval. It writes `result.vtu` (`u`, `p`, `k`, `omega`, `mu_t`),
`forces.csv` (wall zones and resolved force/drag/lift reports), `residuals.csv`,
and `summary.json` (SHA256 inputs, build/device identity and timing). A case
with a force, drag, lift or downforce report stops on its own once its
residual criteria are met and every report has stayed within 0.1 % over 300
iterations (residuals stalled above their criteria are waived once the forces
have held for 600). If the forces instead settle into a bounded cycle that
neither drifts, grows nor decays, as the hex Ahmed body's wake does on finer
meshes at Fluent's momentum relaxation 0.7, the run stops as `statistically_steady`, not converged, and
reports each force's mean and range ([default convergence on fine meshes](validation/results/default_convergence_fine/REPORT.md)).
Such cases run to at most `--iterations` (default 5000); `summary.json`
records `stop_reason` and what was met. Other cases run to `--iterations`
(default 1000). The criteria are editable
([Convergence criteria](docs/app-protocol.md#convergence-criteria-v1-additive)).
Peak device memory is recorded as `null`; the current VTU writer has no
boundary-face wall shear output.

Windows CI produces `kutta-windows-rtx.zip` with HDF5/zlib DLLs, license notices
and [Windows requirements and examples](docs/README-windows.txt).

## Experimental STL meshing

```sh
./build/kutta mesh car.stl --out car-mesh --case --progress jsonl
./build/kutta solve car-mesh/case.kutta.json --backend metal --out results
```

The automatic external-aero mesher castellates a balanced octree around binary
or named-solid ASCII STL, snaps the wall onto the surface and adds prism layers
(snappyHexMesh-style), keeping every cell at orthogonal quality 0.15 or above.
Mesh quality does not by itself establish aerodynamic accuracy or convergence.
See the [mesh setup and CLI contract](docs/app-protocol.md#automatic-stl-meshing-mesh-v1-additive)
and [measured results](validation/results/mesher/REPORT.md).

## Run without a Fluent case

A Kutta case is a mesh plus a small JSON setup that states only what differs
from Fluent's defaults ([docs/kutta-case.md](docs/kutta-case.md)):

```sh
./build/kutta polyhedra car.msh --out car-poly.msh               # tets to polyhedra, prism layers kept (tet meshes only)
./build/kutta new-case car-poly.msh --speed 20 --out car.kutta.json  # zone roles suggested from their names
./build/kutta solve car.kutta.json --backend metal --out results
./build/kutta convert old.cas.h5 --out old.kutta.json            # a Fluent case's settings as a Kutta case
```

Steady solves on tetrahedral regions oscillate instead of converging, so
convert a tet/prism mesh to polyhedra first, as Fluent recommends; hex,
poly-hexcore and polyhedral meshes solve as they are
([car mesh convergence](validation/results/car_mesh_convergence/REPORT.md)).
Conversion leaves about one polyhedron per five tetrahedra. Size the tet region
for the polyhedral cell count: on the Ahmed body, a polyhedral core needs about
as many cells as a hex core to give the same drag
([tet/poly accuracy](validation/results/tet_poly_accuracy/REPORT.md)).

The desktop app does the same with File ▸ New Case from Mesh and File ▸ Save
Case.

## Docs

- [Architecture and acceptance gates](docs/ARCHITECTURE.md)
- [Kutta cases: mesh plus JSON setup](docs/kutta-case.md)
- [Flow solver design](docs/fv-design.md)
- [Research](docs/research/): Fluent numerics and file formats, SST forensics, verification cases, NVIDIA build path

## Third-party data

- `validation/data/nasa_tmr`: NASA Turbulence Modeling Resource, CC0.
- `validation/fluent_examples`: Ansys example data, MIT (see `LICENSE-ansys-example-data`).
- Cavity tables: transcriptions of Ghia et al. (1982) and Erturk et al. (2005), with sources in each `SOURCE.md`.
- Some reference data used locally (ERCOFTAC, UIUC, Southampton thesis) is not redistributed here. Its `SOURCE.md` files list where to get it.
