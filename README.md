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

The CLI prints Fluent-scaled residuals and average seconds per iteration at the
solver's monitor interval. It writes `result.vtu` (`u`, `p`, `k`, `omega`, `mu_t`),
`forces.csv` (wall zones and resolved force/drag/lift reports), `residuals.csv`,
and `summary.json` (SHA256 inputs, build/device identity and timing). The current
adapter exposes no convergence criteria, so runs stop at `--iterations` (default
1000) and record `converged=false`. Peak device memory is recorded as `null`;
the current VTU writer has no boundary-face wall shear output.

Windows CI produces `kutta-windows-rtx.zip` with HDF5/zlib DLLs, license notices
and [Windows requirements and examples](docs/README-windows.txt).

## Docs

- [Architecture and acceptance gates](docs/ARCHITECTURE.md)
- [Flow solver design](docs/fv-design.md)
- [Research](docs/research/): Fluent numerics and file formats, SST forensics, verification cases, NVIDIA build path

## Third-party data

- `validation/data/nasa_tmr`: NASA Turbulence Modeling Resource, CC0.
- `validation/fluent_examples`: Ansys example data, MIT (see `LICENSE-ansys-example-data`).
- Cavity tables: transcriptions of Ghia et al. (1982) and Erturk et al. (2005), with sources in each `SOURCE.md`.
- Some reference data used locally (ERCOFTAC, UIUC, Southampton thesis) is not redistributed here. Its `SOURCE.md` files list where to get it.
