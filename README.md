# Kutta

A local, GPU-native CFD solver for car aerodynamics. Its first goal is to reproduce Ansys Fluent's steady incompressible k-ω SST results on the same mesh, within 1–2% on downforce and drag, much faster on one consumer GPU. The first targets are an NVIDIA RTX GPU on Windows and Apple Silicon on macOS.

**Status: early development.** Done so far:

- device runtime: CPU reference, Metal, CUDA;
- mesh I/O: Fluent `.msh` / `.cas` / `.cas.h5`, NASA Plot3D, VTU;
- linear solvers, including aggregation AMG;
- Fluent-style SST constitutive kernels and smooth-wall closure, verified against
  the saved mixing-elbow fields on CPU and Metal.

The flow solver is in progress; the SST transport/momentum assembly is not yet
integrated. No accuracy or speed claims are made until they are measured; see the gates in [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Build

```
cmake -S . -B build -G Ninja && cmake --build build && ./build/kutta-tests
```

- **CUDA:** add `-DKUTTA_ENABLE_CUDA=ON` (CUDA 13.0, RTX 20-series and newer, driver R580+).
- **HDF5:** Fluent HDF5 import needs HDF5 1.10.2 or newer.

The suite currently has **85 test cases** (all passing on CPU/Metal, including
the external Fluent fixture). The SST defaults include a switchable signed
cross-diffusion cap at `1e6*rho*beta*omega^2`, inferred from Fluent's rpvar name,
and `k_floor=1e-14`. Wall output includes a finite implicit shear coefficient
that **replaces** the wall-face viscous flux, plus a separately documented
wall-face omega value for gradient/reconstruction. See the
[SST model contract and measured errors](src/kutta/fv/turbulence/README.md).

## Docs

- [Architecture and acceptance gates](docs/ARCHITECTURE.md)
- [Flow solver design](docs/fv-design.md)
- [Research](docs/research/): Fluent numerics and file formats, SST forensics, verification cases, NVIDIA build path

## Third-party data

- `validation/data/nasa_tmr`: NASA Turbulence Modeling Resource, CC0.
- `validation/fluent_examples`: Ansys example data, MIT (see `LICENSE-ansys-example-data`).
- Cavity tables: transcriptions of Ghia et al. (1982) and Erturk et al. (2005), with sources in each `SOURCE.md`.
- Some reference data used locally (ERCOFTAC, UIUC, Southampton thesis) is not redistributed here. Its `SOURCE.md` files list where to get it.
