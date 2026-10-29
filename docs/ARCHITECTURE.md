# Kutta architecture

Kutta is a local, GPU-native CFD solver. Its first job is FSAE external aero: reproduce Ansys Fluent's steady incompressible k-ω SST results on the **same mesh** within 1–2% on downforce and drag, and do it much faster on one consumer GPU.

## Acceptance gates

These are targets until measured. Every claim in a report must name the machine, commit, mesh and case it was measured on.

| Gate | Pass condition |
| --- | --- |
| G0 Runtime | Every kernel matches the CPU reference on Metal and CUDA; STREAM-style bandwidth ≥ 80% of the device's peak |
| G1 Laminar | Manufactured solutions show 2nd-order convergence on hex, tet and polyhedral meshes; lid-driven cavity matches Ghia et al. |
| G2 SST | NASA TMR flat plate, bump and NACA 0012 SST results match CFL3D/FUN3D within the TMR grid-convergence band |
| G3 Fluent same-mesh | On a Fluent case's own mesh and settings, ClA and CdA within 1–2% and aero balance within 0.5 percentage points of Fluent's converged result |
| G4 Speed | ≥ 5× faster wall-clock than the team's Fluent CPU run, on the team's own Windows/RTX machine |

## Principles

1. **Same discrete problem as Fluent first.** Read the team's mesh and settings, then match Fluent's discretization choices. That way any difference is our error, and we can measure it.
2. **Memory bandwidth is the budget.** Steady RANS is bandwidth-bound. Keep the whole solve on the GPU, minimize bytes per cell per iteration, and measure against the bandwidth benchmark.
3. **The CPU backend is the oracle.** It is serial and accumulates in double. Every GPU kernel is tested against it on the same inputs.
4. **Deterministic.** Kernels gather and never scatter, and there are no atomics, so a run is bitwise reproducible on a given device.
5. **One kernel source.** Kernels are written once in a restricted C++ dialect (`src/kutta/kernels/prelude.h`) and compiled as C++, CUDA and Metal. `tools/kernelgen.py` generates each backend's entry points and typed host launchers (`kutta::kern::name(dev, n, ...)`), so arguments can never drift out of sync.

## Precision rules

Metal has no FP64, and consumer NVIDIA FP64 runs at 1/64 rate. So:

- Field storage and arithmetic are `float`.
- Geometry is computed on the host in `double` (`MeshGeometry`). The device stores only *differences*: face-to-cell vectors, cell-to-cell vectors, face area vectors and volumes. It never stores absolute coordinates, so float keeps full relative precision even far from the origin.
- Pressure is stored as a gauge value relative to a reference, never as absolute operating pressure.
- Reductions (dot products, norms, force integrals) accumulate in double on CPU/CUDA and in compensated float pairs on Metal. They always return `double`.
- The outer SIMPLE iteration is a defect correction: residuals are recomputed from the current state each iteration, so inner linear solves can stay in float without limiting final accuracy.

## Layers

```
app/       kutta CLI (later: Python API, UI process)
fv/        finite-volume discretization, SIMPLEC, k-ω SST, forces, monitors
linalg/    CSR on device, Krylov solvers, aggregation AMG
mesh/ io/  PolyMesh + MeshGeometry (host, double); Fluent msh / CFF h5, Plot3D readers; VTU writer
runtime/   Device: buffers, kernel launch, reductions. Backends: cpu (oracle), metal, cuda
kernels/   *.kern.h single-source kernels + prelude
```

Contracts live in the headers: `runtime/device.hpp`, `kernels/prelude.h`, `mesh/mesh.hpp`, `io/mesh_io.hpp`, `linalg/csr.hpp` and `linalg/solver.hpp`. Change a contract only deliberately, and update every implementation and this document in the same change.

### Mesh layout

The mesh is face-based: owner/neighbour, with internal faces first and boundary faces grouped by patch (see `mesh.hpp` for the invariants). 2D meshes are extruded one cell thick, with Symmetry patches `front`/`back`. Assembly is per cell over a cell→face CSR list (gather), never per face with atomic scatter.

### Linear solvers

One CSR pattern per mesh, uploaded once. Only the values change each outer iteration. Pressure uses aggregation AMG (pairwise/plain aggregation, structure built once on the host, coarse values refreshed on the device through precomputed fine→coarse entry maps) preconditioning CG. Momentum and turbulence use BiCGStab or Jacobi sweeps.

### Planned flow solver (fv/, wave 2)

The flow solver is a cell-centred finite-volume incompressible solver using pressure-based SIMPLEC with Rhie–Chow face fluxes. Default schemes:

- least-squares cell-based gradients;
- second-order upwind by gradient reconstruction;
- non-orthogonal diffusion correction;
- k-ω SST as Fluent implements it, with Fluent's ω-based wall treatment.

The exact defaults to match are in `docs/research/fluent-numerics-and-files.md`.

## Build and test

```
cmake -S . -B build -G Ninja && cmake --build build && ./build/kutta-tests
```

Tests run every backend that is compiled in and present on the machine (`KUTTA_FOR_EACH_BACKEND`). CUDA builds use `-DKUTTA_ENABLE_CUDA=ON`, with fatbins for sm_75/86/89/120.

## Rules for agents working in this repo

- Stay inside the files your lane owns. To change a contract header, report the change you need instead of making it.
- Every new kernel gets a test that compares it against the CPU backend.
- Don't add GPL/AGPL dependencies. Permissive licenses (MIT/BSD/Apache, HDF5's BSD-style) are fine.
- Don't push or publish. Leave work committed on your branch with a clear message.
- Report measured numbers with machine, commit and command. Never round a failure into a pass.
