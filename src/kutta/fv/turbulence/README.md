# Fluent SST

`sst.hpp` contains the caller sketch. Construct `SstModel(device, mesh, options)`
once, then call `update(inputs)` each outer iteration. Inputs borrow device SoA
buffers, including unlimited gradients of velocity, k and omega. The extra scalar
gradients are necessary for F1 and cross diffusion. Updates allocate nothing,
download nothing, and use seven deterministic gather/pointwise kernels.

Sources and sinks are **per volume**, with `S = source - sink * phi`. Add
`V * source` to the transport RHS and `V * sink` to its diagonal. Destruction and
negative cross diffusion are implicit; positive cross diffusion is explicit.
The caller enforces zero normal k flux on walls and replaces wall-cell omega
rows by identity using `fixed_omega_mask` and `fixed_omega`. The model does not
modify the caller's k/omega buffers. Its documented floors apply when evaluating
coefficients; the transport solver must enforce positivity after solving too.

`wall_shear` is stress **on the fluid**, antiparallel to relative tangential
velocity. Multiply it by face area to compare with `SV_WALL_SHEAR` (force/face).
`wall_tau` is its positive scalar magnitude. `wall_u_tau = sqrt(tau/rho)` is
the friction velocity, distinct from `wall_u_star` of TG 4-418.

Cell diffusivities blend **inverse** sigmas with F1. Internal face diffusivities
linearly interpolate the complete cell coefficients with `DeviceMesh::weight`.
Smooth walls use molecular viscosity; other boundaries extrapolate the owner
coefficient. Transport boundary conditions still belong to the assembler.
For several wall faces in one cell, select the smallest normal distance, breaking
ties by face id. Fluent's rule for this case remains unknown.

## Formula sources

The local [numerics report](../../../../docs/research/fluent-numerics-and-files.md)
§A1–A2 is overridden by [saved-field forensics](../../../../docs/research/fluent-sst-forensics.md)
§7 wherever they differ.

| Formula in `sst.kern.h` | Source |
| --- | --- |
| Strain `sqrt(2 Sij Sij)` and rotation `sqrt(2 Wij Wij)` | TG §4.4.1; forensics §1–2 |
| F1 and signed CDkw; floor only inside F1 | TG 4-107..109, 4-117 |
| F2; strain-limited mu_t; viscosity ratio cap | TG 4-118..120; forensics §2, §6–7 |
| Inverse sigma blending, Gamma = mu + mu_t/sigma | TG 4-103..106 |
| Gk = mu_t S², optional S*Omega or Omega²; production cap | TG 4-81, 4-460..461; numerics §A1 option table |
| Gomega, computed 1994 gamma or 2003 (5/9, .44) | TG 4-110..113; numerics §A1 |
| beta/beta-star destruction, optional low-Re alpha-star/alpha/beta-star | TG 4-114..116 and §4.4.1 low-Re equations; alpha-star(0) = .072/3 |
| u-star, blended u-plus, stress, both y-plus definitions | TG 4-413..418; forensics §3.1–3.3, §7 |
| Harmonic wall production blend | TG 4-420..422 (analytic log-law derivative `1/(kappa*y+)`) |
| Fixed wall-cell omega, power blend | TG 4-423..426; forensics §3.4, §7 |
| Wall k zero flux | TG §4.18.3; forensics §3.6 |

SST constants: a1=.31, beta-star=.09, beta1/2=.075/.0828,
sigma-k1/2=1.176/1, sigma-omega1/2=2/1.168, kappa=.41 for gamma.
**Wall-only** constants: beta=.072, kappa=.4187, E=9.793,
sqrt(Cmu)=.3, Ccalib=1/3, Cexp=1.3. Gamma is computed from these SST
constants, rather than the rounded .5532/.4403 quoted in the report.

The official [SST TG](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_th/flu_th_sec_turb_kw_sst.html)
and [BSL TG](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_th/flu_th_sec_turb_kw_bsl.html)
give the constitutive equations; the
[standard k-omega TG](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_th/flu_th_sec_turb_kw_std.html)
provides the optional low-Re coefficients.

## Verification

Run `cmake -S . -B build -G Ninja && cmake --build build`, then
`./build/kutta-tests --test-case='*SST*,compiled GPU*'`.
The saved-field test uses the exact external `.cas.h5`/`.dat.h5` paths requested
under `/Users/ezzy/.Others/Kutta/validation/fluent_examples/mixing_elbow/` and
skips with a message if either file or HDF5 support is absent. Its private HDF5
C reader checks section ids, overlap, field shapes, coverage and boundary-owner
mapping. It does not depend on another lane's solution reader or `upload_mesh()`.

Measured on Apple M4 Pro, CPU and Metal, 2026-10-29, 17,822 cells / 2,436 wall
faces, Fluent 22.2. Relative absolute errors below use prediction/saved - 1:

| Quantity | Median | p95 | Maximum | Asserted limits (median / p95 / max) |
| --- | --- | --- | --- | --- |
| mu_t inactive identity, zero gradient | 3.73e-8 | 1.06e-7 | 1.85e-7 | 2e-7 / 3e-7 / 5e-7 |
| mu_t, 6,756 limited cells, own LSQ | .04597 | .14284 | .30976 | .055 / .16 / 2 |
| Wall-cell omega | 5.81e-6 | 4.25e-5 | 2.932e-4 | 1e-5 / 5e-5 / 3.2e-4 |
| y-plus (u-star) | 5.03e-6 | 2.014e-4 | .001759 | 1e-5 / 2.5e-4 / .002 |
| y-plus (u-tau from our shear) | 4.72e-6 | 1.855e-4 | .001982 | 1.5e-5 / 2.5e-4 / .0021 |
| Wall shear force magnitude | 9.36e-6 | 3.702e-4 | .003961 | 1.5e-5 / 4e-4 / .004 |

Inverse-square LSQ includes saved boundary values at face centroids. Limiter
classification agrees in 90.30% of cells. In the 2,555 discriminating cells,
median strain/F2 limiter ratio is .999075; the rotation hypothesis gives .240746.
The finite LSQ scatter is explicitly allowed by forensics §2, not evidence of
1% transport/force agreement. Full shear-vector error has max .005435 due to
additional saved directional scatter (max `1+cos` = 1.12e-5); forensics reports
only its median direction cosine, not a full-vector error bound. Kernel vectors
agree with the independent double wall-law oracle to 3.35e-7 maximum relative
error. Analytic limits, all switches, blends, diffusion, wall selection, repeat
determinism and all output fields are compared on CPU and Metal (1e-5 relative
backend tolerance).

Repository-wide `ctest --test-dir build --output-on-failure` currently passes
69/70 cases. The unchanged `real Fluent external mesh acceptance` inventory
test rejects the extra external `LICENSE-ansys-example-data` file; no expected
fixture is missing. Excluding that one case passes the remaining 69. A separate
`KUTTA_ENABLE_HDF5=OFF` build passes the SST unit cases and prints the explicit
saved-field skip. Neither the inventory test nor external fixtures are owned
by Lane F.

## Integration dependency and unverified behavior

**Geometry dependency (resolved):** wall distances, y+ and wall omega match Fluent only with the face-pyramid centroid convention. Production `compute_geometry()` now uses that convention, and the saved-field test runs on production geometry: wall-omega max error 2.9e-4. With the previous warped-face tetrahedral moments the max error was 1.187%.

- Production limiter activity and whether it caps omega production: never
  active here (max production/cap .3595). Switches and synthetic tests exist.
- Wall-cell production blend and its ordering with the cap: not saved by
  Fluent. Here blend first, then apply the cap; disabling wall production
  retains volumetric production, while still returning the wall override.
- y-plus > 8: analytic limits tested, no saved-field validation. Fluent 2026 R1
  must also re-confirm the constants established on this 22.2 fixture.
- Tabulated wall omega: throws explicitly because its table is unpublished.
  Curvature and compressibility switches also throw; their required physics
  is outside this incompressible interface.
- Low-Re mode, multiple-wall-face rule, wall-limiter exemption, wall-face omega
  definition, undocumented cross-diffusion clipping and u-star limiter:
  not established by this fixture. No undocumented clip/limiter is applied.
- Fluent's exact LSQ/WFGC and `kw-sst-strain? #t` semantics remain unknown.
  The explicit strain/vorticity enum does not pretend to decode that flag.
- No converged flow solve, force-integral acceptance, CUDA runtime or performance
  claim is established by these constitutive saved-field tests.
