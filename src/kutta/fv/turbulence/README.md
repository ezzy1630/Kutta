# Fluent SST

Experimental SA and realizable k-epsilon live in separate files; see
[their contracts and unresolved Fluent wall gates](car_models.md). They remain
rejected by the Fluent case adapter. The SST contract below is unchanged.

`sst.hpp` contains the caller sketch. Construct `SstModel(device, mesh, options)`
once, then call `update(inputs)` each outer iteration. Inputs borrow device SoA
buffers, including unlimited gradients of velocity, k and omega. The extra scalar
gradients are necessary for F1 and cross diffusion. Updates allocate nothing,
download nothing, and use eight deterministic gather/pointwise kernels (wall
omega is separate to stay within Metal's buffer-binding limit).

Sources and sinks are **per volume**, with `S = source - sink * phi`. Add
`V * source` to the transport RHS and `V * sink` to its diagonal. Destruction and
negative cross diffusion are implicit; positive cross diffusion is explicit.
By default `limit_cross_diffusion=true` clips **both signs** of the blended
transport term to `|D_omega| <= cross_diffusion_clip*rho*beta*omega^2`, with
`cross_diffusion_clip=1e6`. This formula is inferred from the case rpvar
`sst-limit-cross-diff-to-n-times-w-destruct 1000000.`, not established by saved
source fields. It does not clip F1's CDkw. Disabling it retains signed transport;
a representational saturation of ratios at +/-1e38 still prevents float overflow
in CDkw diagnostics and implicit sink coefficients at tiny omega.
Wall-face k diffusivity is exactly zero; the caller replaces wall-cell omega
rows by identity using `fixed_omega_mask` and `fixed_omega`. The model does not
modify the caller's k/omega buffers. Its documented floors apply when evaluating
coefficients (`k_floor=1e-14`, `omega_floor=1e-20`, matching Fluent min-limit
rpvars); the transport solver must enforce positivity after solving too.

`wall_shear` is stress **on the fluid**, antiparallel to relative tangential
velocity. Multiply it by face area to compare with `SV_WALL_SHEAR` (force/face).
The wall-function stress **replaces the wall-face viscous flux**.
`wall_shear_coefficient = tau_w/|U_t|` is evaluated without dividing by speed,
so it stays finite and positive at zero speed. Freeze this coefficient for each
momentum solve; add `area*C*(I-n*n^T)` to the cell's momentum block and
`area*C*(I-n*n^T)*U_wall` to its RHS. This acts on tangential velocity; applying
`area*C` indiscriminately to all three diagonal entries would damp normal motion.
`wall_tau` is its positive scalar magnitude. `wall_u_tau = sqrt(tau/rho)` is
the friction velocity, distinct from `wall_u_star` of TG 4-418.

Cell diffusivities blend **inverse** sigmas with F1. Internal face diffusivities
linearly interpolate the complete cell coefficients with `DeviceMesh::weight`.
Smooth walls use zero k diffusivity and molecular omega diffusivity; other
boundaries extrapolate the owner coefficient. Transport boundary conditions
still belong to the assembler.
For several wall faces in one cell, select the smallest normal distance, breaking
ties by face id. Fluent's rule for this case remains unknown.

`wall_face_omega` is a **gradient/reconstruction value**, distinct from
`fixed_omega` (the wall-cell row constraint). Following the likely form in
forensics §3.5, use a power blend of `6*nu/(.072*dy^2)` and
`2*u_star/(.3*.4187*dy)`, exponent 1.3. These are the documented fixed constants,
not the free fitted constants; the reported ~0.42% offset is retained. Its role
in Fluent remains inconclusive. It must not become a transport Dirichlet BC.
The value depends on velocity/k and wall geometry, not gradients: a preliminary
update can bootstrap the first omega reconstruction, then refresh it from each
nonlinear iterate. Non-wall entries are zero and must use their actual BC values.

## Formula sources

The local [numerics report](../../../../docs/research/fluent-numerics-and-files.md)
§A1–A2 is overridden by [saved-field forensics](../../../../docs/research/fluent-sst-forensics.md)
§7 wherever they differ.

| Formula in `sst.kern.h` | Source |
| --- | --- |
| Strain `sqrt(2 Sij Sij)` and rotation `sqrt(2 Wij Wij)` | TG §4.4.1; forensics §1–2 |
| F1 and signed CDkw; floor only inside F1 | TG 4-107..109, 4-117 |
| Signed transport cross-diffusion clip | Inferred from `sst-limit-cross-diff-to-n-times-w-destruct` rpvar |
| F2; strain-limited mu_t; viscosity ratio cap | TG 4-118..120; forensics §2, §6–7 |
| Inverse sigma blending, Gamma = mu + mu_t/sigma | TG 4-103..106 |
| Gk = mu_t S², optional S*Omega or Omega²; production cap | TG 4-81, 4-460..461; numerics §A1 option table |
| Gomega, computed 1994 gamma or 2003 (5/9, .44) | TG 4-110..113; numerics §A1 |
| beta/beta-star destruction, optional low-Re alpha-star/alpha/beta-star | TG 4-114..116 and §4.4.1 low-Re equations; alpha-star(0) = .072/3 |
| u-star, blended u-plus, stress, both y-plus definitions | TG 4-413..418; forensics §3.1–3.3, §7 |
| Harmonic wall production blend | TG 4-420..422 (analytic log-law derivative `1/(kappa*y+)`) |
| Fixed wall-cell omega, power blend | TG 4-423..426; forensics §3.4, §7 |
| Wall-face reconstruction omega, fixed-constant choice | Forensics §3.5; likely form, role inconclusive |
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
under `<checkout>/validation/fluent_examples/mixing_elbow/` and
skips with a message if either file or HDF5 support is absent. Its private HDF5
C reader checks section ids, overlap, field shapes, coverage and boundary-owner
mapping. It does not depend on another lane's solution reader or `upload_mesh()`.

Measured on Apple M4 Pro, CPU and Metal, 2026-10-29, 17,822 cells / 2,436 wall
faces, Fluent 22.2. Relative absolute errors below use prediction/saved - 1:

| Quantity | Median | p95 | Maximum | Asserted limits (median / p95 / max) |
| --- | --- | --- | --- | --- |
| mu_t inactive identity, zero gradient | 3.73e-8 | 1.06e-7 | 1.85e-7 | 2e-7 / 3e-7 / 5e-7 |
| mu_t, 6,756 limited cells, own LSQ | .04597 | .14284 | .30976 | .055 / .16 / .32 |
| mu_t, 11,066 inactive cells, own LSQ | 4.10e-8 | .03942 | .26827 | 5e-7 / .05 / unasserted |
| Wall-face reconstruction omega (chosen fixed constants) | .004187 | .004282 | .004319 | .0043 / .0044 / .0045 |
| Wall-cell omega | 5.81e-6 | 4.25e-5 | 2.932e-4 | 1e-5 / 5e-5 / 3.2e-4 |
| y-plus (u-star) | 5.03e-6 | 2.014e-4 | .001759 | 1e-5 / 2.5e-4 / .002 |
| y-plus (u-tau from our shear) | 4.76e-6 | 1.855e-4 | .001982 | 1.5e-5 / 2.5e-4 / .0021 |
| Wall shear force magnitude | 9.44e-6 | 3.702e-4 | .003961 | 1.5e-5 / 4e-4 / .004 |

Inverse-square LSQ uses the model's `wall_face_omega` at walls and saved boundary
values elsewhere. Compared with the saved-wall-SV_O reconstruction, omega
gradients and CDkw change in all 2,436 wall cells (gradient relative-change p95
.005943, max .006273 over all cells), while mu_t changes by exactly zero and no
limiter classifications change on either backend. F2 and mu_t do not depend on
the omega gradient. No tolerance was relaxed (limited mu_t maximum tightened
from 2 to .32); the ~0.42% face-omega offset is measured separately. Limiter
classification agrees in 90.30% of cells. In the 2,555 discriminating cells,
median strain/F2 limiter ratio is .999075; the rotation hypothesis gives .240746.
The finite LSQ scatter is explicitly allowed by forensics §2, not evidence of
1% transport/force agreement. Full shear-vector error has max .005435 due to
additional saved directional scatter (max `1+cos` = 1.12e-5); forensics reports
only its median direction cosine, not a full-vector error bound. Kernel vectors
agree with the independent double wall-law oracle to 2.98e-7 maximum relative
error. Analytic limits, all switches, blends, diffusion, wall selection, repeat
determinism and all output fields are compared on CPU and Metal (1e-5 relative
backend tolerance).

Repository-wide `cmake -S . -B build -G Ninja && cmake --build build &&
./build/kutta-tests` passes all **85/85 test cases**, including the compiled-GPU
availability gate (CPU and Apple M4 Pro Metal both run). SST coverage is 14 cases
plus that GPU gate: signed physical clipping on/off, the original overflow probe,
zero/floored inputs, low-Re/production switches, a double-precision F1
cross-diffusion-branch and transport oracle, implicit wall stress at zero speed,
wall-face reconstruction omega, structural zero k flux, and backend determinism.
The oracle distinguishes Menter's .856 multiplier from Fluent's rounded 1/1.168.
No CUDA runtime or converged-flow claim follows from these checks.

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
- Low-Re mode, multiple-wall-face rule, wall-limiter exemption and u-star limiter:
  not established by this fixture. The optional cross-diffusion clip above is
  an explicit inference.
- Fluent's exact LSQ/WFGC and `kw-sst-strain? #t` semantics remain unknown.
  The explicit strain/vorticity enum does not pretend to decode that flag.
- No converged flow solve, force-integral acceptance, CUDA runtime or performance
  claim is established by these constitutive saved-field tests.
