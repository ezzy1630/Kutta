# FSAE physics implementation report

> V7.2 follow-up: the unchanged gate now passes on CPU and Metal (0.093%
> poly pressure error), and native rotating walls are enabled. See the
> [operator and regression evidence](../poly_pressure/REPORT.md). V12 remains held.
> The original measurements below describe the source identified in this report.

At its recorded source, this was a **partial implementation with two held acceptance gates**. Translating
walls, superficial porous resistance, absolute-velocity MRF, pressure inlets,
uniform mass-flow inlets and porous jumps reach the strict setup adapter.
Rotating walls and polynomial fans have C++ solver paths, but native rotating
walls remain rejected by V7.2 and native fan import remains rejected by V12.
Pressure far-field is rejected with Fluent's incompressible-applicability reason.

The [formulation, input evidence, unsupported variants and app settings](../../../docs/research/fsae-physics.md)
are part of this report. No Fluent solution agreement is claimed.

## Machine and source identity

- Checkout: `<checkout>`.
- Branch: `t3/fluent-fsae-physics`; base commit
  `a0c5caa6079cdd7ccd299ac093d9c9b18373c989`.
- Apple M4 Pro, macOS arm64, 24 GB; AppleClang 21.0.0.21000101.
- Release build; HDF5 and Metal enabled; CUDA disabled.
- [`metadata.json`](metadata.json) records the exact source-file hashes,
  generated kernels and tested executable SHA256. The local commit containing
  this report binds that source tree; nothing was pushed.

## Measured gates

All numeric rows run CPU and Metal. Values below are the larger observed errors
across those backends/configurations unless separate mesh values are shown.
Velocity errors are m/s. Rotational pressure L2 removes only a constant gauge
and is normalized by the analytic pressure change between radii 1 and 2.

| Gate | Required bound | Measured | Result |
| --- | --- | --- | --- |
| Quick tier | `ctest -LE slow` green, including compiled Metal presence and portability lint | 175 doctest cases; 1,540,118 assertions; CTest 2/2 | PASS |
| V7.1 translating Couette | Velocity L2 < `2e-4`; momentum < `1e-5`; CPU/Metal field difference < `2e-5` | Hex `6.87511e-8`; Voronoi poly `1.23545e-5`; momentum < `2.75e-8` | PASS |
| V7.2 rotating Taylor–Couette | Velocity L2 < `.01`; pressure relative L2 < `.01`; momentum < `1e-5` | Hex velocity `.00452144`, pressure `.00309582`; refined poly velocity `.000324131`, pressure **`.0169127`**; poly momentum `9.26867e-8` | **FAIL pressure, CPU and Metal** |
| V7.3 SST wall/frame translation | Velocity difference < `2e-4`; k difference < `2e-7`; relative omega difference < `2e-4` | `3.27826e-7`, `3.63798e-12`, `1.21774e-7` | PASS |
| V8 Darcy–Forchheimer | Pressure Linf < `.007 Pa`; velocity L2 < `1e-4`; momentum < `1e-5` | Expected drop `7 Pa`; pressure error `1.90735e-5 Pa`; velocity `1.36624e-7` | PASS |
| V9 absolute MRF | Velocity L2 < `.01`; pressure relative L2 < `.01`; frame velocity difference < `.003`; momentum < `1e-5` | Solid-body velocity `.00010236`, pressure `.000730474`; Taylor–Couette velocity `.0045225`, pressure `.00349889`; frame difference `2.5034e-6` | PASS |
| V9.2 conformal rotating/stationary interface | Poly Taylor–Couette velocity L2 < `.01`; momentum < `1e-5` | `.00523233`; momentum `4.56267e-8` | PASS, **velocity only** |
| V9.3 rotating + translating frame | Velocity L2 < `.01`; pressure relative L2 < `.01`; momentum < `1e-5` | `.00452067`; `.00303022` | PASS |
| V10 pressure/mass inlets | Velocity L2 < `1e-4`; pressure Linf < `1e-4 Pa`; corrected imbalance < `1e-6 kg/s` | Mass inlet `2.58141e-7`; pressure inlet `1.92284e-6` | PASS |
| V10.2 SST inlet equivalence | Relative k and omega differences from velocity-inlet transport < `1e-4` | `1.48651e-6`, `1.11959e-6` | PASS |
| V11 fan/porous-jump discontinuity | Pressure Linf < `.003 Pa`; velocity L2 < `1e-4`; momentum < `1e-5`; CPU/Metal pressure difference < `1e-4` | Signed steps `+3`, `-3`, `-2.8`, `+2.8 Pa`; max pressure error `9.01222e-6 Pa`; velocity `2.31354e-7` | PASS, C++ fan API |
| V12 native polynomial fan input | Saved native fan case and strict import acceptance | `KUTTA_FSAE_FAN_CASE` absent; native polynomial encoding not implemented | **FAIL, held** |

V7.2 uses 512 hex cells and 6,144 poly cells, 1,200 outer iterations and an
analytic seed. The unchanged residual check confirms relaxation to the discrete
fixed point rather than merely preserving that seed. The 384-cell poly
diagnostic had 9.08% pressure error; refinement reduces it to 1.69%, but does
not meet the 1% acceptance gate. No gate was relaxed to accept this result.

V8 uses 24 x 3 x 1 cells, flow speed 1, rho=2, mu=.1, length=1,
principal inverse permeabilities `(30,200,400)`, C2 `(4,20,40)`, aligned and
0.6-radian directions, and porosities .35 and 1. The second supplied direction
is intentionally nonorthogonal, exercising Fluent's orthogonalization rule.

V9 compares frame angular speeds 0, .4 and 1 rad/s on the same 512-cell hex
annulus. V9.2 couples an inner rotating zone to an outer stationary zone on a
conformal poly annulus; its pressure accuracy is **not** accepted by that gate.
V9.3 adds frame translation `(.7,-.2,0)` with angular speed .4 rad/s.

The explicit CPU/Metal source oracle covers anisotropic off-axis tensors,
relative porous resistance and simultaneous MRF translation/rotation. Rendering
Fluent TG Eq. 2-10 exposed the required `Omega cross (U-V_translation)` term.
The corrected oracle first produced 12 failing assertions with the omitted
translation term, then passed after correction; see
[`mrf_translation_before.log`](mrf_translation_before.log) and
[`mrf_translation_after.log`](mrf_translation_after.log).
The source equation SVG SHA256 was
`82f337498d85d802ba8b1a3581fab9a5a845c9576873e521ab1c6dc13785bf26`.

## Reproduction

```sh
cmake -S . -B build -G Ninja
cmake --build build -j 6
KUTTA_FLUENT_EXAMPLES=<checkout>/validation/fluent_examples \
  ctest --test-dir build -LE slow -V
./build/kutta-tests --test-case='FSAE*' --test-suite=slow
```

The quick result is in [`quick_final.log`](quick_final.log). At the originally recorded source, the scoped slow
command intentionally returned nonzero for V7.2 and V12; its complete output is
[`slow_acceptance.log`](slow_acceptance.log): **11 cases, 9 passed, 2 failed;
7,207 assertions, 7,204 passed, 3 failed**. The failures are exactly the two
polyhedral pressure assertions and V12's absent native fan fixture.
V1–V6 were not rerun or changed.
Earlier coarse rotational diagnostics are preserved in
[`rotation_pressure.log`](rotation_pressure.log).

## Remaining work and integration boundary

1. **Resolved in the follow-up:** V7.2 passes unchanged on CPU and Metal and
   native rotating walls are enabled. The linked reconstruction report records
   the before/after errors, operator orders and remaining validation limits.
2. Obtain a saved Fluent case containing a **polynomial fan thread**, establish
   its active encoding and sign conventions, implement strict adaptation, and
   pass V12. An MRF rotor case is not a substitute. The mathematical fan path
   already passes V11 and bounded-curve CPU/Metal checks.
3. Native mass-flow/porous-jump mappings have synthetic parser/adapter proof and
   an independent parser reference; they still lack a saved-case field comparison.
   Nested frame composition has unit proof, not independent Fluent numerical proof.
4. No CUDA GPU was available. No new feature was compared with Fluent's converged
   fields, forces or FSAE-car answers. G0 CUDA, G2 full SST acceptance, G3 same-mesh
   Fluent agreement and G4 speed remain unverified by this lane.
5. Keep the additive host-mesh `internal_face_zones` contract when integrating.
   Source normals are tested through import and RCM. No app, linalg, runtime or
   turbulence-model implementation was edited. The app-editable settings are
   listed in the linked formulation document.
