# V2 pressure-order diagnosis

## Current conclusion

The original investigation below is historical. The exact-face Gauss probe
measures first-order local cell-centred pressure-gradient error on Freudenthal
tets (finest order 0.999409); opposite-face cancellation gives 1.99791 on hexes.
The actual LSQ/Second Order momentum pressure operator has tet order 1.01406.
Those local operator measurements do not bound the global pressure-solution
order. The earlier claim that second-order pressure was unattainable was not
established by that experiment.

The 2026-11-03 investigation recovers the unchanged V2 tet load gates using
quadratic interior fluxes and face means, cubic Dirichlet-boundary gradients,
and normal-foot pressure observations. CPU and Metal give finest pressure/
viscous/total force orders 3.28712/2.14839/1.99352. These are measured slopes
on the frozen mesh sequence, not formal orders of the method. The momentum
pressure operator is unchanged, yet the solved pressure L2 order becomes
1.91159. Exact-field wall controls and separate diffusion/predictor ablations
identify defects beyond that operator. A full cubic interior passed V2 but
became nonfinite on a smooth tet channel; it is not retained.

See the [dated report](../../../validation/results/tet_prism_convergence/REPORT.md#2026-11-03--v2-tet-loads-and-the-remaining-native-mode)
for the fresh baseline, controls and full-tier qualification. Native SST
convergence and second-order car loads remain unestablished.

## Historical bounded investigation

The tetrahedral pressure gate remains failed, with its assertion `order >= 1.8`
unchanged. The user authorized continuing to V3/V4 after this bounded diagnosis.
These measurements use Apple M4 Pro, CPU and Metal, Release, source `dcf0e42`
plus diagnostic output in `tests/fv/test_navier_stokes.cpp` and
`tests/fv/test_pressure_gradient.cpp`. Both backends agree to printed precision.

| Study | Refinements | Final velocity order | Final pressure order |
| --- | --- | ---: | ---: |
| Extruded hex | 8, 16, 32 | 1.91567 | 1.82218 |
| Full 3D Freudenthal tet | 4, 8, 16 | 1.98262 | 0.942386 |

Velocity error includes all three components. Pressure errors are volume-weighted
and aligned by a constant gauge offset. The finest hex took 840 outer iterations;
the finest tet took 460. All scaled momentum defects were below 1e-7. The tet
provisional continuity L1 defect was 2.34695e-6 and corrected imbalance 1.27294e-6.

Removing boundary-adjacent cells does not restore second-order pressure. Interior
tet pressure L2 errors at 4/8/16 are 0.221108, 0.121516, 0.0670725. Full-domain
errors are 0.275283, 0.132470, 0.0689334. The maximum error is at a wall-adjacent
cell: 0.212151 at (0.984375, 0.515625, 0.03125) on the finest mesh. Excluding only
wall-adjacent cells gives 0.0668975 at that resolution. The extruded hex has no
cells without a symmetry-cap boundary; its wall-excluded error is printed
separately rather than describing an empty interior as zero.

Standalone LSQ gradients of the exact MMS pressure, with exact face values,
give tet L2 errors 0.334045, 0.168973, 0.0848384: final order 0.994008. They are
first-order consistent. The momentum pressure gradient is the Gauss sum of the
specified second-order face pressure; the Rhie–Chow reconstruction uses the LSQ
gradient. At 16^3 their L2 difference is 0.0267736. This is an observed difference,
not proof that it causes the pressure-order failure.

Two controlled experiments were rejected and removed:

- Projecting wall pressure extrapolation into the wall plane (zero normal
  pressure gradient, which the analytic MMS has) gave velocity order 1.83816 and
  pressure order 0.841405. It did not fix the gate and would change general wall
  behaviour.
- Reconstructing Rhie–Chow pressure jumps with the actual momentum Gauss gradient
  gave velocity order 1.80718 and pressure order 0.742126. It did not fix the gate.

Earlier exploratory wide-stencil mass interpolation and a gradient-difference
compensation also failed to recover the required order and are absent from the
retained code. No quadratic/wider diffusion operator or relaxed threshold was
introduced. Root cause is unresolved. A useful next investigation is the local
momentum/diffusion truncation error by tet orientation, separating smooth pressure
error from cell-shape oscillations; this is a hypothesis, not a diagnosis.

The final rerun at `5e527dc` adds a pressure-update cap for pseudo-time stability,
which changes the convergence path only. CPU and Metal again agree: hex velocity
and pressure orders are 1.91521 and 1.82043; tet orders are 1.98241 and 0.942386.
The finest runs take 860 and 520 iterations respectively. Finest tet interior
pressure L2 is 0.0670728, full-domain L2 is 0.0689337, and maximum error is 0.212136
at the symmetry-related wall cell (0.015625, 0.484375, 0.96875). The tetrahedral
failure is unchanged by this convergence-path adjustment.

Reproduce the retained diagnostics with:

```sh
cmake -S . -B build -G Ninja && cmake --build build
./build/kutta-tests --test-case='FV V2*'
./build/kutta-tests --test-case='FV MMS pressure*'
```
