# FSAE wall, porous and reference-frame physics

This lane adds device-native physics and strict setup adaptation. It does not
establish Fluent same-mesh solution agreement (G3), CUDA execution parity, or
full FSAE-car acceptance. Native Fluent **rotating walls and fan polynomial import
remain rejected**: V7.2 retains a polyhedral pressure-accuracy failure, and V12
lacks verified native fan encoding.
The independently callable fan solver is tested against an analytic pressure step.
V12 deliberately fails until a saved native polynomial-fan case is supplied and
the importer implements its verified encoding.

The rotating-wall C++ path passes velocity on hex and poly meshes, and pressure
on hex. On the refined 32-ring/256-sector poly mesh its converged pressure L2
error is 1.69% of the analytic radial pressure range, above the unchanged 1%
limit. The 8-ring/64-sector diagnostic gives 9.08%. These failures are retained;
native rotating-wall support is not enabled on the strength of velocity alone.

## Equations and conventions

The [Fluent User's Guide, walls](https://ansyshelp.ansys.com/public/views/secured/corp/v261/en/flu_ug/flu_ug_bcs_sec_bound_cond.html)
defines translation and rotation relative to either the laboratory or adjacent
zone. Kutta evaluates `U_wall = translation + Omega cross (x-origin)` at each
face, adds adjacent-frame velocity for relative input, and removes only the
normal component of motion relative to that frame. Angular speeds are SI rad/s;
axis vectors are normalized. The existing SST shear closure receives the
resulting absolute wall velocity, so its cell-to-wall tangential difference is
unchanged by a frame translation. Roughness and specified shear remain rejected.

The [porous-medium formulation](https://ansyshelp.ansys.com/public/views/secured/corp/v261/en/flu_ug/flu_ug_bcs_sec_cell_zones.html)
adds `S = -(mu D + rho |v| C/2) v`, where D and C are tensors constructed from
the three principal resistances. Direction 1 is normalized; direction 2 has its
projection onto direction 1 removed and is normalized; direction 3 is their
cross product. Molecular viscosity is used in the Darcy term. Resistance can
use absolute velocity or velocity relative to the porous zone. Porosity is
validated and retained; it does not divide superficial velocity or multiply
steady convection/diffusion. Standard SST transport remains active in porous
zones. Physical-velocity, alternative inertial, power-law, non-unit relative
viscosity, conical/curvilinear, and laminar-zone variants remain rejected.

The [absolute-velocity MRF equations](https://ansyshelp.ansys.com/public/views/secured/corp/v261/en/flu_th/flu_th_sec_rotate_equations.html)
transport absolute momentum with relative mass flux and add
`-rho Omega cross (U_absolute - V_translation)` to the momentum RHS (Eq. 2-10,
checked by rendering its SVG). There is no additional
`-2 rho Omega cross U` term or separately added centrifugal force in this
formulation. Nested frames resolve parent velocity at the child origin, add
translation, and add angular-velocity vectors. Cycles and missing parents are
rejected. The global origin convention is read explicitly for nested motion.
At conformal interfaces, each frame's normal speed must vanish, consistent with
the [MRF mesh restrictions](https://ansyshelp.ansys.com/public/views/secured/corp/v261/en/flu_ug/flu_ug_sec_mzrrf.html).
Absolute velocity and scalar turbulence fields pass across the interface.
Sliding meshes, nonconformal interfaces and the relative-velocity formulation
remain unsupported. Standard initialization honors `init-relative?`.

Pressure inlets use gauge total pressure and the incompressible inlet-plane
Bernoulli relation. On reversed flow they use total pressure as static pressure
and extrapolate transported values. Mass-flow input is converted once from
whole-patch kg/s or uniform kg/(m^2 s) into velocity, accounting for the supplied
Cartesian direction or each face's inward normal. Both use the existing SST
inflow specifications. Nonzero adjacent frame motion at pressure inlets is
rejected; moving-frame pressure-outlet/backflow conventions are also rejected.
Fluent pressure far-field is explicitly rejected because its
Mach/characteristic formulation requires compressible flow.

Interior jumps use distinct pressure traces on the two sides of each face.
Their prescribed discontinuity is removed from pressure-gradient observations
and included in momentum and Rhie-Chow fluxes. The implicit slope changes the
iteration path while retaining the specified pressure-jump fixed point.
The C++ fan curve is `sum(c_j |v_normal|^j)` Pa along the retained source normal,
with reversal and optional polynomial velocity bounds. Porous jumps use
`Delta p = -(mu/alpha + rho C2 |v_normal|/2) v_normal thickness`.
They add momentum resistance and no turbulence-generation term. Native fan
swirl, averaged fan curves, profiles and UDFs remain unsupported.

## Input evidence and strictness

The native MRF and pressure-inlet names were inspected in Ansys' public axial-fan
case, SHA256 `70f7f73bdcf8eb0975b23e607b7eab89b2f4ae2d24bb227301bf8e9adc070444`.
The porous-zone names and active numeric values were inspected in its catalytic
converter case, SHA256
`55ed2521921f8f5b20b9fdfc6801cb7221ce4d72ea09c388f3d3c8a48f4e3117`.
Both come from the MIT-licensed
[Ansys example-data snapshot](https://github.com/ansys/example-data/tree/9674c61c5a3130ffe979561d683348e13575d472).
They are input-convention evidence; neither case was solved or compared with
saved Fluent flow fields in this lane.

The `flow-spec` mass-flow selector (`0` rate, `1` flux, `2` averaged flux),
direction selectors, and porous-jump keys `alpha`, `dm`, `c2` also have an
independent parser reference in
[cffview's boundary definitions](https://github.com/aiwaonly/cffview/blob/main/src/cffview/boundary.py).
These mappings have synthetic parsing and strict-adapter tests, **not** an
independent native mass-flow/porous-jump saved-case comparison. Averaged/profile
mass flux is rejected. Unknown switches remain rejected even when false.

The adapter integration tests reparse modified thread settings in the public
mixing-elbow case; they do not mistake those synthetic modifications for new
Fluent reference solutions. Its existing WFGC opt-in is explicit in that test.
The parser accepts modern vector encodings without demanding inactive scalar
components. Missing active constants or unresolved expressions are rejected.

## Ownership and storage

New implementation lives in `fv/{zone_physics,boundary_physics,interior_jump}.*`,
`setup/fluent_physics.*`, and three new single-source kernel files. Delimited
hooks in `steady_solver.*` handle assembly, flux, initialization and boundaries.
No turbulence-model, linalg, runtime or app implementation was changed.

An additive `PolyMesh::internal_face_zones` field retains special internal zone
IDs and source normals. Fluent import constructs it, and RCM remaps face IDs
and reverses directions when owner/neighbour swap. Ordinary interior faces do
not acquire duplicate metadata. This is a host-mesh contract addition for the
integrating agent to retain.

Resistance tensors are stored once per zone with one cell-to-zone index field;
moving-porous resistance adds a cell reference-velocity field when required.
Sheet coefficients, jumps and slopes scale with active faces, and gathered
updates launch only for affected cells/faces. The host builds injective active
index lists, so these kernels have no atomics or write races. Disabled volume
and jump features allocate no corresponding device fields; boundary override
metadata costs one integer and three floats per boundary face. No field downloads or host field loops were
added to `iterate()`.

## App settings for the app lane

- Wall mode, translation, angular speed, origin, axis and absolute/relative input.
- Porous enable, porosity, three inverse permeabilities, three C2 values,
  two direction vectors and absolute/relative resistance velocity.
- MRF enable, translation, angular speed, origin, axis, parent zone and the
  nested-origin convention; global absolute-velocity formulation and relative
  standard initialization.
- Pressure-inlet gauge total pressure; mass-flow rate versus uniform flux;
  normal/Cartesian direction; the existing SST turbulence specifications.
- Interior porous-jump permeability, thickness and C2. Fan coefficients,
  direction reversal and velocity bounds are C++-only until V12 is satisfied.
- Pressure far-field should display the incompressible rejection reason.

Measured gates and artifact identity are in
[`validation/results/fsae_physics/REPORT.md`](../../validation/results/fsae_physics/REPORT.md).
