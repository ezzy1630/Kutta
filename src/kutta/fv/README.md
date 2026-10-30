# Finite volume

`upload_mesh()` preserves the caller's cell/face numbering. Call
`renumber_cells_rcm()` before computing geometry or fields when loading a case.
The shared `device_mesh.hpp` contract is unchanged. All device geometry is
stored as float differences computed from host double geometry.

Wall distances use a host double BVH over the same arithmetic-center triangle
fans as mesh geometry. The nearest point can lie in a triangle, on an edge or
at a vertex; an infinite wall plane is not substituted for the finite surface.
No wall patches means infinite wall distance. Non-wall normal distances are
zero. Warped polygons follow the mesh's triangle-fan surface convention.

`make_stencil()` builds the sorted cell-adjacency CSR pattern and incidence to
CSR maps once. The LSQ normal matrix uses weights `1/|d|^2`, with neighbour
centroids on internal faces and face centroids on boundaries. Its symmetric
inverse is stored in six SoA planes. `gradient()` takes actual boundary face
values; callers must supply extrapolated values for Neumann conditions.
`limit_gradient()` computes the cell-to-face Barth–Jespersen scalar separately
so diffusion and pressure can use the unlimited gradient.

Transport uses implicit first-order upwind with deferred second-order upwind
and over-relaxed non-orthogonal diffusion. Assembly gathers into each owned
CSR row, including multiple faces between a cell pair. The linear solver
policy is centralized in `transport.hpp`. All operator work remains on Device;
only the linear solver and explicitly requested monitors synchronize with the
host.

The flow driver shares one uploaded CSR matrix workspace across its sequential
equations. Each assembly overwrites the values and keeps each component's RHS
and reconstruction separately. A standalone `TransportEquation` owns its own
workspace unless given a shared matrix with the identical stencil pattern.
Velocity gradients gather all three components in one kernel; scalar transport
retains the scalar gradient entrypoint as an independent oracle.

`SteadySolver` owns an RCM-renumbered mesh and exposes `iterate(n)`, Device
velocity/pressure/mass-flux fields, and host geometry for initialization and
post-processing. Supply exactly one explicit condition per patch. Velocity
inlets and pressure-outlet backflow accept uniform values or host profiles.
Default outlet backflow is normal to the boundary, with its speed obtained from
the computed mass flux and density. Backflow data controls advected momentum;
the outlet flux predictor continues to use the momentum solution. Moving-wall
normal motion is projected away. Symmetry removes normal velocity
and extrapolates tangentially. Pressure outlets impose gauge static pressure.
Wall force is pressure relative to the requested reference, minus the viscous
traction on the fluid, integrated using Device double/compensated reductions.

For extrapolated boundary values, `boundary_least_squares()` precomputes a
separate observation matrix for each equation. Pressure excludes wall, inlet
and symmetry rows; velocity excludes outlet and symmetry rows; k excludes all
but inlet rows; omega additionally includes independent wall reconstruction
values. Excluded rows have zero weighted deltas. A symmetric eigensolve adds
zero-difference constraints only along deficient directions, preserving all
resolved directions. No previous gradient can feed itself through a boundary
value. Standalone operators with supplied independent face data retain the
original full LSQ stencil.

Set `SolverSettings::sst` and supply inlet/outlet `TurbulenceCondition` values
to enable SST. Momentum uses interpolated `mu+mu_t`, explicit transpose stress,
and zero wall-face viscosity. A gathered component diagonal applies
`A*C*(I-n*n^T)` with its cross terms explicit. The rotationally invariant
trace/3 of the three **physical** component diagonals drives Rhie-Chow and its
face interpolation; the mean component SIMPLEC inverse drives correction.
Neither physical diagonal includes relaxation or pseudo-time. Forces use
`-A*wall_shear` on the body. k/omega transport adds volume-weighted implicit
sinks, clips solved fields at configured floors, and imposes wall-cell omega
rows at the original row scale (unscaled identity rows harm inexact-solve
conditioning). Those rows are rebuilt each assembly and excluded from omega
residual scaling. `monitor()` reports both turbulence residuals without changing
the state.

`initialize_fields()` takes original-input cell arrays and preserves their exact
identity through the owned RCM permutation. It primes constitutive coefficients
and mass flux without taking an outer step. `permutation()` exposes both maps
for comparison with source data. `setup/fluent.hpp` supplies the strict case
adapter and source-ID field mapping. It rejects unsupported spatial settings
and initialization methods by name; see V6 for the current elbow blockers.
Total-pressure outlets apply `p=p0-rho*|u_cell|^2/2` only on backflow faces.

```cpp
auto dev = kutta::Device::create(kutta::Backend::Metal);
auto mesh = kutta::make_box_hex(129, 129, 1, {0, 0, 0}, {1, 1, 1./129});
std::vector<kutta::fv::PatchCondition> bc;
for (const auto& patch : mesh.patches) {
  if (patch.name == "zmin" || patch.name == "zmax")
    bc.push_back({patch.name, kutta::fv::Symmetry{}});
  else
    bc.push_back({patch.name, kutta::fv::Wall{{patch.name == "ymax" ? 1. : 0., 0, 0}}});
}
kutta::fv::SolverSettings settings;
settings.viscosity = 1./100; // unit lid speed and cavity length
kutta::fv::SteadySolver solver(*dev, std::move(mesh), std::move(bc), settings);
solver.iterate(100);
auto residuals = solver.monitor();
```

The device must outlive the solver. Initial fields default to zero; `initialize()`
and `set_source()` evaluate supplied host functions once on cell centroids.

The SIMPLEC correction uses the relaxed/pseudo-time consistent inverse
`1/(aP_relaxed + sum(offdiag))`. Rhie–Chow uses `V/aP_physical` without
relaxation or pseudo-time, and momentum-diagonal-weighted face velocities.
The pressure update is capped at `0.7*aP_physical*consistent_inverse` when
pseudo time makes the correction coefficient smaller than the Rhie–Chow
coefficient. This damps pressure overshoot without changing the steady equations.
Pressure face values use the second-order reconstruction in the design spec.
Momentum-weighted face velocities include an affine-exact correction from the
weighted cell interpolation point to the actual face centroid. Symmetry
pressure extrapolation uses only the tangential centroid-to-face displacement.
The undocumented Fluent high-order Rhie–Chow correction and warped-face
gradient correction are not implemented or claimed as matched.

Closed domains pin the cell containing `reference_point` (Fluent's reference
pressure location, default the origin; nearest centroid if no cell contains
it). Both its pressure
correction row and column are eliminated, preserving symmetry for AMG-PCG.
Pressure outlets remove the need for a pin. Disconnected domains are rejected;
one pin would leave additional unconstrained pressure modes. Solver kinds are
centralized in `transport.hpp` and use the current `solver.hpp` factory.

Automatic pseudo-time uses Fluent's conservative length scale,
`.3 * min(cbrt(domain volume), max domain extent) /
max(domain speed, .1*boundary speed, .1*sqrt(outlet pressure range/rho))`, times
the supplied scale. Speed is refreshed at monitor intervals. A fixed positive
time step overrides this. The flow loop has no field downloads. The first five
steps collect continuity normalization; subsequent reductions occur only at
the configured monitor interval or an explicit `monitor()`/`forces()` call.
`report().iteration` is the completed outer count; `residual_iteration` names
the residual snapshot and is empty before the first monitor. A caller can
therefore distinguish current iteration statistics from cached residuals.

Momentum residuals are steady equation L1 defects, reassembled from the current
state, scaled by `sum(aP_physical*|velocity|)`. Continuity is the provisional
mass-creation L1 defect before pressure correction, normalized by its maximum
over iterations 1–5. The corrected mass imbalance is reported separately.
This makes the pressure-correction linear tolerance distinct from outer
convergence. Zero denominators preserve unknown/infinite ratios rather than
reporting a false zero.
