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
normal motion is projected away. Symmetry uses the tangential projection of
cell velocity and exerts no net viscous force. Its normal component coefficient
is implicit, balanced by the same explicit term. Outlets convect cell velocity;
their viscous flux retains the transpose stress. Pressure outlets impose gauge
static pressure.
Wall force is pressure relative to the requested reference, minus the viscous
traction on the fluid, integrated using Device double/compensated reductions.

`boundary_least_squares()` precomputes each field's observation matrix for
laminar and SST flow alike. Pressure includes wall, inlet and symmetry rows
with zero differences (`p_b=p_c`) and prescribed outlet pressure. These LSQ
observations are separate from momentum's extrapolated boundary pressure
`p_c+G·r_b`, including at symmetry. Velocity includes wall, inlet and projected
symmetry data, and omits outlet rows. k includes inlet rows; omega additionally
includes independent wall reconstruction values. Excluded rows have zero
weighted deltas. A symmetric eigensolve constrains only deficient directions,
preserving all resolved directions. No previous gradient feeds itself through
an extrapolated boundary observation. Standalone operators still accept the
caller's independent boundary face data.

SST strain and production use one derivative gradient: plain LSQ `1/d²` without
wall rows plus half a face-jump fit. The fit uses centroid-to-face rows weighted
by `1/|r_f|²`, and jumps between the two unlimited side reconstructions (boundary
value minus the owner reconstruction at boundaries). This is one pass, with its
geometry matrix precomputed. It is a likely Fluent approximation, not the exact
undocumented WFGC algorithm. Convection and diffusion retain plain LSQ; SOU uses
Barth–Jespersen separately per component, with boundary values in the min/max
and evaluations at face centroids.

Set `SolverSettings::sst` and supply inlet/outlet `TurbulenceCondition` values
to enable SST. Momentum uses interpolated `mu+mu_t`, explicit transpose stress,
and zero wall-face viscosity. A gathered component diagonal applies
`A*C*(I-n*n^T)` with its cross terms explicit. The rotationally invariant
trace/3 of the three **physical** component diagonals drives Rhie-Chow and its
face interpolation; the mean component SIMPLEC inverse drives correction.
Neither physical diagonal includes relaxation or pseudo-time. Forces use
`-A*wall_shear` on the body. k/omega transport adds volume-weighted implicit
sinks (Newton `2*rho*beta*omega` for omega with matching
`+rho*beta*omega²` on the RHS), clips solved fields at configured floors, and imposes wall-cell omega
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
Momentum-weighted face velocities average the two limited side reconstructions
with physical a_P weights. Outlet flux predictors use cell velocity. The
Rhie–Chow pressure jump uses momentum's own gradient, Green–Gauss of those same
second-order face pressures. The exact Fluent high-order Rhie–Chow correction
and complete warped-face gradient correction remain open.

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
state. `ResidualValue::scaled` retains the legacy denominator
`sum(aP_physical*|velocity|)` for existing convergence criteria.
`fluent_scaled` uses the per-component `sum(2*(aP_i+rho*V/dt)*|velocity|)`;
`fluent_normalization` exposes that denominator. k/omega use their physical
transport diagonals (including implicit sinks) plus `rho*V/dt`; fixed omega
rows are excluded. The omega sink is Newton-linearized. These diagnostics use
the run's actual pseudo time and flux; comparison with Fluent's saved scale
also depends on its saved pseudo time and convecting flux. No energy equation
is currently solved by this driver.

Continuity is the provisional mass-creation L1 defect before pressure
correction, normalized by its maximum over iterations 1–5. The corrected mass
imbalance is reported separately. This distinguishes linear pressure tolerance
from outer convergence. Zero denominators preserve unknown/infinite ratios
rather than reporting a false zero.
