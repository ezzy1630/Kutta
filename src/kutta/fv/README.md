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
inlets and pressure-outlet backflow accept uniform values or host profiles;
moving-wall normal motion is projected away. Symmetry removes normal velocity
and extrapolates tangentially. Pressure outlets impose gauge static pressure.
Wall force is pressure relative to the requested reference, minus the viscous
traction on the fluid, integrated using Device double/compensated reductions.

The SIMPLEC correction uses the relaxed/pseudo-time consistent inverse
`1/(aP_relaxed + sum(offdiag))`. Rhie–Chow uses `V/aP_physical` without
relaxation or pseudo-time, and momentum-diagonal-weighted face velocities.
Pressure face values use the second-order reconstruction in the design spec.
Momentum-weighted face velocities include an affine-exact correction from the
weighted cell interpolation point to the actual face centroid. Symmetry
pressure extrapolation uses only the tangential centroid-to-face displacement.
The undocumented Fluent high-order Rhie–Chow correction and warped-face
gradient correction are not implemented or claimed as matched.

Closed domains use one reference cell (default RCM cell 0). Both its pressure
correction row and column are eliminated, preserving symmetry for AMG-PCG.
Pressure outlets remove the need for a pin. Disconnected domains are rejected;
one pin would leave additional unconstrained pressure modes. Solver kinds are
centralized in `transport.hpp` and use the current `solver.hpp` factory.

Automatic pseudo-time uses `.3 * min(cbrt(mean volume), max domain extent) /
max(domain speed, .1*boundary speed, .1*sqrt(outlet pressure range/rho))`, times
the supplied scale. Speed is refreshed at monitor intervals. A fixed positive
time step overrides this. The flow loop has no field downloads. The first five
steps collect continuity normalization; subsequent reductions occur only at
the configured monitor interval or an explicit `monitor()`/`forces()` call.

Momentum residuals are steady equation L1 defects, reassembled from the current
state, scaled by `sum(aP_physical*|velocity|)`. Continuity is the provisional
mass-creation L1 defect before pressure correction, normalized by its maximum
over iterations 1–5. The corrected mass imbalance is reported separately.
This makes the pressure-correction linear tolerance distinct from outer
convergence. Zero denominators preserve unknown/infinite ratios rather than
reporting a false zero.
