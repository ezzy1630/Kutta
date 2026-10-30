# Experimental SA and realizable k-epsilon

These models are **not accepted Fluent case models**. The case allowlist still
rejects SA and every k-epsilon model. See the frozen [V7/V8 protocol](../../../../validation/results/car_turbulence/PROTOCOL.md)
and [measurements](../../../../validation/results/car_turbulence/README.md).
Standard/RNG have not been added because SA and realizable acceptance is a
prerequisite. No Fluent executable or converged same-mesh reference was used.

## Caller contract

Choose one of `SolverSettings::sst`, `sa`, or `realizable_ke`. Multiple choices
throw. New solver hooks dispatch into `sa_transport.*` and
`realizable_ke_transport.*`; the shared transport, linear algebra, runtime and
SST kernels are unchanged. All solve work uses Device plus generated kernels.
Host loops only construct options/BCs, initialize fields, or sample monitors.

```cpp
fv::SolverSettings settings;
settings.sa = fv::SaSettings{}; // explicitly resolved smooth walls
// Inlet/backflow: ModifiedViscosity{nu_tilde} or SaViscosityRatio{mu_t/mu}.
// IntensityLength and IntensityDiameter are also accepted for SA.
fv::SteadySolver solver(device, mesh, boundaries, settings);
solver.initialize_sa_fields(velocity, pressure, nu_tilde); // ORIGINAL cell order
solver.iterate(100);
auto residual = solver.monitor().nu_tilde;
```

For realizable k-epsilon, select `realizable_ke=RealizableKeSettings{}` and call
`initialize_ke_fields(velocity,pressure,k,epsilon)`. Use `KeValues{k,epsilon}`,
`IntensityLength`, `IntensityDiameter`, or `IntensityViscosityRatio` at inlets
and pressure-outlet backflow. `monitor().k` and `.epsilon` are separate from
SST's `.omega`; `dissipation()` returns epsilon. Hybrid turbulence initialization
throws for both new models. CLI/app exposure is intentionally absent.

Constitutive-only callers construct `SaModel` or `RealizableKeModel` once per
device/mesh and call `update()` with borrowed Device buffers. Updates allocate
and download nothing. Per-volume sources follow `source-sink*phi`; multiply
both by volume during assembly. Coefficient evaluation clips at model floors,
and transport clips again after solving, as SST does. No host field traversal
occurs during an update or transport solve.

## SA formulation

`SaFormulation::Fluent` implements TG 4-15..29: cb1=.1355, cb2=.622,
sigma=2/3, cv1=7.1, cw2=.3, cw3=2, kappa=.41, and derived cw1. It retains fv1
and fv2 low-Re damping, the gradient-squared source and nonlinear diffusion.
The default production rate is vorticity. `StrainVorticity` instead uses
`Omega+2*min(0,S-Omega)` before the fv2 correction; signed production/destruction
is split into nonnegative explicit source and implicit sink. This is not a
curvature correction or SA-neg. There is no trip source.

Fluent's documented production/destruction has **no ft2**. The explicit
`TmrStandard` variant includes ft2=1.2*exp(-.5*chi^2), its matching destruction
correction, the TMR Note 1(b) safeguard S-tilde>=.3*Omega (applied here to
both production and r), and r<=10. The reference codes' particular safeguarding
choice has not been independently established. V7 selects
this variant rather than presenting a no-ft2 solve as exact standard SA.
Both variants use finite-surface BVH nearest-wall distance, not grid-line distance.

Resolved walls set nu-tilde=0 with molecular face diffusivity mu/sigma and
molecular momentum stress. Internal diffusivity is (mu+rho*nu-tilde)/sigma.
Momentum uses mu+mu_t and transpose stress; wall-force reconstruction uses the
same molecular stress as the existing laminar path. Wall shear is stress on the
fluid, with friction-velocity y+. SA's viscosity ratio BC inverts chi*fv1=ratio;
it does not substitute the undamped nu-tilde/nu ratio. Intensity/length uses
the UG's explicitly undamped estimate Cmu^(1/4)*sqrt(1.5)*U*I*ell, Cmu=.09;
hydraulic diameter uses ell=.07*Dh. Inlets and outlet backflow prescribe the
transported scalar; outflow and symmetry use zero gradient.

Stable fv2 algebra retains its small positive high-chi limit. The r evaluation
is saturated at magnitude 1e4 outside the TMR cap to avoid overflowing its sixth
power; fw is evaluated in scaled form. This is a floating-point extension of
the displayed formula, not evidence for Fluent's singular-rate implementation.

`SaWallTreatment::FluentAutomatic` throws. The current TG states that the blend
is calibrated over the buffer layer but does not publish its complete algebra
or nu-tilde boundary/row treatment. It explicitly distinguishes that blend from
Kader's EWT law. Substituting SST's wall law would define a different model.

## Realizable k-epsilon formulation

The incompressible inertial-frame model implements TG 4-53..61 with variable
Cmu=1/(A0+As*k*Ustar/epsilon), A0=4.04. As uses the full symmetric strain tensor's
third invariant, and Ustar contains both strain and rotation. The extra moving
reference-frame correction is not implemented; inertial frames have zero frame
rotation. C2=1.9, sigma_k=1, sigma_epsilon=1.2. The epsilon production is
rho*C1*S*epsilon, C1=max(.43,eta/(eta+5)), and destruction is
rho*C2*epsilon^2/(k+sqrt(nu*epsilon)). It is not the standard k-epsilon equation.
The default high-Re k production uses (mu_t+mu)*S^2, including the documented
`coeff-ke-mu-prod` coefficient of 1. No production limiter, buoyancy,
compressibility, curvature, roughness, frame rotation or user source is enabled.

Scalable smooth walls use kappa=.4187, E=9.793, wall Cmu=.09 and
max(y-star,11.225) in all wall-function formulas. The effective distance implied
by that clipped y-star enters wall production and epsilon, too. k has zero wall
flux, and wall epsilon rows are constrained at their original row scale. For
multiple wall faces, the nearest normal distance wins, with face-id tie-break,
matching SST's deterministic policy; Fluent's corresponding policy is unverified.
Applying the clipped distance to epsilon/production is an interpretation of the
TG instruction to replace y-star in wall formulas; no saved Fluent fields have
established that interpretation. Fixed epsilon targets honor the configured
epsilon floor, so the row constraint and post-solve field limit agree.
Wall shear replaces wall-face momentum diffusion, using the existing implicit
tangential tensor hook and force integrator. k/epsilon use unlimited LSQ scalar
gradients with independent inlet observations. Velocity rates currently use the
existing plain LSQ momentum gradient, not an exact Fluent WFGC reconstruction.
As in the existing SSTm path, isotropic turbulent stress is omitted. Consequently
no claim about matching Fluent's physical static pressure follows from these tests.

Boundary methods follow the UG: k=1.5*(U*I)^2,
epsilon=.09^(3/4)*k^(3/2)/ell, or epsilon=rho*.09*k^2/(mu*ratio),
with ell=.07*Dh for the diameter method. The empirical inlet Cmu=.09 is distinct
from the variable constitutive Cmu. Direct epsilon and direct omega are distinct
types, so using SST data for k-epsilon throws.

## EWT boundary

`evaluate_ewt()` implements the documented Wolfstein/Chen-Patel algebra:
Re_y=rho*y*sqrt(k)/mu, Re_y*=200, A_mu=70, C_l*=kappa*Cmu^(-3/4),
A_epsilon=2*C_l*, the two damped length scales, algebraic epsilon, and viscosity
blending. Its `EwtBlend{delta_re_y,cmu}` is required explicitly. No Fluent default
is invented for delta_re_y: the TG gives a typical range of 5-20% of 200.
The helper handles zero k and infinite wall distance, and avoids exp subtraction
cancellation near the wall. `rke_ewt_law` implements the Kader velocity/derivative
blend for zero pressure-gradient and zero thermal effects, with a=.01 and b=5.

These helpers **do not implement full EWT transport**. The TG describes epsilon
blending only qualitatively, including switching to algebraic epsilon in a
fully viscosity-affected domain. No undocumented cell/row blending, wall
production integration, pressure-gradient or thermal modification is guessed.
`KeWallTreatment::Enhanced` therefore throws at model construction. EWT and
SA automatic walls require additional implementation/evidence before release.

## Sources and provenance

Equations were read from the 2025 R1 PDF and cross-checked against current 2026 R1
sections. Wall-equation numbering shifted by 19 in 2026 R1 (scalable: 4-376;
two-layer: 4-383..392; Kader: 4-393..395). SA and realizable transport numbering
is unchanged. Numerical source checksums and run identity accompany the results.

- [Fluent SA production](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_th/flu_th_sec_spal_turb_prod.html),
  [SA walls](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_th/flu_th_sec_turb_sa_wallbc.html)
- [Fluent realizable model](https://ansyshelp.ansys.com/public/views/secured/corp/v261/en/flu_th/flu_th_sec_turb_realizke.html),
  [epsilon wall treatments](https://ansyshelp.ansys.com/public/views/secured/corp/v261/en/flu_th/ch04s18s02.html)
- [Fluent boundary specifications](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_bcs_sec_bound_cond.html)
- [NASA SA definitions](https://tmbwg.github.io/turbmodels/spalart.html),
  [plate](https://tmbwg.github.io/turbmodels/flatplate_sa.html),
  [bump](https://tmbwg.github.io/turbmodels/bump_sa.html)
- [MKM channel DNS](https://turbulence.oden.utexas.edu/MKM_1999.html), nominal
  Re_tau=590, actual 587.19 in the supplied file. TMR's Chien k-epsilon is a
  different model and is not used as a realizable-model reference.

## Case parsing and app-lane handoff

The parser records `rp-sa?`, SA production/damping/enhanced-wall switches,
named SA coefficients, `nut/scheme` and raw SA boundary selectors/data. It
records k-epsilon coefficients and existing wall flags. Missing values remain
missing; no option is silently mapped to a different model or wall treatment.
Unknown active SA switches and unknown SA numeric controls remain rejected.
The **model allowlist is unchanged**, and SA's selector mapping and coefficient
key provenance still need a real Fluent case fixture before adapter acceptance.

The app lane would need model-specific data import/initialization (nu-tilde vs
k/epsilon), schemes/relaxation, output fields/residuals, BC-method mapping and
model/wall selection. Do not expose these as verified Fluent models yet. New
library options exist only for controlled verification. No app source was edited.
