# Kutta cases: a mesh and a setup file, no Fluent case

A Kutta case is two files: a mesh, and a small JSON setup (`*.kutta.json`,
version 1) that names it. Everything the setup leaves out is Fluent's default
for a new steady, pressure-based case, so a file states only what differs
from Fluent. `kutta solve`, `kutta inspect`, `kutta export-geometry` and the
desktop app take a Kutta case wherever they take a `.cas`:

```sh
kutta new-case car.msh --speed 20 --out car.kutta.json   # starter from the mesh's zone names
kutta inspect car.kutta.json                             # what it means, issues, zones
kutta solve car.kutta.json --backend metal --out run     # same outputs as a Fluent case
kutta convert old.cas.h5 --out old.kutta.json            # a Fluent case's settings as a Kutta case
```

## How it works

The setup is parsed strictly into typed settings (`src/kutta/setup/kutta_case.hpp`):
unknown properties, wrong types, out-of-range numbers and settings that do not
apply to the chosen model are all reported together with their JSON paths, and
a near-miss key names the property it probably meant. The typed settings are
then *lowered*: Kutta writes the Rampant and Thread variables a Fluent writer
would save for those settings on that mesh, with Fluent's defaults for the
rest, and reads them back with the same reader as a `.cas`. Validation, the
settings inventory and overrides (`--setup`), approximations, the app and every
solve path are shared, so a Kutta case and a Fluent case with the same settings
solve to the same bytes (see [Equivalence](#equivalence)).

The mesh can be a Fluent mesh (`.msh`, `.msh.h5`) or a `.cas`/`.cas.h5`, of
which only the mesh is used. A relative mesh path is relative to the setup
file. Zone types in the mesh file are hints only: a mesher that exports every
zone as a wall is fine, because each boundary zone's condition comes from the
setup. Every boundary zone of the mesh needs an entry under `boundaries`;
interior face zones are interior unless assigned a porous jump. Cell zones
are fluid with the case's default or named material; optional `cell_zones` entries select
constant frame motion and porous resistance by exact mesh zone name.

A Kutta case's identity (`inspect`'s `file.sha256`, `summary.json`'s
`case_sha256`, the app's history) is the SHA-256 of the setup file's and the
mesh file's SHA-256 digests, so changing either is a different case.

## Format

```json
{
  "kutta_case": 1,
  "mesh": "car.msh",
  "models": {"viscous": "k-omega-sst"},
  "material": {"name": "air", "density": 1.225, "viscosity": 1.7894e-05},
  "boundaries": {
    "inlet": {"type": "velocity-inlet", "speed": 20, "direction": [1, 0, 0],
              "turbulence": {"intensity": 0.01, "viscosity_ratio": 10}},
    "outlet": {"type": "pressure-outlet", "gauge_pressure": 0},
    "ground": {"type": "wall", "motion": "moving", "speed": 20, "direction": [1, 0, 0]},
    "symmetry": {"type": "symmetry"},
    "car": {"type": "wall"}
  },
  "reference_values": {"area": 1.1, "length": 1.6, "density": 1.225, "velocity": 20},
  "reports": {
    "drag": {"type": "drag", "zones": ["car"], "direction": [1, 0, 0]},
    "downforce": {"type": "lift", "zones": ["car"], "direction": [0, 0, -1]}
  }
}
```

`kutta_case`, `mesh`, `models.viscous` and `boundaries` are required. Numbers
are SI. Vectors are `[x, y, z]`; on a 2D mesh the z component must be 0.

| Property | Values |
| --- | --- |
| `models.viscous` | `laminar`, `k-omega-sst`, `k-epsilon-standard` |
| `models.transition` | `none`, `intermittency` (SST only) |
| `models.k_omega_sst` | `production_limiter`, `production_limiter_clip_factor`, `low_re_corrections`, `curvature_correction`, `compressibility_effects`, `kato_launder`, `vorticity_based_production`, `cross_diffusion_clip` |
| `models.k_epsilon_standard` | `compressibility_effects`, `source_modification`, `time_marching` |
| `models.intermittency` | configurable transition constants: `ctu1`, `ctu2`, `ctu3`, `cpg1`, `cpg2`, `cpg3`, `cpg1_limit`, `cpg2_limit`, `csep`, `tke_term_coefficient`, `sigma`, `minimum`, `maximum` (requires `transition:"intermittency"`) |
| `material` | default fluid: `name`, `density` (kg/m³), `viscosity` (Pa s) |
| `materials.<name>` | additional constant fluids; `density` and `viscosity` both required |
| `operating_pressure` | Pa |
| `pressure_reference_point` | `[x,y,z]`, Fluent's pressure-reference location (m) |
| `cell_zones.<zone>` | optional `material` name, `frame_motion` and `porous` objects (below) |
| `relative_frame_origins_absolute` | boolean, case-wide interpretation of child-frame origins; default `false` |
| `boundaries.<zone>` | `{"type": ...}` with the type's properties below |
| `reference_values` | `area`, `length`, `density`, `velocity`, `pressure`, `viscosity`, `depth` |
| `methods` | `pressure_velocity_coupling` (`simple`, `simplec`, `coupled`), `pseudo_transient`, `gradient` (`least-squares-cell-based`, `green-gauss-cell-based`), `pressure` (`standard`, `second-order`), `momentum` / `turbulent_kinetic_energy` / `specific_dissipation_rate` / `turbulent_dissipation_rate` / `intermittency` (`first-order-upwind`, `second-order-upwind`), `flux_type` (`auto`, `rhie-chow-momentum`, `rhie-chow-distance`), `gradient_limiter` (`standard`, `differentiable`, `none`), `warped_face_gradient_correction`, `high_order_term_relaxation` |
| `controls` | `flow_courant_number`, `pseudo_time_scale_factor`, `pseudo_time_step` (number or `"automatic"`), `pseudo_relaxation`, `relaxation`, `explicit_relaxation` (each by equation: `pressure`, `momentum`, `density`, `body_force`, `k`, `omega`, `epsilon`, `turbulent_viscosity`, `intermittency`), `limits` (`min_turbulent_kinetic_energy`, `min_specific_dissipation_rate`, `min_turbulent_dissipation_rate`, `max_turbulent_viscosity_ratio`) |
| `initialization` | `method` (`hybrid`, `standard`), `hybrid_iterations`, `velocity`, `gauge_pressure`, `k`, `omega`, `epsilon`, `intermittency` |
| `reports.<name>` | `type` (`force`, `drag`, `lift`), `zones`, `direction`, `coefficient` |
| `monitors.residuals.<residual>` | `check`, `absolute_criterion` |

Boundary types:

| `type` | Properties |
| --- | --- |
| `wall` | `motion` (`stationary`, `moving`, `rotating`); moving: `speed` and `direction`, or `velocity`; rotating: `angular_speed` (rad/s, right-handed), `axis_origin`, `axis_direction`; both: `relative` (to the adjacent cell zone) |
| `velocity-inlet` | `speed` (normal to the boundary), `speed` + `direction`, or `velocity`; `gauge_pressure`, `intermittency`, `turbulence` |
| `pressure-outlet` | `gauge_pressure`, `backflow_pressure` (`total`, `static`), `prevent_reverse_flow`, `intermittency`, `turbulence` (backflow) |
| `pressure-inlet` | `total_pressure` (gauge total pressure, Pa), `gauge_pressure` (supersonic/initial static gauge pressure), optional inward `direction`, `intermittency`, `turbulence` |
| `mass-flow-inlet` | `mass_flow_rate` (kg/s) **or** `mass_flux` (kg/(m² s)), optional inward `direction`, `gauge_pressure`, `intermittency`, `turbulence` |
| `porous-jump` | `permeability` (m², > 0), `thickness` (m, ≥ 0), `inertial_resistance` (1/m, ≥ 0); retained **internal** face zone only |
| `symmetry` | none |
| `interior` | none (interior face zones only) |

`turbulence` is `{"method": ..., ...}`: `intensity-viscosity-ratio` (the
default; `intensity`, `viscosity_ratio`), `intensity-length-scale`,
`intensity-hydraulic-diameter`, `k-omega` (`k`, `omega`) or `k-epsilon`
(`k`, `epsilon`). The method may be omitted when the values name it.

### Wheels, ground and fluid zones

```sh
kutta new-case car.msh --speed 20 \
  --wheel front-left=0.23,0.8,0.6,0.23,0,-1,0 \
  --wheel front-right=0.23,0.8,-0.6,0.23,0,-1,0 --out car.kutta.json
```

Each repeatable `--wheel` is `ZONE=RADIUS,X,Y,Z,AX,AY,AZ`: radius and axis
origin in metres, followed by the signed rotation axis. It targets an exact
wall name (including a `--role ZONE=wall` override). Angular speed is
`freestream speed / radius` in rad/s and motion is absolute (`relative:false`).
For +x flow and +z up, axis −y makes the wheel's bottom move along +x,
consistent with the moving ground. In 2D, use +z for +x motion at the bottom.
The starter requires wheel geometry from the user; names alone do not imply
centres or radii. Wheels without geometry remain stationary with a note.

```json
"cell_zones": {
  "radiator": {
    "porous": {"viscous_resistance": [100, 10000, 10000],
               "inertial_resistance": [5, 100, 100], "porosity": 0.8,
               "direction_1": [1, 0, 0], "direction_2": [0, 1, 0],
               "relative_velocity": true}
  },
  "rotor": {
    "frame_motion": {"angular_speed": 80, "axis_origin": [0.8, 0.6, 0.23],
                     "axis_direction": [0, -1, 0], "velocity": [0, 0, 0]}
  }
}
```

`material` defines the default fluid. Additional `materials` are keyed by name
and require explicit constant `density` and `viscosity`; a cell zone's
`material` selects one. Multiple names are supported only when every active
fluid zone has identical density and viscosity, as the Fluent adapter requires.
Missing names and duplicate definitions of the default material are errors.

Omitting a cell-zone entry leaves it stationary and nonporous. An empty
`frame_motion` object enables zero frame motion; `porous`
requires explicit principal viscous resistances. Viscous and inertial resistance
vectors must be nonnegative. Direction 1 and Direction 2 must be nonzero and
nonparallel; Fluent orthogonalizes Direction 2 against Direction 1. Porosity
is in [0,1]. The implemented formulation transports superficial velocity;
`relative_velocity` chooses resistance relative to the zone frame or absolute
velocity. Both objects may appear on the same zone.

`frame_motion.relative_to` names another cell zone; omitted means the
absolute frame. Missing names and cycles are rejected. Origins of relative
frames are relative to the parent by default; `relative_frame_origins_absolute:true`
uses their absolute positions instead. Translational and angular
velocities compose through the same Fluent adapter as native cases. Only
constant Cartesian frame motion in the absolute-velocity formulation is
supported. A moving zone cannot touch a pressure outlet; at a symmetry plane
its frame velocity must be tangential. These existing solver gates still
apply. On 2D meshes motion stays in x-y and the rotation axis must be along z.

Pressure and mass-flow inlet directions are inward, normalized Cartesian
vectors; omit `direction` for the local inward normal. Mass-flow rate and
mass flux are mutually exclusive and nonnegative. Frame-relative inlet data,
reverse-flow suppression, averaged or spatially varying mass flux remain
unsupported. `gauge_pressure` is preserved as Fluent's initial/static input;
this incompressible solver uses total pressure at a pressure inlet.

Porous jumps use `alpha`, `dm`, `c2` in the native Fluent thread and
`mu/alpha * thickness * velocity + rho/2 * C2 * thickness * |velocity| * velocity`
for the pressure loss. The mesh must retain that internal face-zone identity
and orientation (for example a mesher-exported `porous-jump`/fan zone).
An ordinary `interior` zone whose identity the mesh reader discards cannot be
selected as a jump by inventing a name. No internal faces are inferred from
boundary surfaces.

The exposed non-default model controls are those the strict Fluent adapter
already uses: SST's positive `cross_diffusion_clip`, the intermittency
constants listed above, and k-ε's two solver switches. Other SST/k-ε model
constants remain unsupported. Transition `sigma` must be positive, `csep`
and `tke_term_coefficient` nonnegative, and `0 < minimum < maximum <= 1`.
Turning both k-ε switches off removes the need for `--approximate ke-options`;
leaving either on still requires it. No approximation is silently enabled.

## Defaults

Sources: **UG**/**TG** the Ansys Fluent User's/Theory Guide section
(2026 R1 unless the linked edition says otherwise);
**UG fig** a value shown in a UG dialog figure, recorded in
`docs/research/fluent-numerics-and-files.md`; **O** the value Fluent writes
into new cases where the guides state none, observed identically in the
Fluent 22.2 mixing elbow, the 23.1 `elbow_param` and the 22.2 nozzle examples.

| Setting | Default | Source |
| --- | --- | --- |
| Solver | pressure-based, steady, double precision | UG 9 |
| Operating pressure | 101325 Pa | UG 9.15 |
| Material | `air`, 1.225 kg/m³, 1.7894e-05 Pa s | O |
| SST production limiter | on, clip factor 10 | TG 4.21, UG 16.16.6 |
| SST low-Re, curvature, compressibility, Kato-Launder | off | UG 16.16.11, TG 4.19, UG 16.6.4, TG 4.21 |
| Wall ω treatment | correlation | TG 4.18.3 ("current default") |
| k-ε constants, wall functions | TG 4.3.1.2; standard wall functions | TG, UG 16.16 |
| Turbulence limits | k ≥ 1e-14, ω ≥ 1e-20, ε ≥ 1e-20, μt/μ ≤ 1e5 | O |
| Pressure-velocity coupling | Coupled | UG 37.3 |
| Pseudo-transient (Global Time Step) | on with Coupled | UG 37.14 |
| Pseudo time step | automatic, conservative length scale, factor 1 | UG 37.14.2.2, 37.14.2.3 |
| Flow Courant number | 200 | UG 37.3.1.5 |
| Explicit relaxation (Coupled, pseudo time off) | 0.5 | UG 37.3.1.5 |
| Pseudo-time relaxation | p 0.5, momentum 0.5, density 1, body force 1, k/ω/ε/intermittency 0.75, μt 1 | UG fig, O |
| Under-relaxation (pseudo time off) | p 0.3, momentum 0.7, density 1, body force 1, k/ω/ε/intermittency 0.8, μt 1 | UG fig, O |
| Gradient | least squares cell based, standard limiter | UG 37.3 |
| Pressure interpolation | second order | UG 37.3 |
| Momentum | second-order upwind | UG 37.3 |
| k, ω, ε, intermittency | first-order upwind | UG 37.3 |
| Flux type | Auto Select (momentum-based Rhie-Chow) | UG 37.3.2 |
| Warped-face gradient correction | on for 3D meshes with polyhedral cells, else off | UG 37.3 |
| Initialization | hybrid, 10 iterations | UG 38, O |
| Reference values | area 1, length 1, density 1.225, velocity 1, pressure 0, viscosity 1.7894e-05, depth 1 | UG 42.13 fields, O values |
| Inlet/outlet turbulence | intensity 5 %, viscosity ratio 10 | UG 8.4 |
| Intermittency at inlets/outlets | 1 | O |
| Velocity inlet | magnitude normal to the boundary | UG 8.4 |
| Pressure outlet | gauge 0 Pa, total backflow pressure, normal backflow, reverse flow not prevented | UG 8.4, O |
| Wall | stationary, no slip; moving walls relative to the adjacent zone, direction +x | UG 8.4.15, O |
| Wall/frame rotation | zero angular speed, origin [0,0,0], axis [0,0,1]; frame translation [0,0,0], parent absolute | UG 8.2.1.1.6, 10.2; O for stored axis/zero values |
| Pressure inlet | gauge total and initial/static pressures 0 Pa; inward normal, absolute frame | UG 8.4.3.1.1.3, 8.4.3.3 |
| Mass-flow inlet | mass-flow rate 1 kg/s, inward normal, initial pressure 0 Pa, absolute frame | UG 8.4.5.1, 8.4.5.3 |
| Porous fluid zone | viscous resistance required; inertial resistance [0,0,0], directions +x/+y, porosity 1, superficial velocity, relative resistance on, relative viscosity 1 | UG 8.2.3.6–8.2.3.8; O for inertial zeros and direction axes |
| Porous jump | permeability 1e10 m², thickness 0 m, inertial resistance 0 /m (no loss) | UG 8.4.22.2, fig 8.81 |
| Relative frame origins | relative to parent (`false`) | O (`reference-frame/relative/origin-absolute?`); conversion preserves the saved flag |
| Pressure reference point | [0,0,0] m | O |
| SST negative-cross-diffusion clip | 1e6 | O |
| Intermittency constants | ctu1 100, ctu2 1000, ctu3 1; cpg1 1, cpg2 −0.5, cpg3 0; limits 1.5/3; csep 1, tke coefficient 1, sigma 1; intermittency limits 1e-7/1 | TG 4.8; O for stored limits |
| k-ε source modification / time marching | both on; require named approximation | O (undocumented Fluent solver switches) |
| Residual monitors | every solved equation checked, absolute criterion 1e-3 | UG 37.17 |
| Report coefficients | drag and lift reports are coefficients; force reports are forces | O |
| Lift direction | +y; drag and force direction +x | O |

There is no universal observed viscous-resistance default: the inactive
porous settings in Fluent 22.2 mixing elbow/nozzle store 211100000 /m²,
whereas the 23.1 `elbow_param` stores zero. The guide defines the resistance
but does not assign a numerical default, so JSON requires
`porous.viscous_resistance` explicitly rather than choosing one writer's
value. Conversion always writes the saved three resistances.

The boundary and porous defaults above are documented in the
[Fluent User's Guide, §8.4.3, §8.4.5, §8.4.15 and §8.4.22](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_ug/flu_ug_bcs_sec_bound_cond.html),
including the [porous-jump default dialog](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_ug/graphics/g_flu_ug_BCPorousJump.png),
and [§8.2.1 and §8.2.3](https://ansyshelp.ansys.com/public/views/secured/corp/v261/en/flu_ug/flu_ug_bcs_sec_cell_zones.html).
Frame equations follow [Theory Guide §2.2.1, absolute-velocity formulation](https://ansyshelp.ansys.com/public/views/secured/corp/v261/en/flu_th/flu_th_sec_rotate_equations.html).
Angular speeds in JSON are always SI rad/s, regardless of Fluent's display units.

Fluent's explicit relaxation for Coupled without pseudo time is documented as
0.5 (UG 37.3.1.5), but the writers in the corpus store 0.75; a Kutta case
uses the documented value, and `kutta convert` writes 0.75 explicitly when a
case has it.

## What a version 1 Kutta case cannot express

`kutta convert` lists every setting it cannot write and refuses when one would
change the solution. Not expressible, by design of version 1: energy (with
constant properties it does not change the flow, so convert notes and drops
it), gravity, solid zones, nonuniform material properties, far fields, native fan
polynomial input, wall roughness and specified shear, turbulence models other than laminar, SST
(with or without intermittency) and standard k-ε, non-default SST/k-ε constants other than the exposed controls,
node-based gradients, PISO, surface reports other than force/drag/lift, and
profiles. Expressions with a constant value are written as that value.
Run controls (iteration count, convergence switch) stay command-line options.

The approximations (`--approximate`) a Fluent case needs are still needed by a
Kutta case with the same physics: `wfgc` on a 3D polyhedral mesh, where
Fluent's default turns warped-face gradient correction on; `compressibility`
and `ke-options` for k-ε with compressibility effects and Fluent's default k-ε
solver options. The exception is `reverse-flow`: Fluent 18 and older left
`prevent-reverse-flow?` unwritten, while a Kutta case always states it.

## Equivalence

`tests/data/kutta_cases/` holds the Kutta cases `kutta convert` writes for four
Fluent cases of `validation/fluent_examples` (mesh paths relative to that
root). `tests/app/test_kutta_case_cli.cpp` checks, in the quick tier, that
`kutta convert` still writes exactly those files, and in the slow tier that
each Kutta case, solved on every backend present, writes the same
`fluxes.csv`, `forces.csv`, `loads.json` and `result.vtu` bytes, and the same
`residuals.csv` apart from wall time, as its Fluent case:

| Fluent case | `--approximate` (both) | Result |
| --- | --- | --- |
| `mixing_elbow/mixing_elbow.cas.h5` (V6) | `wfgc` | byte-identical |
| `mixing_elbow/elbow_param.cas.h5` | `compressibility,ke-options` | byte-identical |
| nozzle (`flprj_3891/…/nozzle_para_named_1.cas.h5`) | `wfgc` | byte-identical |
| `adjoint_cylinder_tutorial/cylinder_tutorial.cas` | `reverse-flow` (Fluent case only) | identical apart from the Fluent case's approximation notice |

The FSAE ground-effect wing (`tools/fsae_aero/naca_ground_case.py`) is
generated both ways: `--write cas` edits the mixing elbow's settings into a
Fluent case, `--write msh` and `--write kutta` write the mesh and a Kutta
case from the script's options alone. Both are generated fixtures
(`validation/generated_fixtures.tsv`), and the slow tier checks that the
Kutta case solves to the same bytes as the Fluent case.

The physics extension adds quick-tier typed/default/error-path tests and a
slow real-CLI test using eighteen small native Fluent encodings: rotating and
translating walls (components and direction), pressure inlets (normal and
direction), mass-flow rate and flux, porous resistance, MRF, relative porous
MRF, porous jumps, transition constants, k-ε solver options, named materials,
nested frames with relative and absolute origin coordinates, absolute porous
MRF, and a nonzero pressure-reference point. Every pair
compares all four output files byte for byte on each available backend.
A quick CLI test also writes and solves a starter with a rotating wheel.
These are input-equivalence tests, not new physical-accuracy claims.

The corpus physics test scans all native cases and attempts the strict Fluent
setup first. Runnable originals must convert and match outputs. Held originals
are reported as `KUTTA_CORPUS_HOLD` with their exact existing setup errors;
conversion cannot establish equivalence to an original that cannot run.
See [the extension evidence](../validation/results/mesh_to_run/PHYSICS.md)
for the observed corpus holds and CPU/Metal baseline comparisons.

A small, unrestricted circular-wheel/road mesh is committed at
`tests/data/meshes/fsae_wheel.msh` (216 quadrilateral cells, 2D). Regenerate it
with `uv run tools/fsae_aero/wheel_starter_mesh.py --out PATH.msh`. Its wheel
has radius 0.2 m and centre [0.4,0.25,0]:

```sh
kutta new-case tests/data/meshes/fsae_wheel.msh --speed 20 \
  --wheel wheel_front=0.2,0.4,0.25,0,0,0,1 --out wheel.kutta.json
kutta solve wheel.kutta.json --backend cpu --iterations 40 --out wheel-cpu
kutta solve wheel.kutta.json --backend metal --iterations 40 --out wheel-metal
```

This fixture checks setup, rotation direction, ground consistency and solve
execution; its coarse mesh is not a wheel-aerodynamics accuracy certificate.
