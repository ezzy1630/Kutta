# Kutta case physics extension evidence

Observed on macOS with an Apple M4 Pro, CPU and Metal, from `t3/expand-kutta-case-physics` based on
main `1e5d5ae835c06eb16d6bf6ee1f8d2cd05b95ace5`. Implementation commit
`d8cc3d3917d02d2a2da0392437b41d3d038d94c7`; the acceptance commit containing this report
adds the tested CLI gates and analytic wheel fixture.
No solver equations, `src/kutta/fv`, `src/kutta/kernels`, or Qt `app/` changed.
No restricted fixture data is included.

## Expressible settings

Version 1 now includes pressure inlets, uniform mass-flow rate/flux inlets,
constant Cartesian MRF (including named parent frames and both origin
conventions), anisotropic superficial porous resistance (absolute or relative
velocity), retained internal porous jumps, multiple material names with the
same constant fluid properties, pressure-reference position, SST negative
cross-diffusion clipping, intermittency model constants and k-epsilon solver
switches. Rotating and translating walls already existed; their real native
conversion paths now have byte-equivalence coverage. `new-case --wheel`
uses supplied radius/origin/signed axis and freestream speed to create
absolute rotating walls with angular speed `speed / radius`.

All new inputs use typed settings and strict JSON readers. Unknown keys,
wrong types, invalid ranges, incompatible options, absent zone/material
names, parallel porous axes and cyclic frame references are rejected with
paths. Conversion rejects unresolved active expressions, missing active
native fields and unsupported formulations rather than substituting defaults.
`docs/kutta-case.md` cites the Fluent guides and distinguishes documented
from observed defaults. Viscous resistance is required because the saved
Fluent defaults differ between corpus writers.

Native topology checking now recognizes BC codes 4 (pressure inlet), 20
(mass flow) and 14 (porous jump). A quick test checks each matching code and
rejects a mismatched code. This is file-adapter plumbing, not solver physics.

## Observed equivalence

The slow native test writes actual `.cas` files on the small committed box
mesh: ordinary defaults come from the lowerer, then independent Scheme
edits supply the feature under test. It invokes the real `convert` CLI and
solves the native and converted inputs on CPU and Metal for three iterations.
All eighteen scenarios match `fluxes.csv`, `forces.csv`, `loads.json` and
`result.vtu` byte for byte; residual rows also match after removing wall time:

- rotating wall; translating wall by components; translating wall by direction;
- pressure inlet by normal and Cartesian direction;
- mass-flow rate and mass flux;
- anisotropic porous resistance; constant MRF; relative porous MRF;
- porous jump; non-default transition constants and SST clipping;
- k-epsilon switches off; two named fluids with identical properties;
- nested frames with relative origins and with absolute origins;
- absolute porous MRF; a nonzero pressure-reference point.

These are input-format equivalence tests. The native fixtures are synthetic
writer encodings, not new saved Fluent reference solutions. Existing physical
accuracy gates remain unchanged.

The existing mesh-to-run equivalence tests also pass on both backends for
mixing elbow, elbow_param, nozzle, the Fluent 18 cylinder (its existing
provenance-notice exception), and the generated NACA ground-effect wing.
The NACA `.msh` and `.kutta.json` were generated into a local runs directory
and checked against the existing fixture ledger; its `.cas` was symlinked
from the external runs fixture. No external fixture was written.

An unchanged main executable was built and retained before implementation.
Each of the following was solved with it and with the implementation for
40 iterations, cold-started without a data file:

| Source | Options | CPU | Metal |
| --- | --- | --- | --- |
| `mixing_elbow/mixing_elbow.cas.h5` | `--approximate wfgc` | all four files byte-identical to main | all four files byte-identical to main |
| `mixing_elbow/elbow_param.cas.h5` | `--approximate compressibility,ke-options` | all four files byte-identical to main | all four files byte-identical to main |

This comparison has sixteen file equality checks, each also recorded with a
SHA-256 in the local `build/physics-evidence/baseline-comparison.tsv`.

The analytic circular-wheel/road mesh has 216 quadrilateral cells, radius
0.2 m and centre [0.4,0.25,0]. Its generator reproduces the committed mesh
byte for byte (SHA-256
`421da9f13f3bdcb74f7c310df79cb5ec2755aef1268fcf7bdf2bc383737c0dac`).
The real starter sets speed 20 m/s, axis +z, angular speed 100 rad/s,
absolute wall rotation, and +x ground motion. It solves for 40 iterations on
CPU and Metal with no approximations. Both stop at the iteration limit;
convergence and wheel-aerodynamics accuracy are not claimed.

## Corpus holds, not equivalence claims

The corpus inventory found 74 `.cas`/`.cas.h5` files. Twenty-one use the
relevant boundary/zone features. None can run unchanged through the existing
strict mesh/Fluent adapter, so **native corpus equivalence for the new physics
is blocked**, not passed. The slow corpus test records each original error
as `KUTTA_CORPUS_HOLD`; it compares outputs whenever an original is runnable.
The sweep also considers named fluids, configurable model options and
nonzero pressure-reference coordinates. There is no runnable original for a native-to-Kutta equivalence comparison.

| Cases (paths relative to the corpus) | Count | Original holds (representative, not exhaustive) |
| --- | --- | --- |
| `axial_fan_performance_curve/` root and three project copies | 4 | gravity, PRESTO, motion-spec 1; some unresolved expressions |
| `borealis_fracture/` three saved cases | 3 | gravity, node-based gradients; one alternative porous inertial formulation |
| `catalytic_converter/catalytic_converter_param.cas` | 1 | laminar zone in turbulent flow and unresolved inlet profile |
| `examples/MixingTank-DOE/test-laminar-visc.cas.h5` | 1 | gravity, PRESTO, motion-spec 1 and specified shear |
| `examples/Tyler-Sofrin-Modes-Compressor/` two cases; `other_folders/result_files/fluent-axial_comp/` three states | 5 | adapted/nonconformal/overset mesh interfaces |
| `optislang/nozzle-2D-WB.cas.h5` | 1 | ideal-gas density and unresolved total-pressure expression |
| `other_folders/pytwin/other_files/HX_CFD.cas.h5` | 1 | unresolved mass-flow expression |
| `peridoic_rot/`, `peridoicrot_withpipes/`, `periodic_rot/` | 3 | unsupported periodic boundaries and other original settings |
| `vki_turbine/vki_turbine.cas` | 1 | adapted/hanging-node/nonconformal mesh |
| `wigley_hull/wigley.cas.h5` | 1 | unverified vector/switch encodings and other strict adapter holds |

Gravity, solids, far fields, roughness, specified shear and native polynomial
fan import are still rejected by the solver adapter. Exposing them would
require solver/importer work beyond this task's authorization to expose
existing supported physics. CUDA compilation/execution and Windows MSVC
were not available locally; portability is inferred from standard C++20 and
the unchanged backend interfaces, not claimed as observed runtime proof.

## Reproduction and checks

From `<checkout>`, with `KUTTA_FLUENT_EXAMPLES` and `KUTTA_RUNS` pointing to
the external fixtures:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/kutta-tests --test-suite-exclude=slow
./build/kutta-tests --test-case='Native FSAE Fluent physics*,Corpus FSAE physics*'
./build/kutta-tests --test-case='Kutta cases converted from Fluent*,The NACA ground-effect wing from a bare mesh*'
python3 tools/lint_portability.py
```

Generate missing NACA files into a local runs directory using
`uv run tools/generate_fixtures.py --runs DIR --only wing/naca4412_ground.msh --only wing/naca4412_ground.kutta.json`;
symlink the existing donor `.cas` beside them and set `KUTTA_RUNS=DIR` for the
NACA check. Never generate into the external fixture checkout.

Final quick tier: **276 tests, 7,457,586 assertions, zero failures**.
Final affected slow selection: **4 tests, 564 assertions, zero failures**
(the corpus holds are reported, not counted as output equivalence).
Local evidence is under `build/physics-evidence/`: `quick-final.log`,
`slow-final.log`, `native-final.log`, `corpus-final.log`, wheel setup/run artifacts, and baseline
runs with their comparison TSV. Build output and scratch data remain ignored.
