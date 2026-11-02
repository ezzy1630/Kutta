# Ahmed 35 degree drag diagnosis

The 0.26 drag gate and its 8% tolerance are unchanged. The original finest
full-model result, Cd 0.3137 on 1.33M cells, remains a failure. This investigation
adds a smaller refinement family, separates pressure and friction forces, and
fixes errors in the reporting code. It does not yet establish a mesh-independent
SST drag or a new solver-arithmetic defect. The measured compact-family
refinement response is Cd 0.36207 → 0.33032; it is not a solver recovery.
The original +20.7% result is not fixed, and its numerical/model split remains
unidentified. Native-mesh robustness failures and unsuccessful method controls
are recorded below.

## Source and method

The four hex-mesh runs and the initial numerical controls use the Release build of
`8866cc0e16df5b0642b2243611099ddccd7f3cd3`,
Apple Metal on the shared M4 Pro, macOS 26.6.2 and AppleClang 21.0.0. The benchmark binary SHA-256 is
`9dfa41d6c775d95b8367b5eaf2a54c309230a8cea14be8ad32572158bdce8a31`.
No source under `src/`, solver default, test gate or iteration budget changes.
The benchmark also exposes the existing pressure-interpolation setting, with its
second-order default preserved. The two pressure-scheme controls use a rebuild
of `57b593e9019b6723d0eb0d960ef1100fe7c6e744`. Their receipt records the actual
scheme and executable hash; an unconfirmed explicit scheme is rejected.

The case is the existing full Ahmed body with 35 degree slant, L 1.044 m,
W 0.389 m, H 0.288 m, 0.222 m slant length, R 0.100 m nose radii and 0.050 m
clearance. It omits stilts. The domain, stationary ground, slip far boundaries,
Re_L 2.78e6, Tu 1%, viscosity ratio 1 and reference area 0.112032 m² are unchanged.
The wall layers start at 15 µm, grow by 1.35 and span 17.283 mm in 20 layers.
SST uses its existing defaults. Convection is second-order upwind. SIMPLEC uses
0.7 momentum and 1.0 pressure relaxation; the explicitly named Coupled control
uses its existing 0.5/0.5 relaxation defaults.

Each run retains the existing 5000-iteration cap, residual criterion 1e-3 and
relative lift-and-drag change criterion 5e-4 over 1000 iterations. There is no
wall-clock stop. A run which misses any convergence condition is labelled
unconverged. For such runs, the quoted Cd is the mean of the last 2000 iterations,
as in the earlier reports; it is an iterative statistic, not a physical time
average. Force components and field comparisons use the final iterate, not that
mean. Wall times are omitted because the GPU was shared. Only one heavy solve
ran at a time in this lane.

`MallocSpaceEfficient=1` is set before launching the benchmark. Unlike `kutta`
it does not re-execute itself with that allocator setting. Peak memory is the
macOS `phys_footprint` reported by `/usr/bin/time -l`, including Metal buffers.

## Mesh study

The original mesher refines a large outer volume. Its threshold reaches the
0.4 m far-field spacing at 3 m from the body, and the wake box has a 2 m
transition thickness. The compact family changes those distances to 0.6 and
0.3 m. It keeps the domain, body surface algorithm, wall-normal layers and
wake-box extents. Each level divides the specified body surface spacing by
sqrt(2). These are local refinements, not geometrically similar uniform grids;
cell-count^(-1/3) would give a misleading refinement ratio.

The new options default to the original values. The recorded regression meshes
and `validation/generated_fixtures.tsv` are unchanged. Regenerating the default
L1 15 µm mesh reproduces the ledger SHA-256
`9f39867109f04f2f0cc2bcc1207e195d7029f68f4d094c3676d67911918cb135`
byte for byte. Mesh hashes, metadata and
portable solve commands are in `summary.json`. Rebuilding compact L1 through
the retained CLI reproduced every byte after the mesh's comment line. The
comment differed because the exploratory call used integer `35` and the CLI
used float `35.0`.

| run | cells | iterations | Cd | last-2000 range | max residual | converged | peak GB |
| --- | --- | --- | --- | --- | --- | --- | --- |
| compact_L1 | 327,584 | 5000 | 0.36207 | 0.36197–0.36215 | 9.4e-06 | no | 1.03 |
| compact_L2 | 563,140 | 5000 | 0.34275 | 0.33833–0.34680 | 0.0074 | no | 1.68 |
| compact_L2.75 | 919,620 | 5000 | 0.33032 | 0.32498–0.33650 | 0.0087 | no | 2.66 |
| reference_L1 | 446,364 | 5000 | 0.35448 | 0.35443–0.35454 | 5.1e-06 | no | 1.36 |

![Mesh trend and final-iterate drag components](study.png)

The original-grading L1 control has exactly the same body-wall face centroids
and areas as compact L1. Changing only off-body grading raises the last-2000
mean Cd by **0.00759**, from 0.35448 to 0.36207. Compact and original mesh
levels therefore cannot be pooled into a single refinement series. The current
ledger L1 contains 446,364 cells; the older archived study reported 446,128.
Those archived values are not asserted to be exact same-mesh reruns.

### Richardson extrapolation

Use `Cd(h) = Cd_infinity + C h^p`, with h the specified surface spacing.
A three-point algebraic fit alone is not evidence of the asymptotic range.
The analysis code supports unequal refinement ratios and explicitly rejects
nonmonotonic series and series with no positive-order solution.

The earlier thin-wall study's exact archived means are 0.35451567, 0.33264874,
0.31440100 and 0.30308473 at L1, L2, L3 and L3.5. Fitting L1–L3 gives
p = 0.5221 and **Cd_infinity = 0.2224**. The next triplet, L2/L3/L3.5,
has no positive-order fit. The extrapolated 0.2224 is therefore rejected.
It is not a prediction or a gate pass. Those archived runs predate this branch;
they are retained evidence from the external-3D report, not fresh runs here.

The fresh compact triplet gives a formal p = 0.5042 and Cd_infinity = 0.2416. Using the eight combinations of the observed iterative range endpoints gives 6 positive-order fits, with Cd_infinity 0.2209–0.3362. This endpoint exercise is a sensitivity check, not a confidence interval. Fixed wall-normal/far-field resolution, unconverged snapshots and the exterior-mesh control prevent a mesh-independent model-error attribution. An extrapolation inside tolerance is not treated as a physical-gate pass.

## Force accounting and location of the error

For fluid-outward face normal n, pressure drag is `sum(Cp n_x A)/A_ref` and
friction drag is `-sum(tau_x A)/(0.5 A_ref)`. `walls.csv` stores stress on the
fluid, hence the minus sign. The front includes the nose and its 100 mm fillets;
the slant and base use their CAD orientation. The remaining roof, underside
and sides are also included so that the components close to the force monitor.
No pressure offset is fitted to improve the comparison.

| run | total Cd | pressure Cd | friction Cd | front P/F | slant P/F | base P/F | other P/F |
| --- | --- | --- | --- | --- | --- | --- | --- |
| compact_L1 | 0.36204 | 0.31780 | 0.04424 | -0.00170 / 0.00866 | 0.18796 / -0.00058 | 0.13154 / 0.00000 | 0.00000 / 0.03616 |
| compact_L2 | 0.34636 | 0.30174 | 0.04462 | -0.00367 / 0.00864 | 0.17215 / -0.00050 | 0.13326 / 0.00000 | 0.00000 / 0.03648 |
| compact_L2.75 | 0.33650 | 0.29148 | 0.04502 | -0.00529 / 0.00867 | 0.16075 / -0.00040 | 0.13603 / 0.00000 | 0.00000 / 0.03676 |
| reference_L1 | 0.35450 | 0.31058 | 0.04391 | -0.00385 / 0.00862 | 0.18506 / -0.00055 | 0.12937 / 0.00000 | 0.00000 / 0.03585 |

Pressure drag of an individual open surface depends on the pressure datum.
A constant Cp shift changes the front, slant and base components while cancelling
on the complete closed body. Therefore a negative front component alone does
not prove a front-pressure defect. Liu and Moser explicitly warn about the
reference-pressure dependence of components. Experiment and outlet-gauge CFD
components must be interpreted with that limitation.

Edge bands use Cartesian distances: the first 40 mm in x after the roof edge,
and the outer 30 mm in y/z on the base.

| run | slant first 40 mm | rest of slant | base outer 30 mm | base centre |
| --- | --- | --- | --- | --- |
| compact_L1 | 0.05931 | 0.12865 | 0.07955 | 0.05198 |
| compact_L2 | 0.04882 | 0.12333 | 0.07639 | 0.05687 |
| compact_L2.75 | 0.04266 | 0.11808 | 0.07501 | 0.06101 |

A frozen-field check uses exact mesh wall-owner cells, matched to the exported
cell centroids within 1e-8 m. On compact L1, replacing extrapolated wall pressure
with owner-cell pressure changes slant Cd 0.1879629 → 0.1879338. The leading
40 mm changes 0.0593139 → 0.0592874. This direct load-evaluation effect is tiny
compared with the observed suction error. It is not an ablation of pressure
reconstruction inside the momentum equations. `bench_wall_owners` makes the
owner map without constructing a solver or modifying a field.


The compact-family final-iterate slant pressure drag falls by 0.02722, while
friction changes by only +0.00077. Of the slant change, 0.01665 occurs in the
first 40 mm in x. Base pressure drag rises by 0.00449 over these snapshots.
This is evidence of a mesh-sensitive roof/slant suction field, not evidence that
all of the total-drag excess is mesh error.

## Experiment, CFX and wake profiles

The unchanged Cd 0.26 gate comes from Ahmed et al. at 60 m/s, Re_L 4.29e6.
The pressure/profile comparison uses Lienhart–Becker at 40 m/s, Re_L 2.78e6.
That experimental/setup difference cannot be silently assigned to SST model error.

[ERCOFTAC KBwiki AC1-05, Lienhart–Becker experiment](https://kbwiki.ercoftac.org/w/index.php/Test_Data_AC1-05)
provides rear Cp, hot-wire slant profiles and LDA wake planes,
including x = 80, 200 and 500 mm behind the base. It provides no total drag for
the Lienhart–Becker experiment. Coordinates in the rear Cp and velocity files are numerically in mm,
including hot-wire files whose metadata says m. Velocity columns are selected
by their metadata names: experimental BL files have U/W, whereas CFX BL files
have U/V/W. U and W are divided by the experimental bulk velocity, 40 m/s.

Pressure comparisons use the nearest wall face within 6 mm. Velocity comparisons
use the nearest cell within 6 mm for BL profiles and 30 mm for wake planes.
These are bounded sampling diagnostics, not interpolation or additional gates.
Counts are reported because coarser meshes do not support every near-wall point.
Different coverage must not be mistaken for a mesh-convergence trend in RMS error.
Steady/iterative snapshots are compared to experimental mean velocities; this
also limits the interpretation of oscillatory runs.

Finest compact mesh, final iterate. Counts are matched/eligible.

| reference | region | taps | reference mean Cp | Kutta mean Cp | Cp RMS error |
| --- | --- | --- | --- | --- | --- |
| experiment | slant | 139/139 | -0.2876 | -0.3831 | 0.1960 |
| experiment | base | 98/99 | -0.2203 | -0.2546 | 0.0893 |
| experiment | slant_leading_40mm | 42/42 | -0.2690 | -0.5137 | 0.3371 |
| experiment | slant_remaining | 97/97 | -0.2956 | -0.3265 | 0.0763 |
| CFX-5 SST | slant | 1099/1110 | -0.2334 | -0.3689 | 0.1711 |
| CFX-5 SST | base | 1027/1041 | -0.1734 | -0.2825 | 0.1627 |
| CFX-5 SST | slant_leading_40mm | 358/360 | -0.2330 | -0.4303 | 0.2473 |
| CFX-5 SST | slant_remaining | 741/750 | -0.2336 | -0.3393 | 0.1177 |

Common supported taps across all listed fresh runs. experiment slant: 87, experiment base: 68, CFX-5 SST slant: 562, CFX-5 SST base: 521.

| run | experiment slant RMS | experiment base RMS | CFX slant RMS | CFX base RMS |
| --- | --- | --- | --- | --- |
| compact_L1 | 0.3090 | 0.1556 | 0.3390 | 0.2158 |
| compact_L2 | 0.2258 | 0.1176 | 0.2676 | 0.1742 |
| compact_L2.75 | 0.1932 | 0.0871 | 0.1733 | 0.1580 |
| reference_L1 | 0.3113 | 0.1715 | 0.3173 | 0.1878 |

Finest compact mesh, final iterate. Errors are normalized by bulk velocity.

| reference | station, mm | matched/available | U RMS error / U∞ | W RMS error / U∞ |
| --- | --- | --- | --- | --- |
| experiment | x=80 | 840/840 | 0.1199 | 0.0664 |
| experiment | x=200 | 840/840 | 0.1393 | 0.0511 |
| experiment | x=500 | 840/840 | 0.0840 | 0.0257 |
| experiment | BL_x=-163 | 11/14 | 0.0874 | 0.0581 |
| experiment | BL_x=-123 | 18/20 | 0.1615 | 0.0629 |
| experiment | BL_x=-83 | 21/24 | 0.1438 | 0.0682 |
| experiment | BL_x=-43 | 19/27 | 0.1305 | 0.0774 |
| experiment | BL_x=-3 | 17/30 | 0.1060 | 0.0741 |
| CFX-5 SST | x=80 | 3619/3619 | 0.1285 | 0.0607 |
| CFX-5 SST | x=200 | 3428/3428 | 0.1718 | 0.0798 |
| CFX-5 SST | x=500 | 3336/3336 | 0.0976 | 0.0319 |
| CFX-5 SST | BL_x=-163 | 30/36 | 0.1119 | 0.0659 |
| CFX-5 SST | BL_x=-123 | 35/39 | 0.1554 | 0.0578 |
| CFX-5 SST | BL_x=-83 | 36/41 | 0.1284 | 0.0481 |
| CFX-5 SST | BL_x=-43 | 36/44 | 0.0747 | 0.0384 |
| CFX-5 SST | BL_x=-3 | 26/34 | 0.0875 | 0.0447 |

| run | roof -203 / slant -163,-123,-83,-43,-3 mm | reversed slant area |
| --- | --- | --- |
| compact_L1 | - R R R R R | 86.1% |
| compact_L2 | - R R R R A | 86.6% |
| compact_L2.75 | - R R R R A | 87.4% |
| reference_L1 | - R R R R R | 85.8% |

The finest compact snapshot still has excessive leading-slant suction: mean Cp
-0.514 against -0.269 over the supported experimental taps. The remaining slant
is much closer, -0.327 against -0.296. Common-population errors decrease with
refinement on both slant and base. Inspection of bounded nearest-cell
centreline wake profiles shows too little upper-wake velocity deficit at 80 and
200 mm: forward flow recovers too quickly there. The table gives the full-plane
aggregate errors rather than treating a centreline impression as a wake-length
measurement. The last slant station also reattaches in the finer snapshots,
where the experimental and CFX profiles remain reversed.

## What published SST results establish

[Liu and Moser, *Numerical modeling of airflow over the Ahmed body*, 2003,
Table 1](https://www.jasonyunlong.com/pdf/Vancouver2003_AB.pdf) lists SST Cd
0.241 at 35 degrees. It lists front/base/slant components 0.026/0.107/0.108.
The SST setup is not described sufficiently to make this a matched solver
oracle. Its component sum equals its tabulated total, so the table also cannot
supply an independent friction benchmark. The paper lists LSTM base/slant
components 0.129/0.121, distinct from Ahmed's 0.095/0.090. The earlier report's
0.055 friction estimate is the remainder of Ahmed's total and pressure components,
not a measured Lienhart–Becker friction force.

[Demuren and Trehan, *Computational aerodynamics of the Ahmed model car body*,
FEDSM2006-98285, Table 3](https://citeseerx.ist.psu.edu/document?doi=df6700ff91de91eb85a8c1277510500079798402&repid=rep1&type=pdf)
reports steady Fluent 6.2 SST Cd 0.28 at 35 degrees, against its quoted
experimental 0.257. They use second-order upwind convection and 2.88M tetrahedra
after a grid study with k-epsilon. Their SST solution misses the observed slant
separation and wake structure. Its apparently close drag is not full validation.

The [ERCOFTAC CFX-5 SST fields, Kuntz, Menter and Durand](https://kbwiki.ercoftac.org/w/index.php/CFD_Simulations_AC1-05)
are the closest available Ansys-family field
comparison: 2.5M cells on a half model, y+ around 1, second-order convection,
Tu 1% and viscosity ratio 1. They are not a same-mesh Fluent solution.
In the archived 25-degree cross-check, CFX and Kutta both fail to reproduce the experiment's slant
reattachment. That cross-check supports a model limitation for topology; it does
not measure the SST contribution to Kutta's 35 degree total-drag error.

These results establish that published SST calculations can give much less drag
than the current coarse Kutta meshes. They do not establish a unique expected
SST Cd, a transferable model-error percentage, or permission to change the gate.

## Demonstrated fixes

The reporting code had two errors, reproduced with analytic synthetic surfaces
before checking the corrected behavior:

- `ahmed_breakdown` always used the 35 degree slant mask. A 25 degree slant face
  at x = -0.19 m contributed zero. It now uses the run's slant angle and the
  force monitor's requested reference area. The analytic example gives
  0 → 0.1056545654 with A_ref = 2 m².
- The Cp comparison used a height range to identify slant taps. It included
  49 side-wall taps at 35 degrees and 38 at 25 degrees. It now selects the
  slant plane, allowing 0.5 mm for CFX's exported coordinate offset and rounded
  experimental coordinates. The full 35 degree slant reference population is
  139 taps, mean Cp -0.287554, instead of 188 taps, mean -0.265926. At 25 degrees
  it is 143 taps, mean -0.597972, instead of 181, mean -0.554547. A run's 6 mm
  sampling bound can reduce those populations further.

Both focused checks fail against `8866cc0` reporting code and pass after the
fix. They change the interpretation of earlier Cp tables, not the solver's
fields or drag. No solver default was adjusted to hit a literature value.

## What was tried and what remains

The prior 50 → 15 µm wall-spacing experiment changes Cd by less than 0.4% at
matching mesh levels. The existing poor-cell pressure fix already makes those
wall layers run. Neither can explain the original 20.7% finest-mesh excess.
The earlier report locates that robustness defect in sheared first-layer cells
at sharp edges. It is already fixed on this branch's base, not a new recovery
claimed here.

The compact mesh saves memory, but its measured +0.00759 exterior-grading
effect at L1 prevents treating it as an interchangeable version of the original
family. The memory-scale report already measures 4.50 GB for the original
1.33M-cell path. That is why this study stops its hex family at 919,620 cells
and a measured 2.66 GB, rather than reproducing the original finest mesh.

A half-model generator was attempted to keep the original exterior resolution
at half the cell count. Gmsh either invalidated the discrete layer shell when
meshing its shared rim, or rejected the merged surface with
`PLC Error: A segment and a facet intersect at point`. No valid half mesh or
half-model solve resulted. That attempted code was removed.

The native tet/prism mesh retains 55,371 tetrahedra and 74,960 prisms before
hex subdivision. Its 130,331 cells cost less than 0.5 GB. The following controls
all retain the 5000-iteration cap and original convergence criteria. Failed
runs have no physical Cd result; their final forces are not drag predictions.

| control | iterations | last-2000 mean Cd | max residual | result | peak GB |
| --- | --- | --- | --- | --- | --- |
| compact L2, SIMPLEC, differentiable, second-order pressure | 1500 | invalid | 22.1 | nonfinite solver input after 1500-iteration monitor | 1.684 |
| native, SIMPLEC, Standard limiter, second-order pressure | 250 | invalid | 0.688 | nonfinite solver input after 250-iteration monitor | 0.418 |
| native, identical default after adapter rebuild | 250 | invalid | 0.688 | nonfinite solver input after 250-iteration monitor | 0.417 |
| native, Coupled, differentiable, second-order pressure | 5000 | 0.42867 | 0.0227 | unconverged | 0.480 |
| native, SIMPLEC, Standard limiter, Standard pressure | 5000 | 0.49832 | 0.000317 | unconverged | 0.440 |

Changing only native pressure interpolation from second-order to Standard
prevents the nonfinite failure for all 5000 iterations, but Cd is 0.49832 and
the force criterion is still missed. This isolates pressure interpolation as
a stability-sensitive choice on this mesh. It does not isolate a faulty
pressure term, establish which field is accurate, or justify a default change.

The two native default runs, before and after adding the explicit scheme flag,
have identical values in every monitored column except seconds, through all
five monitors before failure. Thus the new adapter has not changed the default
failure. Late failed fields contain extreme velocities on flat roof and
underside wall layers as well as the base; they are too late to identify the
first failing face. A numeric location of the maximum after blow-up is not a
root cause.

The differentiable limiter is an existing opt-in, with its default coefficient
unchanged. Fluent's [Theory Guide 23.3.4.3, Gradient
Limiters](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_th/flu_th_sec_grad_limit.html)
describes convergence problems with non-differentiable limiters. That supports
the experiment, not a claim that Kutta uses identical Fluent coefficients. On
this target the existing opt-in does not provide an accepted solution. The
Coupled-plus-differentiable run changes two choices and therefore isolates only
the combination, not either choice separately.

Source inspection also rejected a missing face-value bound as a proposed fix:
`fv_moment_momentum` already bounds its extrapolated value to the sample cloud
and applies the minimum cell limiter; the limiter itself is clamped to [0,1].
Adding another clamp would not address a demonstrated missing invariant. No
speculative kernel fix or tuned default is retained.

The measured attribution is therefore limited:

- **Mesh:** compact L1 → L2.75 lowers the last-2000 mean Cd by 0.03175,
  with unchanged solver settings and wall-normal layers. This is a measured
  response to refinement, not the remaining mesh error on the original 1.33M
  mesh. The separate exterior-grading effect is +0.00759 at L1. Neither can be
  converted into a percentage of the original finest-mesh excess.
- **Kutta numerics/robustness:** the native default failure and the hex
  differentiable-limiter failure are reproducible limitations. Direct wall
  pressure load reconstruction contributes only about 0.00003 to coarse slant
  Cd. No controlled experiment yet identifies a faulty arithmetic term that
  explains the hex pressure excess. The reporting defects above are fixed.
- **SST:** the 25-degree topology cross-check is consistent with a model
  limitation. Published 35-degree SST Cd spans at least 0.241–0.28 in the cited
  unmatched setups. A numeric SST share of the current excess remains
  unidentified; assigning the residual error to SST would be unsupported.

A matched, converged Fluent solve on one of these exact meshes is unavailable.
That is the remaining way to distinguish a Kutta discretization discrepancy
from different SST formulations/settings without relying on unrelated literature
meshes. Further refinement also needs a valid half mesh or better local cell
allocation before a mesh-independent Cd can be claimed under 3 GB.

## Verification and reproduction

Release configuration and build succeeded on Metal/CPU. The two analytic
reporting regressions and four analysis checks pass. Portability lint and
`git diff --check` pass. The default mesh regenerates to its ledger hash, and
all fresh runs fit the measured 3 GB budget. Exported integrated forces close
to the monitor within the recorded CSV precision; exact closure errors and
wall-owner identity errors are retained in `summary.json`. Own monitor histories
and portable receipts are in `evidence/`. Failed receipts mark physical fields
null because post-failure exports are not valid solutions. `-` in the topology
table means no supported slant faces at that station; -203 mm is on the roof
for the 35-degree geometry, not an attached-flow observation.

The complete quick tier passes: **277 tests, 7,457,592 assertions, zero
failures**, with `--test-suite-exclude=slow` and the prescribed fixture
environment. No skips or expected-failure markers were added. This validates
the unchanged default path and integration; it does not recover the Ahmed
physical gate.

The new benchmark's invalid pressure-scheme value is rejected before solver
creation. Explicit scheme requests are checked against the benchmark's emitted
configuration before analysis. Neither CUDA compilation nor CUDA runtime was
verified in this lane.

The new code does not change solver arithmetic. A full slow-tier comparison is
therefore not claimed or substituted by the earlier lane's comparison. The main
baseline XML is read-only and was not copied into this report. No slow XML,
restricted curves, generated meshes or cell fields are committed.

From `<checkout>`, with gmsh 4.15.2, numpy, scipy and matplotlib in `build-py`,
and a read-only link to the restricted ERCOFTAC tree at
`validation/data/ercoftac_ahmed_body`:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 2
build-py/bin/python tools/fsae_aero/ahmed.py --slant 35 --level 2 \
  --first 1.5e-5 --refinement-reach .6 --wake-padding .3 \
  --out runs/ahmed_drag/mesh/compact_L2.msh
build-py/bin/python tools/ahmed_drag/run.py \
  --mesh runs/ahmed_drag/mesh/compact_L2.msh --out runs/ahmed_drag/compact_L2
build/bench_wall_owners runs/ahmed_drag/mesh/compact_L2.msh body \
  > runs/ahmed_drag/compact_L2/owners.csv
# Other refinement levels: 1 and 2.75. Run one solve at a time.
# Existing method controls: --limiter differentiable, --coupled,
# --pressure-scheme standard. Native mesh: ahmed.py --cell-type tet-prism.
# Native default is known to fail; this is a diagnostic, not an accepted case.
build-py/bin/python tools/ahmed_drag/analyze.py runs/ahmed_drag/compact_L1 \
  runs/ahmed_drag/compact_L2 runs/ahmed_drag/compact_L2.75 --out runs/ahmed_drag/summary.json
build-py/bin/python -m unittest discover -s tools/fsae_aero -p 'test_*.py'
build-py/bin/python -m unittest discover -s tools/ahmed_drag -p 'test_*.py'
python3 tools/lint_portability.py
KUTTA_FLUENT_EXAMPLES=<reference-checkout>/validation/fluent_examples \
KUTTA_RUNS=<reference-checkout>/runs build/kutta-tests --test-suite-exclude=slow
```
