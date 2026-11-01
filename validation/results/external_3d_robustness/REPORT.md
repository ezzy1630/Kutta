# External 3D robustness: y+ ~ 1 wall layers at sharp edges, and Coupled on the Ahmed body

**Both failures are fixed, and the Ahmed body now runs at y+ ~ 1 (15 um
first layer, body y+ mean 0.8-0.9) under SIMPLEC and Coupled, on four mesh
levels up to 2.1M cells.**

- Thin layers: the second-order face pressure of first-layer hexes sheared
  across sharp edges (Fluent orthogonal quality down to 0.014) amplified
  gradient noise about 20 times. Following Fluent's poor-mesh numerics
  (cell-quality threshold 0.05), those cells now take their face pressures
  from the cell value.
- Coupled: sink cells had `a_P` near zero, so their Rhie-Chow coefficient was
  about 850 times normal in the residual, while the Jacobian floored it. The
  residual now uses the same floor under Coupled. The same fix makes the h/c
  0.224 wing converge from `f770a31`, where it diverged.
- Ahmed 35 deg at y+ ~ 1: Cd 0.3545, 0.3326, 0.3144 and **0.3031 (+16.6%
  against 0.26)** from L1 to L3.5, within 0.4% of the 50 um values at the same
  level, so wall spacing was not the cause of the excess. The Cd is still
  falling and misses the 8% gate. SIMPLEC and Coupled agree to 0.1% where
  SIMPLEC settles. At 25 deg the slant stays separated, as in CFX-5 SST and
  unlike the experiment (Coupled matches CFX at every station, SIMPLEC is
  attached at the last).
- SIMPLE and SIMPLEC are byte-identical on meshes with no flagged cell.
  Coupled cases and the TMR bump move by <= 9e-4 on asserted errors. The slow
  tier has 0 regressions and the same 65 known-failure assertions.
- Not fixed: Coupled's continuity plateau (forces settle, residuals do not),
  and the tet duct and channel under Coupled.

## Reproduction (main `f770a31`)

Meshes are `tools/fsae_aero/ahmed.py` level 1 (35 deg slant, all-hex: tets
and prisms split by `subdivide.py`), run with `bench_fsae_aero` under the
FSAE lane's Ahmed settings (SST, Re_L 2.78e6, Tu 1%, uniform freestream
start, Fluent defaults: SIMPLEC 0.7/1.0, Coupled 0.5/0.5 and CFL 200).

300 iterations on Metal (CPU: 60), main `f770a31` against the candidate.
"Non-finite" is the linear solver rejecting a non-finite input; the
iteration given is the last 50-iteration monitor before it.

| first layer | scheme | main | candidate: Cd at 300 | candidate: max residual, max velocity |
| --- | --- | --- | ---: | --- |
| 15 um | SIMPLEC | non-finite after 50 (Cd 1e7, max velocity 4.6e5 U at 50) | 0.3506 | 1.2e-3, 1.62 U |
| 15 um | SIMPLEC, CPU, 60 iterations | non-finite after 50 (max velocity 1.7e6 U) | 0.4202 (still settling) | 2.7e-2, 1.61 U |
| 25 um | SIMPLEC | non-finite after 200 (max velocity 5.8e5 U at 200) | 0.3498 | 1.3e-3, 1.62 U |
| 50 um | SIMPLEC | runs, Cd 0.3509 | 0.3509 | 1.2e-3, 1.63 U |
| 15 um | Coupled | non-finite before 50 | 0.3551 | 5.1e-3, 1.62 U |
| 25 um | Coupled | non-finite before 50 | 0.3534 | 4.7e-3, 1.63 U |
| 50 um | Coupled | non-finite before 50 (iteration 9-10 with a monitor every iteration) | 0.3539 | 6.1e-3, 1.62 U |

Cells flagged by the fix (below): 2,127 (15 um), 371 (25 um) and 10 (50 um)
of 446,128.

## Diagnosis

Two independent mechanisms. Neither is the turbulence model or its wall
treatment, and neither is the tet boundary-row loop fixed by the tet
robustness lane.

### 1. Thin layers: second-order face pressures of sheared wall cells

**Where it starts.** Cell dumps after 1, 2, 5, 10 and 20 iterations
(`evidence.tar.gz: reproduction/`):

- Iterations 1-2: the impulsive start drives the first-layer cells on the
  flat nose to 38 U (15 um) and 139 U (50 um). Both meshes do this, and the
  50 um mesh recovers, so it is not the failure.
- From iteration 5 the growth sits in the first-layer cells at the base's
  lower outer corners (x = 0, y = +-0.182, z = 0.055: base, side and
  underside meet) and along the lower longitudinal edges; by iteration 20 the
  pressure there is 1e5 and the velocity 6e3 U.

**Which cells.** `orthogonal_quality()` (below) and a scratch probe give, for
the first-layer cell at that corner (`evidence.tar.gz: mesh_quality/`; the last two columns are mesh-wide minimum and count):

| first layer | worst wall-face angle | face skewness (link lengths) | aspect ratio | Fluent orthogonal quality | cells with quality <= 0.05 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 15 um | 89.2 deg | 20.3 | 1,629 | 0.014 | 2,127 |
| 25 um | 88.6 deg | 13.0 | 1,013 | 0.023 | 371 |
| 50 um | 87.1 deg | 6.1 | 489 | 0.045 | 10 |

The layer is marched along vertex normals, which at a sharp edge point along
the bisector, and the prisms are split into hexes. The corner cell's faces are
slightly warped (its 15 um-thick top face varies by 3 um in height), and its
centroid lands 0.33 mm sideways from both of its large faces' centroids while
the cell is 15 um thick. Every quantity that scales with that offset over the
thickness grows as 1/first-height, which is why 50 um ran and 25 um did not.
All 2,127 flagged cells are in the body's wall layers: 1,626 within 20 mm of
a sharp edge (637 upper longitudinal and slant-side, 629 lower longitudinal,
360 base perimeter) and 501 in the first layers of flat faces.

**Which term.** Ablations on the 15 um mesh, 60 iterations, Metal, each one
switch on main (diagnostic build, not committed; `evidence.tar.gz:
diagnosis/`, the switches are in `diag.patch`):

| change | max velocity at 60 (U) |
| --- | ---: |
| none (main; this run stopped at 40) | 7e4 at 40 |
| laminar (no SST, no wall treatment) | 4e8 at 40 |
| hybrid initialization | 5e5 |
| first-order upwind for every equation | 1.4e10 |
| distance (not momentum) face-velocity weights | 1.5e7 at 55 |
| **Fluent standard pressure interpolation** | **1.6** |
| second-order pressure, cell value at the faces of cells with quality <= 0.05 (wall and interior faces) | **1.6** |
| the same at their boundary faces only | 3e4 |
| the same at their interior faces only | 1.6e9 at 35 |

So the amplifier is the second-order face pressure (TG 23-43,
`p_f = (p_0 + grad p_0 . r_0 + p_1 + grad p_1 . r_1)/2`) in those cells. Its
Green-Gauss pressure force divides the difference between the two large
faces by the cell volume. With `r_0` 0.33 mm sideways on both faces, half the
difference between neighbouring cells' tangential gradients times 0.33 mm
lands on a 15 um cell: a gain of about 20 on gradient noise that its own
pressure correction feeds back. The Rhie-Chow jump is not involved: on these
faces its skew part is already faded to the link derivative.

### 2. Coupled: the Rhie-Chow coefficient of sink cells

Coupled on the 50 um mesh fails at iteration 9-10 on main. Its first
blow-up is not in the wall layer: a sign-alternating pressure spike at
x = -0.12 m, y = +-0.10 m, z = 0.315 m, about 2 cm above the slant, in
centimetre-sized tet-derived hexes of quality 0.9.

**Jacobian check.** A diagnostic took the block update `x` at a given
iteration, applied a 5% step, recomputed the steady residual and compared
each cell's true change with the Jacobian's prediction `-0.05 J x`
(`evidence.tar.gz: diagnosis/jacobian_check`):

| iteration | step | momentum: predicted / true change (u) | continuity: predicted / true change |
| --- | --- | --- | --- |
| 2 | u and p | 8.9e-5 / 8.0e-5 | 1.2e-4 / 1.6e-4 |
| 6 | u and p | 8.2e-5 / 3.8e-5 | 7.1e-5 / **5.0e-3** |
| 6 | u only | | 7.1e-5 / 4.9e-3 |
| 6 | p only | | 7.1e-5 / 5.0e-3 |
| 6 | **none** | | 7.1e-5 / **4.9e-3** |

Momentum is predicted adequately. Continuity changes by 70 times the
prediction **even with a zero step**: the mass flux used for the block RHS,
recomputed at the same state, moves by that much. The worst cell, at the
spike, is predicted +2.9e-6 and moves -2.7e-3. There the physical momentum
diagonal is 5.3e-8 against an off-diagonal sum of 4.5e-5: a sink cell, which
the current flux only feeds. Kutta's assembled `a_P` carries the net
outflow (README), so in such cells it is near zero, and the Rhie-Chow
coefficient `V/a_P` was about 850 times its value elsewhere. Every small
change of `a_P` between iterations (the turbulence update changes `mu_t`)
then swings the face flux. The coupled-external lane floored `a_P` at the
off-diagonal sum **inside the Jacobian only**; the residual's Rhie-Chow and
face-velocity weights kept the raw value. SIMPLEC is shielded because its
pressure correction re-satisfies continuity every iteration, so its sink cells
do not persist. A smaller pseudo-time step made Coupled fail sooner, as
expected when momentum is held back and the flux is set by Rhie-Chow alone.

Flooring that `a_P` in the residual too makes the 50 um mesh run (300
iterations, max velocity 1.6 U, Cd 0.354). On the 15 um mesh Coupled still
diverged until both mechanism 1's treatment and a matching Jacobian were in:
the own-cell extrapolation term `Q` of the Coupled pressure-force block
assumed `p + grad p . r` in cells that no longer extrapolate.

Ruled out on the way, with the same 50 um mesh: first-order convection,
standard pressure, distance weights and all three together (diverge by 10-20);
laminar flow; a compact Rhie-Chow jump (worse); recomputing the flux with the
current `a_P` before the block solve (fails at once); a smooth
`a_P - net outflow` without the floor's max (fails at once, float32
cancellation where outflow dwarfs the off-diagonals).

## How Fluent handles it

- **Poor mesh numerics.** Fluent applies "a local numerical correction to the
  transport equation discretization in poor-quality cells" and enables it by
  default for cells of orthogonal quality 0. `solve/set/poor-mesh-numerics/
  cell-quality-based?` extends it to quality <= 0.05, adjustable with
  `set-quality-threshold`. Its default "1st order" correction "uses locally
  low-order discretization and neglects some non-orthogonal gradient
  contributions in the diffusive fluxes"; "2nd order" assembles the bad
  cells' gradients from better neighbours; `orthogonality-enhancing-cell-centroids`
  relocates poor cells' centroids
  ([Fluent User's Guide 36.22.5, Robustness with Meshes of Poor Quality](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/flu_ug/flu_ug_sec_solve_conv_stab.html)).
- **Orthogonal quality.** The minimum over a cell's faces of the cosines
  between the face area vector and the vectors from the cell centroid to the
  face centroid and to the neighbour centroid
  ([User's Guide 6.2.2, Mesh Quality](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_GridQuality.html)).
  `orthogonal_quality()` implements that; boundary faces have only the
  face-centroid cosine.
- **Second-order pressure.** TG 23-43 extrapolates both sides with the cell
  gradients, as Kutta does; the Theory Guide warns that the scheme "may have
  some trouble ... with a bad mesh" (tet robustness report).
- **SST wall treatment.** Fluent's omega-based wall treatment is
  y+-insensitive and Kutta's is the same correlation; the laminar ablation
  diverged the same way, so the wall treatment is not the cause here.
- **Coupled.** Fluent relaxes its Coupled equations with a flow Courant number
  and pseudo-transient explicit factors (coupled-external report). It does not
  document whether its momentum `a_P` carries the net mass outflow (the Theory
  Guide gives only `a_P phi = sum a_nb phi_nb + b`, TG 23.2.1). Patankar's
  standard form uses continuity to drop that term, which gives
  `a_P >= sum a_nb`; Kutta now enforces that bound where Coupled's iterate
  needs it.

**Deviations, and why.** Fluent's default marks only quality-0 cells; that
would mark none of these (minimum 0.014) and changes nothing. Kutta uses the
documented 0.05 option by default, and corrects only the pressure
reconstruction: the convection and diffusion parts of Fluent's "1st order"
were not needed for stability (first-order convection alone still diverged),
and leaving them keeps those terms second order. Fluent's documentation does
not say which terms its correction changes, so this is the closest
documented equivalent, not a reproduction.

## The fix

Commit `9a7885b` (tests and tooling in `d405ab4`):

1. **`orthogonal_quality()`** (`mesh/mesh.cpp`): Fluent's cell metric, above.
2. **Poor-cell pressure faces** (`SolverSettings::poor_cell_quality`, default
   0.05; negative disables). For second-order pressure the solver flags cells
   at or below the threshold when it is built. A new kernel,
   `fv_poor_cell_pressure_faces`, runs after `fv_pressure_faces` and replaces
   a flagged cell's side of each face pressure with its cell pressure,
   boundary faces included (prescribed-pressure faces keep their value). The
   Rhie-Chow jump is left alone. Meshes with no flagged cell allocate and
   launch nothing, so they are bitwise unchanged; the Ahmed meshes pay 4 bytes
   per cell. `SteadySolver::poor_cells()` reports the count and
   `bench_fsae_aero` prints it.
3. **Matching Coupled Jacobian** (`CoupledSolver`, new `low_order_pressure`
   argument): flagged cells contribute no own-cell extrapolation term `Q` to
   the pressure-force block.
4. **Coupled `a_P` floor in the residual** (`coupled_floor_physical`, called
   from `update_momentum_diagonal` only under Coupled): the physical diagonal
   that Rhie-Chow and the momentum-weighted face velocities use is floored at
   its off-diagonal sum, the same floor and arithmetic as
   `coupled_copy_diagonal`. The three momentum rows share one CSR value
   buffer and identical off-diagonals, so the last assembled row serves all
   three. At a fixed point the floor adds at most the cell's continuity
   residual (`a_P - sum|a_nb|` is the net outflow plus boundary and sink
   terms). SIMPLE and SIMPLEC do not run it.
5. Tooling: `run.py --first` and fractional `--level` for Ahmed meshes
   (the mesh name gains `_first15um`; the default 50 um meshes keep their
   names), and the analysis script `analyze.py` beside this report.


## Ahmed 35 deg mesh study at y+ ~ 1

Same geometry, domain and conditions as the FSAE aero report. The difference
is the wall layer: 20 prism layers from a 15 um first height (growth 1.35,
17.3 mm, about the total thickness of the 50 um meshes' 16 layers), all split
to hexes. Level 3.5 is `ahmed.py --level 3.5`, a surface size between L3 and
L4. Each run was capped at 5,000 iterations and 5,400 s on Metal. The L2
pairs ran two at a time on the shared GPU and the rest ran alone, so the
times are indicative only. Cd for a run that did not settle is the mean of
the last 2,000 iterations. "Flagged cells" are the cells with orthogonal
quality <= 0.05.

| level | body triangles | cells | flagged cells | scheme | Cd | Cd range, last 2000 | vs 0.26 | friction / front / slant / base Cd | CL | stop | solve s | peak GB |
| --- | ---: | ---: | ---: | --- | ---: | --- | ---: | --- | ---: | --- | ---: | ---: |
| L1 | 3,748 | 446,128 | 2,127 | SIMPLEC | 0.3545 | 0.3545-0.3546 | +36.4% | 0.044 / -0.004 / 0.185 / 0.129 | 0.220 | residuals 5.1e-6 | 632 | 1.80 |
| | | | | Coupled | 0.3548 | 0.3545-0.3551 | +36.4% | 0.044 / -0.004 / 0.186 / 0.129 | 0.223 | continuity plateau 4.7e-3 | 553 | 2.13 |
| L2 | 5,696 | 757,796 | 1,262 | SIMPLEC | 0.3326 | 0.3272-0.3383 | +27.9% | 0.044 / -0.005 / 0.163 / 0.137 | 0.146 | oscillating, residual 7.5e-3 | 1614 | 3.13 |
| | | | | Coupled | 0.3288 | 0.3285-0.3290 | +26.5% | 0.044 / -0.005 / 0.163 / 0.127 | 0.141 | continuity plateau 5.1e-3 | 1263 | 3.67 |
| L3 | 9,638 | 1,444,908 | 139 | SIMPLEC | 0.3144 | 0.3109-0.3175 | +20.9% | 0.045 / -0.008 / 0.151 / 0.129 | 0.090 | oscillating, residual 6.2e-3 | 1382 | 5.73 |
| | | | | Coupled | 0.3094 | 0.3093-0.3094 | +19.0% | 0.045 / -0.007 / 0.149 / 0.123 | 0.086 | continuity plateau 3.5e-3 | 1490 | 6.78 |
| L3.5 | 13,332 | 2,108,192 | 42 | SIMPLEC | 0.3031 | 0.3028-0.3036 | **+16.6%** | 0.046 / -0.008 / 0.144 / 0.121 | 0.070 | residuals 6.2e-4, Cd drifting 0.1% | 2075 | 7.93 |
| | | | | Coupled | 0.3031 | 0.3030-0.3031 | **+16.6%** | 0.046 / -0.008 / 0.144 / 0.121 | 0.071 | continuity plateau 2.2e-3 | 2461 | 9.36 |

Body y+ is 0.82-0.88 mean and 1.64-1.70 maximum on every level, against
the 50 um meshes' 3 and 7. The peak footprint is 3.8-4.1 KB per cell under
SIMPLEC and 4.4-4.8 KB under Coupled. L3.5, at 2.1M cells and 9.4 GB, is the
largest level that fits the ~12 GB budget. By extrapolation, L4 would be about
3M cells and over 12 GB, so it was not run.

**Against the experiment.** Cd still misses 0.26 by 16.6% at L3.5 and is
still falling: -0.022, -0.018, then -0.011 for half a level. The series is
not in its asymptotic range, so no extrapolated value is claimed and the 8%
gate is not met. **The y+ ~ 1 wall spacing is not the cause of the excess.**
At the same level, the 50 um meshes gave 0.3531, 0.3322 and 0.3137 (FSAE aero
report), within 0.4% of these values. Friction stays at 0.044-0.046, close to
the 0.055 the reference components imply. The excess is still slant pressure
(0.144 against 0.090), which falls at every level, and base pressure (0.121
against 0.095), which falls only from L2 on (0.129, 0.137, 0.129, 0.121). At L3.5 the base Cp mean matches the experiment (-0.216
against -0.220, RMS 0.042 over 99 taps), and the slant suction is still too
strong (-0.300 against -0.266, RMS 0.109 over 188 taps).

**SIMPLEC against Coupled.** At L1 and L3.5, where SIMPLEC settles, the two
agree to 0.1% and 0.05% in Cd. At L2 and L3, SIMPLEC does not settle: its Cd
oscillates by 2-3% with a raw continuity residual of about 2e-2. Coupled holds a
steady value 1.1-1.6% below SIMPLEC's mean, with almost the same components
except base pressure. SIMPLEC oscillated the same way at 50 um on L2 and
L3 (FSAE aero report). Coupled runs to the iteration cap on its continuity
plateau, not on a
diverging residual (open item 1).

**Wake topology**, at the centreline stations x = -163, -123, -83, -43 and
-3 mm (R reversed, A attached; experiment and CFX-5 SST are R R R R R):

| level | SIMPLEC | Coupled | reversed share of slant area |
| --- | --- | --- | --- |
| L1 | R R R R R | R R R R R | 0.86 / 0.85 |
| L2 | R R R R A | R R R R A | 0.87 / 0.87 |
| L3 | R R R R A | R R R R A | 0.88 / 0.88 |
| L3.5 | R R R R A | R R R R A | 0.89 / 0.89 |

This is the 50 um result again: a fully separated slant except the last
station, 3 mm from the base edge.

## Ahmed 25 deg topology vs CFX-5 SST

Level 2, 15 um (762,644 cells, 1,227 flagged). Topology at the same stations,
x = -163 ... -3 mm:

| source | topology | reversed share of slant area | Cd (experiment 0.285) |
| --- | --- | ---: | --- |
| experiment (Lienhart & Becker) | R R A A A (separates, reattaches by -83) | | 0.285 |
| CFX-5 SST (CFD11) | R R R R R | | |
| Kutta SIMPLEC, 50 um (FSAE aero report) | R R R R R | 0.62 | 0.321 |
| **Kutta SIMPLEC, 15 um** | **R R R R A** | 0.64 | 0.3182 (+11.7%), converged at 2,850 iterations, residual 4.1e-6 |
| **Kutta Coupled, 15 um** | **R R R R R** | 0.57 | 0.3221 (+13.0%), continuity plateau 1.5e-2 |

At y+ ~ 1 Kutta still has the CFX-5 SST topology: the slant separates and
does not reattach, unlike the experiment, which most steady RANS models miss
(ERCOFTAC evaluation, FSAE aero report). The one difference is SIMPLEC's last
station, attached 3 mm from the base edge, as at 35 deg. Rear Cp shows the
same picture. The slant mean is -0.460 against the experiment's -0.556 (the
reattaching bubble's suction peak is missing) and CFX's -0.320. The base is
-0.241 against -0.174 (experiment) and -0.175 (CFX). Coupled ends 1.2% above
SIMPLEC because it has not converged continuity (open item 1); the two
schemes' slant and base pressure components differ by 0.003 and 0.000.

## Tet duct and tet channel under Coupled

CPU, 400 iterations (the coupled-external lane's setups: `duct.py tet`, Re
1e5; `bench_convergence --case channel`).

| case | main Coupled | candidate Coupled | candidate SIMPLEC |
| --- | --- | --- | --- |
| tet duct, 9,729 cells | diverges between 100 and 150 (residual 8e24) | bounded, not converged: max residual 0.95, Cd 0.0165 | not converged: max residual 2.4, Cd 0.0163 |
| tet channel | stalls: scaled continuity 5.00, raw 1.45e-2 | stalls: 4.99, raw 1.45e-2 | scaled continuity 0.54, raw 1.1e-3 |

The `a_P` floor alone stops the duct diverging (the diagnosis build's floor
switch did the same), so it is mechanism 2 again; I did not locate its first
cells. Neither
case converges under Coupled. The channel is unchanged: the
coupled-convergence report already records it unconverged in both schemes
with oscillating raw continuity, and the Coupled report traces it to the
Rhie-Chow jump on its skewed tetrahedral wall layers. Neither mesh has a
flagged cell, so the pressure-face change is not involved. The duct does
not converge under SIMPLEC in 400 iterations either, so it is not a
Coupled-only defect there.

## Well-conditioned cases

The two parts of the fix reach different cases. The pressure-face change
runs only where a cell has orthogonal quality <= 0.05, under second-order
pressure. The `a_P` floor runs on every Coupled solve, but changes only cells
whose diagonal is below its off-diagonal sum.

**Byte-for-byte** (`parity.sh`, `kutta solve --iterations 40`, CPU and Metal;
fluxes, forces, loads, VTU and residual history):

| case | scheme | CPU | Metal |
| --- | --- | --- | --- |
| V6 mixing elbow | SIMPLEC | identical | identical |
| elbow_param | SIMPLE | identical | identical |
| V6 mixing elbow | Coupled | differs | differs |
| nozzle | Coupled | differs | differs |

Both Coupled cases have sink cells, so the floor changes them. At 40
iterations, the elbow's scaled continuity goes from 0.1016 to 0.0999 and the
nozzle's from 1.545e-3 to 1.517e-3 (CPU). Their converged effect is
quantified in the slow tier below.

**Wing in ground effect, Coupled** (Zerihan h/c 0.224, level 2, 195,581
cells, no flagged cells; `bench_fsae_aero --re 4.62e5 --downforce
--moving-ground --tu 0.001 --ratio 1 --coupled`, uniform start):

| build | Metal | CPU, 60 iterations |
| --- | --- | --- |
| main `f770a31` | non-finite before iteration 50 (3 runs, 2 option sets) | residual 4.9e16 at 50 |
| candidate | residuals met at 800 iterations, 259 s: lift 1.401506, drag 0.030516 | bounded, still settling |

The coupled-external report gives this point as converged under Coupled in
750 iterations with lift 1.401503 and drag 0.030511, against SIMPLEC
1.401467 and 0.030513. The candidate reproduces those forces to 3e-6 and
5e-6. **Main as committed does not**: the f770a31 sources in the baseline
build match the commit byte for byte, and the mesh has the same cell count.
I did not find why that report's runs converged. With no flagged cells, the
residual-side `a_P` floor is the only change that runs on this mesh, so it is
what makes the case converge here.

**Slow tier** (`slow-main.xml.gz` and `slow-cand.xml.gz`, compared with
`tools/forensics/compare_slow_assertions.py`; `slow-comparison.json.gz`):
32 common cases, **0 regressions, 0 recoveries, 0 missing assertions**, 7 new
assertions in 1 new case (the y+ ~ 1 Ahmed test). Main has 39,538 passed and
65 failed assertions; the candidate has 39,545 and 65. The known failures are
the same assertions on both sides: V6 25, V2 8 (including tet 8->16 velocity
order), V5 plate 6, V5 bump 9, V7 SA 16, V12 polynomial fan 1. V8,
Taylor-Couette, the Fluent coverage cases and the FSAE aero tests pass on
both sides.

73 asserted values moved, none across its bound:

| case | why it moves | largest change |
| --- | --- | --- |
| V6 mixing elbow (Coupled), 35 values | `a_P` floor | profile error 0.00524 -> 0.00518 (bound 0.005, known failure both sides); pressure drop error 0.01969 -> 0.01967; raw continuity at the fixed point 2.080e-3 -> 2.069e-3 |
| nozzle (Coupled), 11 values | `a_P` floor | velocity error 0.016034 -> 0.016034, pressure error 0.0052699 -> 0.0052699 (relative changes ~1e-5) |
| V5 TMR bump SST, level 2, 12 values | flagged cells | Cf error at x = 0.75: 0.15118 -> 0.15128 (bound 0.1, known failure); Cp error changes <= 2.3e-4 absolute |
| V7 SA TMR bump, level 2, 14 values | flagged cells | Cd error 0.16173 -> 0.16144 (CPU, known failure); Cf and Cp error changes <= 9e-4 absolute |
| Ahmed 35 deg L1, 50 um | 10 flagged cells | Cd relative to the 0.3531 guard: 8.7e-5 -> 2.8e-5 |

elbow_param (SIMPLE), the TMR flat plate and every case not listed are
unchanged. **The TMR bump is not bitwise identical, and this was not
expected.** Its level-2 grid (353 x 161) has 2,068 of 56,320 cells at
orthogonal quality 0.036, and level 1 has 1,236 at 0.018. These are the
first-layer cells far upstream and downstream: 1.77 long, 8e-4 thick, and
trapezoidal, because the first-layer height changes by 7% across each cell.
That moves the centroid 0.011 along the wall, against a half-height of
4e-4, which gives cosine 0.036 between the wall-normal face and its
centroid vector. They are poor by Fluent's definition, so Fluent's
cell-quality-based option at 0.05 would flag them too. The effect
is below 1e-3 absolute on every asserted error.

## Tests

- `tests/mesh/test_mesh.cpp`: orthogonal quality of a sheared box equals
  `1/sqrt(1 + t^2)` for shears 0, 1 and 20.
- `tests/fv/test_flow_kernels.cpp`: the poor-cell kernel leaves faces bitwise
  alone without flags, and with a flag takes the flagged cell's pressure at
  interior and wall faces, keeps prescribed-pressure faces, and leaves the
  Rhie-Chow jump unchanged.
- `tests/fv/test_fsae_aero.cpp` (slow): Ahmed 35 deg L1 at 15 um, 2000
  iterations. SIMPLEC must finish without failure with residuals met, y+ mean
  < 1.5, and Cd within 2% of 0.3545 (regression guard: this run's value, not
  the experiment). Coupled must finish without failure, with Cd within 1% and
  lift within 3% of SIMPLEC. On main both solves go non-finite within 50
  iterations. In the slow tier it passed its 7 assertions with two solves of
  about 4.5 minutes each (Cd 0.35452 SIMPLEC, 0.35482 Coupled). The mesh comes from
  `run.py ahmed --points 35 --level 1 --first 1.5e-5` in
  `KUTTA_FSAE_AERO_MESHES`, as the other FSAE meshes do.

**Quick tier** (`kutta-tests --test-suite-exclude=slow`, fixtures present):
main 256 cases, 7,453,834 assertions; candidate 258 cases, 7,453,857
assertions, all passed. `tools/lint_portability.py` is clean.

**Slow tier**: as above, run concurrently for main and candidate
(`-r=xml --success`), about 2 hours each on the shared machine.

## Open items

1. **Coupled's continuity plateau on the Ahmed body.** Coupled no longer
   diverges, but its scaled continuity stops at 2e-3 to 1.5e-2 (raw about
   1e-3), so it never meets the residual criteria. Forces settle anyway. At
   35 deg L1 and L3.5 they match converged SIMPLEC to 0.1%, but at 25 deg L2
   Cd sits 1.2% above SIMPLEC's converged value. What was measured: the
   plateau does not move with the pseudo-time step, it scales with the
   explicit relaxation, and the block solve converges at every iteration.
   That points to a mode the Jacobian under-predicts, as in mechanism 2, but
   I did not find it. Part of the gap is the scaling: Fluent's scaled
   continuity divides by the largest imbalance of the first five iterations,
   which on 35 deg L1 is 124 under SIMPLEC (impulsive start) and 0.33 under
   Coupled. The same criterion therefore asks Coupled for a raw imbalance
   about 375 times smaller. In raw terms Coupled's plateau (1.5e-3) is 4
   times SIMPLEC's converged value (3.6e-4).
2. **Tet duct and tet channel under Coupled** do not converge (above). The
   duct no longer diverges.
3. **Ahmed 35 deg Cd** is +16.6% at the largest level that fits, and still
   falling. It needs a finer mesh than this machine holds, or a local
   refinement of the slant and wake, to say where it lands.
4. **The coupled-external report's wing result** (h/c 0.224 converged under
   Coupled on main) does not reproduce from `f770a31`. This branch converges
   it to the same forces.
5. **Threshold default (owner call).** Fluent's default marks only
   quality-0 cells. Kutta defaults to Fluent's documented 0.05 option, which
   also reaches the TMR bump's thin far-field cells (changes <= 9e-4 on its
   asserted errors). `poor_cell_quality = 0` would restore the bump bit for
   bit, but the Ahmed body at y+ ~ 1 would diverge again by default. I
   recommend keeping 0.05.
6. Fluent's correction is not documented term by term. Kutta changes only
   the pressure reconstruction in flagged cells (see Deviations).
7. CUDA: both new kernels go through `kernelgen.py`. Locally they passed a
   syntax check against CUDA stubs only; CI compiles them.

## Reproduce

Meshes need gmsh in `build-py` (see `tools/fsae_aero/run.py`); `<main
build>` is a Release build of `f770a31`.

```sh
# Ahmed meshes: L1 at 15 and 25 um, and the study levels (run.py makes them on demand)
build-py/bin/python tools/fsae_aero/ahmed.py --slant 35 --level 1 --first 1.5e-5 --out runs/mesh/ahmed35_L1_first15um.msh
# reproduction table (add --coupled; CPU: --backend cpu --iterations 60)
build/bench_fsae_aero --mesh runs/mesh/ahmed35_L1_first15um.msh --re 2.78e6 --length 1.044 --area 0.112 \
  --walls body --lift-z --tu 0.01 --window 1000 --iterations 300 --out runs/repro/cand_simplec_15
# mesh study and tables
python3 tools/fsae_aero/run.py ahmed --points 35 --level 2 --first 1.5e-5 --iterations 5000 --max-seconds 5400 --out runs/study
python3 tools/fsae_aero/run.py ahmed --points 35 --level 2 --first 1.5e-5 --iterations 5000 --max-seconds 5400 \
  --extra=--coupled --tag _coupled --out runs/study
build-py/bin/python validation/results/external_3d_robustness/analyze.py runs/study
# well-conditioned parity
validation/results/external_3d_robustness/parity.sh <main build>/kutta build/kutta runs/parity
# slow tier, per build, then the comparison (addresses and continuity values in INFO are masked)
KUTTA_FSAE_AERO_MESHES=<dir with the FSAE meshes and ahmed35_L1_first15um.msh> \
  ./build/kutta-tests --test-suite=slow -r=xml --success > slow-cand.xml
for t in main cand; do sed -E 's/0x[0-9a-f]+/0xADDR/g; s/continuity=[0-9.eE+-]+/continuity=X/g' slow-$t.xml > m-$t.xml; done
python3 tools/forensics/compare_slow_assertions.py m-main.xml m-cand.xml > slow-comparison.json
```

`evidence.tar.gz` holds the reproduction logs and histories, the cell dumps
and mesh-quality probes, the diagnostic switches (`diag.patch`,
`diag2.patch`, not committed code), the ablation logs and the study's
histories. It contains no meshes, cell fields or reference data.
