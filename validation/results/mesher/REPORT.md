# External-aero mesher: castellate, snap and prism layers

A snappyHexMesh-style pipeline replaces the cut-cell mesher (minimum
orthogonal quality 0.003, diverging solves). Every mesh below has zero invalid
cells and orthogonal quality (OQ) ≥ 0.15, and both coarse cases run 5000 Metal
iterations without diverging. **Neither solve meets the qualification
criteria**: Ahmed stalls on continuity with Cd 0.454, far above experiment;
FSAE forces are steady (CL 1.46, CD 1.33, 23% front) but k runs away in a few
stagnant cells. The default FSAE mesh passes the gates; its layer stage took 61 min.

## Round 3 (2026-11-13): convergence, layer time, default size, contract

**Status against the round 3 goal** (default or coarse `kutta mesh` → default
`kutta solve` converges, fits 12 GB, meshes in minutes): meshing in minutes
and fitting 12 GB are met; **convergence is not**. Only Ahmed coarse converges
(as it did before). The default FSAE mesh used to take 53 min and now takes
under 5 and runs 5,000 iterations in 8.55 GB without diverging (without the
notch removal it diverged by iteration 31), but both FSAE meshes end at the
iteration limit with k stalled (residual ~1–3) and forces cycling. The six app
contract issues are fixed.

Code: commits on `t3/mesh-convergence-memory` from e08ec08 (see `git log`);
every solve below
uses main's starter (SIMPLEC, momentum 0.9, statistically-steady stop, 5000
cap), the solver unchanged. Convergence is reported by the solve's `verdict`.
Machine and contention as in round 2: another lane's app scenario and up to
three of these solves shared the GPU, so times are contended.

### Before / after

"Before" is the round 2 mesher (fb44554) on the same STL; its meshes were
re-solved here with main's starter so only the mesh differs.

| | Ahmed coarse | Ahmed default | FSAE coarse | FSAE default |
|---|---|---|---|---|
| cells, before → after | 144,877 → 144,919 | 148,705 → 122,151 | 1,312,313 → 1,391,056 | 4,417,972 → 4,224,992 |
| mesh time (layer stage), before | 2.0 s | 3.3 s | 64 s (41 s) | 3,197 s (3,097 s) |
| mesh time (layer stage), after | 1.7 s | 3.8 s | 46 s (21 s) | 286 s (189 s) |
| mesher max RSS, before → after | 0.16 → 0.20 GB | 0.17 → 0.14 GB | 1.12 → 1.40 GB | 3.65 → 3.59 GB |
| verdict, before | `converged`, 621 it | `diverged`, iteration 36 | `statistically_steady`, 1,940 it | not solved (5000 iterations on 4.4M cells unaffordable) |
| verdict, after | `converged`, 621 it | `statistically_steady`, 1,600 it | `iteration_limit`, 5,000 it | `iteration_limit`, 5,000 it (5,028 s) |
| final residuals continuity / k / ω, before | 8.6e-4 / 2.3e-5 / 1.6e-5 | (diverged) | 6.1e-2 / 2.1 / 0.93 | |
| final residuals continuity / k / ω, after | 8.6e-4 / 2.5e-5 / 1.6e-5 | 7.5e-3 / 2.9e-4 / 4.5e-5 | 9.2e-2 / 3.3 / 0.50 | 8.4e-2 / 1.0 / 0.50 |
| forces, before | Cd 0.4535, CL +0.448 | | mean CL 1.413 (0.916–1.477), CD 1.332 | |
| forces, after | Cd 0.4517, CL +0.447 | mean Cd 0.4469 (0.4466–0.4471), CL +0.056 | mean CL 1.359 (0.649–1.474), CD 1.280 (0.823–1.721) | last CL 1.219, CD 1.219; window mean CL 1.229 (1.090–1.577), CD 1.232 (1.202–1.505) |
| front-axle share of downforce, before → after | | | 21.5% → 21.8% | 14.6% |
| device memory (B/cell), after | 0.29 GB (2,031) | 0.25 GB (2,063) | 2.91 GB (2,094) | 8.55 GB (2,025) |
| estimated solver memory, after | 0.30 GB | 0.26 GB | 2.92 GB | 8.87 GB |

Before/after FSAE coarse use `fsae-coarse.json` from each round's
`generate_geometry.py` (round 2: 3 mm gap floor, 1 cell across; round 3: 4 mm,
2 cells, which opens front slot 1–2). Ahmed reference: A = 0.112032 m²,
L = 1.044 m, 40 m/s, stationary ground (`{"type":"wall","motion":"stationary"}`,
not a report zone). FSAE: A = 1 m², 20 m/s, moving ground, rotating wheels;
front-axle share from `kutta aero --aero` with the front contact patch at
x = −0.78 m and a 1.73 m wheelbase. Ranges are min–max over the
statistically-steady stop's window. CL for Ahmed is lift (positive up), for
FSAE downforce.

**The round 2 coarse tuple with this round's code diverges** (iteration 3,082;
rear slot 1–2 closed, front slot 2–3 one cell across), which is why the coarse
setup changed.

### Layer collapse: minutes, not an hour

Each collapse round re-ran the whole-mesh displacement smoothing (Laplacian
over a 3-point band, 12 Jacobi sweeps, quality retreat). The smoothing is
linear in the source displacements, so a round now smooths only the change:
an exact delta solution over the 12 band hops reached from the sources that
changed (`Displacement::run`, `displace.cpp`). Scales reset only within 6 hops
of a changed source; when a retreat would halve a source outside that set,
the set grows and the round restarts, otherwise compounding halvings cascade
(an earlier attempt without the restart left 13.5% of wall unlayered).
Settle re-examines only faces on edges whose use changed; its fixed point is
the same (byte-identical mesh in the check run). A unit test checks the
incremental result against a fresh displacement of the same target
(< 1e-12 m). Default-size FSAE: layer stage 3,097 s → 189 s, whole mesh
3,197 s → 286 s (contended); max RSS 3.59 GB against 3.65 GB.

### Default sizing and memory

Omitted `refinement_boxes` now add a `body` box (2 × `max_size`, 0.1 L around
the body, 0.25 L behind) and a `wake` box (4 × `max_size`, 1 L behind), and
the report lists both with their achieved octree cell. `min_size` defaults to
0.008 m (octree 7.8 mm, was 0.006 m) with the gap floor unchanged at
`min_size/4` (1.95 mm), so every 12–20 mm wing slot keeps two or more cells.
Synthetic FSAE car: 4,224,992 cells (2.49M core, 1.74M prism), estimated
8.87 GB for the default solve, under the 10 GB warning.

**A coarser default diverged.** `min_size` 0.016 m with a 2 mm gap floor gave
4.08M cells (only 3% fewer) and **diverged at iteration 31**: k 5.8e8, ω at
its 1e-20 floor and μt at its cap in 2 mm cells where front-wing slot 1–2
(13 mm, 4 cells across) meets the left endplate, around (−1.204, 0.554,
0.159) m; the worst cell is nearly stagnant (|u| 0.4 m/s). It is not the
starter: the same mesh at momentum 0.7 diverges at iteration 32. The
0.008 m mesh (2 mm slot cells next to 7.8 mm, not 15.6 mm, surface cells)
does not. A 4 mm floor with 0.016 m (2.43M cells) ran 60 iterations but
closes rear slot 1–2, so it was not adopted. The 2.5–3.5M-cell target is not
met by a mesh that runs. `estimated_solver_memory` is cells × 2,100 B (SIMPLEC, measured
2,003–2,232 B per cell on Metal for this mesher's meshes; Coupled 2,600,
measured 2,528–2,589 in round 2), with a warning above 10 GB. It is not
measured on CUDA, whose context adds memory.

### The k runaway: what the mesh does and does not explain

- **Where.** On FSAE coarse, k runs to 1e11–1e13 with μt at its cap in a few
  stagnant prism cells where a closed or one-cell slot remnant meets an
  endplate (rear slot 1–2 at about (0.996, 0.47–0.56, 0.958) m; front slot 1–2
  at the endplate). Those cells are walled on three or four sides, with
  near-zero convective flux. Forces cycle by ±0.5% and continuity stalls near
  3e-2. On the 4.2M default-like mesh (all slots open, old starter) k still
  stalled at residual ~1 after 3,000 iterations, with downforce and drag ranges
  of 4% and 2%.
- **Experiment 1, close one-cell gaps** (`round3-thin-gap-experiment.patch`):
  remove every fluid cell with the body on two opposite sides (and, in the
  second rule, on at least three sides), repeatedly. It removed the slot
  remnants but **diverged at iteration 37** under main's starter: k 6.8e13,
  ω down to 1e-20 and μt at its cap in OQ ≈ 0.98 layer cells at front-wing
  element 1's trailing edge. The same setup without it runs to the 5000
  limit. **On the coarse car it is not applied** (below).
- **Experiment 2, gap-limited layers** (cap total layer height at a fraction of
  the local gap): diverged at about iteration 60. Not adopted.
- **At the default size the same removal is what lets the solve run.**
  Without it, the default car diverges by iteration 27–31 with k 3e10 and
  ω at its floor in 2 mm cells at front-wing slot 1–2 and the left endplate,
  around (−1.20, 0.554, 0.169) m (only the left side: the octree lattice is
  not symmetric there). With it (9,318 cells removed, every slot still
  2–4 cells across) the solve runs on (table). So the removal is applied only
  when `proximity_cells` ≥ 3 (the default): then a one-cell-thick fluid cell
  is a notch, not a slot. Below 3 the setup asks for thin slots, and closing
  them diverges on the coarse car. This rests on these two cars, not a
  general argument.
- **Default sizing** keeps every slot open (the round 2 suggestion).

**Candidate solver defect (not changed; the solver belongs to another lane).**
Reproduction from this commit:

```sh
git apply validation/results/mesher/round3-thin-gap-experiment.patch   # removal at any proximity_cells
cmake --build build --target kutta
uv run tools/mesher/generate_geometry.py --out g
./build/kutta mesh g/fsae.stl --setup g/fsae-coarse.json --out m --case --force \
  --speed 20 --wheel wheel-fl=0.23,-0.78,0.72,0.23,0,-1,0 \
  --wheel wheel-fr=0.23,-0.78,-0.72,0.23,0,-1,0 \
  --wheel wheel-rl=0.23,0.95,0.72,0.23,0,-1,0 --wheel wheel-rr=0.23,0.95,-0.72,0.23,0,-1,0
./build/kutta solve m/case.kutta.json --backend metal --iterations 100 --progress jsonl </dev/null
```

Rerun on this commit with the patch: 1,385,567 cells (2,219 removed as thin),
0 cells below OQ 0.15, and the solve stops with `verdict: "diverged"` at
iteration 37 ("nonfinite input vector"; last finite residuals at iteration
20: continuity 0.979, k 0.122, ω 0.005). The diverging cells pass every mesh gate (OQ ≈ 0.98, positive volume);
k-ω should not blow up there, so the next step is in the solver's k-ω
source/limiter treatment, not the mesh.

### Contract fixes for the app (all v1 additive, documented in app-protocol.md)

1. **Cell types vs `kutta inspect`.** `cell_types` counts by origin: octree
   cubes with 2:1 hanging faces are `hex` (inspect: poly, more than six
   faces) and layer cells are `prism` (inspect: hex, six quads). It is kept;
   `cell_types_by_shape` (tet / hex / poly) uses inspect's classifier, six
   quadrilateral faces hex, four triangular faces tet, otherwise poly. On an
   FSAE mesh it equals `kutta inspect`'s counts, and a unit test recounts it.
2. **`proximity_min_size: 0` rejected.** It is now accepted and means
   `min_size/4`, as the omitted default always did.
3. **Progress fraction went backwards** (0.80 then 0.78 when layers restarted
   a stage). Progress is now a running maximum across stages: nondecreasing
   over all 41 events of an FSAE mesh, and in a unit test.
4. **Tyre–ground contacts.** Each `castellation/gaps` entry adds
   `contact: true` when the zones touch (width ≤ 1e-6 of the domain diagonal);
   FSAE reports 4, one per tyre.
5. **Histogram key.** Volume-ratio histograms count faces; they add `faces`
   and keep `cells` as a documented deprecated alias with the same counts.
6. **`skewness_max` 1.36.** The field reported OpenFOAM-style face skewness:
   the face centroid's distance from where the owner–neighbour centroid link
   crosses the face plane, over the link length. That is unbounded; sheared
   thin layer cells exceed 1. It moves to
   `face_skewness_max`; `skewness_max` is now Fluent's equiangle skewness,
   max over faces of max((θmax−90)/90, (90−θmin)/90) on each face's corner
   angles, in 0..1. Hanging nodes on straight edges are not corners (an
   interior angle ≥ 150° is skipped), so 2:1 faces are not flagged.
   `skewness_method` names it and `faces_skewness_above_0_95` counts the
   worst. FSAE coarse: `skewness_max` 1.0 with 47 faces above 0.95
   (default: 137), slivers at snapped trailing edges; `face_skewness_max`
   1.53 (default 2.07).

### Open problems (round 3)

1. **Convergence.** No FSAE mesh converges under the default solve. Coarse
   cycles, downforce coefficient 0.65–1.47 over the last window (k residual
   ~3); default 1.09–1.58 (k residual 1.0), with the front-axle share 14.6%
   against coarse's 21.8%, so the balance is not mesh-independent either. Ahmed default stops statistically steady, not converged.
   The diverging and stalling cells pass every mesh gate, so the next lever is
   the k-ω treatment in stagnant junction cells (solver lane), or a mesher
   change that avoids 2 mm cells meeting 3–4 walls at slot/endplate
   junctions.
2. **Candidate solver defect** above (k-ω blow-up in a good layer cell).
3. **Volume-ratio maxima and snap error** (round 2 item 4): not attempted.
4. **Found by the post lane** (`validation/results/mesher_post/REPORT.md`, not
   fixed here): horizontal ledges on the snapped Ahmed body at octree levels
   (2.1e-4 m² at z = 0.1); a non-manifold fin 8 mm outside the Ahmed side,
   the source of the 7.84 mm snap-error maximum; and about 100 tangled cells
   on FSAE coarse that cannot be sliced (inverted fan tetrahedra).
5. `tools/mesher/qualify.py` still applies round 2's stricter force criterion
   (5e-4 over 50 samples); its Ahmed ground now uses the stationary-ground
   form.

The sections below are round 2 (38b1389), kept for reference.

## Reproduce

```sh
cmake -S . -B build -G Ninja && cmake --build build
uv run tools/mesher/generate_geometry.py --out scratch-geometry
KUTTA_MEMORY_PHASES=1 uv run tools/mesher/qualify.py --kutta build/kutta \
  --geometry scratch-geometry --out scratch-qualification --iterations 5000
./build/kutta mesh scratch-geometry/fsae.stl --out scratch-default --case --json
```

`qualify.py` runs `mesh`, `new-case` (FSAE wheels rotating, from
`fsae-wheels.json`) and `solve` (Coupled, SST, second order; WFGC off, as it is
not implemented) and exits 1 unless the solve converges: scaled residuals
≤ 1e-3 and drag and downforce changing < 5e-4 over 50 samples after iteration
500, cap 5000. Ahmed: 40 m/s, A = 0.112032 m², stationary ground; FSAE: 20 m/s,
A = 1 m², moving ground. Apple M4 Pro (24 GB), Release/Metal, configured with
`-DCMAKE_OSX_SYSROOT=$(xcrun --sdk macosx26.5 --show-sdk-path)` (not committed).
GPU and CPU were shared, so **all times are contended**; no CUDA machine was
available. Hashes, settings and results: [receipt.json](receipt.json).

## Method

**Castellate** removes octree leaves inside or crossed by the STL (2:1 hanging
faces). **Snap** moves wall points to the least-squares intersection of their
faces' tangent planes (implicit feature capture), decaying over three point
layers, and retreats wherever a cell would fall below OQ 0.15, so geometric
error is reported, not hidden. **Layers** push the snapped mesh off the walls
and fill the gap with stacks that collapse locally where quality fails;
terminations split so unlayered walls stay on the surface, or stay lifted.

## Mesh quality

| | Ahmed coarse | Ahmed default | FSAE coarse | FSAE default |
|---|---|---|---|---|
| cells | 144,877 | 148,705 | 1,312,313 | 4,417,972 |
| hex / poly / prism | 93,313 / 474 / 51,090 | 91,061 / 3,312 / 54,332 | 524,228 / 6,069 / 782,016 | 2,407,850 / 31,750 / 1,978,372 |
| hex share of core | 99.5% | 96.5% | 98.9% | 98.7% |
| invalid cells; OQ min | 0; 0.154 | 0; 0.161 | 0; 0.150 | 0; 0.150 |
| cells below OQ 0.15 / 0.2 | 0 / 7 | 0 / 37 | 0 / 93 | 0 / 689 |
| core volume ratio p99 / p999 / max | 9.6 / 9.6 / 23 | 10.1 / 36.6 / 151 | 9.5 / 19.5 / 124 | 9.6 / 19.2 / 159 |
| core faces with ratio > 20 | 1 of 279,523 | 380 of 288,662 | 1,226 of 1.53M | 6,679 of 7.21M |
| layer-face ratio p999 / max | 6.7 / 51 | 17.1 / 295 | 22.1 / 378 | 18.9 / 383 |
| snap error p99 / max | 0.1 / 7.8 mm | 5.3 / 12.3 mm | 2.8 / 13.5 mm | 2.7 / 15.4 mm |
| mesh time; peak RSS | 2.0 s; 108 MB | 6.6 s; 173 MB | 90 s; 1.0 GB | 3,768 s; 3.1 GB |

OQ histogram (0.15–0.2 / –0.3 / –0.5 / –0.7 / –0.9 / –1; none below 0.15): Ahmed
coarse 7 / 70 / 651 / 1,441 / 23,822 / 118,886; FSAE coarse 93 / 1,240 /
11,545 / 31,796 / 163,783 / 1,103,856. All four are in the receipt.

- **Volume ratio** outside layers is within ~20 at p999 except Ahmed default
  (36.6); the maxima (23–159) are not, where snapping squeezes a fine cell
  beside a 2:1 coarser one. A snap-time neighbour-size check cut the Ahmed
  default maximum to 37 but broke exact cube-corner capture, so it was dropped.
- **Feature edges**: cube corners on and off the lattice are exact to 1e-9 m; on
  the cars the tolerance is the snap error (p99 ≤ 5.3 mm, max 15.4 mm), at sharp
  trailing edges and concave junctions where full snapping would break OQ 0.15.
- **Slots** (12–20 mm): FSAE default has 3–4 cells across each; the coarse setup
  asks for 1, and its **rear-wing slot 1–2 (12.8 mm) closes**.

## Prism layers

`layers` takes `count`, `first_height`, `growth` and `zones` (default 4 layers,
2 mm, 1.3: 12.4 mm). Coverage is a fraction of the selected wall area.

| | full / thinned / none | first height mean (min) | split / lifted points |
|---|---|---|---|
| Ahmed coarse (3 layers) | 94.9 / 1.5 / 3.5% | 2.0 mm (0.5) | 30 / 844 |
| Ahmed default | 75.7 / 12.6 / 11.7% | 1.7 mm (0.5) | 364 / 2,476 |
| FSAE coarse | 76.2 / 16.9 / 6.9% | 1.25 mm (0.5) | 4,587 / 11,340 |
| FSAE default | 58.0 / 24.0 / 18.0% | 1.0 mm (0.5) | 21,817 / 37,119 |

**FSAE y+ regime: wall functions, y+ ≈ 50.** A 2 mm first layer gives an
estimated y+ p50 of 52 at 20 m/s, 83–84% of wall area in 30–300 and none below
5 (flat-plate estimate), for Kutta's y+-insensitive SST wall treatment. y+ ≈ 1
needs 20–40 µm first cells and ~25 layers: over twice the cells, past 12 GB.

## Solves (Metal, coarse meshes)

| | Ahmed coarse | FSAE coarse |
|---|---|---|
| converged | no: continuity 1.7e-2 (others ≤ 6e-4) | no: k 0.94, continuity 3.6e-2, ω 0.13 |
| iterations; solve time | 5000; 974 s | 5000; 4,106 s |
| final coefficients | Cd 0.454 (pressure 0.412, friction 0.042), CL +0.451 | CL 1.460, CD 1.330, front axle 23.3% |
| range, last 1000 iterations | Cd 0.4531–0.4544 | CL 1.445–1.465, CD 1.325–1.334 |
| device memory | 366 MB, 2,528 B/cell | 3.40 GB, 2,589 B/cell |

**Ahmed**: Cd 0.454 is 75% above experiment (0.26) and above Kutta's
body-fitted study (0.353 / 0.332 / 0.314 on 0.40M / 0.69M / 1.33M cells, y+ ≈ 3);
with ~13 mm wall cells and y+ ≈ 108 this is a pipeline check, not accuracy.

**FSAE**: forces settle within ±0.7% but the run is unconverged, so not
accepted. A = 1 m²; balance is the front-axle share of downforce from the
pitching moment about the rear contact patch (axles at x = −0.78, 0.95 m). The
first unstable field is **k**: 1.7e11, with μt at its 1.79 Pa·s cap, in stagnant
prism cells (|u| < 1 m/s) at (0.93, 0.545, 0.95) m, where the closed rear-wing
slot 1–2 meets the left endplate. Before the tyre fix it sat where the front
wing entered a tyre (CD spiked to 16.6). Velocity, pressure and forces stay
bounded: a local turbulence stall, not a divergence. The solver is unchanged.

## Default FSAE mesh and memory

Defaults give 4,417,972 cells (2.44M core, 1.98M prism). Octree, castellate and
snap take 92 s; **layer collapse took 229 rounds and 3,660 s**: a collapse
front moving along long edges re-runs the whole-mesh displacement each round.
A 10-iteration Metal solve used **11.42 GB of device memory, 2,584 B per cell**
(13.8 GB footprint; phase marks, not a guaranteed peak); the coarse meshes used
2,528–2,589 B per cell and the default without layers 2,656.

**12 GB RTX: probably not.** 11.4 GB leaves no room for a CUDA context; at this
cost about 4.3M cells is the limit (estimated from Metal, not measured on
CUDA). The default without layers fits (2.44M cells, 6.5 GB measured).

**Polyhedra and hanging faces cost little**: storage is face-based, so memory
follows faces and cells, not cell type. FSAE coarse has 3.17 faces per cell
(about 3.0 conformal, +6%), 97.3% quad faces (4.03 nodes per face) and 95.7%
six-faced cells; the ~3% with 9 or more are coarse cells with hanging sides.

## Fixes and root causes in this lane

The mesher is f51b228. Later commits stop layer collapse spinning when only
distant stacks move a cell, raise its round guard from 50 to 1000, and stop
counting merge-collapsed cells as hexes (aa02e5f, 20c4fb3, 38b1389). FSAE
divergences fixed: layer "chimney" stacks walled on 3–4 sides; the tangent
wheel–ground wedge (10 mm plinths, rotating wheels); disconnected pockets; and
(20c4fb3) a synthetic car whose wing, floor and sidepods passed through the tyres.

## Open problems

1. **FSAE k runaway** in stagnant layer cells at closed slots and junctions.
   Next: the default mesh with open slots (5000 iterations on 4.4M cells were
   not affordable here) and k-ω source linearization at tiny convective flux.
2. **Ahmed** continuity stall; Cd needs a mesh study on this mesher.
3. **Layer-collapse time** (61 min): recompute the displacement only near
   dropped stacks, or taper terminations so retreat stops advancing per round.
4. **Volume-ratio maxima** beside 2:1 transitions (a fix must keep exact
   corner capture); snap error up to 15 mm where retreat protects quality.
5. **Coarse FSAE closes rear-wing slot 1–2**: use `proximity_cells` ≥ 2.
6. The app on main expects a tunnel zone `sides`; the mesher writes
   `side-low` and `side-high` (app unchanged).

## Verification

Quick tier (`ctest -LE slow`, fixtures linked locally) passed at 38b1389: 326
cases, 7,996,651 assertions, and the portability lint. Focused mesher tests
cover castellation on and off the lattice, corner capture, layer volume,
gap sizing, curved-body snap and the CLI. FSAE runs and default meshes used the
38b1389 binary; the Ahmed solve used one built from the aa02e5f tree, whose
mesh is byte-identical at 38b1389 (solver sources unchanged in this lane).
