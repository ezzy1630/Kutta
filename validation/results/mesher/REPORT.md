# External-aero mesher: castellate, snap and prism layers

A snappyHexMesh-style pipeline replaces the cut-cell mesher (minimum
orthogonal quality 0.003, diverging solves). Every mesh below has zero invalid
cells and orthogonal quality (OQ) ≥ 0.15, and both coarse cases run 5000 Metal
iterations without diverging. **Neither solve meets the qualification
criteria**: Ahmed stalls on continuity with Cd 0.454, far above experiment;
FSAE forces are steady (CL 1.46, CD 1.33, 23% front) but k runs away in a few
stagnant cells. The default FSAE mesh passes the gates; its layer stage took 61 min.

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
