# External-aero mesher: partial implementation, qualification failed

The STL → mesh → case → solver path exists. Both generated geometries mesh with
zero invalid cells, and a tiny generated case solves for three CPU iterations
in the quick tier. **The requested aerodynamic acceptance is not achieved:**
Ahmed and FSAE Metal runs diverge; there are no accepted Cd, CL or balance values.
Prism layers, solver-quality cut-cell conditioning and the 12 GB target remain.

## Reproduce

```sh
cmake -S . -B build -G Ninja
cmake --build build
uv run tools/mesher/generate_geometry.py --out scratch-geometry
KUTTA_MEMORY_PHASES=1 uv run tools/mesher/qualify.py --kutta build/kutta \
  --geometry scratch-geometry --out scratch-qualification --iterations 5000
./build/kutta mesh scratch-geometry/fsae.stl --out scratch-default --case --progress jsonl
```

The qualification script deliberately exits **1** when the solve fails or does
not converge. It calls `mesh`, the separate `new-case` command and `solve`, and
keeps failure receipts. It sets Ahmed to 40 m/s, area 0.112032 m², length 1.044 m
and stationary ground; FSAE uses 20 m/s, moving ground and reference area 1 m².
Qualification wheels remain stationary. Both explicitly disable WFGC because generated faces are planar. Otherwise they
keep starter SST, Coupled and second-order transport defaults. Criteria are all
scaled residuals ≤1e-3 and drag/downforce relative changes <5e-4, over 50 monitor
samples after iteration 500 (20 iterations per sample), with a 5000-iteration cap.

Measured on Apple M4 Pro, 24 GB unified memory, macOS 27.2.1, AppleClang 21,
Release/Metal. The current machine's default SDK/linker combination failed even
the CMake C probe; local configuration used
`-DCMAKE_OSX_SYSROOT=$(xcrun --sdk macosx26.5 --show-sdk-path)`.
No machine-specific setting was committed. Exact executable/input hashes,
source commits, settings and results are in [receipt.json](receipt.json).
GPU/CPU were shared with other lanes; all times below are contended wall times,
**not performance claims**. No CUDA/RTX or Windows/Linux runtime was available.

## What works

- Binary/ASCII STL, multiple files/named solids, SI/mm coordinates, named walls.
  Vertex welding (1e-7 body diagonal), duplicate/degenerate facet removal and
  winding repair are counted. Open/nonmanifold edges are diagnosed, not hidden.
- Three-ray majority classification of connected-shell unions, without flood
  fill. A small-leak regression stays inside; normals supplied by STL are ignored.
- Tunnel, z=0 contact, y symmetry, refinement boxes, sharp-edge and separate-shell
  proximity sizing, and enforced face/edge/corner 2:1 octree balance.
- Shared piecewise-linear distance-field clipping, convex local agglomeration,
  and bounded 1% isosurface snapping. Uncut neighbours retain planar polygons.
- Closure, directed manifold edges, positive volume, owner orientation, convexity
  and an orthogonal-quality validity floor of 1e-8 are enforced before output.
- Fluent ASCII writer round-trips through the existing reader. `--case` shares
  `fsae_starter`; wheel syntax, moving ground and exact body roles are supported.
- Strict path-qualified mesh JSON, JSONL stages/errors/done, input SHA256,
  process memory, timing, quality statistics, y+ estimates, overwrite protection.
- Reproducible rounded 35° Ahmed geometry reuses the existing surface generator.
  Synthetic FSAE includes three-element wings, endplates, floor/diffuser,
  monocoque, nose, sidepods and four wheels flattened 3 mm onto the ground.

## Measured mesh results

| | Ahmed coarse | FSAE coarse | FSAE default diagnostic |
|---|---:|---:|---:|
| Cells | 24,353 | 208,815 | 3,319,824 |
| Faces | 91,365 | 1,013,514 | 19,946,816 |
| Nodes | 33,449 | 347,873 | 7,169,281 |
| Geometric hexes | 20,058 | 81,745 | 1,224,945 |
| Negative / nonconvex cells | 0 / 0 | 0 / 0 | 0 / 0 |
| Missing body zones | 0 | 0 | 0 |
| Minimum orthogonal quality | 0.00876 | 0.00379 | 0.00292 |
| Maximum face-centre skew / centre distance | 38.28 | 44.75 | 88.39 |
| Maximum adjacent volume ratio | 125,394 | 296,296 | 728,853 |
| Mesh seconds | 0.70 | 12.07 | 427.25 |
| Mesh peak footprint | 17.7 MB | 176 MB | 3.80 GB |

The default diagnostic predates the final distance-query optimization and
additional report fields; its geometry algorithm is the same. It is **not** a
5–15M-cell or GPU-fit qualification. Only 37% of default FSAE cells are geometric
hexes; a hex-dominant result by cell count is not established. All are stored as
Fluent polyhedra, including geometric cubes with subdivided planar faces.

## Solve results and failure boundary

The recorded coarse qualification fails on both cases with
`solver.cpp:117: solve: nonfinite input vector`. Ahmed's last completed monitor
is iteration 160; FSAE's is iteration 5. The solve-process wall times were 3.6 s
and 10.7 s, respectively, including setup. These are failed runs, not solve times
to a converged result. Cd, CL, balance and settled-force timings are **missing**.

For context, the existing body-fitted Ahmed report records Cd 0.3137 on 1.33M
cells against experiment 0.26 (itself a failed drag gate). The present mesher has
no converged Cd to compare to either. Domain, wall resolution and initialization
also differ; this is not a same-mesh solver comparison.

The meshes are topologically valid but contain very poor cut-cell shapes and
extreme adjacent volume ratios. Snapping/agglomeration removed the first
FSAE version's ~1e-17 m³ slivers and ~5.9e9 volume ratio, but did not recover
convergence. SIMPLEC, standard-pressure and reduced-Courant controls also failed.
With Courant 1, time-scale factor 0.1 and 0.2 flow relaxation, FSAE's omega
residual reached 85 at monitor 20 before failure. These controls do **not** prove
a solver defect. The relative contributions of cell geometry, reconstruction
and SST stiffness remain unresolved; no solver numerics or gates were changed.
A synchronizing field diagnostic also found a distinct startup failure:
`initialize_hybrid` → potential gradient → `fv_hybrid_velocity` produces
|w|=321,588 m/s **before iteration 1**, at (-0.993678,-0.661715,0.312276) m,
volume 1.3993e-10 m³ and orthogonal quality 0.6603. This is not a failed SST
iteration or simply the minimum-orthogonal cell. Standard freestream
initialization avoids that starting field but still diverges after monitor 20.
The receipt preserves both observations; neither is a converged solution.

Observed device allocations at phase marks were 65.4 MB / 2,685 bytes per cell
for Ahmed and 694 MB / 3,324 bytes per cell for FSAE (Metal, Coupled cached).
Observed process peaks were 198 MB and 912 MB. Phase marks are lower bounds on
true peak memory. Extrapolating the FSAE observation gives ~16.6 GB at 5M cells,
before host/setup headroom; this is an inference, not a CUDA measurement or a
12 GB fit claim. No settled FSAE solve on the M4 Pro is established either.

## Verification and remaining work

Build and portability lint pass locally. Full quick tier passes with fixtures:
323 cases / 7,996,603 assertions on the final source; the focused selection
passes 7 cases / 93 assertions (including one existing AMG case). Focused mesher tests
cover parsing/repair/leaks, inside queries, refinement and balance, convex cuts,
writer round-trip, JSON errors/overwrite protection and CLI mesh+solve.

1. Condition cut cells beyond a validity floor: remove extreme skew/volume
   ratios while preserving convexity, zones and geometric error. Local convex
   merging alone is insufficient. Isolate the first unstable stencil/field
   before attributing a solver-side cause.
2. Resolve close surfaces in the same connected shell, validate slot/thin-sheet
   preservation quantitatively, and qualify arbitrary defective CAD exports.
   Current proximity refinement uses separate connected shells; larger cracks
   remain unrepaired and sub-grid features can disappear.
3. Add real boundary layers with collision/contact collapse. Nonzero layer
   requests currently fail explicitly; no prism-layer mesh is claimed.
4. Repeat Ahmed/FSAE convergence and force acceptance, refinement studies and
   CUDA memory measurements. Default cell-count/hex dominance and 12 GB targets
   are unmet. Run the portability CI matrix; lint is not MSVC/GCC build proof.

Large generated meshes/STLs/runs are removed after recording the small receipt.
