# Car meshes: a converged default and a convergence stop you can trust

**Every steady car case tested now converges under the defaults and stops
itself.**
- The starter template (`kutta new-case`) uses SIMPLEC.
- A tet/prism mesh is converted to polyhedra first (`kutta polyhedra`).
- Runs with force reports stop when their forces have settled.
- Before this change, no 3D case converged under the defaults, and every run
  went to its iteration limit.

The tet/prism Ahmed mesh still gives forces about 32% above the hex control at
level 1, both before and after conversion. That gap belongs to the native mesh,
not to the conversion or the stop (see "Not solved").

## What changed

1. **Forces-settled stop** (`src/kutta/setup/convergence.*`, `settings.cpp`,
   `app/solve.cpp`).
   - Checking is on by default for any case, Kutta or Fluent, that has a
     force, drag, lift or downforce report. Cases without one are unchanged.
   - Conditions: the case's scaled residual criteria (Fluent's) **and** every
     force coefficient's range over the last 300 iterations ≤ 1e-3. The range
     is relative to the largest |coefficient|, so a near-zero side force
     cannot block the stop.
   - A settled sample is confirmed one iteration later, so a cycle whose
     period divides the monitor interval cannot fake a stop.
   - Residuals can stall above their criteria while the forces hold, as k, ω
     and continuity do on external aero. If that lasts twice the window (600
     iterations), the residual condition is waived and recorded as
     `residuals_waived: true`.
   - A run is never called converged while its forces drift or oscillate.
   - The contract is v1 additive (`docs/app-protocol.md`, "Convergence
     criteria"):
     - setup `convergence.forces`;
     - an iteration condition `{"kind":"forces",…}`;
     - `residuals_waived` in `done` and `summary.json`;
     - inventory entries `monitors/forces-settled/{active,tolerance,window}`.
2. **`kutta polyhedra MESH --out POLY.msh`** (`src/kutta/mesh/polyhedral.*`,
   `io/fluent_msh_write.cpp`).
   - Fluent's tet-to-polyhedra conversion: the median dual of the tet region.
   - Prisms, hexes and pyramids are kept. Shared triangles are split into
     vertex quads, so there are no hanging nodes.
   - Volume, patches and the surface are preserved exactly.
   - At sharp convex edges a dual cell's centroid would fall behind a face.
     Those vertices keep hex quarters instead: 7, 11 and 4 vertices on Ahmed
     L1, L2 and L3.
   - Ahmed L1 goes from 130,331 to 86,215 cells (55,371 tets → 11,149
     polyhedra) in 0.85 s. L2 goes from 217,977 to 133,976 cells, and L3 from
     409,113 to 232,438.
3. **Starter template** (`fsae_starter`).
   - It writes `pressure_velocity_coupling: simplec`.
   - On a 3D polyhedral mesh it also writes
     `warped_face_gradient_correction: false`, so no `--approximate` is needed.
   - On a mesh with tets, `new-case` notes that it should be converted first.
   - A case that states its own scheme, Kutta or Fluent, keeps it.
4. **Octree test meshes** (`make_octree_box`, `io/octree.cpp`). These are
   2:1-balanced meshes whose transition cells are polyhedra with hanging face
   nodes, the kind of transition the mesher lane will produce.

Tests: `tests/mesh/test_polyhedral.cpp` (conversion, `.msh` round trip,
octree balance, uniform flow on CPU and Metal, including the committed 11-tet
mesh `tests/data/meshes/native_farfield_tets.msh`), 4 new cases in
`tests/setup/test_convergence.cpp`, and the app protocol checks.
Gates:
- Quick tier: 324 of 324 cases, 7,996,896 assertions.
- Slow cases 1–4 (the Kutta-case CLI tests) against main `0747475`: 0
  regressions and 0 changed values.
- The only identity changes are the dropped "checking is off" log line and a
  test-line shift.
- The native-physics equivalence test pins Fluent's Coupled default, keeping
  its coverage.

## Why conversion and SIMPLEC, not a tet numerics fix

- **Raw tets** don't converge under either scheme:
  - Coupled: drag range 1.1% and lift 4.6% over the last 300 iterations.
  - SIMPLEC with bench relaxation: still unsettled at 5000 iterations
    (continuity 8e-3, Cl 0.55–0.61).
  - The 11-tet reproducer shows the tet discretization itself fails under
    default numerics, and its dual does not. The earlier limiter-path tricks
    were not repeated.
- **Coupled** limit-cycles on every 3D external case: hex lift 2.3%, poly
  drag 0.2%, still creeping at 3000. This is the Barth-Jespersen cycle in
  `validation/results/coupled_plateau/REPORT.md`, still unfixed.
- **SIMPLEC with Fluent's default relaxation** settles all of them.
- **Pressure relaxation 1.0** instead of 0.3 settles poly L1 in 681 iterations
  rather than 861. That is not worth leaving Fluent's default for.

## Results under the new default (Metal, M4 Pro)

| case | cells | iterations | solve | stop | forces |
| --- | --- | --- | --- | --- | --- |
| Ahmed hex control L1 | 446,364 | 921 | 96 s | residuals + forces | Cd 0.3545, Cl 0.2199 |
| Ahmed native tet/prism L1 → poly | 86,215 | 861 | 27 s | forces (residuals waived) | Cd 0.4692, Cl 0.5251 |
| Ahmed native L2 → poly | 133,976 | 661 | 31 s | residuals + forces | Cd 0.4388, Cl 0.4293 |
| Ahmed native L3 → poly | 232,438 | 961 | 68 s | forces (residuals waived) | Cd 0.3943, Cl 0.3115 |
| NACA 4412 in ground effect, SST + transition (the case states Coupled) | 53,881 | 2601 | 62 s | forces (residuals waived) | Cd 0.0232, downforce 1.2136 |
| octree backward step, 2:1 transitions | 9,404 | 621 | 10 s | residuals + forces | floor −4.0793, step 0.8292 |
| same step, uniform fine control | 69,632 | 621 | 16 s | residuals + forces | floor −4.1389, step 0.8368 |

- **Hex control:** the converged answer equals the long-run hex reference
  (Cd 0.3545, `validation/results/ahmed_drag`).
- **Octree step:** forces through the transitions are within 1.4% (floor) and
  0.9% (step) of a mesh with 7.4 times as many cells.
- **Octree wall-mounted cube** (9,048 and 94,392 cells): the wake is unsteady
  under both schemes and both meshes, so no steady run settles. Only 4.6% of
  the per-iteration velocity change was in the transition cells. It stays in
  the tests as uniform flow only.

## Time to a converged answer (paired, alternating, Metal)

- Two reps alternated the base `0747475` build with this branch.
- The GPU was shared with other agent lanes throughout.
- Results were bit-identical between reps; only wall time varied (rep 1 /
  rep 2).
- "Before" is the old default: Coupled with checking off, so every
  before-run ends at the 3000 cap.

| case | before | after | mean |
| --- | --- | --- | --- |
| Ahmed tet/prism L1 | raw tets: 80 / 99 s, unsettled (drag 1.1%, lift 4.6%) | convert 0.9 s + 33 / 33 s, converged at 861 | no answer → 34 s |
| same mesh as poly, old template | `--approximate wfgc`: 85 / 152 s, unsettled (drag 0.21%, lift 0.31%) | as above | no answer → 34 s |
| Ahmed hex L1 | 391 / 437 s, unsettled (lift 2.3%) | 110 / 105 s, converged at 921 | no answer → 108 s |
| wing in ground effect | 80 / 89 s, drag range 0.11% | 64 / 74 s, converged at 2601 | 84 → 69 s, same forces |
| octree step | `--approximate wfgc`: 26 / 28 s, settled but not stopped | 11 / 11 s, converged at 621 | 27 → 11 s, same forces |

## Not solved

1. **Native tet/prism Ahmed forces are about 32% above the hex control at
   L1.**
   - Native Cd: 0.469, 0.439, 0.394 at L1, L2, L3.
   - Hex Cd: 0.3545, 0.3326, 0.3144.
   - The raw tet run wanders over the same range (Cd 0.43–0.47), so the
     conversion does not cause the gap.
   - The gap is 32% at L1 and L2 and narrows to 25% at L3.
   - I infer, but have not shown, that the cause is the native mesh's
     resolution near the body; the hex L1 mesh has 3.4 times the cells. That
     belongs to the mesher lane.
   - The Ahmed 0.26 drag gate is unchanged and still fails for both families.
   - Follow-up ([tet/poly accuracy](../tet_poly_accuracy/REPORT.md)): the cause
     is the core cell count, not the solver. The hex control is the native
     mesh subdivided, so it has about 20 times the native poly core's cells.
     A polyhedral core of matched cell count comes within 2–3% of hex.
2. **Coupled still limit-cycles on 3D external aero.** The template avoids it.
   - A Fluent case saved with Coupled still runs Coupled, and the forces
     criterion will correctly refuse to call it converged.
   - Whether to switch such a case automatically is an owner call. I
     recommend not: keep the case faithful.
3. **SIMPLEC with intermittency transition oscillates on the wing**
   (downforce 1.4%, drag 8%).
   - The wing case states Coupled, which settles.
   - The template does not enable transition. SST alone converges under
     SIMPLEC in 1001 iterations (20 s).
4. **CUDA was not compiled here** (no toolkit on this Mac). The new code is
   host-only. CI builds CUDA.

## Reproduce

The Ahmed cases are the existing Ahmed 35° bench setup (unit speed and
density, Re_L 2.78e6, area 0.112032) plus the two methods the starter writes.
The tet/prism mesh is the native mesher's L1 output (sha256 `35cee20e…`).
Restricted inputs are not committed. Per-run rows are in
[evidence/runs.csv](evidence/runs.csv).

```sh
kutta polyhedra native_L1.msh --out native_poly.msh
kutta solve ahmed_poly.kutta.json --backend metal --out run   # converges at 861
```
