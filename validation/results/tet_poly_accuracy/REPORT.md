# Tet/prism and polyhedral Ahmed drag: a mesh-resolution gap, not a solver defect

**Cause: (a), mesh resolution, specifically the volume ("core") cell count. No
solver change.**
- The "hex control" is not an independent mesh. `tools/fsae_aero/ahmed.py`
  builds it by splitting the same native tet/prism mesh: 4 hexes per tet and
  3 per prism.
- After `kutta polyhedra`, the native core has about 1/20 as many cells as
  the hex control's core (11,255 against 221,484 at L1). The wall layers are
  identical in height, layer count and y+.
- With the hex control's exact wall layer and the native polyhedral core, Cd
  stays at the native value. With a polyhedral core of about the hex core's
  cell count, Cd comes within 2–3% of hex.

## Decisive experiments (Ahmed 35°, L1 surface unless noted, Metal)

Same setup as the car-mesh lane: unit speed and density, Re_L 2.78e6, Tu 1%,
viscosity ratio 1, A_ref 0.112032, SIMPLEC. Rows are `bench_fsae_aero` runs
(for `walls.csv`), pressure relaxation 0.3 (1.0 for hex and native L1). Where
both were run they match `kutta solve` within 0.0002. All rows converged
(residuals and forces) except where marked.

| mesh | wall layer | core cells | total cells | Cd | Cl |
| --- | --- | ---: | ---: | ---: | ---: |
| hex control L1 | 3 hexes per prism | 221,484 hex | 446,364 | 0.3545 | 0.2199 |
| hex L1 written as polyhedra (every cell type 7) | same | same | 446,364 | 0.3545 | 0.2199 |
| native L1 → poly (the reported case) | prisms | 11,255 poly | 86,215 | 0.4692 | 0.5252 |
| **hybrid: hex control's wall layer + native poly core** | 3 hexes per prism | 11,255 poly | 236,135 | **0.4616** | **0.5138** |
| native, core at level 3 → poly | prisms | 36,331 poly | 111,291 | 0.4125 | 0.3613 |
| native, core at level 4 → poly | prisms | 82,044 poly | 157,004 | 0.3872 | 0.3047 |
| native, core at level 5 → poly | prisms | 192,275 poly | 267,235 | **0.3663** | **0.2546** |
| hybrid: hex wall layer + level-5 poly core | 3 hexes per prism | 192,275 poly | 417,155 | 0.362* | 0.213* |
| native, surface 0.0198 m, L1 core → poly | prisms (9,652 tris) | 14,949 poly | 207,989 | 0.4246 | 0.4010 |
| native L1 uniformly refined (8 tets/tet, 4 prisms/prism) → poly | prisms (14,992 tris) | 82,436 poly | 382,276 | 0.3824 | 0.2894 |

\* Not settled at 2,500 iterations: drag drifting 0.3639 → 0.3624 and lift
0.2015 → 0.2126 over the last 1,500, residual plateau 3.7e-3. This is the only
unconverged row.

What the rows rule out:
- **Cell representation:** hex-as-poly is bit-identical to hex. The solver is
  face-based; the cell type only sets the warped-face gradient flag, which
  Kutta does not implement either way.
- **Wall treatment, prism layers, y+:** swapping in the hex control's exact
  224,880-cell layer moves Cd only 0.4692 → 0.4616. Every mesh has 20 layers
  from 15 µm (17.28 mm total), with body y+ 0.85–0.90 mean and ≤ 1.7 max.
- **Core resolution:** refining only the core takes Cd 0.469 → 0.413 → 0.387
  → 0.366. At a matched core cell count (192k poly against 221k hex) the
  native mesh is within 3.3% of hex, and within ~2% with the hex wall layer.
  That is smaller than the hex family's own L1 → L2 change (7%), so there is
  no poly-specific excess to attribute to the scheme.

## Where the drag differs (pressure Cd by region, bench `walls.csv`)

| region | hex L1 | native L1 poly | hybrid L1 | level-5 core | refined ×2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| nose (front 100 mm) | −0.0039 | 0.0167 | 0.0179 | 0.0023 | 0.0066 |
| slant | 0.1851 | 0.2594 | 0.2640 | 0.1772 | 0.1976 |
| slant, first 40 mm | 0.0585 | 0.0897 | 0.0865 | 0.0642 | 0.0576 |
| base | 0.1293 | 0.1465 | 0.1339 | 0.1408 | 0.1325 |
| friction (all) | 0.0439 | 0.0466 | 0.0458 | 0.0460 | 0.0456 |

Of the L1 pressure-drag excess (0.112), the slant gives 0.074 (0.031 of it in
the first 40 mm behind the roof edge), the nose 0.021 and the base 0.017. All are
set by the core cells above the 17 mm prism stack. Friction differs by < 0.003.

## Discretization screen on the native L1 poly mesh (bench, one switch each)

Default Cd / Cl 0.4692 / 0.5252.
- Differentiable momentum limiter: 0.5048 / 0.6236.
- Standard pressure interpolation: 0.5313 / 0.6619.
- Normal-foot pressure boundary gradient: 0.4692 / 0.5250.
- First-order k/ω convection: 0.4698 / 0.5105.
- Unlimited, or poor-cell treatment off: nonfinite in iteration 1.

No Fluent-documented alternative moves the poly result toward hex.

## Per level, default `kutta solve` (no code change: before = after)

| level | hex | hex as poly | native tet/prism → poly |
| --- | --- | --- | --- |
| L1 | 0.3545 / 0.2199 (446k) | 0.3545 / 0.2199 | 0.4692 / 0.5252 (86k) |
| L2 | 0.3274 / 0.1590 at the 1000 cap; 0.3292 / 0.1566 at 4000 (758k, unconverged) | 0.3274 / 0.1590 at 1000 (identical) | 0.4388 / 0.4293 (134k) |
| L3 | 0.3183 / 0.0938 at 1000; 0.3144 / 0.1044 at 4000 (1.44M, unconverged) | not run (identical by construction) | 0.3944 / 0.3118 (232k) |
| L1 surface, level-5 core | | | 0.3663 / 0.2546 (267k) |

Values are Cd / Cl (total cells). Every native poly row converged by
residuals and forces; the hex L2 and L3 rows did not (see open problems). The
0.26 experimental gate is unchanged and still fails on every mesh. That gap
is the open SST/resolution question in [Ahmed drag](../ahmed_drag/REPORT.md),
not a mesh-type question.

## Open problems (not addressed here)

- **Default SIMPLEC does not converge on the finer hex Ahmed meshes.** L2 and
  L3 stall at continuity 2.9e-3 and 4.4e-3 (k about 4e-3), and their forces
  never hold within 0.1% over 300 iterations by 4000. The car-mesh lane
  checked hex L1 only, which converges.
- **`--unlimited` and `--poor-cell-quality 0` fail on the native poly mesh**
  with a nonfinite linear-solver input in the first iteration. Both are
  off-default and were not investigated.

## What this means for teams bringing tet/prism meshes

- `kutta polyhedra`, a median-dual conversion like Fluent's, turns about five
  tets into one polyhedron (55,371 → 11,149 at L1). A tet core sized for a
  hex-like cell count is about five times coarser once converted.
- At a matched polyhedral cell count Kutta gives the hex answer. The README
  now says to size the tet region for the polyhedral cell count.
- Nothing here suggests that a poly mesh of normal density from ANSA or
  Fluent Meshing would read 25–32% high.

## Reproduce

Meshes are regenerated, not committed (native L1 sha256 `35cee20e…`, as in
the car-mesh lane).

```sh
uv run --with gmsh==4.15.2 --with numpy python tools/fsae_aero/ahmed.py --slant 35 --first 1.5e-05 --cell-type tet-prism --level 5 \
  --body-surface-size 0.03959797974644666 --out nC5.msh      # L1 surface + level-5 core
kutta polyhedra nC5.msh --out nC5_poly.msh
kutta solve ahmed.kutta.json --backend metal --out run    # mesh: nC5_poly.msh
bench_fsae_aero --mesh nC5_poly.msh --re 2.78e6 --length 1.044 --area 0.112032 --walls body \
  --lift-z --tu 0.01 --urf-p 0.3 --window 300 --iterations 2500 --cells --out run_b  # walls.csv
```

The hybrid meshes came from a scratch script that was not kept. It replaced
each prism in the `kutta polyhedra` output with `subdivide.py`'s three hexes.
The conversion already splits each shell triangle into the same three vertex
quads, so the interface is conformal.
