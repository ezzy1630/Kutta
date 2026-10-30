# Flow solver design (fv/)

What the wave-2 lanes build, and why. Fluent behaviour referenced here is documented in `docs/research/fluent-numerics-and-files.md` (sections A1–A4); unknowns that the saved Fluent fields can settle are in `docs/research/fluent-sst-forensics.md`.

## Goal and the one rule

Same mesh, same settings → Fluent's converged answer within ~1% on forces. A steady converged solution depends only on the **discrete equations**, not on the path the solver takes to reach it. So:

- **Match exactly:** gradients, face reconstruction, limiter, Rhie–Chow face flux, diffusion and non-orthogonal correction, SST source terms and the wall treatment, boundary conditions, and force integration. Each comes from the case file, never from our own defaults.
- **Free to differ:** pressure–velocity coupling, relaxation, pseudo-time and linear solvers. These only change how fast we reach the solution. Pick what is fastest on a GPU.

One exception: the Rhie–Chow coefficient `d_f` depends on `a_P`, which includes pseudo-time and relaxation terms. That makes the converged face flux depend weakly on the path, as an O(h²) pressure-smoothing term. Use Majumdar's relaxation-independent form so that our converged flux matches any consistent Fluent run. Then measure the residual effect on the elbow case.

## Algorithm

1. **Coupling.** Use SIMPLEC with pseudo-transient relaxation (`ρV/Δτ` added to `a_P`, Fluent's automatic conservative Δτ), segregated. It is simpler and needs less memory than a block-coupled solve. All per-cell work is embarrassingly parallel, and the pressure solve reuses `AmgPcg`. Fluent's block-coupled default is also implemented (`src/kutta/fv/coupled_solver.hpp`). It solves the same residual, so both share the fixed point; measurements are in `validation/results/coupled/REPORT.md`.
2. **One outer iteration:**
   - gradients (LSQ) for u, v, w, p, k, ω, with the limiter;
   - assemble and relax momentum, then a few Jacobi or BiCGStab sweeps;
   - Rhie–Chow face flux;
   - pressure correction (AMG-PCG to rel 1e-2…1e-3);
   - correct flux, velocity and pressure;
   - k and ω assembly and solve (positivity-preserving: implicit sinks, clip at floors);
   - μt update;
   - residuals and monitors (forces) every N iterations.
3. **Convergence.** Judge convergence by force plateaus and scaled residuals. Fluent users judge by CL/CD history, and so do we.

## Discretization (Fluent defaults, each switchable from `FluentCase`)

| Item | Implementation |
| --- | --- |
| Gradient | Least squares, cell-based. The 3×3 normal matrix uses inverse-distance² weights, whatever reproduces Fluent's "Gram–Schmidt weights per face". The inverse is stored per cell (6 floats). Boundary faces contribute their face values. |
| Limiter | Barth–Jespersen ("standard"), applied cell-to-face, on the reconstruction gradient only |
| Convection | Upwind, first or second order. SOU: φ_f = φ_up + ψ∇φ_up·r_up,f. Deferred correction: the implicit part is first-order upwind and the explicit part is the high-order correction |
| Pressure at faces | "Second order": p_f = ½(p0+p1) + ½(∇p0·r0 + ∇p1·r1). On boundaries, extrapolate with the gradient |
| Mass flux | Rhie–Chow momentum-based: J_f = ρ_f ū_f·S + d_f [(p0+∇p0·r0) − (p1+∇p1·r1)], with d_f from face-averaged V/a_P |
| Diffusion | Over-relaxed non-orthogonal correction. The implicit part is Γ |S|²/(S·d) (φ1−φ0); the explicit part is Γ ∇φ_f·(S − d |S|²/(S·d)) |
| Warped-face gradient correction | Undocumented in Fluent. Measure the sensitivity on poly meshes before chasing it |
| Turbulence | Fluent SST: research A1, with the corrections established in `fluent-sst-forensics.md` §7. SST-2003m structure with 1994 γ; μt limiter on **strain S** (established); production limiter 10 (not exercised by any data yet, so keep "limiter on ω production" switchable); k has a zero-gradient wall BC; **ω is fixed in wall cells** (not at wall faces) from Eq. 4-423 with β_i = 0.072 and κ = 0.4187; y⁺ uses the centroid-to-wall-face-plane normal distance; blended wall-cell P_k (A2, not yet verified) |
| Wall distance | Exact distance from the cell centroid to the nearest wall surface, computed on the host in double with a BVH, stored as float. This matches `SV_WALL_DIST` (established) and feeds F1/F2. Wall functions use each wall face's own normal distance |

## Device data layout (SoA, all float unless noted)

- **Per face** (internal faces first): `owner`, `neighbour` (int), `S` (3), `d_PN` (3, cell-to-cell), `r_P` and `r_N` (3 each, cell-to-face), `w` (linear interpolation weight). Stored as differences only (see the precision rules).
- **Per cell:** `V`; the cell→face CSR (`cf_offsets`, `cf_faces`, `cf_sign`) for gather assembly; the LSQ inverse (6); wall distance; fields `u, v, w, p, k, ω, μt`; gradients (3 per field).
- **Per boundary face:** the BC kind code and BC values. Patch ranges are contiguous, so the BC kind is uniform per patch.
- **Matrices:** one CSR pattern (cell adjacency) shared by every equation. A face→CSR-entry map gives the assembly target.
- **Cell order:** run RCM renumbering on load. Face arrays follow the owner order.

**Traffic budget:** target ≤ ~2 KB moved per cell per outer iteration, excluding the AMG pressure solve. Measure it per kernel against `bench_bandwidth`. This turns "how fast" into an engineering number before any comparison with Fluent.

## Verification ladder

| Step | Test | Pass |
| --- | --- | --- |
| V1 | Scalar convection–diffusion, manufactured solution on hex, tet and poly-dual meshes | 2nd order (L2 slope ≥ 1.9) diffusion; SOU convection ≥ 1.8 |
| V2 | Laminar NS manufactured solution (steady), plus independent pressure-operator and wall-force refinement | Velocity ≥ 1.9 on both meshes; hex pressure ≥ 1.8; tet pressure ≥ measured momentum pressure-operator order − 0.1 (operator ≥ 0.9). Exact-face Gauss: tet order 0.9–1.1, hex ≥ 1.9. Wall pressure, viscous and total force orders ≥ 1.8; retain failures |
| V3 | Re=1000 cavity at 65²/129²/257² against Erturk 601² tables | Both centreline max errors decrease from 129² to 257², finest-pair order ≥ 1.6; report Richardson primary-vortex values against Erturk/Botella. Re=100/400 Ghia deviations are reported |
| V4 | Fluent cavity mesh (`validation/fluent_examples/lid_driven_cavity`) runs end-to-end | Converges; centreline velocities match a generated mesh with the same grid lines within 2e-6 |
| V5 | SST flat plate, bump and NACA 0012 (`validation/data/nasa_tmr`) | cf and CL/CD within the TMR grid-family band, after the documented compressibility offset |
| V6 | Mixing elbow (`validation/fluent_examples/mixing_elbow`), same mesh as Fluent | Cell fields u, p, k, ω and wall shear compared against `.dat.h5`; report the distributions |

Every step runs on CPU and Metal. CUDA joins once we have hardware.

V2's former second-order tet pressure gate is incompatible with the deliberately
matched Gauss momentum pressure force: midpoint-face quadrature errors leave a
first-order pointwise gradient on Freudenthal tets even with exact face pressure.
Opposite-face cancellation gives second order on hexes. This is an operator
limit, not a reason to change the discretization. The consistency gate is paired
with the independent ≥1.8 integrated-force gate; it does not declare first-order
forces acceptable. See [measured verification](../src/kutta/fv/VERIFICATION.md).

The `slow` doctest suite contains long refinements and cavities; omit it with
`--test-suite-exclude=slow` for CI/quick loops. Running `./build/kutta-tests`
without filters exercises the full ladder implemented so far, including retained
failing accuracy gates.
