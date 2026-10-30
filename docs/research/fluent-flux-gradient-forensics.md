# Fluent face mass flux and cell velocity gradient: forensics from saved fields

Date: 2026-10-29. Companion to `fluent-sst-forensics.md` (the "SST report") and `fluent-numerics-and-files.md` (the "main report"). Same data: Fluent's mixing-elbow example, `validation/fluent_examples/mixing_elbow/mixing_elbow.cas.h5` + `.dat.h5`, written by Fluent 22.2 in double precision (`rp-double? #t`).

**Implementation qualification, 2026-10-29:** the side-mean inference below is
incomplete. It fits this saved flux substantially better, but its implementation
loses hex wall-load order and adds tet solution-order failures. Affine exactness
does not resolve the wall-layer curvature defect. The solver retains the
second-order point formulation; see the later interpolation-order section of
[`VERIFICATION.md`](../../src/kutta/fv/VERIFICATION.md). Sections 3 and 5 of
`face_flux.py` now compare the rejected corrections and retained boundary policy.
"Kutta today" in the original tables means the historical pre-round-1 solver.

**Question.** Kutta's SST solver converges about 3.9 % (velocity), 3.1 % (pressure) and 15.6 % (k) away from Fluent on this mesh. At Fluent's saved fields its scaled residuals are about 1000× Fluent's, continuity included. Which discrete formulas differ?

**Tags** (as in the SST report)

- **established**: reproduced to near round-off, or an unambiguous margin over every alternative.
- **likely**: a regression coefficient or structure test points to it, but the reconstruction is not exact.
- **open**: not settled.

**Scripts** (parse-only; usage as in `tools/forensics/cff.py`, e.g. `python -I tools/forensics/face_flux.py <case_dir> <stem> <cache_dir>`)

- `tools/forensics/fvops.py`: sparse gradient, face-reconstruction and limiter operators.
- `tools/forensics/strain_gradient.py`: §1.
- `tools/forensics/face_flux.py`: §2–§3.

Each runs in about 1 s on an M4 Pro once `cff.py` has cached the geometry.

Every table number below is printed by one of these two scripts. Numbers explicitly marked as scratch probes (rejected variants in §1.2, the checks in §1.3, the wall-row scan in §2.5) came from throwaway probes that were not kept; they can be rebuilt from `fvops.py` in a few lines.

---

## 0. Results

**Bottom line:**

- Neither target was reached exactly.
- Both probes give strong structural answers, and Kutta's current flux formula is clearly not Fluent's.
- At Fluent's saved state, Kutta's flux formula has a per-face error of 8.1e-3·|F|. Its divergence is 2.04e-3 kg/s, which reproduces Kutta's measured initial continuity of 2.07e-3 kg/s.
- The reconstructed Fluent formula cuts the per-face error 14× (5.9e-4) and the divergence 2.7× (7.4e-4 kg/s).
- Fluent's own final continuity residual is 2.5e-6 kg/s, so about 300× of the divergence is still unexplained. Most of it sits in the two wall-adjacent cell layers.

### Face mass flux (goal 2)

1. **Boundary faces: established.** `SV_FLUX` = ρ U_face·S to 4e-15 relative on inlets and the outlet. Walls and the symmetry plane are exactly 0.
2. **Interior velocity interpolant: established.** It is the mean of the two side reconstructions u_c + G_c·r_c, with r_c = x_f − x_c.
   - Errors are 6.2e-4·|F| for this interpolant.
   - The alternatives are 3.5e-3 for "mean value + mean gradient·(x_f − x_mid)", 3.5e-3 for the distance-weighted point-corrected form, and 4.6e-3 for no gradient.
   - Kutta's `fv_momentum_faces` uses the point-corrected form.
3. **Reconstruction gradient: established for its weighting.** It is LSQ cell-based with 1/d² weights. The flux error has a sharp minimum at exactly p = 2 (§2.3).
   - **Likely:** it is Barth-Jespersen *limited*. The limited-minus-unlimited term has regression coefficient 1.03, and it is the largest single improvement.
   - Kutta passes unlimited gradients.
4. **Side weighting by momentum a_P: likely.** The fit is (a_P0 u_0r + a_P1 u_1r)/(a_P0 + a_P1), with coefficient 1.02–1.07 across three a_P proxies. This matches the documented a_P weighting [S4] applied to the *reconstructed* side velocities.
5. **Rhie-Chow pressure term: likely.** It is d_f·[(p0 + G_m0·r0) − (p1 + G_m1·r1)], where G_m is the momentum equation's own pressure gradient: Green-Gauss of the second-order face pressures.
   - Coefficient: **0.998** with G_m, versus 0.70 with the LSQ pressure gradient Kutta uses, and 0.33 with cell-based Green-Gauss.
   - d_f = ρ V̄/ā_P·|A|²/(A·d). a_P is convection + diffusion; adding pseudo-time ρV/Δτ fits no better.
6. **Pressure outlet: established that the flux uses the *cell* velocity, not a reconstruction.**
   - Error with the cell velocity is 1.4e-4·|F|, versus 4.3e-3 with reconstruction, which is what Kutta does.
   - **Likely:** the RC term is d_b[(p_c + G_m·r_b) − p_b], coefficient 1.0–1.24.
7. **Open: the residual.**
   - 3.3e-4·|F| on deep faces with no limiting. This is about 6 % of the reconstruction term, i.e. a few-percent gradient mismatch.
   - Most of the remaining divergence is in wall-adjacent cells and their neighbours (§2.5).

### Cell velocity gradient (goal 1)

8. **Fluent's strain rate S is not computed from plain 1/d² LSQ: established.** With exact saved inputs, it misses Fluent's S by a median of 5.0 % (4.0 % on boundary-free cells). Not one cell matches to better than 1.5e-5. Every plain LSQ weighting, Green-Gauss (cell or node) and vertex-stencil LSQ is worse.
9. **Structure: likely.**
   - G_S = G_LSQ(1/d²) + ½·C(G_LSQ). C is a per-cell least-squares fit, with rows from the centroid to each face centroid weighted 1/|r_f|², of the face reconstruction jumps J_f = (φ_nb + G_nb·r_nb) − (φ_c + G_c·r_c).
   - Solving per cell for the coefficient that reproduces Fluent's S exactly gives a median of **0.499** (IQR 0.42–0.58). Alternative correction operators scatter 2–4× wider.
   - This matches the undocumented rpvars `recon/cell-lsf? #t` and `recon/wfgc-jump-enhanced? #t` (warped-face gradient correction, WFGC).
   - **Inferred:** this is WFGC, or close to it.
10. **Near walls: likely.** Building the base LSQ without wall-face rows reduces the layer-2 error from 9.3 % to 2.3 %.
11. **Effect on the turbulence model.**
    - μt relative L2 at the saved state drops from **3.41 %** (Kutta today; VERIFICATION.md reports 3.40 %) to **1.03 %**.
    - Limiter on/off agreement rises from 90.3 % to 97.1 %.
    - The median S error falls from 5.0 % to 1.07 %, and to 0.74 % on boundary-free cells.
12. **The derivative gradient and the reconstruction gradient differ: likely.**
    - The flux gives the jump correction a coefficient of 0.03, so it is absent there.
    - The flux prefers plain LSQ 1/d² (plus the limiter).
    - The S data rejects plain LSQ.

---

## 1. Velocity gradient behind the SST strain rate

### 1.1 Method

- In cells where the SST μt limiter is active, S_F = ρ a1 k / (`SV_MU_T`·F2) is Fluent's own strain rate. Here a1 = 0.31 and F2 uses `SV_WALL_DIST`. The SST report §2 establishes this form.
- That gives 6,756 cells.
- "Clean" means at least 3 face hops from any boundary: 2,276 cells.
- F2 is not the cause of the mismatch. Cells with F2 = 1 to 1e-12 show the same errors.

### 1.2 Candidates, |S/S_F − 1|

| Gradient operator | All limited: median | p90 | Clean: median | p90 | μt rel. L2 (all cells) | Limiter on/off agreement |
|---|---|---|---|---|---|---|
| LSQ 1/d⁰ | 1.34e-1 | 2.88e-1 | 7.42e-2 | 3.72e-1 | 4.36 % | 0.805 |
| LSQ 1/d¹ | 8.09e-2 | 1.90e-1 | 6.43e-2 | 2.52e-1 | 3.91 % | 0.861 |
| **LSQ 1/d² (Kutta today; `lscb-wf 2`)** | 4.98e-2 | 1.25e-1 | 4.03e-2 | 1.19e-1 | 3.41 % | 0.903 |
| LSQ 1/d³ | 5.02e-2 | 1.53e-1 | 5.85e-2 | 1.66e-1 | 4.19 % | 0.833 |
| Green-Gauss cell, arithmetic faces | 6.43e-2 | 1.39e-1 | 6.57e-2 | 1.65e-1 | 3.91 % | 0.884 |
| Green-Gauss of second-order faces (LSQ2) | 1.81e-2 | 1.13e-1 | 1.11e-2 | 4.51e-2 | 2.15 % | 0.884 |
| LSQ2 + ½ jump fit | 1.49e-2 | 1.12e-1 | 7.39e-3 | 3.12e-2 | 2.11 % | 0.920 |
| **LSQ2 + ½ jump fit, wall rows dropped from the base LSQ** | **1.07e-2** | **5.78e-2** | **7.39e-3** | **3.12e-2** | **1.03 %** | **0.971** |

Also tested and rejected in scratch probes (all ≥ 4.9 % median):

- Boundary rows: none, normal-projected.
- Vertex-neighbour LSQ, p = 0–3.
- Green-Gauss node-based with simple, 1/d and 1/d² nodal averages.
- Area-, normal-distance- and face-distance-weighted LSQ.
- Iterated Green-Gauss.
- Quadratic two-ring LSQ (≈10 %).
- S from the limited gradient.
- Invariants other than √(2SᵢⱼSᵢⱼ): a free fit returns 0.97·2SᵢⱼSᵢⱼ + 0.09·2WᵢⱼWᵢⱼ.

**On the quadratic fit.** It is 10 % from Fluent and 11–14 % from the best candidate. So the candidate is close to Fluent's specific operator, not merely "more accurate".

### 1.3 Structure test

For each clean cell, find t with S(G_LSQ2 + t·C) = S_F:

| C | t* p10 / p25 / **p50** / p75 / p90 | IQR |
|---|---|---|
| **jump fit, face rows, w = 1/\|r_f\|²** | 0.264 / 0.419 / **0.499** / 0.585 / 0.756 | **0.166** |
| jump fit, w = 1/\|r_f\|⁰ | 0.157 / 0.389 / 0.507 / 0.725 / 1.223 | 0.336 |
| jump fit, w = 1/\|r_f\|⁴ | 0.175 / 0.327 / 0.461 / 0.549 / 0.735 | 0.222 |
| Green-Gauss(second-order faces) − LSQ2 | 0.294 / 0.659 / 0.898 / 1.051 / 1.35 | 0.392 |
| LSQ 1/d³ − LSQ2 | −0.277 / 0.05 / 0.31 / 0.628 / 1.597 | 0.578 |

**Scratch-only checks** (not in the script):

- Per-face factors instead of ½ (distance, normal-distance and volume interpolation weights) are all wider.
- Reusing the centroid-row LSQ weights for the jump fit gives IQR 0.41.
- Computing the jumps from limited gradients helps only marginally (clean median 7.4e-3 → 6.6e-3).
- Iterating G_{k+1} = G_LSQ + ½C(G_k) oscillates. The single pass is best.

**What is established vs inferred**

- **Established:** the jump-correction direction with a constant ½ fits Fluent far better than anything else tried. Its median per-cell coefficient is 0.499.
- **Likely:** this is Fluent's "LSF"/WFGC correction.
- **Open:** what produces the remaining ±17 % spread in t*, worth 0.7 % in S on clean cells. Candidates are the exact jump definition, a geometric weight, and Fluent's own face centroids on warped faces.
  - Face warp on this mesh is small (median 0.7 % of √A), and the error does not correlate with warp (Spearman −0.08).
  - Area-weighted face centroids beat node-average centroids.

### 1.4 Near walls

Median |S/S_F − 1| by wall layer (layer 1 is wall-adjacent and never limited):

| | L2 | L3 | L4 | L5 | L6 | L7 | L8 | deeper |
|---|---|---|---|---|---|---|---|---|
| LSQ 1/d² | 6.4e-2 | 5.9e-2 | 6.0e-2 | 3.1e-2 | 5.0e-2 | 4.1e-2 | 2.1e-2 | 1.1e-2 |
| + ½ jump fit | 9.3e-2 | 1.4e-2 | 1.1e-2 | 9.4e-3 | 8.7e-3 | 4.0e-3 | 2.3e-3 | 1.9e-3 |
| + ½ jump fit, no wall rows in base | **2.3e-2** | 1.4e-2 | 1.1e-2 | 9.4e-3 | 8.7e-3 | 4.0e-3 | 2.3e-3 | 1.9e-3 |

**Reading the table**

- The layer-2 gain comes only through the wall cell's gradient, which enters layer-2 face reconstructions. Wall cells are never limited, so their own S cannot be tested.
- "No wall row" is the best simple variant, but it is a proxy. The wall-face rpvars `recon/cell-lsq/cortype 1` and `recon/bc-minmax-id-new 3` are unexplained.
- The flux test in §2.5 is ambiguous on the same question.

---

## 2. Face mass flux

### 2.1 Data facts

- `SV_FLUX` (all 91,581 faces) is positive c0 → c1. Its sum|div| is 7.7e-8 kg/s, so it is the post-solve, conservative flux.
- Fluent's final continuity residual is 9.47e-4 × 2.614e-3 = 2.48e-6 kg/s. That is the imbalance of its flux formula recomputed at (nearly) the saved fields, so the exact formula should miss `SV_FLUX` by ≈1e-5 of |F| per face.
- `SV_PSEUDO_DT` is all zero. The global pseudo time step is stored as the rpvar `pseudo-auto-time-step 0.04265362362279325`.
- Interior `SV_P` is the arithmetic mean of the two cells (your probe), so it says nothing about the scheme.

### 2.2 Documented form

Fluent 12 TG eq. 18.4-6 [S4] and the main report's 23-46..49:

J_f = ρ_f (a_P,c0 v_n,c0 + a_P,c1 v_n,c1)/(a_P,c0 + a_P,c1) + d_f[(p_c0 + ∇p_c0·r_c0) − (p_c1 + ∇p_c1·r_c1)] (+ J_ho)

- d_f is a function of ā_P.
- The high-order term J_ho is undocumented (`pressure/ho-rc-scheme 2` in this case).
- The 2025 R1 beta doc [S5] calls the standard flux "high-order, momentum-weighted interpolation".

### 2.3 Regressions on deep faces

The base is F − ½(ρS·u_0r + ρS·u_1r), using unlimited LSQ 1/d² reconstructions.

- **Quiet faces:** 15,324 faces at least 3 hops from boundaries with both cells unlimited. Base residual 6.2e-4·|F|.
- **Velocity interpolant:** see §0 item 2. The side-reconstruction mean wins by 6×.
- **Reconstruction weight** (after RC), median |err|/|F|:

| Weight | 1/d¹ | 1/d^1.5 | **1/d²** | 1/d^2.5 | 1/d³ |
|---|---|---|---|---|---|
| \|err\|/\|F\| | 5.4e-4 | 4.1e-4 | **3.3e-4** | 4.9e-4 | 6.6e-4 |

- **Terms present** (least-squares coefficient; 1 means full strength):

| Fit (faces) | Coefficients | R² |
|---|---|---|
| limited − unlimited reconstruction, + RC (all 26,748 deep faces) | **limited 1.03**, RC 0.68 (LSQ grad) | 0.65 |
| RC jump with momentum gradient G_m = GG(second-order p_f), + a_P-weighting term | **RC 0.998, a_P 1.07** | 0.61 |
| RC jump with LSQ 1/d² gradient | RC 0.70, a_P 1.04 | 0.55 |
| RC jump with Green-Gauss cell gradient | RC 0.33, a_P 0.95 | 0.30 |
| RC averaged-gradient form, (p0 − p1) + ½(G0 + G1)·d, momentum gradient | RC 0.93, a_P 0.94 | 0.54 |
| + the S-derived jump correction in the reconstruction | **jump-fit 0.03** | 0.61 |
| a_P proxy = convection only / + diffusion / + ρV/Δτ | RC 0.93 / 1.00 / 1.08; a_P 1.02 / 1.07 / 1.06 | 0.60 / 0.61 / 0.61 |

**Notes on the table**

- The a_P-weighting term is ½·(a0 − a1)/(a0 + a1)·(ρS·u_0r − ρS·u_1r). That is the difference between a_P-weighted and arithmetic averaging of the reconstructed side velocities.
- a_P is a proxy built from the saved fluxes: Σ max(F_out, 0) + μ_eff diffusion. It is not Fluent's coefficient. The R² values are capped by the unexplained remainder (§2.5).

### 2.4 Assembled formula vs Kutta, at Fluent's saved fields

Saved boundary fluxes are kept, so the divergence measures the interior formula.

| Interior flux formula | \|err\|/\|F\| median | p90 | sum\|div\| (kg/s) | of which layer 1 / 2 / 3 / deep |
|---|---|---|---|---|---|
| a_P-weighted limited LSQ2 reconstructions + RC(G_m) | 6.2e-4 | 1.5e-2 | 8.79e-4 | 2.9e-4 / 2.6e-4 / 6.5e-5 / 1.7e-5 |
| **same, wall faces dropped from the wall cells' LSQ** | **5.9e-4** | **1.2e-2** | **7.44e-4** | 2.0e-4 / 2.1e-4 / 6.5e-5 / 1.7e-5 |
| same, unlimited velocity gradient | 9.5e-4 | 2.0e-2 | 1.29e-3 | |
| same, RC jump with LSQ pressure gradient | 7.9e-4 | 2.0e-2 | 1.07e-3 | |
| **Kutta kernels today** (`fv_momentum_faces` + `fv_rhie_chow`, same a_P proxy) | **8.1e-3** | **1.3e-1** | **2.04e-3** | 4.4e-4 / 4.3e-4 / 2.5e-4 / 4.6e-5 |
| Fluent's own formula (from its residual) | ~1e-5 | | 2.48e-6 | |

Kutta's 2.04e-3 kg/s reproduces its measured initial continuity at the saved fields: 2.068e-3 kg/s, VERIFICATION.md V6.

### 2.5 What remains open

1. **Wall-adjacent cells.**
   - Wall-cell gradient variants pull in opposite directions (scratch probe scanning the wall-row LSQ weight from 0 to 1).
     - Lateral faces between two wall cells prefer no wall row: 1.3e-3·medF.
     - Faces from a wall cell to layer 2 prefer about half the wall-row weight.
     - Zeroing the wall-cell gradient is worse everywhere.
   - Fluent probably treats the wall-cell velocity gradient or reconstruction in a way not captured here.
   - Layers 1 and 2 hold 55 % of the remaining divergence.
2. **Interior.**
   - 3.3e-4·|F| on quiet deep faces. That is about 6 % of the reconstruction term (median 5.4e-3·|F| there), consistent with the few-percent gradient mismatch S shows.
   - The S-derived correction does not fix it (coefficient 0.03). So either the reconstruction gradient has its own deviation from LSQ2, or the limiter stencil differs.
3. **Limiter details.**
   - My Barth-Jespersen implementation uses the face-neighbour min/max including boundary face values, cell-to-face, scalar per component. That is what `fv_limiter` does.
   - The vector form (minimum ψ over components) is clearly worse.
   - The rest of Fluent's stencil is not pinned. It has `recon/limit/epsilon 0.01`, `recon/lim/method 6` and `recon/bc-minmax-id-new 3`.
4. **d_f.** Its exact a_P (diagonal per component vs trace, relaxation, pseudo-time) cannot be separated with the proxy. The data says "about convection + diffusion, without ρV/Δτ".

### 2.6 Pressure outlet (155 faces, p_b = 0, no backflow)

| Face velocity | RC jump gradient | No RC | RC ×1 | Fitted RC coefficient |
|---|---|---|---|---|
| **cell velocity** | momentum G_m | 1.5e-3 | **1.4e-4** | 1.24 |
| **cell velocity** | LSQ | 1.5e-3 | 2.8e-4 | 1.01 |
| reconstructed (Kutta) | momentum | 3.8e-3 | 4.3e-3 | −0.44 |
| reconstructed (Kutta) | LSQ | 3.8e-3 | 4.3e-3 | −0.36 |

---

## 3. Other reasons the Fluent state is not a fixed point of Kutta's equations

These are measured at Fluent's saved state. Each moves Kutta's equations, not just its convergence path.

1. **Continuity.** All of Kutta's initial continuity defect at the saved state is explained by its flux formula (2.04e-3 reproduced vs 2.07e-3 measured).
   - Kutta has four structural differences: interpolant form, unlimited gradients, LSQ instead of momentum pressure gradient in the RC jump, and reconstructed outlet velocity.
   - Fixing these cuts the defect 2.7× at the saved state, and the per-face flux error 14×.
2. **Eddy viscosity.** μt L2 3.41 % → 1.03 % with the corrected derivative gradient (§1). This feeds the k/ω residuals of 1.2e-2/3.1e-2 that VERIFICATION.md reports.
3. **Gradient roles must be split.**
   - The flux reconstruction wants plain LSQ 1/d² (+ limiter).
   - The S derivative wants the jump-corrected gradient.
   - Kutta uses one LSQ gradient for both.
   - Which gradient Fluent uses for momentum diffusion, k/ω production and cross-diffusion is **open**. The saved data does not separate these.
4. **Not investigated.** The momentum equation's own residual at Fluent's state: SOU limiter, diffusion non-orthogonal correction, wall stress. That is the next probe, and `fvops.py` has the pieces.

---

## 4. Changes Kutta's kernels would need

Anchors are function names in the worktree `t3-265ec8e3`. Line numbers moved during this task; the worktree is still being edited.

1. **`fv_momentum_faces`** (`src/kutta/kernels/fv_flow.kern.h`)
   - Interior faces: replace the point-corrected mean-gradient form with the a_P-weighted mean of limited side reconstructions:
     `u_f = w (u_P + ψ_P∘G_P · r_P) + (1 − w)(u_N + ψ_N∘G_N · r_N)`, `w = a_P/(a_P + a_N)`.
     - ψ is the per-component Barth-Jespersen factor, the same as `fv_limiter`, applied to each velocity component.
     - G is the plain LSQ 1/d² gradient.
   - **Outlet faces** (kind 2/4 with outflow): use the cell velocity u_P for the mass flux, not `u_P + G·r`.
     - The reconstructed value may still be what is advected; that was not tested.
2. **`update_flux`** (`src/kutta/fv/steady_solver.cpp`): compute and pass the limited velocity gradient (or the limiter factors) into `fv_momentum_faces`.
3. **`fv_pressure_faces` / `fv_rhie_chow`**: build the RC jump from the momentum-equation gradient G_m = (1/V) Σ_f p_f A_f, where p_f = ½(pp + pn) is the face pressure `fv_momentum_source` already integrates.
   - The jump is `(p_P + G_m,P·r_P) − (p_N + G_m,N·r_N)`. Kutta currently uses `grad_p`, which is LSQ.
   - Keep `ρ·V̄/ā·|A|²/(A·d)` with the physical (no pseudo-time) diagonal; the data is consistent with that.
   - This needs one extra Green-Gauss gradient pass over `face_pressure`.
   - Do the same at outlet faces: `(p_P + G_m·r_b) − p_b`.
4. **SST strain rate** (`sst.cpp`, the `vector_gradient` call before `f.strain`): use G_S = G_LSQ + ½·C(G_LSQ). C is the per-cell LSQ fit, face-centroid rows weighted 1/|r_f|², of the jumps J_f. The base G_LSQ omits wall-face rows.
   - Cost: one extra face loop and a precomputable 3×3 per cell.
   - Keep the plain LSQ gradient for flux reconstruction.
   - This is **likely, not exact**: expect μt within about 1 % L2 of Fluent, not round-off.
5. **Regression tests these support.**
   - On any Fluent case/data pair, the flux formula's sum|div| at the saved fields should fall well below Kutta's current value. On this case: 2.0e-3 → 7.4e-4 kg/s.
   - μt L2 at the saved fields should fall below 1.1 %.
   - Neither is a round-off gate. Both are ratchets until the open items close.

## 5. Sources

- [S1] PyFluent settings API, `least_squares.weighting_factor`: "0=unweighted, 1=inverse distance weighted, 2=inverse distance squared weighted". https://fluent.docs.pyansys.com/version/dev/api/solver/_autosummary/settings/least_squares.html
- [S2] Fluent 2024 R2 Theory Guide 23.3.3, "Evaluation of Gradients and Derivatives": LSQ cell-based via Gram-Schmidt weights per face, Green-Gauss node-based face average. https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/flu_th_sec_eval_derivatives.html
- [S3] Fluent 2024 R2 User's Guide, convergence and stability: WFGC "improve[s] gradient accuracy for all gradient methods"; it targets high-aspect-ratio cells, non-flat boundary-layer faces and centroids outside the cell. Fast mode costs about 22 % more memory. No algorithm is given. https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_sec_solve_conv_stab.html
- [S4] Fluent 12.0 Theory Guide §18.4.2 (ENEA mirror), eq. 18.4-6: momentum-weighted (a_P) face-flux averaging. d_f is a function of the averaged a_P. https://www.afs.enea.it/project/neptunius/docs/fluent/html/th/node372.htm
- [S5] Fluent 2025 R1 beta features 19.5, "Modified Momentum-Based Rhie-Chow": the standard flux uses "high-order, momentum-weighted interpolation". https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/fluent_beta_doc/flu_beta_solver_rc_mod.html
- [S6] Fluent 2024 R2 Theory Guide 23.3.4, gradient limiters: standard = Barth-Jespersen, cell-to-face default. Formulas are images only. https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/flu_th_sec_grad_limit.html
- Case rpvars used as clues: `recon/lscb-wf 2`, `recon/cell-lsf? #t`, `recon/cell-lsf-weighted? #t`, `recon/wfgc-jump-enhanced? #t`, `recon/cell-lsq/cortype 1`, `recon/pressure-recon-grad-ggcb? #t`, `pressure/ho-rc-scheme 2`, `pbs/flux-index 1`, `pseudo-auto-time-step 0.04265…`.
  - `pressure-recon-grad-ggcb? #t` fits the finding that the RC jump uses a Green-Gauss pressure gradient. It is not identical: plain cell-based Green-Gauss scores 0.33, while Green-Gauss of second-order faces scores 1.0.
