# Fluent SST forensics from saved fields: what Fluent actually computes

Date: 2026-10-29. Companion to `fluent-numerics-and-files.md` (cited as "the main report"). This document settles its §0 item 10 and Open items 1–3 as far as one real solved case allows.

**Tags**

- **established**: reproduced from Fluent's own saved numbers to near round-off, or with an unambiguous margin.
- **likely**: the best explanation, but it rests on an inference or a fit.
- **inconclusive**: the data cannot decide it.

**Data**

- `validation/fluent_examples/mixing_elbow/mixing_elbow.cas.h5` + `.dat.h5` (ansys/example-data, MIT; checksums in `SOURCE.md`).
- Written by **ANSYS_FLUENT 22.2 Build 10212**. Kutta targets 2026 R1; see the version caveat in §8.

**Scripts:** `tools/forensics/*.py`. Run order and usage are in the docstring of `tools/forensics/cff.py`. Every number below is printed by those scripts.

---

## 0. Results

1. **The μt limiter uses strain S, not vorticity Ω, with `kw-sst-strain? #f`.** Established.
   - `SV_MU_T` equals ρk/ω to ≤ 2.2e-16 in 62.1 % of cells and is smaller everywhere else.
   - In the limited cells, the S hypothesis reproduces Fluent's limiter value to a median of 0.1 % in the cells where S and Ω differ by more than 20 %. Ω is off by a factor of about 4 there.
2. **Wall-adjacent cell ω is fixed (a cell value, not a face Dirichlet value) by TG Eq. 4-423..4-426.** Established, max error 2.9e-4 over all 2,436 wall cells. Two constants differ from what the main report assumed:
   - **β_i = 0.072**, the standard k-ω `kw-beta-i`, not SST's β_i,1 = 0.075. Using 0.075 is wrong by 1.6–4 %.
   - **κ = 0.4187** (`sa-kappa`), not `sst-kappa` 0.41.
3. **Wall-function length and velocity scales.** Established, both to ≤ 4e-3.
   - Δy = wall-normal distance from the cell centroid to the wall face's own plane.
   - u* = √((μ|U_t|/Δy + C_μ^½ ρk)/ρ), Eq. 4-418.
   - τ_w follows the 4-413..4-417 blend exactly.
   - `SV_WALL_YPLUS` uses u*. `SV_WALL_YPLUS_UTAU` uses u_τ = √(τ_w/ρ).
4. **`SV_WALL_SHEAR` is a per-face force (τ_w·|A_f|), not a stress.** It points antiparallel to the near-wall tangential velocity. Established. This corrects the main report §B6.
5. **`SV_WALL_DIST` (cell field) = exact distance from the cell centroid to the faceted wall surface.** Established, median relative difference 4e-6.
6. **k at walls is zero-gradient**: wall-face k = cell k, median ratio 1.000000. Established within 0.5 %.
7. **Production limiter, wall-cell production, the S/Ω meaning of `kw-sst-strain? #t`, and the wall law at y⁺ > 8**: inconclusive. This case never activates the production limiter and has y⁺ ∈ [0.3, 8.0]. §7 lists what would settle each.

---

## 1. Case and method

**Case**

| Item | Value |
|---|---|
| Mesh | 17,822 cells (17,093 polyhedra, 729 hex); 91,581 faces; 66,417 nodes |
| Boundary zones | walls `wall-elbow` and `wall-inlet` (2,436 faces); `outlet`; `cold-inlet`, `hot-inlet`; `symmetry-xyplane` |
| Fluid | air, ρ = 1.225, μ = 1.7894e-5 (constant); energy on |
| Solver | Coupled (`flow/scheme 24`), pseudo-transient automatic Δτ |
| Discretisation | SOU for momentum, k and ω (`mom/scheme 1`, `k/scheme 1`, `omega/scheme 1`); Second Order pressure (12) |
| Gradients | LSQ cell-based, WFGC on (`recon/cell-lsf? #t`) |
| Convergence | 294 iterations; stopped on the default 1e-3 criteria (§5) |
| Wall-adjacent y⁺ | 0.30 … 8.00, median 4.74 |

**Geometry (own computation)**

- **Faces.** Each polygon is fanned into triangles about its node average. Area vector = Σ triangle areas. Centroid = area-weighted triangle centroids.
- **Cells.** Pyramids from each face to the mean of the cell's face centroids. Closure max |Σ A_out| / V^(2/3) = 1.4e-15.
- **Orientation.** The right-hand-rule normal points into c0 for 100 % of faces, matching the main report. Mass imbalance per cell is ≤ 5.6e-7 of the max face flux.

**Gradients (own computation)**

- LSQ cell-based with neighbour weights 1/dᵖ, p ∈ {0, 1, 1.5, 2, 2.5, 3, 4}, plus cell-based Green-Gauss.
- Boundary faces act as neighbours at the face centroid, carrying the saved boundary-face `SV_U/V/W`. Walls are 0, verified.
- S = √(2 S_ij S_ij), Ω = √(2 W_ij W_ij).

**Wall distance (own computation).** Three candidates (§4).

**Method principle.** Test algebraic identities between saved fields, so solver convergence and transport discretisation drop out.

- In cells where Fluent's μt limiter is active, ρ·a1·k / `SV_MU_T` recovers Fluent's internal X·F2 exactly.
- Only the velocity gradient (needed for S and Ω) is not saved. That gradient is the sole source of scatter in §2.

---

## 2. Eddy viscosity and the S-vs-Ω limiter (main report Open item 1)

### 2.1 Gradient-free facts — established

- R = `SV_MU_T`·ω/(ρk) has min 0.3020 and max 1 + 2.2e-16.
- In 11,066 cells (62.1 %), |R − 1| ≤ 2.2e-16: μt = ρk/ω exactly. So α* = 1, and μt is computed from the saved (final) k and ω, not lagged.
- In the other 6,756 cells (37.9 %), the a1 limiter is active.
- The limiter is active in no wall-adjacent cell (0/2,436). It is active in 1,840/2,436 second-layer cells, 1,136 third-layer, 795 fourth-layer, 506 fifth-layer, and 2,479/5,870 cells deeper in.
- Max μt/μ = 20.7, so `turb-viscosity/max-limit 1e5` is inactive.

### 2.2 Hypothesis test

Error metric: |μt,pred / `SV_MU_T` − 1|, where μt,pred = ρ a1 k / max(a1 ω, X F2) and y = `SV_WALL_DIST`. Using my own exact wall distance gives identical numbers.

"Agreement" is the fraction of cells where the hypothesis turns the limiter on or off exactly as Fluent did. Errors in the inactive cells come only from the hypothesis wrongly switching the limiter on.

| Gradient | Hypothesis | Agreement | Active cells (n = 6,756): median, p95 | Inactive cells (n = 11,066): median, p95 (p99) |
|---|---|---|---|---|
| none | **no limiter** (μt = ρk/ω) | — | 0.171, 1.588 | 0, 0 (0) |
| LSQ unweighted | **S** | 0.805 | 0.075, 0.258 | 0, 2e-16 (9.1e-3) |
| LSQ unweighted | Ω | 0.668 | 0.145, 1.588 | 0, 2e-16 (8.1e-2) |
| LSQ 1/d | **S** | 0.861 | 0.063, 0.217 | 0, 2e-16 (6.0e-2) |
| LSQ 1/d | Ω | 0.723 | 0.121, 1.588 | 0, 2e-16 (9.2e-2) |
| LSQ 1/d² | **S** | **0.903** | **0.046, 0.143** | 0, 0.039 (0.126) |
| LSQ 1/d² | Ω | 0.752 | 0.091, 1.588 | 0, 0.069 (0.145) |
| Green-Gauss | **S** | 0.884 | 0.060, 0.161 | 0, 0.078 (0.181) |
| Green-Gauss | Ω | 0.737 | 0.100, 1.586 | 0, 0.107 (0.225) |

**Discriminating subset.** These are limited cells where |S − Ω|/max(S, Ω) > 20 % (2,555 cells, LSQ 1/d²). The ratio r = (X F2)_own / (X F2)_Fluent:

| X | Median | IQR | p5–p95 |
|---|---|---|---|
| **S** | **0.9991** | [0.960, 1.021] | [0.886, 1.092] |
| Ω | 0.241 | [0.091, 0.556] | [0.036, 0.773] |

**Conclusion — established.**

- Fluent's μt limiter uses the strain-rate magnitude S, exactly as TG Eq. 4-118 says. `kw-sst-strain? #f` does **not** select vorticity.
- Ω's p95 error of 1.59 is gradient-independent. Ω·F2 < a1ω in cells where Fluent cut μt to 0.3–0.6 of ρk/ω, so no gradient error can rescue Ω.
- What `kw-sst-strain? #t` would do is **inconclusive**; it needs a case saved with that flag set.

### 2.3 Residual scatter is our gradient, not the formula — likely

**Weight-exponent scan (limited cells, X = S):**

| Weight | Median \|r − 1\| | p95 | Signed median |
|---|---|---|---|
| 1/d⁰ | 0.134 | 0.380 | −0.129 |
| 1/d¹ | 0.081 | 0.245 | −0.072 |
| 1/d² | 0.050 | 0.152 | −0.031 |
| 1/d^2.5 | 0.046 | 0.159 | −0.007 |
| 1/d³ | 0.050 | 0.197 | +0.001 |
| 1/d⁴ | 0.055 | 0.329 | +0.015 |

**By wall layer (1/d²):**

| Layer | Median r |
|---|---|
| 2 | 0.944 |
| 3 | 0.979 |
| 4 | 0.945 |
| 5 | 0.989 |
| deeper | 0.996 |

**Conclusions**

- **Likely:** the scatter sits in the high-aspect prism layers, where unweighted LSQ is known to be poor. So the remaining 5 % is gradient reconstruction, not the model.
- **Likely:** Fluent's cell LSQ is distance-weighted, roughly inverse-square. This is consistent with the unexplained rpvar `recon/lscb-wf 2`.
- **Inconclusive:** the exact weighting and the effect of WFGC (`recon/cell-lsf? #t`, `recon/cell-lsf-weighted? #t`).

### 2.4 F2 — likely as documented

- F2 = tanh(Φ2²) with y = `SV_WALL_DIST` is consistent with the data.
- Binned by F2, r (1/d²) is 0.913 for [0.5, 0.8) (110 cells), 0.915 for [0.8, 0.95), 0.955 for [0.95, 0.999) and 0.977 for ≥ 0.999.
- The mild trend is within gradient scatter and the sample at low F2 is small. F2 is not independently pinned down.

---

## 3. Near-wall treatment (main report §A2, Open item 3)

All 2,436 wall faces have exactly one wall-adjacent cell, and each such cell has exactly one wall face. Notation:

- Δy_n = (x_P − x_f)·n̂_f, the normal distance from the cell centroid to the wall face's plane.
- Δy_c = |x_P − x_f|.
- U_t = tangential cell velocity relative to the wall (the walls are stationary here).

### 3.1 `SV_WALL_SHEAR` — established

- It is tangential: |n̂·τ|/|τ| ≤ 1.8e-14.
- Its direction is antiparallel to U_t: cos = −1.000000.
- Its magnitude is τ_w·|A_f|. Dividing by face area gives a stress that matches the wall law in §3.3 to a median of 9e-6. Taken as a stress, it would be 5e-5 of μU_t/Δy.
- For body forces, the viscous force on the wall is presumably −Σ_f `SV_WALL_SHEAR` (likely, from the sign). Check it against a `.out` report file when one exists.

### 3.2 y⁺ definitions — established

Ratio of the saved field to ρ·u·d/μ, by distance and velocity scale:

| Saved field | Δy_n, u* (4-418) | Δy_n, u_τ = √(τ_w/ρ) | Δy_n, C_μ^¼ k^½ | `SV_WALL_DIST`, u* | Δy_c, u* |
|---|---|---|---|---|---|
| `SV_WALL_YPLUS` | **median 1.000002, max \|−1\| 1.8e-3** | 1.086 | 2.455 | 1.000125, max 1.6e-2 | 0.99991, max 0.196 |
| `SV_WALL_YPLUS_UTAU` | 0.921 | **median 0.999999, max \|−1\| 1.6e-4** | 2.261 | — | — |

- So **y⁺ = ρ u* Δy_n / μ, with u* = √((μ|U_t|/Δy_n + C_μ^½ ρ k_P)/ρ) and C_μ^½ = 0.3** (TG 4-417/4-418).
- The distance is the wall-normal distance to the cell's own wall face, **not** `SV_WALL_DIST` and not centroid-to-face-centroid.
- The k-based u* = C_μ^¼ k^½ alone is clearly not what Fluent uses.

### 3.3 Wall shear law — established

- τ_w = ρ u* |U_t| / u⁺, with u⁺ = [(y⁺)^-4 + (u⁺_turb)^-4]^(-1/4) and u⁺_turb = ln(E·max(y⁺, 0.2))/κ.
- With the documented κ = 0.4187, E = 9.793: τ_w/pred has median 1.000004, |−1| median 9e-6, p95 3.7e-4, max 3.9e-3.
- A free fit gives κ = 0.4206, E = 9.97. Its rms is no better (2.6e-4), so the documented constants stand.
- In this y⁺ range, the blend raises τ_w about 1.7 % (median) above the laminar μ|U_t|/Δy_n.

### 3.4 ω in wall-adjacent cells — established

Ratio of saved cell `SV_O` to each candidate:

| Candidate | Median ratio | \|ratio − 1\| p95 | Max |
|---|---|---|---|
| **Eq. 4-423, β_i = 0.072, C_calib = 1/3, C_exp = 1.3, κ = 0.4187, β* = 0.09** | **1.000002** | **4.3e-5** | **2.9e-4** |
| same, β_i = 0.075 (SST β_i,1) | 1.0163 | 2.7e-2 | 4.0e-2 |
| same, β_i = 0.0828 (SST β_i,2) | 1.0546 | 9.4e-2 | 0.14 |
| same, y = `SV_WALL_DIST` instead of Δy_n | 0.99966 | 1.1e-2 | 6.0e-2 |
| same, u* = C_μ^¼ k^½ | 1.4875 | 0.57 | 0.83 |
| viscous 6ν/(0.075 Δy²) | 0.700 | 0.52 | 0.64 |
| TMR 60ν/(0.075 Δy²) | 0.070 | 0.95 | 0.96 |
| log u*/(√β* κ Δy) | 1.484 | 1.31 | 10.9 |
| Menter automatic √(ω_vis² + ω_log²) | 0.633 | 0.53 | 0.64 |

**Free fit** with C_calib = 1/3 fixed:

- β_i = 0.072000
- C_exp = 1.29995
- κ√β* = 0.125613 (0.4187 × 0.3 = 0.125610)
- rms log error 2.1e-5

The log branch dominates for y⁺ ≳ 4 (ω⁺_turb/ω⁺_lam = 2.3 at y⁺ = 8), so κ√β* is genuinely pinned, not just the viscous branch.

**Conclusions**

- The saved wall-cell ω equals the formula evaluated with that cell's own final k and U_t. So Fluent **fixes ω in the wall-adjacent cell** at this value every iteration, rather than solving the cell's ω equation against a face value. This resolves main report §A2 [U] and Open item 3.
- The β_i in ω⁺_lam = 6 C_calib/(β_i y⁺²) is the standard k-ω `kw-beta-i = 0.072`. With C_calib = 1/3, ω⁺_lam = 27.78/y⁺², and as y⁺ → 0, ω_P → 27.78 ν/Δy_n². The main report's "2μ/(ρβ_iΔy²)" holds with β_i = 0.072.
- κ is 0.4187 (`sa-kappa` / wall-law κ), not `sst-kappa = 0.41`. Using 0.41 would move κ√β* by 2 %, outside the fit.

### 3.5 Wall-face `SV_O` — form likely, role inconclusive

- Face ω ≈ 2.43 × cell ω (p5–p95 2.31–2.68), falling with y⁺.
- A free fit of ω⁺ = blend(a/y⁺², b/y⁺, c) gives a = 83.26 (6/0.072 = 83.33), b = 15.91 (2/(κ√β*) = 15.92), c = 1.288, rms log error 3.7e-5. In other words, C_calib = 1 and twice the log term.
- With constants fixed at (6/0.072, 2/(κ√β*), 1.3) there is a constant 0.42 % offset (p95 0.43 %).
- Because the cell value is fixed, this face value can only enter through the wall cell's ω gradient (SOU reconstruction into layer 2, CD_kω in the wall cell). Its exact definition is inconclusive.

### 3.6 k at walls — established within 0.5 %

- Wall-face `SV_K`/cell `SV_K`: median 1.000000, |−1| p95 6.0e-4, max 4.6e-3. That is a zero-normal-gradient k at walls, as the TG says.
- The small residual is unexplained. A tangential-gradient reconstruction does not reproduce it.

### 3.7 Not testable here

- **Wall-cell production (4-419..4-422):** not saved, and it cannot be isolated from fields alone. Inconclusive.
- **Whether wall cells are exempt from the μt limiter:** none is limited, but my S·F2/(a1ω) in wall cells ranges up to 0.67 (unweighted) or 1.04 (1/d², 55 cells > 1). The data cannot distinguish "exempt" from "never reaches the threshold". Inconclusive.
- **`turb-wf-ustar-limiter? #t` (coef 0.1):** the u* formula matches without any limiter, so it was inactive or irrelevant here. Inconclusive.

---

## 4. Wall distance (main report Open item 4, partly)

Comparison with the saved cell field `SV_WALL_DIST` (walls only; the symmetry plane is not a wall). Values are |d/`SV_WALL_DIST` − 1|:

| Own definition | All cells: median | p95 | max | Wall-adjacent: median | p95 |
|---|---|---|---|---|---|
| **Exact distance from cell centroid to nearest wall polygon** (fan triangulation about the node average) | **4.0e-6** | **4.5e-5** | 3.0e-2 | 6.7e-6 | 3.6e-5 |
| Normal distance to the plane of that nearest face | 5.3e-4 | 5.4e-3 | 0.79 | 5.6e-4 | 6.5e-3 |
| Distance to nearest wall-face centroid | 1.5e-3 | 5.2e-2 | 1.06 | 7.4e-4 | 2.4e-2 |

- **Established:** `wall-distance/method 2` (with `quad-correction? #t`) gives the true minimum Euclidean distance from the cell centroid to the faceted wall surface. Only 17 cells differ by more than 1e-3 (max 3 %). These are presumably where Fluent's quad-correction or face triangulation differs from my fan.
- **Likely:** F1/F2 use this distance. The μt test cannot separate it from Δy_n, since F2 ≈ 1 near walls. Results with own vs saved distance are identical.
- **Established:** the wall functions (§3) use a different distance, the per-face normal distance Δy_n.

---

## 5. Storage and file semantics (main report Open item 6, and §B6)

- **Cell fields are centroid values — established for this mesh.** My pyramid-decomposition centroids reproduce `SV_WALL_DIST` to a median of 4e-6. They reproduce the y⁺ and ω relations (which go as Δy_n⁻¹ and Δy_n⁻²) to ~1e-5 median.
  - The obvious alternative, the mean of the cell's face centroids, shifts Δy_n by a median of 2.5e-3 (p95 6.5e-3).
  - That would show up as ~5e-3 in ω, so Kutta must use volume-weighted (pyramid) centroids.
- **`SV_MU_T` is consistent with the saved `SV_K`, `SV_O` and velocity** (§2.1), so saved turbulence fields form one consistent state.
- **The case saves `SV_WALL_DIST` (cells) and boundary-face `SV_K`, `SV_O`, `SV_U/V/W`** for every boundary face (ids 1–4,799), not only walls. `SV_WALL_YPLUS_UTAU` also exists (§3.2).
- **Residual arrays `/results/residuals/phase-1/<eq>/data` (n, 4) — established for this file.**
  - Column 0 is the raw residual. Column 1 is the normalisation, so the scaled residual = col0/col1.
  - For continuity, col1 is constant and equals the max of col0 over iterations 1–5.
  - For the other equations, col1 varies per iteration, presumably the global Σ|a_P φ_P| scaling.
  - The run stopped at iteration 294, which is exactly the first iteration where col0/col1 < 1e-3 for continuity, momentum, k and ω (energy is 1.3e-8 < 1e-6). Continuity went from 1.016e-3 to 9.47e-4.
  - Columns 2–3 are all zero here. Their meaning is still inconclusive.

---

## 6. Settings found (main report §A1 option table) — established for this file

Relevant rpvars (full list: `tools/forensics/settings.py`):

| rpvar | Value | Note |
|---|---|---|
| `kw-sst-on?`, `kw-bsl-on?`, `kw-std-on?`, `kw-geko-on?`, `rp-trans-sst?` | #t, #f, #f, #f, #f | plain SST |
| `kw-low-re-mod?` | #f | α* = 1 (confirmed, §2.1) |
| `kw-not-trans-sst-production-limiter?`, `sst-pk-factor` | #t, 10. | production limiter on, C_lim 10 |
| `kw-not-trans-sst-kato-launder?`, `kw-vorticity-pk?` | #f, #f | P_k from S |
| `kw-sst-strain?` | #f | limiter uses S (§2) |
| `curvature-correction?` | #f | |
| `sst-a1`, `sst-beta-i1/i2`, `sst-sig-k1/k2`, `sst-sig-w1/w2`, `sst-kappa` | 0.31, 0.075/0.0828, 1.176/1., 2./1.168, 0.41 | SST constants, stored separately from `kw-*` |
| `kw-beta-i`, `kw-beta-star-inf`, `kw-alpha-inf`, `kw-sig-k/w` | 0.072, 0.09, 0.52, 2./2. | standard k-ω constants; **`kw-beta-i` is used by the SST wall ω (§3.4)** |
| `sa-kappa` | 0.4187 | the κ the wall treatment uses |
| `sst-limit-cross-diff-to-n-times-w-destruct`, `sst-linearize-neg-cross-diff?` | 1e6, #t | undocumented CD_kω clip/linearisation, effect untested |
| `sst-pk-lam-factor`, `sst-po-lam-factor`, `sst-v61?` | 0., 0., #f | undocumented |
| `kw-wall-omega-treatment-r13?`, `kw-set-wall-w?`, `kw-shear-mod?` | #t, #f, #t | 22.2 has no `turbulence/wall-omega-treatment` rpvar |
| `turb-wf-ustar-limiter?`, `turb-wf-ustar-limiter-coef` | #t, 0.1 | undocumented |
| `turb-viscosity/max-limit` | 1e5 | viscosity-ratio cap |
| `recon/cell-lsq?`, `recon/cell-lsf?`, `recon/cell-lsf-weighted?`, `recon/lscb-wf` | #t, #t, #t, 2 | LSQ + WFGC; weight function 2 (§2.3) |
| `wall-distance/method`, `wall-distance/quad-correction?` | 2, #t | §4 |
| `k/scheme`, `omega/scheme`, `mom/scheme`, `pressure/scheme`, `flow/scheme` | 1, 1, 1, 12, 24 | SOU turbulence, coupled |

Inlets have `ke-spec 3` (intensity + hydraulic diameter under the main report's code mapping, I = 5 %). The outlet has `ke-spec 2` with `turb-viscosity-ratio 4`. These matter only for reproducing the case, not for the forensics above.

---

## 7. What Kutta should implement

Constants:

- From the SST model: a1 = 0.31, β* = 0.09, β_i,1/β_i,2 = 0.075/0.0828, σ as in the main report §A1.
- **For the wall treatment only:** β_i,wall = 0.072, κ = 0.4187, E = 9.793, C_μ^½ = 0.3, C_calib = 1/3, C_exp = 1.3.

**Eddy viscosity (every cell, including wall-adjacent):**

```
μ_t = ρ k / max( ω , S F2 / a1 )            S = sqrt(2 S_ij S_ij)  (strain, NOT vorticity)
F2  = tanh(Φ2²),  Φ2 = max( 2√k / (β* ω y) , 500 μ / (ρ y² ω) )
y   = exact min distance from cell centroid to the (faceted) wall surface
μ_t ← min(μ_t, 1e5 μ)                         (turb-viscosity/max-limit; inactive in practice)
```

**Wall-adjacent cell P, with wall face f (per face; one wall face per cell in the reference case):**

```
Δy   = (x_P − x_f) · n̂_f                      wall-normal centroid distance to the face plane
U_t  = (U_P − U_wall) − [(U_P − U_wall)·n̂_f] n̂_f
u*   = sqrt( ( μ |U_t| / Δy + 0.3 ρ k_P ) / ρ )
y⁺   = ρ u* Δy / μ
u⁺   = [ (y⁺)^-4 + ( ln(9.793 · max(y⁺, 0.2)) / 0.4187 )^-4 ]^(-1/4)
τ_w  = ρ u* |U_t| / u⁺                         stress; face force = −τ_w |A_f| Û_t  (= SV_WALL_SHEAR)
ω⁺_lam  = 6·(1/3) / (0.072 · y⁺²)
ω⁺_turb = 1 / (sqrt(0.09) · 0.4187 · y⁺)
ω_P  = (ρ u*² / μ) · ω⁺_lam · [ 1 + (ω⁺_turb/ω⁺_lam)^1.3 ]^(1/1.3)     FIXED cell value, recomputed each iteration
k    : zero normal flux at the wall face
reported y⁺        (SV_WALL_YPLUS)      = ρ u* Δy / μ
reported y⁺(u_τ)   (SV_WALL_YPLUS_UTAU) = ρ sqrt(τ_w/ρ) Δy / μ
```

- Implement "fixed cell value" the usual way: set the wall cell's ω row to identity with RHS ω_P, or use a large-coefficient source.
- When the wall cell has several wall faces, the face rule is **inconclusive**. A sensible choice is the face nearest the centroid; flag it in the code.

**Still per the documentation (unverified by data, keep switchable):**

- Production limiter min(G_k, 10 ρ β* k ω) and its use in G_ω.
- Blended wall-cell P_k (4-419..4-422).
- CD_kω floor 1e-10 and the `sst-limit-cross-diff…` clip.

**Regression tests Kutta can run against any Fluent `.cas.h5`/`.dat.h5` pair** (all rebuild from saved fields, independent of Kutta's convergence):

1. μt = ρk/ω exactly wherever the limiter is inactive.
2. Wall-cell ω, τ_w and both y⁺ fields reproduce to ≤ 4e-3.
3. Own wall distance matches `SV_WALL_DIST` to ≤ 1e-4 median.

`tools/forensics/` already does all three.

**Gradient.**

- To match S where the limiter is active, use distance-weighted LSQ. Inverse-square is the best of the plain weightings tried: median 5 %, p95 15 % error in X·F2, residual bias −3 %.
- Exact Fluent agreement needs its weight function and the WFGC algorithm, both undocumented. Expect this to be the dominant remaining μt mismatch in prism layers.

---

## 8. Still open, and the cheapest data that would settle each

| Open item | Why this case can't settle it | What settles it |
|---|---|---|
| Production limiter form and whether it caps G_ω | Never active: max μt S²/(10ρβ*kω) = 0.36 | Any case with a stagnation region (the FSAE car nose and wings will have one). G_k is not saved, so the test needs a k-equation budget or a saved production field. Ask whether the team can add "Production of k" to the data file quantities. |
| Exact LSQ weights and WFGC | Velocity gradient not saved | Ask for one data file with velocity-gradient (or strain-rate) cell quantities added to the saved data; then S can be compared directly. |
| Wall law at y⁺ 30–80 (FSAE bodywork) and the 2026 R1 `correlation`/`tabulated` choice | This case has y⁺ ≤ 8 and is Fluent 22.2 (`kw-wall-omega-treatment-r13?`; no `turbulence/wall-omega-treatment`) | Rerun `tools/forensics/wall_treatment.py` on the team's 2026 R1 case. If it reports ~1e-5 errors there too, the formula carries over. If not, the residual pattern versus y⁺ shows the `tabulated` shape directly. |
| `kw-sst-strain? #t` semantics | Flag is #f here | A case saved with it #t |
| Wall-cell production, μt-limiter exemption in wall cells, `turb-wf-ustar-limiter` | Not saved / not discriminating | A case with a saved production field, or a Fluent run |
| Meaning of wall-face `SV_O` (§3.5) and residual columns 2–3 | Form fitted, role unknown | Low priority |

**Version caveat.** Every "established" result above holds for Fluent 22.2 on this case. 2026 R1 has a different wall-ω rpvar name, so the wall-ω constants must be re-confirmed on a 2026 R1 file before calling them final. The scripts take the case path as an argument for exactly this purpose.

**Main-report corrections to apply** (not edited here; this task's scope was one new file):

- §A2: wall ω uses β_i = 0.072 and κ = 0.4187, and is imposed as a fixed wall-cell value. Δy is the normal distance to the face plane.
- §A1: `kw-sst-strain? #f` means S.
- §B6:
  - `SV_WALL_SHEAR` is a force per face, antiparallel to the near-wall velocity.
  - Residual col1 is the normalisation.
  - SST data files also carry `SV_WALL_DIST` and `SV_WALL_YPLUS_UTAU`.
