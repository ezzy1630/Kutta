# Fluent numerics and file formats: what Kutta must match and read

Research date: 2026-10-29. Target: Ansys Fluent 2026 R1, steady incompressible RANS k-ω SST, FSAE external aero, same mesh, |ΔCL|, |ΔCD| ≲ 1–2 %.

**Evidence tags**

- **[D]** Documented. A source is given, usually the 2026 R1 Theory Guide (TG) or User's Guide (UG).
- **[O]** Observed directly in a real Fluent file: Ansys' `elbow1.cas.h5`/`.dat.h5` (written by Fluent 24.1 Build 82), or the Fluent 25R2 settings dumps in the cffview repo.
- **[I]** Inferred: my reasoning, or common practice.
- **[U]** Unverified or open.

**Doc access.** The 2026 R1 manuals are public HTML with no login: `https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/<book>/<page>.html`. In the citations below, `TG:` = `…/flu_th/` and `UG:` = `…/flu_ug/`. The equations are SVG images, so I rendered them to read them.

---

## 0. Decision-relevant findings

1. **Fluent's "SST" matches no single NASA TMR variant.** [D]
   - It is closest to **SST-2003m**: strain-based eddy-viscosity limiter, Pk = μt S², production limiter factor 10, CDkω floor 1e-10.
   - Its ω-production coefficients are the **1994 computed** values (α∞,1 = 0.5532, α∞,2 = 0.4403), not 5/9 and 0.44.
   - The wall treatment is not TMR's at all. Fluent uses a blended "y+-insensitive" ω wall value, **zero-flux k at the wall** (not k = 0), and a blended wall-cell production term.
2. **Default convection order for k and ω is first-order upwind.** [D] Many teams switch it to second order. The 25R2 sample case did: `k/scheme 1`, `omega/scheme 1` [O]. We must read this from the case file, never assume it.
3. **Steady single-phase default is the pressure-based Coupled solver** with pseudo-transient Global Time Step (automatic, *conservative* length scale, time scale factor 1), LSQ cell-based gradients, Second Order pressure, SOU momentum, **Rhie-Chow momentum-based** flux (Auto Select), and the Standard (Barth–Jespersen) limiter, cell-to-face. [D]
4. **Warped-Face Gradient Correction (WFGC) is ON by default for 3D meshes containing polyhedra**, which covers poly-hexcore FSAE meshes. Its algorithm is **not documented**. The same goes for the "tabulated" ω-wall table and the wall-distance method. These three are the black boxes most likely to cost the last 1 %. [D/U]
5. **`.cas.h5` + `.dat.h5` hold everything we need**, without changing the team member's workflow. [O]
   - Full face-based mesh (polyhedra included).
   - All settings as Scheme text: models, BCs, discretization, URFs, reference values, report definitions.
   - Cell fields `SV_*`.
   - **Face fields including `SV_P`, `SV_WALL_SHEAR` (n×3), `SV_WALL_YPLUS`, `SV_FLUX`**.
   - Residual history.
6. **Orientation and storage conventions, verified on a real file.** [O]
   - Face node order follows the right-hand rule, and that normal points **into c0**.
   - `SV_FLUX` is **positive out of c0**.
   - Boundary faces come first and have no c1 entry.
   - Vector fields are stored (n, 3) row-major. VTK's CFF reader indexes multi-component data component-major, which looks wrong for (n, 3) data.
7. **Residual history in `.dat.h5` is decimated after 1000 points by default** [D]. Per-iteration force history lives only in the report files (`*.out`). Ask for those.
8. **Fluent has its own native multi-GPU solver** that supports SST, coupled pseudo-transient, LSQ, moving/rotating walls, MRF and porous zones. [D] (UG: `flu_ug_sec_gpu_solver_supported_features.html`) This is the incumbent "fast" baseline Kutta will be compared against, not just CPU Fluent.
9. **Hardware evidence.** Fluent prints a parallel process table (ID / Hostname / Core / O.S. / PID / Vendor = CPU model) and an optional performance timer (average wall-clock per iteration) [D]. Neither is stored in the case or data file, which carry only the version string [O]. So the transcript or batch log is the only hardware evidence.
10. **Cheap way to settle the unknowns.** Saved `SV_MU_T`, `SV_K`, `SV_O`, `SV_U/V/W` and the wall fields let us test several unknowns offline against Fluent's own numbers, without Fluent:
    - which μt limiter is used (the rpvar `kw-sst-strain? #f` exists and is unexplained);
    - the wall-ω and wall-shear formulas;
    - the force integration.

---

## A1. k-ω SST as implemented in Fluent

Sources: [D] TG `flu_th_sec_turb_kw_std.html` (4.4.1), `flu_th_sec_turb_kw_bsl.html` (4.4.2), `flu_th_sec_turb_kw_sst.html` (4.4.3), `flu_th_sec_prod_limit_two_eq_mod.html` (4.21), `flu_th_sec_curv_corr.html` (4.19); UG `flu_ug_sec_turb_setup_options.html` (16.16). Fluent cites Menter 1994, AIAA J. 32(8) as its SST reference [432] (TG `flu_th_biblio.html`).

### Transport equations [D]

```
∂(ρk)/∂t + ∂(ρ k u_i)/∂x_i = ∂/∂x_j(Γ_k ∂k/∂x_j) + G_k − Y_k            (4-71)
∂(ρω)/∂t + ∂(ρ ω u_i)/∂x_i = ∂/∂x_j(Γ_ω ∂ω/∂x_j) + G_ω − Y_ω + D_ω      (BSL/SST)
Γ_k = μ + μ_t/σ_k ,  Γ_ω = μ + μ_t/σ_ω
σ_k = 1/(F1/σ_k,1 + (1−F1)/σ_k,2) ,  σ_ω = 1/(F1/σ_ω,1 + (1−F1)/σ_ω,2)      (4-105, 4-106)
μ_t = (ρk/ω) · 1 / max[1/α*, S F2/(a1 ω)]                                   (4-118)
      α* = 1 by default (Low-Re Corrections OFF); then μ_t = ρ a1 k / max(a1 ω, S F2)
S   = sqrt(2 S_ij S_ij)  ("strain rate magnitude")
F2  = tanh(Φ2²),  Φ2 = max[ 2√k/(0.09 ω y), 500 μ/(ρ y² ω) ]               (4-119, 4-120)
F1  = tanh(Φ1⁴),  Φ1 = min[ max(√k/(0.09 ω y), 500 μ/(ρ y² ω)), 4 ρ k/(σ_ω,2 D_ω⁺ y²) ]   (4-107, 4-108)
D_ω⁺ = max[ 2 ρ (1/σ_ω,2)(1/ω) ∂k/∂x_j ∂ω/∂x_j , 1e-10 ]                     (4-109)
D_ω  = 2 (1−F1) ρ (1/(ω σ_ω,2)) ∂k/∂x_j ∂ω/∂x_j    (not clipped; can be negative)  (4-117)
G_k = μ_t S²                                                                (4-81, 4-461)
G_k ← min(G_k, C_lim ρ ε),  C_lim = 10, ON by default for all ω-equation models    (4-460)
G_ω = (α α*/ν_t) G_k                                                        (4-110)
α   = α∞ = F1 α∞,1 + (1−F1) α∞,2 ,  α∞,i = β_i,i/β*∞ − κ²/(σ_ω,i √β*∞),  κ = 0.41   (4-111..4-113)
Y_k = ρ β* k ω       (BSL/SST: f_β* = 1)                                    (4-114)
Y_ω = ρ β ω² ,  β_i = F1 β_i,1 + (1−F1) β_i,2   (BSL/SST: f_β = 1)          (4-115, 4-116)
y = "distance to the next surface" (wall distance)
```

**Notes on these equations**

- [I] For ω models I read ε in the limiter as β\* k ω, so the limiter is min(G_k, 10 ρ β\* k ω). This is the standard Menter form, but Fluent writes it generically as `C_lim ρ ε`.
- [I] Because G_ω is built from G_k, the k-limiter probably caps ω production too (as in SST-2003). The ν_t in G_ω is presumably the limited μ_t/ρ, which with α\* = 1 gives G_ω = α ρ S² when unclipped.
- [D] Without the low-Re correction, β\* = β\*∞ = 0.09 and β = β_i for incompressible flow, since the compressibility function F(M_t) = 0.

### Constants (SST) [D] and the derived values I computed [I]

| Symbol | Value | Note |
|---|---|---|
| σ_k,1, σ_k,2 | 1.176, 1.0 | Fluent divides by σ: 1/0.85 and 1/1.0 |
| σ_ω,1, σ_ω,2 | 2.0, 1.168 | 1/0.5 and 1/0.856 |
| a1 | 0.31 | |
| β_i,1, β_i,2 | 0.075, 0.0828 | |
| β\*∞ | 0.09 | rpvar `kw-beta-star-inf 0.09` [O] |
| α\*∞, α0 | 1, 1/9 | only used if Low-Re is on; rpvar `kw-alpha-0 0.11111` [O] |
| R_k, R_ω, R_β | 6, 2.95, 8 | low-Re only [D/O] |
| ζ\*, M_t0 | 1.5, 0.25 | compressibility only |
| κ | 0.41 (in α∞) | the wall law uses 0.4187 (§A2) |
| α∞,1 / α∞,2 | **0.5532 / 0.4403** [I computed] | TMR SST-2003 uses 5/9 = 0.5556 and 0.44 |
| C_lim | 10 | "Production Limiter Clip Factor" [D] |

### Option defaults (SST) [D unless tagged]

| Option | Default | Source / rpvar [O = 25R2 sample] |
|---|---|---|
| Low-Re corrections (α\* damping) | OFF | UG 16.16.11; `kw-low-re-mod? #f` |
| Production limiter | ON, C_lim = 10 | TG 4.21, UG 16.16.6; `kw-not-trans-sst-production-limiter? #t` |
| Kato-Launder (G_k = μt S Ω) | OFF (ON only for Transition SST) | TG 4.21; `kw-not-trans-sst-kato-launder? #f` |
| Curvature correction (Smirnov–Menter f_r1 ∈ [0, 1.25], C_curv = 1) | OFF | TG 4.19; `curvature-correction? #f` |
| Compressibility effects | OFF | UG 16.6.4 (option list) [I default] |
| Near-wall ω treatment | `correlation` | TG 4.18.3 ("current default"); `turbulence/wall-omega-treatment correlation` [O] |
| Vorticity-based production | OFF | `kw-vorticity-pk? #f` [O] |
| `kw-sst-strain?` | `#f` [O] | **[U] Undocumented.** It may switch S↔Ω in the μ_t limiter, which would contradict Eq. 4-118's "S is the strain rate magnitude". Resolve with the μ_t test in §0, item 10. |
| ω floor | `omega/min-limit 1e-20` [O] | |

### Fluent SST vs Menter 2003 and NASA TMR variants

TMR definitions are from `https://tmbwg.github.io/turbmodels/sst.html` (the TMR has moved from turbmodels.larc.nasa.gov). [D]

| Item | **Fluent SST (2026 R1 defaults)** | TMR "SST" (Menter 1994) | TMR "SST-V" | TMR "SST-2003" (Menter, Kuntz & Langtry 2003) | TMR "SST-2003m" |
|---|---|---|---|---|---|
| μ_t limiter quantity | **S** (strain) [D]; see the `kw-sst-strain?` caveat | Ω (vorticity) | Ω | S | S |
| k production P | μ_t S² | exact τ_ij ∂u_i/∂x_j, including the −⅔ρkδ_ij term | μ_t Ω² − ⅔ρk ∂u_k/∂x_k | exact | μ_t S² |
| Production limiter | min(P, 10 ρβ\*kω) [I form], default ON | min(P, **20** β\*ρωk), k-equation only | 20, k only | **10**, both equations | 10, both |
| ω production | (αα\*/ν_t)·G_k [D]; G_k is probably the limited value [I] | (γ/ν_t)·P | (γ/ν_t)·μ_tΩ² | (γ/ν_t)·P̃ (the paper's typo αρS² is corrected on TMR) | same as 2003 |
| γ1, γ2 (= α∞,1, α∞,2) | **computed: 0.5532, 0.4403** [I] | computed | computed | 5/9, 0.44 | 5/9, 0.44 |
| CDkω floor | **1e-10** | 1e-20 | 1e-20 | 1e-10 | 1e-10 |
| F1, F2, a1, β, σ | same (σ stored inverted) | same | same | same | same |
| −⅔ρkδ_ij in the momentum stress | [U] not stated. At walls k → 0, so wall forces are barely affected [I] | included | included | included | dropped |
| Wall ω BC | blended y+-insensitive value, §A2 | 10·6ν/(β1 Δd1²) | same | same | same |
| Wall k BC | **zero flux (Neumann)**, even at y+ ≈ 1 [D] | k = 0 | k = 0 | k = 0 | k = 0 |
| Wall-cell production | blended P_k,wall (4-419..422) [D] | none | none | none | none |
| Low-Re damping | optional, OFF | none | none | none | none |
| Freestream sustaining terms | none documented [U] | none (that is the "SST-sust" variant) | none | none | none |

**Bottom line** [I]: implement "SST-2003m with 1994 γ values" plus Fluent's wall model, and make the γ set, the S/Ω limiter and the limiter-on-ω switch configurable. Then calibrate the switches against saved Fluent fields.

---

## A2. Near-wall treatment for ω models ("y+-insensitive") [D]

Source: TG `section_kbf_dl4_vvb.html` (4.18.3); roughness TG `x1-5720008.164.html` (4.18.5).

### Momentum (wall shear)

```
τ_w = ρ u_τ u*                                                       (4-413)
u_τ = U_tan / [ (u⁺_lam)^-4 + (u⁺_turb)^-4 ]^(-1/4)                  (4-414)
u⁺_lam = y⁺ ,  u⁺_turb = (1/κ) ln(E · max(y⁺, 0.2)),  κ = 0.4187, E = 9.793   (4-415, 4-416)
y⁺ = Δy ρ u*/μ   (Δy = wall-normal distance of the wall-cell centroid)   (4-417)
u* = [ (μ |ΔU_tan/Δy| + C_μ^(1/2) ρ k) / ρ ]^(1/2)                    (4-418)
```

- [I] When y⁺ is small, u⁺ = y⁺ and τ_w reduces to μ U_tan/Δy, the exact laminar limit.

### Turbulence

- k gets a zero-flux wall BC at every y⁺, even when the mesh resolves the viscous sublayer. [D]
- Wall-cell production is blended between the two layers:

```
P_k = μ_t (τ_w/μ · ∂u⁺/∂y⁺)²                                          (4-419)
P_k,lam = (ρk/ω)(τ_w/μ)² ,  P_k,turb = (τ_w²/μ) ∂u⁺_turb/∂y⁺ (analytic, from the log law)   (4-420, 4-421)
P_k,wall = (P_k,lam⁻¹ + P_k,turb⁻¹)⁻¹                                  (4-422)
```

- ω at the wall is prescribed:

```
ω_wall = ρ (u*)²/μ · ω⁺_wall                                          (4-423)
ω⁺_wall = ω⁺_lam [ 1 + (ω⁺_turb/ω⁺_lam)^C_exp ]^(1/C_exp)               (4-424)
ω⁺_lam = 6 C_calib/(β_i (y⁺)²),  ω⁺_turb = (1/√β*∞) du⁺_turb/dy⁺        (4-425, 4-426)
C_calib = 1/3,  C_exp = 1.3   (calibrated on plane Couette flow, Re = 1e6)
```

- [I] With C_calib = 1/3 and y⁺ → 0, ω_wall → 2μ/(ρ β_i Δy²), taken with the wall-cell centroid distance. Compare TMR's 60ν/(β1 Δd1²) at the first node. Whether Fluent imposes this as a wall-face Dirichlet value or as a fixed wall-cell value is not stated [U]. Test it against saved face `SV_O`, if present, or against cell values.

### "tabulated" option [D]

- An alternative look-up table for ω⁺_wall, with separate tables for BSL, SST and GEKO.
- Table 4.6 reports flat-plate drag error vs a y⁺ ≈ 0.001 reference:

| Mesh | correlation | tabulated |
|---|---|---|
| y⁺ ≈ 0.1 | 0.47 % | 0.05 % |
| y⁺ ≈ 22 | **5.86 %** | 0.02 % |

- **The table itself is not published** [U]. If a case uses `tabulated` and has y⁺ ≫ 1, we cannot match it from documentation. Its shape could be fitted from Fluent Couette runs, if we ever have Fluent access.
- [I] FSAE y⁺ of 30–80 (§A5) is exactly where the two options disagree. Read the `turbulence/wall-omega-treatment` rpvar from every case.

### Roughness [D]

- Uses the Cebeci–Bradshaw ΔB(K_s⁺) in three regimes: smooth, transitional, fully rough.
- Default for ω models is the virtual wall shift by 50 % of K_s.
- Inputs per wall: `roughness-height` (default 0) and `roughness-const` (default 0.5) [O].
- FSAE bodywork is normally smooth, with K_s = 0 [I].

### Energy

Kader blending. Not needed for us.

---

## A3. Pressure-based solver: discretization and controls

Sources: TG 23.3–23.6 (`flu_th_SolveDiscretization_spatial.html`, `flu_th_sec_eval_derivatives.html`, `flu_th_sec_grad_limit.html`, `flu_th_sec_solver_pb_mom_eq.html`, `flu_th_sec_solver_pb_cont_eq.html`, `flu_th_sec_uns_solve_pv.html`, `flu_th_sec_time_theory_seg.html`, `flu_th_sec_underrelax_pseudo.html`, `flu_th_sec_pseudo_auto.html`); UG 37.x (`flu_ug_sec_solve_choose_disc.html`, `flu_ug_sec_uns_solve_pvel_1.html`, `flu_ug_use_grad_limit.html`, `flu_ug_sec_solve_pseudo.html`, `flu_ug_sec_solve_monitor.html`, `flu_ug_sec_solve_conv_stab.html`). rpvar names come from cffview `reader.py` and the 25R2 sample `examples/general.scm` [O].

### Default settings

| Setting | Fluent default (steady, single-phase) | Formula / notes | rpvar [O] |
|---|---|---|---|
| P-V coupling | **Coupled** + Pseudo Time **Global Time Step** [D] | Implicit p-gradient in momentum, implicit Rhie-Chow face flux, block coupled AMG. Without pseudo time: Flow Courant Number 200, explicit relax p/mom 0.5 [D]. | `flow/scheme` (20 SIMPLE, 21 SIMPLEC, 22 PISO, 24 Coupled); `pseudo-time-method/coupled-pbns/dt-method` (0 = off) |
| Gradient | **Least Squares Cell-Based** [D] | Gram–Schmidt weights per face [D] | `recon/cell-lsq? #t`, `recon/node-lsq?` |
| Convection, momentum | **Second Order Upwind** [D] | φ_f = φ_up + ∇φ_up·r⃗ (Barth–Jespersen linear reconstruction), limited gradient | `mom/scheme 1` (0 FOU, 1 SOU, 4 QUICK, 6 MUSCL) |
| Convection, k and ω | **First Order Upwind** [D] | The FSAE sample used SOU [O] | `k/scheme`, `omega/scheme` |
| Pressure interpolation | **Second Order** [D] | p_f = ½(p0+p1) + ½(∇p0·r0 + ∇p1·r1) (23-43). Other options: Standard (a_P-weighted), PRESTO! (staggered CV), Linear, Body Force Weighted. | `pressure/scheme` (10 Standard, 11 Linear, 12 Second Order, 13 BFW, 14 PRESTO!) |
| Mass flux | **Auto Select → "Rhie-Chow: momentum based"** for most cases; distance-based only for transient SRS + compressible + adapted mesh [D] | J_f = Ĵ_f + d_f[(p0+∇p0·r0) − (p1+∇p1·r1)]; Ĵ_f = ρ_f (a_P0 v0 + a_P1 v1)/(a_P0+a_P1)·n̂ + J_ho,aP (23-46..49). d_f is a function of the face-averaged a_P. J_ho is an **undocumented** high-order non-orthogonal correction [U]. Majumdar [421] cited for under-relaxation consistency. | `pbs/flux-auto-select? #t`, `pbs/flux-index 1` (0 distance, 1 momentum) |
| Gradient limiter | **Standard** (Barth–Jespersen min/max), **cell-to-face** [D] | Scalar limiter on the reconstruction gradient. Also available: multidimensional (Kim–Caraeni–Makarov), differentiable (modified Venkatakrishnan), limiter filter (off). | `recon/limit? #t`, `recon/lim/type 0`, `recon/lim/dx-to-face? #t`, `recon/lim/filter/on? #f` |
| Warped-Face Gradient Correction | **ON automatically for 3D meshes with polyhedra**; fast mode costs ~22 % more memory [D] | Algorithm **not documented** [U] | `recon/cell-lsf? #t` (cffview maps this to WFGC) |
| High Order Term Relaxation | Unchecked in the 2026 R1 dialog figure [D-fig]. When on: factor 0.25 steady / 0.75 transient; type Standard; Flow Variables Only [D] | φ_HO terms ← relaxed (23-15) | `recon/relax/relax?`, `recon/relax/steady-urf`, `recon/relax-limit-mode` |
| Momentum discretization form | "Correction form" ON by default [D] | Affects robustness and round-off, not the converged answer [I] | — |
| Pseudo time, global | Automatic Δτ, **Length Scale Method = Conservative**, Time Scale Factor 1. Δτ_u = 0.3 L/max(0.1 U_bc, U_domain); Δτ_p = 0.3 L/U_ΔP with U_ΔP = 0.1 √((p_max − p_min)_bc/ρ̄); Δτ_rot = 0.1/ω_max; take the min. L_vol = V^(1/3), L_ext = max(Lx, Ly, Lz); conservative L = min(L_vol, L_ext), aggressive L = max [D] | Discrete form ρ V (φ−φ_old)/Δτ + a_P φ = Σ a_nb φ_nb + b (23-101). UG: factor 3–10 "can often increase the convergence rate". | `pseudo-auto-time-step-scale-factor`, `…/length-scale-method` |
| Pseudo explicit relaxation (Coupled) | p 0.5, momentum 0.5, density 1, body forces 1, k 0.75, ω 0.75, μ_t 1, energy 0.75 [D-fig] | φ = φ_old + αΔφ (23-65) | `pressure-coupled/{pressure,mom}/pseudo-explicit-relax`, `{k,omega,turb-viscosity,…}/pseudo-relax` |
| Segregated URFs (SIMPLE) | p 0.3, mom 0.7, k 0.8, ω 0.8, μ_t 1, density 1, body forces 1 [D-fig] | Implicit: a_P φ/α = Σ… + (1−α)/α a_P φ_old (23-66); for the coupled solver (1−α)/α = 1/CFL (23-67) | `<eq>/relax` |
| Precision | Double or single, read from the case | Data arrays may be float32 [O VTK checks precision] | case-config `rp-double?` |

**Typical FSAE choices** [I, from forum and guide practice, the cffview 25R2 sample and the UG advice]:

- Coupled + pseudo-transient.
- LSQ gradients, Second Order pressure, SOU for all equations.
- WFGC on (poly-hexcore).
- Time-scale factor 1–10.
- Sometimes HOTR, or first-to-second-order switching.
- Run 1000–3000 iterations and judge convergence by the CL/CD plateau, not residuals.

### Residuals Fluent reports (pressure-based) [D]

Source: UG 37.17.1.1 (`flu_ug_sec_solve_monitor.html`).

- **Unscaled residual.** R^φ = Σ_P |Σ_nb a_nb φ_nb + b − a_P φ_P| (37-27).
- **Default display is globally scaled.** R^φ = Σ_P |…| / Σ_P |a_P φ_P| (37-28). For momentum, φ_P in the denominator is replaced by the cell velocity magnitude.
- **Continuity, global scaling.** R^c_iter N / R^c_max over the first 5 iterations (37-32). The unscaled value R^c is the sum of |mass creation| (37-30).
- **Local scaling** (optional): RMS of (residual/a_P) over (φ_max − φ_min) (37-29). The enhanced continuity form (37-31) is the default when local scaling is on.
- Normalization by iteration M (default 5) is off by default.
- Default convergence criteria: 1e-3 (1e-6 for energy). With local scaling: 1e-5 steady.
- Residual history goes into the data file automatically. **The store limit is 1000 points; beyond it every other point is discarded** [D].
- [I] For the coupled pseudo-transient solver, a_P presumably includes the ρV/Δτ term, so Kutta's residuals will only be comparable in trend, not in absolute value. Compare forces and fields, not residuals.

---

## A4. Boundary conditions in FSAE runs and how Fluent applies them

Sources: UG `flu_ug_bcs_sec_bound_cond.html` (8.4), `flu_ug_bcs_sec_cell_zones.html` (8.2.3), `flu_ug_sec_mzrrf.html` (11.3), `flu_ug_ReferenceValues.html` (42.13), `flu_ug_sec_op_press.html` (9.15); TG `flu_th_sec_report_force_moment.html` (25.2.1). Keys are from cffview's 25R2 Thread Variables dump [O].

| BC | Fluent behaviour | Thread-variable keys [O] |
|---|---|---|
| **Velocity inlet** | Default turbulence spec **Intensity + Viscosity Ratio**: I = 5 %, μ_t/μ = 10 [D]. k = 3/2 (u_avg I)² (8-75); ω = ρk/μ · (μ_t/μ)⁻¹ (8-80) [D]. Length-scale option: ω = k^½/(C_μ^? ℓ). The 2026 R1 image renders it as C_μ ℓ; the standard form has C_μ^¼ [U]. Inflow mass flux uses only the normal velocity component. If flow exits through the inlet, convection is upwinded from the interior [D]. | `velocity-spec` (0 mag+dir, 1 components, 2 normal), `vmag`, `ke-spec` (1 I+L, 2 I+VR, 3 I+D_h), `turb-intensity`, `turb-viscosity-ratio`, `turb-length-scale`, `frame-of-reference` |
| **Pressure outlet** | Gauge static pressure imposed (default 0). Everything else is extrapolated from the interior. Backflow uses I = 5 %, VR = 10, direction Normal to Boundary [D]. Backflow pressure is Total or Static [D]; the 25R2 sample has `press-spec-gen 0` = Total Pressure [O; default unconfirmed, U]. Average-pressure (strong averaging) is optional. | `p`, `press-spec-gen`, backflow turbulence keys, prevent-reverse-flow |
| **Moving ground** | Wall, "Moving Wall", Translational speed + direction. **Normal components of wall motion are discarded** [D]. Default "Relative to Adjacent Cell Zone"; identical to absolute when the zone is not moving [D]. | `motion-bc` (0 stationary, 1 moving), `relative?`, `vmag`, `wall-translation-vector`, `velocity-components` |
| **Rotating wheels** | Either a rotating wall (speed, axis origin and direction), correct only for surfaces of revolution [D], or an MRF cell zone (frozen-rotor frame motion) for spoked rims [D, UG 11.3.1]. Typical practice: rotating wall on the tyre, MRF only for the rim region [I]. Tyre contact patches and deformed treads have normal velocity components that Fluent drops [D/I]. | `rotating?`, `omega`, `axis-origin-component`, `axis-direction-component` on walls; `mrf-omega` + reference-frame keys on fluid zones |
| **Symmetry (half car)** | Zero normal velocity and zero normal gradients of all variables, i.e. a slip wall [D]. With zero-shear walls, the turbulence wall treatment is removed by default (Neumann) [D]. | zone type `symmetry` (bc-type 7) |
| **Porous radiator** | S_i = −(μ/α · v_i + C2 · ½ρ\|v\| v_i) per direction (Darcy + inertial; 8-1/8-2), with 1/α and C2 given per direction. Superficial velocity formulation by default [D]. Turbulence is solved in the zone unless "laminar zone" is set [D]. Force reports can include porous zones [D]. | fluid-zone keys (porous?, viscous/inertial resistance vectors) [I names] |
| **Operating pressure** | Default 101325 Pa. **All input and reported pressures are gauge.** Not used for constant density [D]. | `operating-pressure` |

### Forces and coefficients [D]

- **Force along a vector a:** F_a = a·F_p + a·F_v.
- **Force coefficient:** C = F_a / (½ ρ_ref v_ref² A_ref). The moment coefficient also uses L_ref.
- **Pressure force uses the reference pressure.** So F_p = Σ_f (p_f − p_ref) A⃗_f [I exact form].
- **Reference values are post-processing only:** area, density, length, pressure, velocity, viscosity [D]. Typical defaults are A = 1, ρ = 1.225, v = 1, L = 1, p = 0 [I]. Read them from the `reference-values` rpvars [O].
- **Force report definitions** (Drag / Lift / Force / Moment): force vector (default drag x̂, lift ŷ), wall zones, per-zone option. Each creates a report file and plot by default [D].
- [I] For half-car runs, check whether A_ref is the half area or the forces get doubled. This is team-specific; read it from the case and report definitions.
- [I] For an exact replication we can recompute forces from the data file: face `SV_P` and `SV_WALL_SHEAR` (per wall face), with face areas from the mesh. This gives a self-check against their `.out` values.

---

## A5. Published FSAE setups

Fluent-specific data is sparse. Most recent public FSAE theses use STAR-CCM+; I note the solver for each.

| Source | Solver | Car | Cells | Wall / y⁺ | Other | Runtime / HW |
|---|---|---|---|---|---|---|
| Monash Motorsport (Ansys blog, 2024-03-04) — https://ansys.synopsys.com/blog/monash-motorsport-tops-fsae-ranks-with-ansys | Fluent | Full car; cornering by re-orienting the mesh (steer, pitch, roll) | "up to 33 million" | — | k-ε, rotating wheels, ground effect | "around an 11-hour run and post-processing" on university HPC |
| ECU-R guest blog (LEAP Australia) — https://www.leapaust.com.au/blog/cfd/ecu-r_analysing-fsae-flow-fields-using-cfd/ | Fluent Meshing + Ansys CFD | straight-line | ~25 M, poly-hexcore, 4 mm refinement | **8 inflation layers, y⁺ 30–80** including the ground | Cp0 / helicity post-processing | — |
| KTH, Dahlberg 2014 — https://www.diva-portal.org/smash/get/diva2:737287/fulltext01.pdf | Fluent (Workbench / ICEM) | 3D vehicle | ~6 M (limited by time) | y⁺ adaptation in Fluent | 2D validation vs a wing-in-ground experiment | — |
| TU Wien, Kirchberger 2023 — https://repositum.tuwien.at/bitstream/20.500.12708/188917/1/Kirchberger%20Michael%20-%202023%20-%20CFD%20Simulation%20and%20Validation%20of%20a%20Formula%20Student...pdf | STAR-CCM+ | **Half car**, symmetry plane, porous radiator, 15 m/s | 10 → 163 M study; design runs at 34 M (33 mm base) | high-y⁺ (4 layers, y⁺ 30–100) on bodywork, low-y⁺ on aero | 1800 iterations (10000 checked: CDA changed 3.9 %, CLA 0 %); CLA ≈ 4.7–5.2 m²; validated with pressure taps | 18 mm mesh is "4 times" the runtime of 33 mm |
| KTH DeV18 2023 — https://www.diva-portal.org/smash/get/diva2:1778751/FULLTEXT01.pdf | STAR-CCM+ 17.06 | Half car, symmetry walls/roof | — | 10 prism layers, stretching 1.4 | rotating wheels 480 rpm at 40 km/h | — |
| ASME JFE 148(3) 031202, cornering FS car — https://asmedigitalcollection.asme.org/fluidsengineering/article/148/3/031202/1222813/Aerodynamics-of-a-Cornering-Formula-Student-Car | STAR-CCM+ | cornering | — | — | steady RANS SST, tyres as rotating walls | — |

[I] Synthesis for planning:

- Half car, 15–20 m/s, moving ground, rotating-wall tyres, porous radiator.
- 20–40 M poly-hexcore cells, 5–10 prism layers, y⁺ either ~1 on wings or 30–80 on bodywork.
- 1–3 k iterations; about half a day on a university cluster or a 16–64-core workstation.
- The friend's actual setup should be read from his case file, not taken from this table.

---

## B6. Fluent Common Fluids Format (HDF5): `.cas.h5`, `.dat.h5`, `.msh.h5`

Sources:

- [D] UG Appendix B.1 `flu_ug_app_cff.html`. It gives the group tree only: settings are "formatted strings", heavy data are 1D/2D datasets, and light data are attributes.
- [D-code] VTK `IO/FLUENTCFF/vtkFLUENTCFFReader.cxx` and `vtkFLUENTCFFInternal.cxx` (https://github.com/Kitware/VTK/tree/master/IO/FLUENTCFF).
- [D-code] pyfluent `src/ansys/fluent/core/offline/{case_file,data_file,_lispy}.py` (https://github.com/ansys/pyfluent).
- [D-code] cffview (https://github.com/aiwaonly/cffview, BSD-3, "only tested with Fluent 25R2").
- [O] `h5ls` / `h5dump` on Ansys' `elbow1.cas.h5`/`.dat.h5` from https://github.com/ansys/example-data/tree/master/pyfluent/file_session (Fluent 24.1, 3D hex, 13,852 cells, 44,500 faces, 16,968 nodes).

All IDs are **1-based**. Arrays are chunked and deflate-compressed (level 1 observed). Every section dataset carries `minId`/`maxId` (and `chunkDim`) attributes. A field can have several numbered sections `1, 2, …`, so always place data by `minId..maxId`, never by position. [O]

### Case file layout [O, with types]

```
/meshes/1                      attrs: dimension(i32)=3, version(i32)=2,
                                      nodeCount/nodeOffset, cellCount/cellOffset, faceCount/faceOffset (u64)
  /nodes/zoneTopology          nZones; datasets id(i32), minId, maxId(u64), dimension
  /nodes/coords/<zoneId>       f64 (n,3), attrs minId,maxId
  /cells/zoneTopology          id, minId, maxId, cellType(i32), name ("a;b;c" one string), dimension, childZoneId, zoneType, fields
  /cells/ctype/<k>             attrs elementType, minId, maxId, [dataset cell-types(i16) when zone is mixed (cellType 0)]
  /cells/partition/<k>         partition ids (ignore)
  /cells/tree, /faces/tree     only with hanging-node adaption: nkids, kids per parent [D-code]
  /faces/zoneTopology          id, minId, maxId, zoneType(i32), faceType(i32), name ("wall;symmetry;…"),
                               c0, c1 (= adjacent CELL-ZONE id; 0 = none), childZoneId, shadowZoneId, flags
  /faces/nodes/<k>             attrs nSections/minId/maxId; nnodes(i16, per face); nodes(u32, flat 1-based node ids)
  /faces/c0/<k>                u32 owner cell per face, attrs minId=1,maxId=nFaces
  /faces/c1/<k>                u32 neighbour cell; present ONLY for faces with two cells (here minId=5889..maxId=44500)
  /faces/interface             non-conformal interfaces (pf0/pf1) [D-code]
/settings                      one string dataset each (all Scheme text):
  "Rampant Variables"          the rpvar alist "(37 ( (name value) … ))": models, numerics, URFs, reference values,
                               named expressions, report definitions, monitors, case-config
  "Thread Variables"           per-zone BCs: "(39 (zone-id zone-type zone-name domain-id)( (key . value) … ))"
  "Domain Variables", "Cortex Variables" (units, graphics), "TGrid Variables", "Origin", "Solver", "Version"
/special/…                     e.g. boundary_layer_flags (cells/faces), Overset_DCI
```

Observed elbow values [O]:

- `faces/zoneTopology`: id = [3,4,5,6,7,9]; zoneType = [3,7,5,10,10,2]; faceType = 4 (quad); names "wall;symmetry;pressure-outlet-7;velocity-inlet-6;velocity-inlet-5;default-interior"; c0 = 2 for all; c1 = [0,0,0,0,0,2].
- **Boundary faces occupy IDs 1..5888 and the interior zone is last.**
- Cell zone 2: `cellType` 4 (hex), name "fluid".
- The zone, cell and face type codes are **the same integers as the legacy bc-type, element-type and face-type tables** in §B7. So one decoder serves both formats.

**Polyhedra and polygons**

- A cell has no node list. It is defined by the faces whose c0 or c1 points to it. Polyhedral zones have `cellType` 7. Mixed zones have type 0 and per-cell codes in `cells/ctype/<k>/cell-types`. [D-code]
- Polygon faces: `faceType` 5, or 0 for mixed. Per-face sizes are in `nnodes`.
- The UG figure also shows an optional `cells/nodes` group [D]. It is absent in the elbow file.

**Orientation** [O, tested on the elbow file]

- The RH-rule normal (n0→n1→n2) points **into c0** for 61/61 sampled boundary faces and 200/200 interior faces. This matches the legacy spec ("thumb toward c0").
- `SV_FLUX` is **positive leaving c0**:

| Zone | Flux sum |
|---|---|
| pressure-outlet | +1.919 |
| inlets | −0.302 and −1.618 |
| walls | 0 |

- So reverse the node order (or negate the normal) to get an outward-from-c0 area vector consistent with Fluent's flux sign.

**Hanging nodes** [D/I]

- Hexcore or poly-hexcore meshes refined by hex subdivision can have hanging nodes (UG 37.3.2).
- Parent faces then sit in "parent" zones (bc-type 31) and inactive parent cells in type-32 zones, linked by the `tree` groups.
- The solver uses the children. Skip the parents.

### Settings parsing [O/D-code]

- Strings are Scheme s-expressions. pyfluent parses them with its own `_lispy.parse`. cffview uses `sexpdata` plus regexes.
- `Rampant Variables` starts with `(37 (` [O].
- Model flags live in the `case-config` sub-alist as dotted pairs: `(rp-3d? . #t)`, `(rp-kw? . #t)`, `(rp-double? . #t)`, `(rp-seg? . #t)`, `(rp-unsteady? . #f)`.
  - The top-level list can hold same-named entries with different values (elbow: top-level `(rp-ke? #f)` vs case-config `(rp-ke? . #t)`).
  - **Use case-config for model flags**, as pyfluent does for `rp-3d?`.
- Some entries appear twice, e.g. `(flow/scheme 20)` and `(flow/scheme 21 integer #f 21 20)`. Take the two-element top-level entry [O/I].
- BC values look like `(vmag (constant . 15) (profile "" ""))`. Expressions appear as `…-sel "constant"` / `…-expr "1"` pairs [O].
- Useful rpvars:

| Group | rpvars |
|---|---|
| Iteration count | `number-of-iterations` (pyfluent `iter_count`) |
| Scheme selections | `flow/scheme`, `mom/scheme`, `k/scheme`, `omega/scheme`, `pressure/scheme` |
| Gradients and limiter | `recon/cell-lsq?`, `recon/cell-lsf?`, `recon/limit?`, `recon/lim/type` |
| Flux | `pbs/flux-auto-select?`, `pbs/flux-index` |
| Pseudo time and relaxation | `pseudo-time-method/coupled-pbns/dt-method`, `pseudo-auto-time-step-scale-factor`, `pressure-coupled/*/pseudo-explicit-relax`, `*/pseudo-relax` |
| SST | `kw-sst-on?`, `kw-bsl-on?`, `kw-geko-on?`, `kw-low-re-mod?`, `kw-not-trans-sst-production-limiter?`, `kw-not-trans-sst-kato-launder?`, `curvature-correction?`, `turbulence/wall-omega-treatment`, `kw-sst-strain?`, `kw-vorticity-pk?` |
| Wall distance | `wall-distance/method` (=2, undocumented [U]), `wall-distance/quad-correction? #t` |
| Operating and reference | `operating-pressure`, `reference-values` |

- Fluent version: `/settings/Version` ("24.1"), `/settings/Origin` ("ANSYS_FLUENT 24.1 Build 82") [O].

### Data file layout [O]

```
/results/1/phase-1/cells/fields      "SV_P;SV_U;…" list
/results/1/phase-1/cells/SV_<name>/<k>    f64 (or f32) (n) or (n,ncomp), attrs minId,maxId
/results/1/phase-1/faces/fields
/results/1/phase-1/faces/SV_<name>/<k>    face data; sections cover only the faces where the field exists
/results/residuals [attr version=3]/phase-1/<eq>/{iterations (n), data (n,4)}   eq = continuity, x-velocity, …, k, omega
/settings/"Case File", "Data Variables" (flow-time, time-step, …), "Origin", "Solver", "Version"
```

**Cell fields** (elbow, k-ε case) [O]: `SV_P` (gauge static pressure), `SV_U`, `SV_V`, `SV_W`, `SV_K`, `SV_D` (ε), `SV_MU_T`, `SV_MU_LAM`, `SV_DENSITY`, `SV_T`, `SV_H`, `SV_BF_V`.

- For SST expect `SV_O` = ω. VTK's name map lists "O" as "specific dissipation rate" [D-code].
- The VTK map also lists `_G` gradient, `_M1`/`_M2` previous-time-level, `_MEAN`/`_RMS`, `WALL_DIST`, `DUDX…` names that may appear depending on settings.

**Face fields** [O]: `SV_P`, `SV_U/V/W`, `SV_K`, `SV_FLUX` (mass flux, all faces), `SV_DENSITY`, plus wall-only fields:

| Field | Shape | Faces |
|---|---|---|
| `SV_WALL_SHEAR` | (3630, 3), row-major | wall faces 1..3630 |
| `SV_WALL_YPLUS` | (3630) | wall faces |
| `SV_WALL_YPLUS_UTAU` | — | wall faces |
| `SV_WALL_V`, `SV_WALL_VV`, `SV_WALL_T_INNER`, `SV_WALL_PRORUS_ZONE_FORCE_MEAN` [sic] | — | wall faces |

The elbow file also holds boundary-coefficient arrays `SV_F_SC_*`, the heat fluxes and `SV_ARTIFICIAL_WALL_FLAG`.

- **Vector layout check** [O]: `SV_WALL_SHEAR` rows are exactly tangent to their faces (n̂·τ/|τ| = 0.000 on 10/10 faces). The component-major interpretation gives 0.14–0.75. So vectors are (n, 3) row-major.
- VTK's `GetData` indexes multi-component arrays as `data[k*maxId + (j-1)]`, which is component-major. I think that is a VTK bug for (n, 3) data [I].
- **Residual data columns:** 4 per iteration. Column 0 is the raw residual (the continuity scale equals the iteration-1 raw value). Column 1 is a constant scale factor. Columns 2–3 are unknown [U]. Verify against the printed scaled residuals in a transcript before using them.

### Open-source readers

| Reader | Handles | Notes |
|---|---|---|
| VTK `vtkFLUENTCFFReader` (ParaView "Fluent CFF Case Reader") | `.cas.h5` + `.dat.h5`: nodes, faces, cells (polyhedra), cell data; face data opt-in | Requires `/meshes` and `/settings`. Reads face nodes as u32 and c0/c1 by sections. Its face-node pointer starts at `minId` of each section [I: check multi-section files]. Multi-component indexing is suspect (above). |
| pyfluent `CaseFile` / `DataFile` / `FileSession` | rpvars via `_lispy`, mesh surfaces, face fields per surface, `iter_count`, `precision`, input/output parameters | Reads `settings["Rampant Variables"][0]`. `.flprj` projects resolve into a `<name>.cffdb` folder. Comment says residuals are not yet exposed. |
| cffview (PyPI 0.6.1) | Dumps Rampant / Thread / Cortex variables to `.scm`, parses solver, BC, material and report settings to JSON | Best map of rpvar names → GUI settings; BC enum codes are in `boundary.py`. Mesh reading reuses VTK logic. |
| meshio `ansys` | **Legacy** `.msh` (ASCII + partial binary) | No mixed binary cells or faces; no CFF. |
| VTK `vtkFLUENTReader` | Legacy `.cas`/`.dat` | Has the legacy sub-section-id → name table (1 PRESSURE, 2 MOMENTUM, 5 TKE, 19 WALL_SHEAR, 63 SDR = ω, 103 MU_TURB, 111–113 X/Y/Z_VELOCITY, …). |
| OpenFOAM `fluent3DMeshToFoam` | Legacy `.msh` | Not inspected [U]. |

[I] Recommendation: write our own reader (HDF5 C API or h5py) against the layout above, and cross-check it against VTK and pyfluent on Ansys' example files. Neither library is complete or bug-free for our needs (face fields, vectors, polyhedra, multi-section).

---

## B7. Legacy `.msh` / `.cas` / `.dat` format [D]

Source: UG B.2 `flu_ug_app_form_case_data.html`.

**General rules**

- Sections are `(index (header…)(body…))`.
- Grid-section header indices are **hexadecimal**; zone IDs in zone sections are decimal.
- Binary sections prefix the index with 20 (single) or 30 (double), e.g. 2010/3010, and end with `End of Binary Section 3010)`.

| Index | Section | Header / content |
|---|---|---|
| 0 | comment | `(0 "text")` |
| 1 | header | program that wrote the file |
| 2 | dimensions | `(2 3)` |
| 10 | nodes | `(10 (zone-id first last type ND)(x y z …))`; zone-id 0 = declaration (total count) |
| 12 | cells | `(12 (zone-id first last type element-type))`; element-type 0 mixed, 1 tri, 2 tet, 3 quad, 4 hex, 5 pyramid, 6 wedge, **7 polyhedral**; type 1 active, 32 inactive (hanging-node parents) |
| 13 | faces | `(13 (zone-id first last bc-type face-type))`; body per face `n0 n1 [n2 …] c0 c1`; mixed (0) or polygonal (5) faces prefix the node count `x`. **RH-rule normal points toward c0**; c = 0 means no cell. bc-type: 2 interior, 3 wall, 4 pressure-inlet, 5 pressure-outlet, 7 symmetry, 8 periodic-shadow, 9 pressure-far-field, 10 velocity-inlet, 12 periodic, 14 fan / porous-jump / radiator, 20 mass-flow, 24 interface, 31 parent, 36 outflow, 37 axis; +1000 for non-conformal intersection zones. face-type: 0 mixed, 2 linear, 3 tri, 4 quad, 5 polygon. |
| 18 | periodic shadow faces | pairs |
| 58 / 59 | cell tree / face tree | hanging-node hierarchy |
| 61 | interface face parents | non-conformal |
| 39 / 45 | zone | `(39 (zone-id zone-type zone-name domain-id)((cond . value)…))`; 45 = name only |
| 37 / 38 | rpvars / domain variables (case) | same Scheme alists as CFF `Rampant` / `Domain Variables` [O: pyfluent splits legacy `.cas` at `"(37 ("`] |
| 40 | partitions | |
| 33 | grid size (`.dat`) | cells, faces, nodes check |
| 300 (2300/3300) | data field | `(300 (sub-id zone-id size n-time-levels n-phases first last)(…))`; sub-ids are defined in the install's `xfile.h`, VTK copies them |
| 302 | residuals | `(302 (n residual-id size domain-id)(…))` |

---

## B8. Transcript (`.trn`) and report files (`.out`)

### Transcript [D]

UG 3.13 `flu_ug_TranscriptFile.html`.

- Records all console input and output, including session start, end and total time and **the wall-clock time per iteration**.
- Off unless the user starts it, or sets Preferences → General → **Automatic Transcript**. Automatic files are named `fluent-YYYYMMDD-HHMMSS-<pid>.trn` in the working directory (or the temp directory if that is not writable). The auto-transcript default is not stated [U].
- In batch or cluster runs the console goes to the scheduler's stdout file (e.g. `slurm-*.out`, `*.log`). That file is the transcript equivalent [I].

**What to parse** [D unless tagged]

1. **Iteration table** [I format; seen in user logs, e.g. https://www.cfd-online.com/Forums/fluent/204063-nita-fractional-step-residuals-convergence-print.html]:
   - A header line `iter  continuity  x-velocity  y-velocity  z-velocity  k  omega  [report defs…]  time/iter`.
   - Per-iteration rows of scaled residuals in `%.4e`.
   - Report definitions with "print" enabled add columns.
   - The last column is the estimated remaining time `h:mm:ss` followed by the remaining iteration count. It is not the per-iteration time [I, forum consensus].
   - The header is reprinted periodically; transient runs add `Flow time = …s, time step = n` lines.
2. **Hardware.** The parallel table `ID  Hostname  Core  O.S.  PID  Vendor`, rows `n0…nN-1` plus `host`. Core is printed as `k/<cores on machine>` and Vendor is the CPU model string. Documented example: UG 44.7 `flu_ug_sec_parallel_connect.html`. Fluent prints this at parallel start-up [I, widely seen in logs] and on demand via Parallel → Network → Connectivity [D]. Rows = number of solver processes; distinct Hostnames = machines. GPU-solver sessions also print GPU info [U].
3. **Timing.** `parallel/timer/usage` output [D, UG 44.8 `flu_ug_parallel_monitor_optimize.html`]: "Performance Timer for N iterations on M compute nodes", "Average wall-clock time per iteration: x sec", "Total wall-clock time", "Simulation wall-clock time for N iterations". Present only if the user or journal asks for it; otherwise derive timing from the transcript or file timestamps [I].

### Report files [D/O]

UG 42.2.2 `flu_ug_reporting_sec_monitoring_solution.html`.

- Every report definition (Drag / Lift / Force / Moment …) gets its own report file and plot **by default** (Preferences → Automatically Plot and File).
- The file is written in the working directory. The "Full File Name" adds a run index (e.g. `report-file-0.out`).
- Format [O, cffview `examples/mass-rfile.out`]:

```
"mass-rfile"
"Iteration" "mass"
("Iteration" "mass(inlet)" "mass(outlet-oil)")
2948 0.02622239847736852 -0.02622218077450807
```

- Line 1 is the file title, line 2 the definition names, line 3 the column labels (per-zone columns are named `def(zone)`). Then one whitespace-separated row per write.
- Transient files add a time-step column and a `flow-time` column [I].
- **When a run is continued, iteration numbers continue and earlier rows remain**, so take the last occurrence of each iteration [O: the sample starts at 2948].
- Precision is full double [O].

---

## C. What to ask the team member for (no extra work for him)

**Ask for these from one finished run (same save):**

1. **`<case>.cas.h5` and `<case>.dat.h5`.**
   - Together they give us the mesh, every model, BC, numerics, URF and reference-value setting, and the report definitions.
   - They also give the final cell and face fields: p, u, k, ω, μ_t, wall shear, y⁺ and face pressure.
   - Plus residual history (decimated after 1000 points).
   - If he has several, ask for the final pair, not autosaves. Autosaves are named `<root>-<iteration>.cas.h5/.dat.h5` [I].
2. **All `*.out` report files from the run directory.** These hold the per-iteration CL, CD and force histories, written automatically for each report definition [D].
3. **Whatever console record exists:** a `.trn` (auto-transcript, `fluent-YYYYMMDD-HHMMSS-pid.trn`), or the batch or cluster stdout log.
   - This gives hardware (process table and CPU model), core count, wall clock and printed residuals.
   - If none exists, ask one question instead: machine, cores and roughly how long the run took.
4. **If he works in Workbench:** the system's Fluent folder (`<project>_files/dp0/FFF/Fluent/` [I]) contains all of the above. Zipping it is enough.

**Not needed:** the `.msh.h5` (the case file already contains the mesh), screenshots, and the force-report printout (we recompute it from the data file and cross-check against the `.out` files).

---

## Open items

1. **`kw-sst-strain? #f`.** Does Fluent's μ_t limiter use S or Ω by default? Test by recomputing μ_t from saved `SV_K`, `SV_O` and velocity gradients against `SV_MU_T`.
2. **Production limiter.** Is ε = β\*kω, and does the cap apply to G_ω too? Test the same way, in stagnation cells.
3. **Wall ω.** Imposed at the wall face or in the wall cell? Exact wall-cell production? Compare against saved `SV_WALL_SHEAR` and `SV_WALL_YPLUS`.
4. **WFGC, `J_ho` and the wall-distance method (=2).** Undocumented. Measure the sensitivity: run Fluent with WFGC off if the team can; otherwise quantify on our side.
5. **Inlet length-scale formula exponent** (C_μ vs C_μ^¼). Irrelevant if the viscosity-ratio spec is used, which is the default.
6. **Columns 2–3 of the `/results/residuals/*/data` (n,4) array.**
7. **Defaults not stated in the UG text:** Backflow Pressure Specification, Automatic Transcript, HOTR (seen only in a figure).
8. **Whether the GPU solver's discretization differs from the CPU solver's.** It has its own rpvars, e.g. `gpuapp/use-omega-wall-extrpln-bc? #t` [O]. This matters if the team runs Fluent GPU.

## Sources

- Ansys Fluent 2026 R1 Theory Guide: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_th/ (pages named inline).
- Ansys Fluent 2026 R1 User's Guide: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_ug/ (pages named inline).
- NASA/TMBWG Turbulence Modeling Resource, SST variants: https://tmbwg.github.io/turbmodels/sst.html
- Menter, F. R., AIAA J. 32(8), 1598–1605, 1994 (Fluent [432]); Menter, Kuntz & Langtry 2003 (via TMR).
- VTK Fluent CFF reader: https://github.com/Kitware/VTK/tree/master/IO/FLUENTCFF ; legacy reader: https://github.com/Kitware/VTK/blob/master/IO/Geometry/vtkFLUENTReader.cxx
- pyfluent offline readers: https://github.com/ansys/pyfluent/tree/main/src/ansys/fluent/core/offline
- cffview: https://github.com/aiwaonly/cffview (examples/general.scm, boundary.scm, mass-rfile.out; src/cffview/reader.py, boundary.py)
- meshio Fluent reader: https://github.com/nschloe/meshio/blob/main/src/meshio/ansys/_ansys.py
- Ansys example files (Fluent 24.1): https://github.com/ansys/example-data/tree/master/pyfluent/file_session
- FSAE sources: listed in §A5.
