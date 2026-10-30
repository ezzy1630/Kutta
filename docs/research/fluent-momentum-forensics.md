# Fluent momentum, energy, k and ω discretisation: forensics from the saved state

Date: 2026-10-29. Round 2 after `fluent-flux-gradient-forensics.md` ("round 1"). Companions: `fluent-sst-forensics.md` ("SST report") and `fluent-numerics-and-files.md` (its corrections box is authoritative).

**Implementation qualification, 2026-10-29:** the saved-state rankings do not
establish a consistent general boundary/flux scheme. V2 rejects the independently
reconstructed face mean, stored zero-difference pressure observations at skewed
face centroids, and an unrestricted cancellation of the symmetry normal
coefficient. Projected symmetry values must be associated with their normal
foot, and symmetry must still constrain normal velocity. See the later
interpolation-order investigation in
[`VERIFICATION.md`](../../src/kutta/fv/VERIFICATION.md) for the retained
formulation and the unresolved conflict with this saved-state fit.

**Data.** The same case as before: `validation/fluent_examples/mixing_elbow/mixing_elbow.{cas,dat}.h5`. It is Fluent 22.2, double precision, pressure-based Coupled, pseudo-transient, k-ω SST, LSQ 1/d² gradients, WFGC on, Second Order pressure, SOU for momentum/energy/k/ω, standard (Barth-Jespersen) limiter.

**Question.** Under which discrete momentum (then k, ω) equations is Fluent's saved state a fixed point? The method evaluates each candidate's steady residual per cell at the saved fields, with the saved `SV_FLUX` as the convecting flux. A candidate is scored by Fluent's own scaling (Σ|R| / Σ a_P|U|) against Fluent's final scaled residuals.

**Tags** (as in the earlier reports)

- **established**: reproduced to near round-off, or an unambiguous margin over every alternative.
- **likely**: the best of the tested alternatives by a clear margin, but on top of an unexplained floor.
- **open**: not settled.

**Scripts** (parse-only; usage as in `tools/forensics/cff.py`; each runs in 1–4 s once `cff.py` has cached the geometry)

- `tools/forensics/momentum.py`: building blocks (residual terms, limiter variants, scalar transport, Fluent's residual columns).
- `tools/forensics/momentum_residual.py`: §1, §2, §6.
- `tools/forensics/energy_residual.py`: §3, §4.
- `tools/forensics/turbulence_residual.py`: §1.3, §5.

Every number below is printed by one of them, except the items listed as scratch probes in §8.

---

## 0. Results

**Bottom line.** No tested discretisation makes the saved state a fixed point.

- **The floor is the same in every equation.** The best momentum operator leaves a scaled residual of about 2e-4 against Fluent's 1.4e-6 (99× in Σ|R|).
  - The scalar energy equation, which has no pressure, Rhie-Chow, transpose or wall-function terms, leaves 1.98 W against Fluent's 3.1e-3 W (630×).
  - k and ω are 220–1100× off.
  - Round 1's flux floor (3.3e-4·|F| on quiet faces) is the same phenomenon.
- **Every individual operator choice still has a sharp, unambiguous optimum**, and it is always at Fluent's documented setting. That makes the ranking usable even though the absolute target was not reached.
- **The floor is about 6–7 % of the SOU reconstruction term** in momentum and energy. It sits almost entirely on polyhedral cells and in steep or limited regions.
- **It is orthogonal to every operator perturbation tried.** That covers gradient weights, the WFGC jump correction, limiter forms (including a free per-cell limiter factor), reconstruction point, geometry, the convecting flux and diffusion scale.
- **Its cause is open.** §4 gives the evidence and the one data item that would settle it.

**Established**

1. **Momentum a_P, read from Fluent's own residual normalisation** (col1 = Σ a_P,i |U_P|, one value per component):
   a_P,i = 2 × [Σ max(F_out, 0) + Σ_f μ_eff α_f + ρV/Δτ + Σ_wall c_w A (1 − n_i²) + Σ_sym μ_eff A/Δn · n_i² + Σ_inlet μ_eff A/Δn]
   - α_f = |A|²/(A·d), c_w = τ_w/|U_t|, Δτ = rpvar `pseudo-auto-time-step`.
   - The x−y difference, which only walls produce, gives a factor of exactly **2.00000** on the per-component wall term c_w A (1 − n_i²).
   - The z−x difference, after the wall part, equals μ_eff A/Δn on the symmetry-normal component (ratio **1.00004**).
   - The totals match to 0.14 % (x, y) and 0.6 % (z). Without ρV/Δτ they are 9.5 % off.
   - The factor 2 is structural. Its cause is inferred (an implicit 1/0.5 on the pseudo-time momentum system, `mom/pseudo-relax 0.5`).
2. **Energy a_P = cp Σ max(F_out, 0) + Σ k_eff α + ρ V cp/Δτ**, with Σ a_P T over *absolute* T, matches Fluent's normalisation to 0.09 %. There is no factor 2 here.
3. **k a_P = Σ max(F_out, 0) + Σ Γ_k α + ρV/Δτ + ρβ*ωV** (implicit dissipation), to 0.26 %.
4. **Rhie-Chow d_f uses the physical a_P without ρV/Δτ.** The RC coefficient is 0.995 / 1.001 (two averaging forms) against 1.093 with ρV/Δτ and 2.19 with the normalisation's 2×a_P. This replaces round 1's proxy.
5. **SOU face value** = φ_up + ψ_up G_up·(x_f − x_up), with:
   - G the LSQ 1/d² gradient;
   - ψ the Barth-Jespersen factor, per component, over all face neighbours plus the cell's own boundary-face values, evaluated to face centroids;
   - x_f the area-weighted face centroid.

   Each choice is a sharp minimum in both momentum and energy (§2.2, §3.2).
6. **Momentum convection, interior and outlet:**
   - Outlet faces convect the **cell** velocity (2.4e-5 on outlet cells), not u_P + G·r_b (4.3e-5) or the saved face velocity (2.9e-5).
   - `SV_WALL_SHEAR` enters as the force on the fluid. The opposite sign gives 3.7×.
7. **Wall-cell k production** = harmonic combination p_lam·p_log/(p_lam + p_log), with:
   - p_lam = μ_t (τ_w/μ)²
   - p_log = τ_w²/(μ κ y⁺)
   - y⁺ = ρ u* Δy/μ
   - κ = 0.4187

   This gives a wall-cell k residual of **1.6e-6**, against 3.5e-5 for μ_t S² and ≥ 8e-5 for every other form tested (§5.1). It is what Kutta's `sst_wall` already computes.
8. **Boundary-face storage:** saved wall, symmetry and inlet `SV_P` equal the cell pressure exactly. Symmetry `SV_U/V/W` is the tangential projection of the cell velocity. The outlet face velocity is the cell velocity with a flux-consistent normal correction (median 1.5e-3 |U|).

**Likely** (clear margins over alternatives, on top of the floor)

9. **Boundary face pressure in the momentum equation is extrapolated, p_c + G·r_b,** at walls, symmetry and inlets. The saved p_b = p_c is a storage value. Wall cells 1.83 → 1.43e-4; symmetry cells 1.81 → 1.48e-4. Kutta's tangential-only symmetry form gives 1.65e-4.
10. **The pressure LSQ includes wall rows** with the boundary value p_b = p_c. Dropping them raises the wall-cell residual from 1.43e-4 to 2.35e-4. Inlet rows are mildly favoured: inlet cells go from 8.1e-6 to 9.2e-6 without them. Symmetry rows are neutral: 1.495 vs 1.491 e-4 on symmetry cells.
11. **The velocity LSQ includes symmetry rows** with the projected value. Dropping them, as Kutta does, raises the symmetry-cell residual from 1.81e-4 to 2.33e-4.
12. **The outlet face carries the transpose viscous flux μ_eff (∇u_P)ᵀ·A.** Outlet cells 2.43e-5 → 6.9e-6, a factor of 3.5.
13. **The symmetry plane exerts no net viscous force.** The μ_eff A/Δn term is in a_P (item 1) but is cancelled by a matching explicit term. Symmetry cells 1.81 → 1.61e-4; doubling it gives 2.04e-4.
14. **Production uses the same strain rate as the μ_t limiter.** Substituting Fluent's own S_F = ρa₁k/(μ_t F2) in limiter-active cells cuts their k residual 3.4× (3.2e-5 → 9.5e-6) and their ω residual 1.7×.
    - The best computed gradient for S is LSQ 1/d² without wall rows + ½ jump correction (round 1). It lowers the total k residual 1.6×.
15. **G_ω = α ρ P_k/μ_t** (= αρS² below the production limiter), not α(ω/k)P_k. The ω residual is 2.14 vs 2.48.
16. **The saved wall-face ω enters the ω LSQ rows.** Dropping the rows gives 2.31; using the cell value gives 2.55; the saved value gives 2.14.
17. **ω sink linearisation is Newton (2ρβω V)**, matching `kw-dissipation-newton? #t`. The normalisation matches to 1.2 % with fixed wall cells excluded, against 11.7 % with ρβωV.

**Open**

18. **The floor (§4).** It is common to all equations, sits on polyhedral cells, and is not explained by anything tested.
19. **Which gradient the momentum diffusion uses.** The data cannot tell: LSQ 1/d² and the WFGC candidate differ by < 0.5 %.
20. **Face μ_eff interpolation:** harmonic is 7 % better on x-momentum, neutral in energy.
21. **The ω residual in wall layers 2–4** stays at 0.55 even with Fluent's own S. Cross-diffusion is ambiguous: removing it improves layer 2.

---

## 1. Residual normalisation: what Fluent's a_P contains

`/results/residuals/phase-1/<eq>/data` stores the raw residual (col0) and the normalisation (col1) per iteration (SST report §5). At iteration 294:

| Equation | col0 (raw Σ\|R\|) | col1 | scaled |
|---|---|---|---|
| x / y / z momentum (N) | 2.357e-6 / 6.556e-6 / 1.701e-6 | 1.68637 / 1.68650 / 1.69876 | 1.40e-6 / 3.89e-6 / 1.00e-6 |
| energy (W) | 3.127e-3 | 231924 | 1.35e-8 |
| k | 2.496e-7 | 4.8375e-3 | 5.16e-5 |
| ω | 1.920e-3 | 91.093 | 2.11e-5 |
| continuity (kg/s) | 2.476e-6 | 2.614e-3 (max of iterations 1–5) | 9.47e-4 |

The run has three phases: iterations 1–100, 101–200 and 201–294. col1 jumps between them, so these are separate settings or runs. Only the last phase matches the saved fields.

### 1.1 Momentum (`momentum_residual.py` §1)

Pieces, each Σ(·)|U_P|:

| Piece | Value |
|---|---|
| outflow Σ max(F_out, 0) | 0.706397 |
| interior Σ μ_eff α | 0.061922 |
| ρV/Δτ | 0.074535 |
| walls c_w A (1 − n_i²), x / y / z | 8.674e-4 / 9.308e-4 / 5.496e-4 |
| symmetry μ_eff A/Δn n_i², x / y / z | 0 / 0 / 6.512e-3 |
| inlets μ_eff A/Δn | 6.70e-4 |

**Two exact checks**

- (col1_x − col1_y) / (wall_x − wall_y) = **2.00000**.
  - Walls are the only x/y-asymmetric term.
  - This fixes both the factor 2 and the per-component wall coefficient c_w A (1 − n_i²), with c_w = |SV_WALL_SHEAR|/(A |U_t|).
- (col1_z − col1_x)/2 − (wall_z − wall_x) = 6.512676e-3, against Σ μ_eff A/Δn n_z² |U| = 6.512385e-3 (ratio **1.00004**).
  - μ_eff is the *cell* μ + μ_t.
  - Δn is the normal distance from the centroid to the face plane.

**Totals**

| Form | col1 / formula, x / y / z |
|---|---|
| 2 × (all pieces) | 0.9986 / 0.9986 / 1.0063 |
| same without ρV/Δτ | 1.095 |
| 1 × (all pieces) | 1.997 |

The 0.14 % (x, y) remainder is consistent with the μ_eff interpolation choice. It is not resolved.

### 1.2 Energy (`energy_residual.py` §1)

- cp Σ outflow T = 189512; Σ k_eff α T = 21361; Σ ρ V cp/Δτ T = 21250.
- The sum, 232123, is 0.99914 of col1.
- Without the pseudo-time term the ratio is 1.0998.
- T is absolute, which is why Fluent's energy residual is so small when scaled.

### 1.3 k and ω (`turbulence_residual.py` §1)

| a_P candidate | k: col1/Σ | ω: col1/Σ, fixed wall cells excluded |
|---|---|---|
| out + dif + ρV/Δτ | 1.113 | 1.246 |
| + implicit sink ρβ*ωV (k) / ρβωV (ω) | **1.0026** | 1.117 |
| + Newton sink 2ρβωV (ω) | — | **1.012** |
| k sink without ρV/Δτ | 1.100 | — |

### 1.4 Rhie-Chow a_P (`momentum_residual.py` §8)

On round 1's quiet deep faces, using its flux model: a_P-weighted limited reconstructions + RC with G_m = Green-Gauss of the second-order face pressures.

| a_P in d_f | RC coefficient | median \|err\|/\|F\| |
|---|---|---|
| **outflow + μ_eff α** | **0.995** (1.001 with mean(V/a)) | 3.306e-4 |
| + ρV/Δτ | 1.093 | 3.356e-4 |
| 2 × (… + ρV/Δτ) | 2.19 | 4.26e-4 |

**Conclusion.** d_f uses the physical a_P, without pseudo-time and without the normalisation's factor 2, exactly as Kutta's `fv_rhie_chow` comment intends. Adding wall and symmetry terms to the trace changes nothing (scratch).

---

## 2. Momentum residual

### 2.1 Baseline and where it sits (`momentum_residual.py` §2)

**Baseline operator**

- SOU with BJ-limited LSQ 1/d².
- Second Order pressure with LSQ 1/d².
- Arithmetic μ_eff.
- Over-relaxed non-orthogonal correction with the averaged cell gradient.
- Transpose term.
- Saved wall shear as the wall force.
- Symmetry μ_eff A/Δn on u_n.

**Term sizes (scaled x)**

| Term | Scaled x |
|---|---|
| convection | 4.26e-3 |
| pressure | 4.20e-3 |
| interior viscous | 5.8e-4 |
| boundary viscous | 3.1e-4 |
| **total** | **2.18e-4 / 2.15e-4 / 1.87e-4** (x / y / z) |

Σ|R| is 1.05e-3 N, against Fluent's 1.06e-5 N.

**Distribution**

| Region | Share of Σ\|R\| | Median \|R\|/(outflow·\|U\|) |
|---|---|---|
| wall layer 1 | 17.5 % | 2.8e-2 |
| wall layer 2 | 11.3 % | 3.5e-3 |
| ≥ 3 hops from any boundary | 29.6 % | 5.9e-4 |
| … of which limiter inactive | 19.0 % | 5.7e-4 |
| limiter active (any component) | 50.0 % | 2.1e-3 |
| cells with 6 faces (772) | 0.9 % | **6.3e-5** |
| cells with > 6 faces | 99.1 % | **1.37e-3** |

**Single-term scale scans** (deep Σ|R| for R + (c − 1)·T, c = 0.9 / 0.95 / 1.0 / 1.05 / 1.1):

- SOU correction: 3.71 / 2.72 / 3.10 / 4.91 / 7.05 e-4.
- Pressure second-order correction: 3.17 / 3.09 / 3.10 / 3.20 / 3.35 e-4 (flat near 1).

Both minima are within 5 % of c = 1, and neither removes more than 12 %, so the remainder is not a mis-scaled term.

### 2.2 Convective reconstruction (`momentum_residual.py` §3)

Deep Σ|R| (scaled x):

| Gradient / limiter | deep Σ\|R\| | scaled x |
|---|---|---|
| LSQ 1/d^1 / 1.5 / **2** / 2.5 / 3, BJ | 5.22 / 3.49 / **3.10** / 4.68 / 6.59 e-4 | 3.19 / 2.40 / **2.18** / 3.14 / 4.58 e-4 |
| LSQ2 + t·J (jump correction), t = −0.25 / 0.25 / 0.5 | 4.38 / 4.38 / 6.82 e-4 | |
| LSQ2 unlimited | 4.70e-4 | 3.43e-4 |
| BJ without boundary values in min/max | 3.10e-4 | 2.42e-4 |
| BJ with one factor for all components | 6.63e-4 | |
| smooth Venkatakrishnan function | 3.40e-4 | |
| Green-Gauss (second-order faces / cell) + BJ | 7.39 / 9.46 e-4 | |
| velocity LSQ without wall rows | 3.10e-4 | 2.13e-4 |

**Conclusion.** The reconstruction gradient is plain LSQ 1/d², not the jump-corrected gradient. The limiter is per-component Barth-Jespersen with boundary values in the min/max. This confirms round 1's flux result on an independent equation.

### 2.3 Pressure (`momentum_residual.py` §4)

| Face pressure gradient | scaled x | deep |
|---|---|---|
| LSQ 1/d^1 / **2** / 3 | 2.21 / **2.18** / 2.25 e-4 | 3.11 / 3.10 / 3.26 e-4 |
| LSQ2, BJ-limited | 2.13e-4 | 3.06e-4 |
| Green-Gauss of second-order faces (round 1's RC gradient) | 2.55e-4 | 3.64e-4 |
| Green-Gauss cell-based | 3.18e-4 | 5.89e-4 |
| LSQ2 + ½ jump correction | 2.45e-4 | 3.42e-4 |
| no reconstruction (arithmetic) | 4.22e-4 | 7.90e-4 |

**Conclusions**

- The momentum pressure term is Σ p_f A_f with second-order faces built from the LSQ 1/d² pressure gradient.
- The RC jump uses Green-Gauss of those same faces, i.e. (Σ p_f A)/V (round 1). That is the standard consistent pairing.
- Limiting the pressure gradient is worth 1 %. That is inside the noise of the floor and left open.

**Boundary pressure** (scaled z 1.87e-4 → 1.53e-4). Symmetry cells:

| Symmetry face pressure | Symmetry-cell Σ\|R\| |
|---|---|
| p_c + G·r_b | 1.48e-4 |
| p_c + G·t (Kutta) | 1.65e-4 |
| saved p_c | 1.81e-4 |

### 2.4 Viscous terms and the WFGC gradient (`momentum_residual.py` §5)

All variants lie within ±7 % of the baseline:

| Variant | scaled x |
|---|---|
| baseline | 2.18e-4 |
| no transpose | 2.18e-4 (z 1.87 → 1.78) |
| − ⅔ μ ∇·u | 2.18e-4 |
| harmonic μ_eff | **2.03e-4** |
| normal-distance-weighted μ_eff | 2.14e-4 |
| distance-weighted face gradient | 2.17e-4 |
| no non-orthogonal correction | 2.27e-4 |
| minimum-correction split | 2.18e-4 |
| WFGC gradient in Laplacian and transpose | 2.18e-4 (deep 3.108 vs 3.102 e-4) |

**Role of the jump-corrected (WFGC) gradient**

- It is rejected for convective reconstruction in momentum (§2.2), energy (§3) and flux (round 1).
- It is neutral in diffusion.
- It is the best computed gradient for the strain rate S (round 1, §5.2).

### 2.5 Boundary treatments (`momentum_residual.py` §6)

Σ|R| (x + y + z) on the cells next to each boundary type:

| Candidate | wall cells | layer 2 | sym cells | inlet cells | outlet cells |
|---|---|---|---|---|---|
| baseline | 1.830e-4 | 1.188e-4 | 1.811e-4 | 7.7e-6 | 2.43e-5 |
| **boundary p = p_c + G·r_b** | **1.429e-4** | 1.169e-4 | **1.476e-4** | 6.8e-6 | 2.40e-5 |
| … + pressure LSQ without wall rows | 2.351e-4 | 1.368e-4 | 1.495e-4 | 8.0e-6 | 2.37e-5 |
| … + pressure LSQ without wall/symmetry rows | 2.373e-4 | 1.395e-4 | 1.491e-4 | 8.1e-6 | 2.37e-5 |
| … + pressure LSQ outlet rows only (Kutta) | 2.376e-4 | 1.395e-4 | 1.493e-4 | 9.2e-6 | 2.37e-5 |
| velocity LSQ without wall rows | 1.678e-4 | 1.090e-4 | 1.776e-4 | 1.75e-5 | 2.40e-5 |
| velocity LSQ without symmetry rows | 1.884e-4 | 1.266e-4 | 2.334e-4 | 7.7e-6 | 2.44e-5 |
| wall normal viscous stress added | 2.291e-4 | 1.188e-4 | 1.824e-4 | 9.6e-6 | 2.48e-5 |
| **no symmetry normal force** | 1.832e-4 | 1.184e-4 | **1.605e-4** | 7.8e-6 | 2.41e-5 |
| symmetry 2 μ_eff A/Δn | 1.829e-4 | 1.194e-4 | 2.040e-4 | 7.7e-6 | 2.46e-5 |
| outlet convects saved face velocity | 1.837e-4 | 1.195e-4 | 1.814e-4 | 7.7e-6 | 2.90e-5 |
| outlet convects u_P + G·r_b | 1.858e-4 | 1.224e-4 | 1.814e-4 | 7.7e-6 | 4.34e-5 |
| **outlet carries μ(∇u)ᵀ·A** | 1.807e-4 | 1.159e-4 | 1.801e-4 | 7.7e-6 | **6.9e-6** |
| inlet orthogonal diffusion only / none | 1.829 / 1.884 e-4 | | | 7.7e-6 / 1.13e-5 | |

**Wall-cell velocity gradient** (round 1 open item). Dropping wall-face rows improves wall cells (−8 %) and layer 2 (−8 %), but doubles inlet-cell residuals, because inlet cells touch `wall-inlet`. It is not decisive. The established part: the wall-cell gradient is not zero (round 1), and the wall face carries only the wall-function traction, with no normal stress.

### 2.6 Assembled (`momentum_residual.py` §7)

| Assembly | scaled x / y / z | deep Σ\|R\| |
|---|---|---|
| best (extrapolated boundary p, outlet transpose, no symmetry force) | 2.04 / 2.13 / 1.36 e-4 | 3.10e-4 |
| Kutta-like* | 2.42 / 2.55 / 1.88 e-4 | 3.07e-4 |
| Fluent | 1.40 / 3.89 / 1.00 e-6 | |

\* Kutta-like: velocity LSQ rows only at walls and inlets; pressure LSQ rows only at outlets; normal-distance-weighted μ_eff and face gradient; reconstructed outlet velocity; full outlet viscous flux.

The interior is the same for both; Kutta's differences are all at boundaries. The wall shear is taken from Fluent here; Kutta computes its own.

---

## 3. Energy: the clean test bed (`energy_residual.py`)

The energy equation shares only the scalar SOU machinery with momentum. Its fixed point is a direct test of reconstruction, limiter and diffusion.

### 3.1 Data facts

- `SV_H` = cp (T − 298.15) exactly.
- Walls are adiabatic: T_b = T_c to 6e-14.
- Inlets are at 293.15 / 313.15 K.
- `SV_HEAT_FLUX` at inlets and the outlet is not a diffusive flux (scratch: 6.5e5 × k_eff A ΔT/Δn), so it is unused.
- Global balance closes: signed Σ R = −2.4e-3 W against an enthalpy throughput of 18.75 W.

### 3.2 Operator candidates

Σ|R| in W; Fluent 3.1e-3 W:

| Candidate | Σ\|R\| (W) | deep |
|---|---|---|
| first-order upwind | 29.56 | 7.75 |
| SOU, LSQ2 unlimited | 5.61 | 1.72 |
| **SOU, LSQ2 + BJ** | **1.984** | 0.560 |
| LSQ 1/d^1 / 1.5 / 2.5 / 3 + BJ | 2.51 / 2.13 / 2.35 / 2.97 | |
| LSQ2 + 0.25 / 0.5 jump correction | 3.66 / 5.94 | |
| harmonic k_eff / no non-orthogonal correction | 1.995 / 2.009 | |
| BJ without boundary values in min/max | 1.984 | |
| smooth Venkatakrishnan function | 2.011 | |
| two-cell bounds per face (`allneighbors? #f` reading) | 23.5 | |
| BJ to the neighbour centroid (cell-to-cell) | 7.52 | |
| BJ bounds widened by 1 % / 5 % of the local range | 2.34 / 3.76 | |
| reconstruct to the centroid-line/face-plane intersection, or to the midpoint | 6.05 / 6.95 | |
| diffusion × 0.8 / 0.9 / **1.0** / 1.1 / 1.25 | 2.37 / 2.14 / **1.98** / 1.93 / 1.96 | |

The same operator wins as in momentum. The remainder is 7.2 % of the SOU correction (27.6 W).

### 3.3 Where it sits

- 1.80 of the 1.98 W is in cells with |∇T| > 100 K/m, the mixing layer. 32 % of those cells are limited, against 56–66 % in the quieter bins.
- Cells whose own and upwind reconstructions are *all* unlimited carry 0.23 W. There the per-cell ratio R / SOU correction is scattered around zero: p10 / p50 / p90 = −0.081 / +0.007 / +0.096.

### 3.4 Neither the gradient nor the limiter explains it

- **Gradient-operator basis fit** on the unlimited cells, where the residual is linear in the gradient. Directions: jump correction, LSQ1 − LSQ2, LSQ3 − LSQ2, GG(second-order) − LSQ2, GG cell − LSQ2.
  - Every coefficient is 0.025–0.036.
  - The joint fit removes 4 % (0.233 → 0.222).
- **Bounded per-cell limiter fit.** A free ψ_c ∈ [0, 1] in each of the 9,497 cells BJ limits, chosen to minimise the residual.
  - It reaches only 1.28 W, still 400× Fluent.
  - The fitted ψ correlates 0.88 with BJ's (median 0.30 vs 0.27).
  - No limiter can account for the floor, and BJ is close to the best possible one.

---

## 4. The floor

### What it is not

Tested, each with ≤ 7 % effect:

- **Gradient:** LSQ weight power; weight definitions (1/|r_f|², 1/(d·n)², area-weighted, cos²-weighted, doubled boundary rows); WFGC jump correction at any strength; Green-Gauss variants; linear combinations of all of these.
- **Limiter:** Barth-Jespersen variants, Venkatakrishnan, two-cell bounds, cell-to-cell, ε-widened bounds, face-wise clipping, and a free bounded ψ per cell.
- **Reconstruction point:** intersection point, midpoint, normal projection.
- **Scheme:** central, averaged-gradient SOU and blended variants are all 1.2–3.5× worse.
- **Geometry:**
  - face centroids by projected-area fan, |area| fan, iterated fan, first-node fan, projection onto the mean plane;
  - cell-centroid noise of 1e-5 to 1e-4·h changes Σ|R| by < 0.05 %;
  - cells close to 1e-15, Σ x_f·A = 3V to 1e-13, LSQ exact on linear fields to 1e-12;
  - all inputs are float64.
- **Convecting flux:** round 1's reconstructed RC flux instead of `SV_FLUX` gives 2.02 W (non-conservative) / 5.99 W (conservative).
- **Diffusion:** total, laminar and turbulent scales, k_eff interpolation.

### What it is

- About 6–7 % of the SOU correction in momentum and energy.
- On polyhedral cells: 6-face cells have a 22× smaller relative momentum residual.
- Concentrated where gradients are steep or the limiter acts.
- Neighbouring residuals are positively correlated (+0.21).
- The implied state error δT = (first-order + pseudo-time operator)⁻¹ R is smooth (neighbour correlation +0.64) and localised: max 0.55 K, median 1.3e-5 K.

### Two hypotheses

**H1. An undocumented reconstruction ingredient on polyhedral (warped-face) cells.** This would be part of WFGC (`recon/cell-lsf? #t`, `recon/wfgc-jump-enhanced? #t`, `recon/cell-lsq/cortype 1`) that changes face values in a way no linear combination of the tested operators reproduces.

- For: 6-face cells are 22× better.
- For: round 1's S analysis independently found the derivative gradient differs from LSQ 1/d² by a jump-type correction whose per-cell coefficient scatters ±17 %.

**H2. Fluent's reported residual is not the imbalance at the saved fields** (for example, it is evaluated after the linear solve, with the explicit relaxation 0.5 / 0.75 applied afterwards). The saved state would then still drift by about 1e-4 relative per iteration.

- For: every knob sits at its documented optimum.
- For: the floor is orthogonal to all of them, has a similar relative size in every equation, and implies a smooth state error.
- Against: the residual history behaves like a nonlinear imbalance. In phase 1 the k and ω residuals chatter with period 2 (2.390e-3 / 2.416e-3, 6.758e-3 / 6.559e-3; `momentum_residual.py` header).
- Against: Fluent's continuity residual is computed from the fields. If the state drifted at the H2 rate, the next iteration's value would jump about 300×.

The data favours H1 but cannot exclude H2.

### What would settle it

1. **Two consecutive Fluent data files** (iterations 294 and 295, or 294 and 300) for this case.
   - The field change per iteration measures H2 directly.
   - Under H1 the change is ~1e-6 relative; under H2 it is ~1e-4 and shaped like δT above.
2. **Fluent's per-cell residual fields**, written as a data quantity with the expert option that keeps them. These locate the floor cell by cell.
3. **A mixing-elbow re-run with WFGC off** (`recon/cell-lsf? #f`). If the floor collapses to Fluent's level on that file, H1 is confirmed and the jump-correction family is the place to look.

---

## 5. k and ω (`turbulence_residual.py`)

### 5.1 Wall-cell k production: established

Wall-cell k Σ|R| over 2,436 cells. Fluent's whole-domain k residual is 2.5e-7.

| Production in wall cells | Σ\|R\| |
|---|---|
| μ_t S² (LSQ2), no wall treatment | 3.49e-5 |
| none | 2.27e-4 |
| τ_w²/μ | 1.0e-3 |
| τ_w u*/(κΔy) / τ_w u_τ/(κΔy) | 4.2e-4 / 3.6e-4 |
| μ_t(τ_w/μ)² alone / τ_w²/(μκy⁺) alone | 2.9e-4 / 3.0e-4 |
| **harmonic(μ_t(τ_w/μ)², τ_w²/(μκy⁺)), κ 0.4187, y⁺ from u*** | **1.60e-6** |
| same, κ 0.41 | 3.3e-6 |
| same, y⁺ from u_τ | 1.5e-5 |
| power-2 / power-4 blend / min | 8.3e-5 / 1.25e-4 / 1.46e-4 |

In wall cells the k equation is almost at Fluent's fixed point. That is consistent with the floor living in the SOU reconstruction, which matters little in these cells.

### 5.2 Strain rate and production

Production with LSQ2 S and the harmonic wall production:

| S gradient | k Σ\|R\| (scaled) | ω Σ\|R\| (scaled) |
|---|---|---|
| LSQ2 | 5.63e-5 (1.16e-2) | 2.14 (2.35e-2) |
| LSQ2 + ½J | 5.15e-5 (1.07e-2) | 2.16 (2.37e-2) |
| LSQ2 without wall rows + ½J | **3.47e-5 (7.2e-3)** | **1.86 (2.04e-2)** |
| Fluent | 2.50e-7 (5.2e-5) | 1.92e-3 (2.1e-5) |

- The k/ω residual is 53–59 % of the SOU correction, against 7 % for momentum and energy. The turbulence equations have their own deficit on top of the floor.
- By wall layer (LSQ2):
  - k: 1.6e-6 / 1.8e-5 / 1.1e-5 / 7.2e-6 / 3.3e-6 for layers 1 / 2 / 3 / 4 / deep.
  - ω: wall cells are fixed; 0.89 / 0.45 / 0.35 / 0.015 for layers 2 / 3 / 4 / deep.
- **Fluent's own S.** Using S_F = ρa₁k/(μ_t F2) where the μ_t limiter is active (6,756 cells):
  - The k residual of those cells falls 3.2e-5 → 9.5e-6 (LSQ2 base) and 1.8e-5 → 9.5e-6 (nowall + ½J base).
  - Their ω residual falls 0.94 → 0.55.
  - So production and the μ_t limiter use one strain rate, and S accuracy controls much of the k residual.

### 5.3 ω near-wall ingredients

S from LSQ2:

| Variant | ω Σ\|R\| | L2 | L3 | L4 |
|---|---|---|---|---|
| saved wall-face ω in the LSQ rows | 2.14 | 0.89 | 0.45 | 0.35 |
| ω LSQ without wall rows | 2.31 | 1.06 | 0.45 | 0.35 |
| wall-face ω = cell ω | 2.55 | 1.30 | 0.45 | 0.35 |
| harmonic Γ_ω | 2.13 | 0.89 | 0.44 | 0.34 |
| no cross-diffusion | 1.99 | 0.75 | 0.45 | 0.33 |
| S from nowall + ½J | 1.86 | 0.75 | 0.34 | 0.34 |
| G_ω = α(ω/k)P_k | 2.48 | 1.14 | 0.51 | 0.34 |

**Not tested:** the production limiter, which is inactive in this case (SST report).

---

## 6. Items asked to pin down

- **Limiter stencil: established.** Barth-Jespersen, scalar per component, cell-to-face (face centroids), min/max over the cell, all face neighbours, and the cell's own boundary-face values.
  - Boundary values matter at boundaries: 2.18e-4 with them vs 2.42e-4 without.
  - Vector, two-cell, cell-to-cell, ε-widened and smooth variants are all worse.
  - `recon/lim/allneighbors? #f` therefore means "face neighbours, not node neighbours".
  - `recon/limit/epsilon 0.01` has no measurable role here.
- **Wall-cell velocity gradient: open (weakly likely: no wall-face rows).** The wall face carries only the wall-function traction (§2.5).
- **a_P in Rhie-Chow: established** (§1.4). It is the physical convection + diffusion diagonal, without ρV/Δτ or the factor 2.

---

## 7. Kernel changes Kutta needs

Anchors are function names in the worktree `/Users/ezzy/.t3/worktrees/Kutta/t3-265ec8e3`. It was read, not modified; another agent is editing it.

These sit on top of round 1's §4, which still stands. That includes the flux interpolant, limited reconstructions in `fv_momentum_faces`, the RC jump with (Σ p_f A)/V, outlet flux from cell velocity, and the SST strain gradient.

1. **Outlet convected momentum** (`fv_velocity_boundary`, kinds 2/4 with outflow): set the boundary velocity to u_P, not u_P + G·r_b. **Established** (§2.5).
   - Today this face value is both convected and used in the gradient rows.
   - Keep the face reconstruction out of convection.
2. **Boundary face pressure** (`fv_pressure_boundary`): use p_c + G·r_b at symmetry faces too, i.e. drop the tangential projection there. **Likely.**
3. **Pressure LSQ stencil** (`pressure_lsq`, `independent_pressure` in `steady_solver.cpp`): include wall faces, and preferably inlet faces, as rows, with the value p_b = p_c (zero-gradient data). Keep the extrapolation for the face pressure itself.
   - **Likely** for wall rows: wall cells 1.43e-4 with them, 2.38e-4 with Kutta's outlet-only rows.
   - Inlet rows: inlet cells 8.1e-6 → 9.2e-6 without them.
   - Symmetry rows are neutral (1 %).
4. **Velocity LSQ stencil** (`velocity_lsq`, `independent_velocity`): include symmetry faces with the projected cell velocity u_P − (u_P·n)n as the row value. **Likely** (1.81 vs 2.33 e-4). Outlet rows are neutral.
5. **Symmetry momentum** (`fv_velocity_boundary` type 3 + `fv_transport_assemble` + `fv_transpose_stress`):
   - Face velocity = projected **cell** velocity, without the G·t reconstruction. **Established** as Fluent's stored value.
   - **Likely:** no net viscous force, i.e. zero tangential diffusive flux, no transpose flux and no net normal force.
   - Keep μ_eff A/Δn n_i² on the diagonal for the residual scale and stability, balanced by the same explicit term.
   - Kutta today applies D(u_b − u_P) with the reconstructed tangential value plus a transpose flux.
6. **Outlet viscous flux** (`fv_transpose_stress` boundary faces): keep μ_eff(∇u_P)ᵀ·A at outlets. **Likely**, factor 3.5. Kutta already does this.
7. **Wall shear** (`fv_wall_momentum`, `fv_sst_force`, `sst_wall`): no change.
   - Fluent's structure c_w A(1 − n_i²) in a_P per component, with the wall-function traction replacing both stress parts, is what Kutta does. **Established.**
   - No wall-normal viscous stress.
8. **Wall-cell k production** (`sst_wall`): no change. Kutta's harmonic blend is Fluent's. **Established.**
9. **Strain rate for production and the μ_t limiter** (`sst_rates` input in `SstModel::update`): one S for both, from LSQ 1/d² without wall-face rows + ½ jump correction (round 1 item 4). **Likely**: k residual 1.6× lower.
   - Keep plain LSQ 1/d² for SOU reconstruction everywhere.
10. **G_ω** (`sst_sources`): keep the α ρ P_k/μ_t form. **Likely.**
11. **ω sink linearisation** (`sst_sources` → `omega_sink`): use the Newton form 2ρβω for the implicit coefficient, with RHS +ρβω². This changes convergence, not the fixed point. **Established** from the normalisation, 1.2 %.
12. **Residual monitors** (`fv_monitor_cell`, `fv_scalar_monitor`). Needed if Kutta's printed scaled residuals are to be comparable with Fluent's:
    - momentum denominator = Σ 2 (a_P,i + ρV/Δτ)|U| per component, not trace/3 without pseudo-time;
    - energy denominator = Σ (a_P + ρVcp/Δτ) T_abs;
    - k denominator = Σ (a_P + ρV/Δτ + ρβ*ωV)|k|;
    - ω denominator = Σ (a_P + ρV/Δτ + 2ρβωV)|ω|, excluding fixed wall cells.
    - **Established** from §1.
13. **Not supported by this data:**
    - the jump-corrected gradient in convective reconstruction or diffusion (rejected or neutral);
    - harmonic μ_eff (open);
    - any limiter change beyond the Barth-Jespersen definition above.

**Regression gates these support** (on any Fluent case/data pair; all are ratchets, not round-off gates, until §4 closes):

- momentum scaled residual at the saved state ≤ 2.1e-4 (best assembly);
- energy ≤ 2.0 W;
- wall-cell k residual ≤ 2e-6;
- the a_P normalisation identities of §1 to ≤ 0.2 % (momentum x/y) and ≤ 0.1 % (energy).

---

## 8. Scratch probes (not in the kept scripts)

Each is reproducible in a few lines from `momentum.py`:

- face-centroid variants (|area| fan, iterated fan, node average: momentum 4.2e-4) and first-node fan / mean-plane projection (energy 2.005 / 1.984 W);
- cell-centroid noise;
- extra LSQ weight definitions;
- face-wise clipping (energy 5.42 W);
- the exactly determined free-ψ solve (residual 7e-3 W, but ψ from −19 to +22: meaningless);
- round 1's flux as the convecting flux;
- per-face geometry correlations (all |ρ| < 0.16);
- `SV_HEAT_FLUX` ≠ diffusive flux;
- RC coefficient with the wall and symmetry trace in a_P (unchanged, 0.995).

## 9. Sources and clues

- Fluent Theory Guide 2024 R2 §23.3.3.4: LSQ cell-based via Gram-Schmidt weights. It gives no weighting or WFGC formula. https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/flu_th_sec_eval_derivatives.html
- Fluent Theory Guide 2024 R2 §23.3.4, gradient limiters (formulas are images). https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_th/flu_th_sec_grad_limit.html
- Fluent User's Guide 2024 R2 §36.15: unscaled residual = summed imbalance; momentum scaled by Σ a_P|U|. https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_sec_solve_monitor.html
- Case rpvars used as clues:
  - `mom/pseudo-relax 0.5`, `pseudo-auto-time-step 0.04265…`, `kw-dissipation-newton? #t`
  - `recon/lim/allneighbors? #f`, `recon/limit/epsilon 0.01`, `recon/lim/dx-to-face? #t`
  - `recon/cell-lsf? #t`, `recon/wfgc-jump-enhanced? #t`, `recon/cell-lsq/cortype 1`
  - `mom/symmetry-option 2`, `mom/shear-option 2`, `poutlet/full-viscous-force? #f`
  - `energy/inlet-diffusion? #t`, `keprt 0.85`
