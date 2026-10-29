# Kutta verification and validation cases: reference data

Compiled 2026-10-29. Every number below comes from a file in `validation/data/` or from a cited source. Each data directory has a `SOURCE.md` with URLs, the fetch date, license terms and a sha256 for every file. All downloaded files were treated as untrusted data: they were parsed, never executed.

Kutta solves incompressible RANS with k-omega SST. The ladder below goes from code verification (does the discretization converge at the design order?), to model verification (is SST implemented as the TMR defines it?), to validation against experiment, to cross-checks against Ansys Fluent. The end target is Fluent agreement on FSAE car aero.

## Summary

| # | Case | Verifies | Data on disk | Status |
|---|---|---|---|---|
| 1 | Lid-driven cavity, Re 100/400/1000 | Laminar NS discretization, order of accuracy | `ghia1982_cavity/`, `erturk2005_cavity/` | Downloaded and cross-checked |
| 2 | Blasius laminar flat plate | Laminar BL, wall shear, leading-edge handling | formulas below | Formula only |
| 3a | TMR 2DZP flat plate, SST | SST implementation, Cf/CD grid convergence, incompressible offset | `nasa_tmr/2DZP_flatplate/` | Downloaded |
| 3b | TMR 2DBC bump-in-channel, SST/SST-V/SST-2003 | SST with pressure gradient and curvature, CL/CDp/CDv | `nasa_tmr/2DBC_bump/` | Downloaded |
| 3c | TMR 2DN00 NACA 0012, SST/SST-V | Lifting airfoil, farfield, against experiment | `nasa_tmr/2DN00_naca0012/` | Downloaded |
| 4 | Selig S1223, Re 1e5-3e5 | Low-Re high-lift airfoil against experiment | `uiuc_lsat_s1223/` | Downloaded |
| 5 | Zerihan single-element wing in ground effect | Ground effect, moving ground, 2D and 3D | `zerihan_tyrrell026_ground_effect/` | Derived from thesis |
| 6 | Ahmed body 25/35 deg | 3D bluff-body separation, drag, Fluent cross-check | `ercoftac_ahmed_body/` | Downloaded |
| 7 | FSAE car stand-ins | Full-car workflow | none (licensing) | Referenced only |

## Compressibility offset (applies to every TMR case)

All TMR reference results come from compressible codes run at M = 0.15-0.2. TMR says so explicitly: incompressible results "may be close - but not quite the same - as the grid is refined". Expected offsets:

- **Skin friction (measured, 2DZP, M = 0.2):** TMR supplies results from incompressible codes in `2DZP_flatplate/sst/cf_incomp_results_sstv.dat`. Cf at x = 0.97 on the 545x385 grid is **0.0027174** (FUN3D incompressible) and **0.0027212** (SC/Tetra), against **0.0026909/0.0026905** from compressible CFL3D/FUN3D. Incompressible Cf is therefore **about +1.0% higher**. Compare Kutta to the incompressible references first.
- **Pressure-driven quantities (CL, Cp, CDp):** Prandtl-Glauert scaling 1/sqrt(1-M^2) predicts the compressible values are **about 1.1% higher at M = 0.15** (NACA 0012) and **about 2.1% higher at M = 0.2** (bump). Expect Kutta's incompressible CL to be lower by roughly that much. For NACA 0012 at alpha = 10 deg, that gives CL about 1.066 instead of 1.0778. This is an inviscid-linear estimate, so treat the real offset as 0.5-1.5% (inferred, not measured).
- **Freestream turbulence in incompressible units.** TMR sets k/a_inf^2 = 9e-9 and omega*mu/(rho*a_inf^2) = 1e-6. Converting with U = M*a_inf gives:

| case | M | Re (per unit length) | k/U^2 | omega*L/U | Tu | mu_t/mu |
|---|---|---|---|---|---|---|
| 2DZP | 0.2 | 5e6 | 2.25e-7 | 125 | 0.039% | 0.009 |
| 2DBC | 0.2 | 3e6 | 2.25e-7 | 75 | 0.039% | 0.009 |
| 2DN00 | 0.15 | 6e6 | 4.0e-7 | 266.7 | 0.052% | 0.009 |

  The Tu and mu_t/mu columns reproduce the values TMR states on its pages. For 2DZP, TMR also shows that changing Tu from 0.01% to 1% and mu_t/mu from 0.001 to 10 moves Cd and Cf by less than 0.1% (table on `flatplate_sst.html`), so these settings are not sensitive.

## 1. Lid-driven cavity (Ghia, Ghia & Shin 1982; Erturk et al. 2005; Botella & Peyret 1998)

**Verifies:** laminar incompressible NS discretization, the pressure-velocity coupling, and the observed order of accuracy. This is code verification, not turbulence.

**Conditions:** unit square, lid at y = 1 moving at u = 1, no-slip walls, Re = U*L/nu = 100, 400, 1000. Steady.

**Reference data:**
- `validation/data/ghia1982_cavity/ghia1982_u_vertical_centerline.csv`: u(y) at x = 0.5, 17 points for each Re.
- `validation/data/ghia1982_cavity/ghia1982_v_horizontal_centerline.csv`: v(x) at y = 0.5, 17 points for each Re.
- Ghia's own solution is 2nd order on a 129x129 grid. The coordinates `*_exact = (i-1)/128` are the true sample locations.
- The CSVs are a majority vote over six independent transcriptions, because no machine-readable original exists. Three disagreements were resolved this way; details are in `SOURCE.md`. **v(Re=400, x=0.9063) = -0.23827 is a suspected typo in the original paper. Exclude it from pass/fail.**
- `validation/data/erturk2005_cavity/`: Re = 1000 centerline u and v from a 601x601 grid (23 points each), plus primary-vortex data. **Botella & Peyret spectral: |psi|min = 0.1189366, omega_c = 2.067753, centre (0.5308, 0.5652).** Ghia gives 0.117929 at (0.5313, 0.5625), which is 0.85% low. Erturk's 601^2 2nd-order solution gives 0.118781 (0.13% low).

Selected values:

| quantity | Re=100 | Re=400 | Re=1000 |
|---|---|---|---|
| u(y=0.5) Ghia | -0.20581 | -0.11477 | -0.06080 (Erturk -0.0620) |
| min u on the vertical centerline (Ghia) | -0.21090 at y=0.4531 | -0.32726 at y=0.2813 | -0.38289 at y=0.1719 |
| v(x=0.5) Ghia | 0.05454 | 0.05186 | 0.02526 (Erturk 0.0258) |

**Recommended acceptance:**
- **Order of accuracy (the real test):** at Re = 1000 on grids of 64^2, 128^2, 256^2 and 512^2, compute |psi|min. Its error against Botella's 0.1189366 should fall with an observed order of **2.0 +- 0.3** for a 2nd-order scheme. Why: Botella is spectral, so its error is about 1e-7, far below ours, which makes it a valid exact solution. Erturk's sequence (401 -> 513 -> 601: 0.30% -> 0.18% -> 0.13% low) is a check on the expected magnitude. Also run the method of manufactured solutions for formal order. That needs no external data and isolates the scheme from the lid-corner singularity, which can pull the observed order on the cavity below 2.
- **Profile agreement:** on 256^2 or finer, max |u - u_Ghia| and max |v - v_Ghia| should be **<= 0.01** for Re = 100 and 400, and **<= 0.02** for Re = 1000. Against Erturk at Re = 1000, use **<= 0.005**. Why: Ghia's own 129^2 discretization error is about 1% of U at Re = 1000 (the psi_min gap above), so demanding tighter agreement with Ghia would penalise a more accurate solver.

## 2. Laminar flat plate (Blasius)

**Verifies:** laminar boundary-layer resolution, wall shear evaluation, and inflow and leading-edge handling. No download is needed.

**Conditions:** zero-pressure-gradient plate starting at x = 0, uniform U, Re_x = U*x/nu. Use about Re_L = 1e4-1e5 so the flow stays laminar and the BL is thick enough to resolve, with a slip or symmetry boundary ahead of the leading edge, as in the TMR 2DZP layout.

**Reference (Blasius similarity solution, f''(0) = 0.332057):**
- Cf(x) = 0.664 / sqrt(Re_x)
- One-sided CD over length L: 1.328 / sqrt(Re_L)
- delta* = 1.7208 x / sqrt(Re_x); theta = 0.664 x / sqrt(Re_x); H = 2.59
- delta_99 = 4.91 x / sqrt(Re_x)
- u/U profile: f'(eta) with eta = y*sqrt(U/(nu*x)). Solve f''' + 0.5*f*f'' = 0 (f(0) = f'(0) = 0, f'(inf) = 1) numerically in the test harness instead of storing a table.

**Recommended acceptance:** Cf within **1%** of Blasius for 0.1 < x/L < 0.9 on a grid with at least 20 cells across delta_99. Profile u/U within 0.5% of U. Observed order about 2 under refinement. Exclude x/L < 0.1: Blasius has a leading-edge singularity and boundary-layer theory has O(Re_x^-1/2) corrections there that the NS solution correctly captures.

## 3. NASA TMR SST verification and validation

Source: https://tmbwg.github.io/turbmodels/, which moved from turbmodels.larc.nasa.gov; GitHub `TMBWG/turbmodels`, CC0. The SST variants differ:

- **SSTm:** standard SST with production approximated by mu_t*S^2 (exact for incompressible flow).
- **SST-Vm:** vorticity-based production.
- **SST-2003m:** the 2003 constants and limiter.

On 2DZP all three give nearly identical results; they differ slightly on the bump. Pick the variant Fluent uses and verify it against the matching TMR files. My inference is that Fluent's default SST, which applies a production limiter and uses the strain rate in the eddy-viscosity limiter, is closest to SST-2003m. That is unverified: confirm it in the Fluent Theory Guide before choosing the reference files. Equation definitions are saved in `nasa_tmr/model_docs/sst.html`.

All TMR grid families are nested: every coarser grid takes every other point of the next finer one. All grids are Plot3D, formatted, double precision, 2D (nbl = 1).

### 3a. 2DZP zero-pressure-gradient flat plate

**Verifies:** SST implementation (wall omega BC, F1/F2 blending, eddy viscosity limiter), Cf grid convergence, and the compressible/incompressible offset.

**Conditions:** M = 0.2 (compressible references); Re = 5e6 per unit length; plate from x = 0 to 2, so Re_x = 1e7 at the end; reference length 2 (CD uses area 2 per unit span); domain x in [-0.333, 2], y in [0, 1]; symmetry boundary ahead of the plate. Fully turbulent: SST "activates" on its own within x < 0.01.

**Grids on disk:** `nasa_tmr/2DZP_flatplate/grids/`: 35x25 (816 cells), 69x49 (3,264), 137x97 (13,056), 273x193 (52,224, gz). The finest TMR level, 545x385 (208,896 cells), is in the nasa.gov zip but was not downloaded; the reference values below include it. Wall y+ is about 1.7 on 35x25 and about 0.1 on 545x385.

**Reference values (SST-Vm, `2DZP_flatplate/sst/`):**

| N cells | Cf(x=0.97) CFL3D | FUN3D | **FUN3D incompressible** | SC/Tetra incompressible | GHOST incompressible | CD CFL3D | CD FUN3D |
|---|---|---|---|---|---|---|---|
| 208,896 | 0.00269085 | 0.00269055 | **0.00271736** | 0.00272119 | 0.00271085 | 0.00285332 | 0.00284417 |
| 52,224 | 0.00268299 | 0.00268166 | **0.00270981** | 0.00271731 | 0.00270016 | 0.00284557 | 0.00282131 |
| 13,056 | 0.00266477 | 0.00265845 | **0.00269119** | 0.00270149 | 0.00267572 | 0.00282597 | 0.00277329 |
| 3,264 | 0.00262625 | 0.00260951 | **0.00265305** | 0.00267663 | 0.00261878 | 0.00278507 | 0.00267868 |
| 816 | 0.00255183 | 0.00251562 | **0.00257930** | 0.00264699 | 0.00250388 | 0.00270623 | 0.00251199 |

TMR's GCI analysis on the finest three grids gives an apparent order for Cf(0.97) of 1.21 (CFL3D) and 1.39 (FUN3D), with a fine-grid GCI of about 0.26-0.28%. For CD the apparent orders are 1.34 and 1.07. Order 2 is not expected here: the leading-edge singularity, together with the 1st-order upwinding these codes use for the turbulence advection terms, keeps the observed order between 1 and 2. Other files give Cf(x) along the whole plate, u+ vs y+ at x = 0.97 and 1.90, and profiles of mu_t, k, omega and F1/F2 at x = 0.97 for profile-level debugging.

**Recommended acceptance (Kutta is incompressible):**
- Cf(0.97) on each grid should fall inside the spread of the three incompressible codes. On 273x193 that band is 0.0027002-0.0027173, roughly +-0.4% about 0.00271. Pass if Kutta is within **+-1%** of the FUN3D-incompressible value on the 137x97 and 273x193 grids.
- Cf(0.97) should converge monotonically, with an observed order between 1 and 2 across 69x49 -> 137x97 -> 273x193. The extrapolated value should be within **0.5%** of the incompressible fine-grid value, about 0.00272.
- Profiles at x = 0.97: u+(y+) on the log law (kappa = 0.41, B = 5.0) within 2%, and the mu_t/mu peak within 3% of `mut_0.97_sstv.dat` after allowing for the 1% offset.
- Why these tolerances: the incompressible codes disagree by about 0.4% on the same grid, and the compressibility offset is about 1%. Matching the compressible numbers to better than about 1.5% would require modelling something Kutta deliberately omits.

### 3b. 2DBC bump-in-channel

**Verifies:** SST with adverse and favourable pressure gradients and wall curvature, with integrated pressure and viscous forces.

**Conditions:** M = 0.2; Re = 3e6 per unit length; viscous wall for x in [0, 1.5]; bump y = 0.05*sin^4(pi*x/0.9 - pi/3) for 0.3 <= x <= 1.2, otherwise y = 0; symmetry wall from x = -25 to 0 and from 1.5 to 26.5; top boundary at y = 5 (symmetry); reference length 1.5.

**Grids on disk:** `nasa_tmr/2DBC_bump/grids/`: 89x41 (3,520 cells), 177x81 (14,080), 353x161 (56,320, gz). The finer levels, 705x321 and 1409x641, are referenced only. Wall y+ is about 0.95 on 89x41.

**Reference values, SSTm (`2DBC_bump/sst/force_convergence_sst.dat`, `cf_convergence_sst.dat`):**

| N cells | CL CFL3D | CL FUN3D | CD CFL3D | CD FUN3D | CDp CFL3D | CDv CFL3D | Cf(0.75) CFL3D | Cf(0.75) FUN3D |
|---|---|---|---|---|---|---|---|---|
| 901,120 | 0.025047 | 0.025080 | 0.0036045 | 0.0035926 | 0.00042098 | 0.0031835 | 0.0058482 | 0.0058538 |
| 225,280 | 0.024974 | 0.025085 | 0.0036027 | 0.0035906 | 0.00042018 | 0.0031825 | 0.0058240 | 0.0058294 |
| 56,320 | 0.024827 | 0.025098 | 0.0036071 | 0.0035734 | 0.00043165 | 0.0031755 | 0.0057689 | 0.0057645 |
| 14,080 | 0.024507 | 0.025119 | 0.0037071 | 0.0036106 | 0.00055434 | 0.0031528 | 0.0056222 | 0.0056126 |
| 3,520 | 0.023519 | 0.025143 | 0.0045619 | 0.0040563 | 0.00147877 | 0.0030831 | 0.0051639 | 0.0052198 |

Cf convergence at x = 0.6321975 and x = 0.8678025 is in the same file; finest-grid CFL3D values are 0.0050623 and 0.0024894. SST-Vm finest CFL3D: CL = 0.025015, CD = 0.0035993. SST-2003m finest CFL3D: CL = 0.025042, CD = 0.0035821. CDp and CL converge slowly on coarse grids, and differently for the two codes, which makes them sensitive diagnostics of how the wall pressure and the bump-edge geometry are handled.

**Recommended acceptance:**
- Cf at the three stations on 353x161: within **+-2%** of the CFL3D/FUN3D mean on the same grid, after first raising the reference by about 1% for incompressibility.
- CL on 353x161: within **3%** of the reference grid value, after dividing the reference by about 1.02 for compressibility. CDv within 2%.
- CDp: grid-converging and below 5e-4 by the 353x161 grid. Don't apply a percentage tolerance, because it is a small difference of large numbers.
- Why: at 56k cells the two reference codes already differ by 1.1% in CL and 0.9% in CD.

### 3c. 2DN00 NACA 0012 validation

**Verifies:** a lifting airfoil with farfield treatment against experiment, and an SST code-to-code comparison on the same grid.

**Conditions:** M = 0.15 (TMR recommends this for compressible codes; Kutta runs incompressible); Re = 6e6 per chord; fully turbulent; alpha = 0, 10, 15 deg. The airfoil is the sharp-TE scaled NACA 0012, y = +-0.594689181*(0.298222773*sqrt(x) - 0.127125232*x - 0.357907906*x^2 + 0.291984971*x^3 - 0.105174606*x^4), with a maximum thickness of 11.894% of chord. The C-grid farfield is about 500c away. If you use a smaller domain, apply the point-vortex farfield correction (Thomas & Salas 1986).

**Grids on disk:** `nasa_tmr/2DN00_naca0012/grids/`: 113x33 (3,584 cells; 65 points on the airfoil), 225x65 (14,336), 449x129 (57,344, gz) and **897x257 (229,376, gz)**. The 897x257 grid is the one every reference result below was computed on, which is why it was added. The 1793x513 level is referenced only. The surface definition is `n0012points_superbig_clust_fix.dat`.

**Reference values, CFD on 897x257 (TMR summary tables; CFL3D files in `sst/`):**

| code / model | CL(10) | CL(15) | CD(0) | CD(10) | CD(15) |
|---|---|---|---|---|---|
| CFL3D SST | 1.0778 | 1.5068 | 0.00809 | 0.01236 | 0.02219 |
| FUN3D SST | 1.0840 | 1.5109 | 0.00808 | 0.01253 | 0.02275 |
| NTS SST | 1.0765 | 1.5100 | 0.00809 | 0.01251 | 0.02187 |
| CFL3D SST-V | 1.0778 | 1.5060 | 0.00828 | 0.01245 | 0.02224 |

**Experiment (`experiment/CLCD_Ladson_expdata.dat`):** Ladson (NASA TM 4074), Re = 6e6, M = 0.15, transition tripped. TMR considers this the most appropriate dataset for fully turbulent CFD forces. With 180 grit: alpha = 10.18 gives CL = 1.0809, CD = 0.01165; alpha = -0.03 gives CD = 0.00803; alpha = 15.24 gives CL = 1.5169, CD = 0.01870. Gregory & O'Reilly Cp (Re = 3e6, tripped) is TMR's preferred pressure dataset. No skin-friction data exists.

**Recommended acceptance:**
- On 897x257: CL(10) within **1.5%** of 1.066 (the CFL3D value with the about 1.1% compressibility correction removed). CD(0) within **5%** of 0.0081. CD(10, 15) within **5%** of the code mean. TMR reports the codes themselves differ by up to 1% in CL and 5% in CD on this grid.
- Against Ladson: CL slope and CL(10) within **2%**. CD(0) within **5%** of 0.0080-0.0081.
- On 449x129: the same quantities should already be within 2x those tolerances, showing convergence toward the 897x257 values.
- TMR's grid-convergence study for this airfoil (the "numerics" page, 2DANA) has **SA results only**, no SST, so no SST GCI table exists for 2DN00.

## 4. Selig S1223 high-lift low-Re airfoil (UIUC LSATs)

**Verifies:** validation against experiment for the kind of profile used on FSAE wings. It also exercises SST at transitional Reynolds numbers, where a fully turbulent model is known to be imperfect. This case is mainly a **Fluent cross-check**: same mesh, same model, compare Kutta with Fluent.

**Conditions:** 2D, incompressible, free transition in the experiment, UIUC 3x4 ft tunnel, chord 12 in. Typical FSAE wing Reynolds numbers are 2e5-6e5 (c = 0.25-0.35 m at 10-25 m/s).

**Data:** `validation/data/uiuc_lsat_s1223/`
- `s1223_vol1_drag_polars.csv`: CL, CD vs alpha, clean, at Re about 1.0e5, 1.5e5, 2.0e5, 2.5e5 and 3.0e5. Also a Gurney flap (1.04% c) case and a VG case at 2e5.
- `s1223_vol1_lift_curves.csv`: lift curves through stall at Re 61k-302k. Clean CLmax is **2.188 at 16.9 deg (Re 1.98e5)**, **2.199 at 15.9 deg (2.52e5)** and **2.212 at 15.6 deg (3.02e5)**.
- `s1223_vol5_lift_moment.csv`: the 2012 re-test with CL and Cm, clean at 80k-250k, plus Gurney flaps of 1.04, 1.56, 3.12 and 4.17% c and trip strips.
- Coordinates: as designed (`raw/s1223.dat` = `lsat_vol1_COORD01_s1223.dat.txt`) and **as tested** (`lsat_vol1_COORD01_s1223.dap.txt`). Use the as-tested coordinates for experimental comparison and the design coordinates for the Fluent cross-check. Coordinates are not GPL-licensed; keep them internal.

Selected clean polar points:

| Re | alpha | CL | CD |
|---|---|---|---|
| 2.00e5 | -0.00 | 1.083 | 0.0191 |
| 2.00e5 | 4.76 | 1.510 | 0.0245 |
| 2.00e5 | 7.93 | 1.753 | 0.0314 |
| 1.99e5 | 13.17 | 2.094 | 0.0413 |
| 2.99e5 | 0.48 | 1.131 | 0.0164 |
| 2.99e5 | 8.21 | 1.800 | 0.0257 |
| 2.99e5 | 14.74 | 2.175 | 0.0371 |

**Published RANS on S1223** (referenced, not downloaded):
- Jha et al., "CFD predictions of the aerodynamic characteristics of S1223 airfoil at low Reynolds' number", AIAA Aviation 2024. SST, SA, SST-gamma and LES at Re = 2e5. RANS shows large unrecovered separation at 14 deg where LES stays attached.
- Boughou et al. (2022-2023) and an IEEE paper (doc 9836229): SST-gamma beats fully turbulent SST on drag by about 8.6% or more.
- Forum reports (CFD Online threads 129467 and 247008) describe fully turbulent SST in Fluent underpredicting CL.

I found no paper with a mesh-independence study and tabulated SST values that would serve as a hard reference.

**Recommended acceptance:**
- **Against Fluent (same mesh, same SST settings, fully turbulent):** CL within **1%** and CD within **3%** for alpha up to 10 deg. This is the gate that matters for the end goal.
- **Against experiment:** CL within **+-0.10 (about 7%)** in the linear range at Re 2e5-3e5, CD within **+-25%**, and CLmax/stall angle as a reported trend only. Why: a laminar separation bubble dominates at these Re. Fully turbulent SST does not model it, so tight agreement would be luck. Adding transition (gamma or gamma-Re_theta) is the route to tighten this later.
- **Fallback multi-element data:** the AIAA HiLiftPW 30P30N was not fetched. The TMR 2D multi-element case (`multielementverif.html`, a CRM-HL cut at M = 0.2, Re = 5e6, alpha = 16) has **SA-neg results only and no experiment**.

## 5. Single-element wing in ground effect (Zerihan & Zhang 2000; Zerihan PhD 2001)

**Verifies:** validation of ground effect (moving ground, ride-height sensitivity, the force-reduction phenomenon near the ground) on a race-car wing profile. This is the most FSAE-relevant validation case available with open geometry.

**Conditions:** Tyrrell 026 front-wing main element at 80% scale. c = 223.4 mm, span 1100 mm (AR 4.92), TE thickness 1.65 mm, endplates 250x100x4 mm. U about 30 m/s, so Re_c = 0.43-0.46e6. Moving ground (rolling road). Free transition. **Reference incidence alpha = 1 deg is a true chord incidence of 3.45 deg**: quoted angles are measured from a line 2.45 deg off the chord. Ride height h is the distance from the ground to the lowest point of the suction surface, measured at zero incidence. Freestream runs were at h/c = 3.36.

**Data:** `validation/data/zerihan_tyrrell026_ground_effect/`
- `tyrrell026_single_element_coords_ref_incidence_1deg.csv`: thesis Table 2, 44 points per surface, already rotated to the 3.45 deg true incidence (verified from the TE height).
- `fig7a_3d_CL_vs_hc_alpha1deg_digitized.csv`: 3D downforce coefficient vs h/c at alpha_ref = 1 deg, 16 points, digitized automatically at 300 dpi, CL +-0.005. Freestream CL = 0.69. Highlights: **0.835 at h/c = 0.672, 0.945 at 0.448, 1.084 at 0.313, 1.270 at 0.224, 1.420 at 0.179, 1.615 at 0.133, a maximum of 1.722 at 0.082, then 1.610 at 0.043.**
- `fig7b_3d_CD_vs_hc_alpha1deg_digitized.csv`: 3D CD vs h/c, from 0.0280 at h/c = 0.672 to 0.0893 at 0.044.
- `table15_2d_sectional_Cl_vs_hc.csv`: **sectional** Cl at the wing centre, from integrated surface pressures, against the thesis's own 2D CFL3D SA and k-omega results (20-31k-cell grids). Experimental Cl: 0.766 (freestream), 0.902 (0.671), 1.009 (0.448), 1.145 (0.313), 1.286 (0.224), 1.362 (0.179), 1.385 (0.134), 1.371 (0.090). The 2D k-omega results ran 14% high at large heights and matched within 0-5% below h/c = 0.22. The thesis calls the pressure-integrated experimental Cl "only an estimate" because 45 taps miss the LE suction spike.
- `text_reference_values.csv`: freestream CLmax = 1.35 at alpha 11.3 deg; lift slope 4.57/rad; CLmax in ground effect = 2.26 at h/c = 0.134, alpha 12.3 deg.

**Later RANS work on this geometry** (referenced): Mahon & Zhang (2005, 2006, J. Fluids Eng.) ran 2D Fluent RANS of this profile in ground effect, comparing turbulence models, transition and wake. Zhang & Zerihan cover the double-element version (J. Aircraft 2003; open AIAA paper at eprints.soton.ac.uk/22605).

**Recommended acceptance:**
- **3D, h/c >= 0.18 (attached, force-enhancement region):** CL within **5%** of the digitized curve and CD within **10%**. Why: the flow is attached, so the main errors are digitizing (under 0.5%), free versus fully turbulent transition, and tunnel and belt effects. The thesis 2D k-omega results fall within 5% in this range.
- **3D, h/c < 0.13 (force-reduction region with TE separation):** reproduce the trend, i.e. a CL maximum between h/c = 0.07 and 0.12 followed by a drop. Accept a CLmax within 10%. Steady RANS separation prediction is the limit here.
- **2D sectional Cl:** compare trend and code-to-code. Match the thesis CFL3D k-omega within **5%** at the same h/c. Do not hold Kutta to the pressure-integrated experiment, which the thesis itself says underestimates Cl.

## 6. Ahmed body, 25 and 35 deg slant

**Verifies:** 3D bluff-body separation, base pressure and drag. It is a known weak spot of steady RANS, which makes it a strong **Fluent cross-check** and only a cautious experimental validation.

**Geometry** (Ahmed et al. SAE 840300; ERCOFTAC): L = 1.044 m, W = 0.389 m, H = 0.288 m, 100 mm front radii, slant length 0.222 m, sharp rear edges, 50 mm ground clearance on four stilts of 30 mm diameter. Frontal area 0.112 m^2 without stilts. No CAD is hosted; build it from the dimensions.

**Experiments:**
- **Ahmed 1984 (DLR):** 60 m/s, Re_L = 4.29e6, balance drag. Cd at **25 deg = 0.285**: this is the value cited by CFD Support and by Hachimy et al. 2022. SimScale reads 0.2875 off Ahmed's Fig. 4. Cd at **35 deg is about 0.26**: Liu & Moser (ETH) quote 0.260; 0.257 is also widely cited. Component breakdown, front/base/slant/friction:
  - 25 deg (as cited by Hachimy 2022): Ck = 0.0222, Cb = 0.0726, Cs = 0.136, Cr = 0.052.
  - 35 deg (as cited by Liu & Moser): Ck = 0.020, Cb = 0.095, Cs = 0.090.
  - **I could not check these against SAE 840300 itself (paywalled).** Treat them as secondary-sourced: the totals are consistent across sources to about 1%, but the component values are less certain.
- **Lienhart & Becker 2003 (LSTM, MOVA):** 40 m/s, Re_L = 2.78e6, Tu = 0.25%. LDA planes, hot-wire BL profiles over the slant, rear Cp and the inlet profile. **No drag.** On disk: `ercoftac_ahmed_body/exp_lienhart_becker/` (52 files).

**CFD references:**
- `ercoftac_ahmed_body/cfd11_cfx5_sst/` (53 files): CFX-5 SST from Menter's group (now Ansys). 2.5M cells, half model, y+ about 1, 2nd order, same stations as the experiment.
- **ERCOFTAC evaluation:** at 35 deg, all models reproduce the massively separated wake and the right Cd trend, but not the exact level. At 25 deg, "the great majority" of steady RANS runs miss the reattaching 3D wake and predict massive separation instead. Only LES and a low-Re RSM got it right. Menter et al. (2024) add that even a 560M-cell LES stayed attached and only octree meshes reproduced separation and reattachment.
- **Fluent SST, 25 deg (Hachimy et al., CFD Letters 14(5), 2022):** 40 m/s, half model, y+ about 50 (wall functions, not wall-resolved). Meshes of **0.38M / 0.57M / 1.19M cells gave Cd = 0.295 / 0.289 / 0.282**. Components on the fine mesh: Cr 0.0374, Ck 0.030, Cb 0.082, Cs 0.0984. The total matches within 1%, but the components are off by 13-35%.
- **OpenFOAM SST, 25 deg (SimScale validation):** 60 m/s, y+ about 1. Meshes of **3.75M / 6.92M / 10.1M cells gave Cd = 0.2985 / 0.2915 / 0.2835** against 0.2875.
- **Transient RANS, 35 deg (Liu & Moser, ETH):** SST gave Cd = 0.241 (-7.3%), k-eps 0.242, k-eps-v2 0.264 and RSM 0.282, against 0.260.

**Recommended acceptance:**
- **35 deg against experiment:** Cd within **+-8%** of 0.26 (the observed RANS spread is -7% to +9%). Rear Cp and the LDA planes at x = 0, 80 and 200 mm should be qualitatively right: a separated slant and the C-pillar vortex pair.
- **25 deg:** use it **only against Fluent**, with identical mesh, y+ strategy and SST variant. Require Cd within **3%** and the same wake topology (attached/reattaching vs separated slant). Do not gate on the experimental 0.285. An SST result can land near 0.285 for the wrong reason, as the component errors in Hachimy show. Report the topology explicitly.
- Run a mesh study with at least three levels on every gate. The cited studies show Cd moving 3-5% from coarse to fine.

## 7. FSAE car stand-ins

Nothing suitable was found that is both openly licensed and has published, reproducible forces:

- **github.com/TheMahaanOne/fsae-aero-2026** (Team Fateh Racing, 2026). This has a full car, front, rear and side wings, and the bluff body as **STL and Parasolid** (`cad/`, about 9 MB). Forces were computed with STAR-CCM+ and Fluent using steady incompressible k-omega SST, a trimmed-hex mesh, y+ 1-3, a moving ground and rotating wheels at 10 m/s. Published results: **downforce 397.3 N, drag 110.9 N, L/D 3.58**; rear wing 160 N, front wing plus whiskers 90 N, side wings 70 N; full-car downforce 270-351 N across ride heights. **The repository has no license, so all rights are reserved.** Nothing was copied into `validation/`. Ask the author for permission (or a license) before using it. Even then, it is single-team CFD rather than measurement, so use it as a workflow test, not truth.
- **ZHAW bachelor thesis, "Aerodynamic Package of a Formula Student Racing Car"** (digitalcollection.zhaw.ch): wind-tunnel scale-model forces at 17-23 m/s against CFD. This is a possible experimental stand-in, but I did not check whether it includes geometry.
- **Gupta & Saxena, "Aerodynamics analysis of a Formula SAE car"** (repository.up.ac.za, 2017) and similar theses report Fluent full-car CL and CD without releasing geometry.
- **Larger automotive datasets** under CC-BY-SA: DrivAerML (500 DrivAer variants, high-fidelity CFD) and WindsorML (355 Windsor-body variants, WMLES with more than 280M cells). These are road-car rather than FSAE and very large, so I did not download them. The Windsor body with force tables is the most practical next step after Ahmed.

**Recommendation:** treat cases 5 and 6 as the FSAE stand-ins until a team case arrives. The ground-effect wing covers the dominant downforce physics and Ahmed covers bluff-body drag. For a full-car cross-check, get the team's geometry and their Fluent case (mesh plus setup), then compare on identical meshes first.

## Not obtained

- Ghia 1982 and Ahmed 1984 originals (both paywalled). Their values come from cross-checked secondary transcriptions, documented above.
- Zerihan & Zhang 2000 J. Aircraft PDF (paywalled). The same data was taken from the open thesis.
- TMR finest grids (2DZP 545x385, 2DBC 705x321 and 1409x641, 2DN00 1793x513), TMR Tecplot contour zips and the HiLiftPW 30P30N grids and data. All are referenced but not needed for the first gates.
- An SST GCI study for NACA 0012 (TMR only has one for SA).
- S1223 experimental data at Re = 5e5 (UIUC tops out at 3e5) and any S1223 SST mesh-independence study with tables.
- Any openly licensed FSAE geometry with published forces.
