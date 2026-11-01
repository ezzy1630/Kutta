# 2D laminar cylinder: a converged legacy Fluent reference

Branch `t3/dbfd9a6f`, based on main `ccd9ed4`. Apple M4 Pro, CPU reference backend and Metal; CUDA was not run (CI only). The Fluent corpus is read through `KUTTA_FLUENT_EXAMPLES` and never committed.

## Result

The Ansys adjoint-tutorial cylinder (Fluent 18.0, 2D, laminar Re 40) now runs through `kutta solve` with one named approximation, `reverse-flow`. Fluent saved this solution at iteration 66 with scaled continuity 9.32e-6; the case's own criterion is 1e-5. Kutta's cold start, 1000 iterations, against Fluent's saved fields:

| | CPU | Metal |
| --- | ---: | ---: |
| Velocity relative L2 | 0.0645% | 0.0645% |
| Pressure relative L2 | 0.517% | 0.517% |
| Inlet − outlet pressure drop | 0.509% | 0.509% |
| Wall force (pressure + shear) | 0.0705% | 0.0705% |
| Drag report `cd-1`: Kutta / Fluent | 1270.86 / 1271.74 (0.069%) | same |
| Lift report `cl-1`: Kutta / Fluent | 0.332 / 0.488 (\|ΔCl\|/Cd = 1.23e-4) | same |
| Pressure-outlet faces with inflow: Kutta / Fluent | 0 / 0 of 30 | same |

CPU and Metal agree to the last bit: the two `result.vtu` files are byte-identical (SHA-256 `3be71422…`). Logs, residuals, forces and comparisons are in [`runs/cpu`](runs/cpu) and [`runs/metal`](runs/metal).

These runs used this branch's working tree before commit. Their `summary.json` therefore records the base commit `ccd9ed4` as `kutta_commit`. With reference values ρ = 2, U = 1 and A = 1, the case's coefficients equal the force in newtons. The physical drag coefficient is 1271 / (½ · 1 · 40² · 1) = 1.59.

This pair is now a slow-tier gate, at about 1.5x these values. The 3D cases are unchanged: V6, `elbow_param` and the nozzle stay bitwise identical to main on CPU and Metal (see "3D parity").

To get here:

1. **2D cases** are admitted as a one-cell extrusion with symmetry caps.
2. **Fluent's Green-Gauss cell-based gradient** was added, and is selected only by cases that store it.
3. **A legacy (non-HDF5) `.dat` reader** was added. The legacy mesh reader now skips cell-partition sections and rejects nonconformal meshes by name.
4. **Pressure-outlet keys:** Fluent 18.0 does not write `prevent-reverse-flow?` or `p-profile-multiplier`. The first is the `reverse-flow` approximation. The second is a convergence note: it multiplies a constant zero gauge pressure, so it has no effect.

**The corpus data file is damaged upstream.** It lost 15 bytes to a text-mode transfer (see "The data file"). Kutta rejects it with a message naming the cause. The comparisons use a checked restoration. The owner should add the restored file to the corpus.

## Implementation and why it is faithful

### 2D representation

Fluent solves a 2D case on a plane, per unit `reference-depth`. Kutta reads the 2D mesh as one cell through the thickness (the existing `extrude_2d`, thickness 1), with two zone-less caps, `front` and `back`. Setup admits a 2D case only if all of the following hold (`planar_extrusion` in `src/kutta/setup/fluent.cpp`):

- exactly those two caps exist;
- each has one face per cell;
- neither has a Fluent zone of the same name;
- the mesh's z extent equals the case's `reference-depth` to 1e-12.

The caps become symmetry planes. Why this is the planar operator, not an approximation:

- **Convection and diffusion.** On a symmetry cap the normal velocity is zero, and so is the normal gradient of every other quantity. The caps therefore carry no convective flux and no tangential stress.
- **Gradients.** A cap's face value equals its cell's value. The cap pair's area vectors are equal and opposite, so the caps add nothing to a Green-Gauss sum. The side faces' area divided by cell volume is the edge length divided by cell area: the 2D quantity.
- **Forces.** Areas carry the depth, so a force with depth 1 is Fluent's force per unit depth. That is why the depth must equal `reference-depth`.

`tests/fv/test_green_gauss.cpp` checks the gradient claim on a mixed quad/tri extrusion of thickness 0.37. The 3D Green-Gauss gradient equals the planar evaluation to float precision, and its z component stays below 1e-12 relative.

The solved w stays at float noise: a residual of about 1e-11, and a wall force z of 4e-8 N against 1271 N in x. In a 2D case, setup requires Fluent's data to have no `SV_W`, and compares w against zero.

### Green-Gauss cell-based gradient

The case selects the Green-Gauss cell-based gradient (`gradient-option 0`). Kutta implements Theory Guide 2024 R2:

- §23.3.3.1 Eq. 23-30: ∇φ = (1/V) Σ φ_f A_f;
- §23.3.3.2 Eq. 23-31, the face value as the arithmetic mean of the two cell values: φ_f = (φ_c0 + φ_c1)/2.

The face mean makes the gradient linear in the cell values. Since Σ A_f = 0 over a closed cell, Eq. 23-30 can be written as differences: (1/V) Σ (A_f/2)(φ_nb − φ_c) over internal faces plus (1/V) Σ A_f(φ_b − φ_c) over boundary faces.

That is the same row-times-difference form as Kutta's least-squares fit. `green_gauss_cell` (`src/kutta/fv/operators.cpp`) therefore returns a fit with:

- rows A_f/(2V) on internal faces;
- rows A_f/V on boundary faces;
- an identity inverse.

The existing `fv_gradient` and `fv_vector_gradient` kernels evaluate it unchanged on CPU, CUDA and Metal, so no new kernel was needed.

The test checks this against a literal evaluation of Eq. 23-30/23-31 on three meshes: a hex box, a tet box and the mixed planar extrusion. Errors are below 2e-5 of the gradient scale (float rows and values).

How it is selected:

- `SolverSettings::gradient` defaults to `LeastSquares`.
- Setup selects Green-Gauss only when the case stores it. Every LSQ case keeps LSQ: the quick test checks this for all six corpus pairs, and the parity runs below are bitwise identical.
- Laminar momentum and pressure use it, and so do realizable and standard k-epsilon transport.
- Green-Gauss with SST or SA is rejected as unsupported.

Kutta's least-squares path adds a quadratic projected face reconstruction. With Green-Gauss the face value uses only the gradient, Fluent's documented linear reconstruction, so the projection is not built.

### Pressure-outlet keys

Fluent 17.2 and 18.0 cases lack `prevent-reverse-flow?` and `p-profile-multiplier`. Cases from 21.2 onward store them, as `#f` and `1.`. Kutta treats them as follows:

<a id="reverse-flow"></a>
- **`prevent-reverse-flow?` absent → `--approximate reverse-flow`.** When this key is absent, Fluent's default and its old behaviour are undocumented. The global `reversed-flow/pressure-bc? #f`, present in every case, is also undocumented. Kutta never suppresses reverse flow, and an 18.0-or-earlier case without the key now needs the opt-in.
  - The effect can only appear on outlet faces where flow enters the domain. The cylinder outlet has 30 faces. None has inflow in Fluent's saved converged `SV_FLUX`, or in Kutta's converged CPU and Metal solutions (the `backflow` count in the slow gate). The measured effect on this case is therefore zero.
  - Other cases need their own check. A case that stores the key keeps the existing rule: it must be `#f`.
- **`p-profile-multiplier` absent.** When the gauge pressure is the constant 0, any multiplier gives 0, so setup records a convergence note. Otherwise the multiplier must be 1, as before.

### Remaining settings on this case

- **Pressure:** standard pressure interpolation (already supported).
- **Momentum:** second-order upwind.
- **Coupling:** Fluent ran Coupled without pseudo time. Kutta runs SIMPLEC under-relaxation, with the existing convergence-path note. This changes the iteration path, not the fixed point.
- **Flux:** the case predates `pbs/flux-index`, so the momentum-weighted Rhie-Chow flux (Fluent 12.0 TG 18.4.1) applies, with the existing note.
- **Initialization:** the case's hybrid initialization.

### Legacy `.dat` reader

`src/kutta/case/fluent_data.cpp` (`LegacyData`) reads Fluent's text and binary section format:

- **Sections read.**
  - `(4 (60 …))` gives the byte order: 60 = little-endian.
  - 33 gives the grid size, checked against the case. A periodic-shadow zone (type 8/1008) counts in the case but not the data.
  - 37 gives the variables and the iteration count.
  - 300/2300/3300 are fields: header `(id zone size ntime nphase first last)`, values row-major.
  - 302/2302/3302 are residual histories: rows of (iteration, raw, scale).
- **Field names** come only from the preceding `(0 "SV_P, domain 1, …")` comment.
- **Other sections** (for example 3314 integer data, and time levels, phases or domains other than the mixture) are skipped. Each skip is recorded as an informational issue.
- **Body closing.** Field bodies close after the trailer (`…\nEnd of Binary Section 3300))`), residual bodies before it (`…)\n\nEnd of Binary Section 3302)`). Every field and residual body in the cylinder (restored), the three Borealis files and VKI has exactly one of these two separators. The reader matches them exactly, with no whitespace skipping, so a body one byte short cannot pass by consuming the separator's LF.
- **Residual equations.** Equation IDs are named in the order of the case's `residuals/settings` when the counts agree, and the mapping is logged as informational. The legacy file has no columns 2/3, so they are NaN.
- **Checks.** Values must be finite and fall inside a known zone. Each field must have the right component count and cover every fluid cell. One ID must not carry two names, and chunks must not overlap.

The shared section scanner moved to `src/kutta/case/legacy.cpp`, with no behaviour change for case files. `read_data` dispatches by extension: `.h5` goes to CFF, anything else to legacy.

The legacy mesh reader (`src/kutta/io/fluent_msh.cpp`) now:

- skips binary section 40 (cell partition IDs);
- rejects face zones with type codes of 1000 or more (nonconformal parents and shadows) by name. Before, it failed with "invalid face node count".

## The data file

`adjoint_cylinder_tutorial/cylinder_tutorial.dat` in the corpus came from ansys/example-data. That repository's `.gitattributes` is `* text=auto`, and the file arrived in commit `f09067df48` (#395, 2025-09-11). The evidence that a text-mode transfer stripped the CR from every CR LF pair:

- 4,390 lone CR bytes survive, but the file has no CR LF pair at all.
- Five binary bodies are shorter than their headers, by 15 bytes in total.
- The `.cas`, the VKI pair (gzipped upstream) and the Borealis files (Dataverse) are intact.

| Section | Field | Bytes lost |
| --- | --- | ---: |
| `3300 (1 11 …)` | `SV_P`, cells | 3 |
| `3300 (1 8 …)` | `SV_P`, interior faces | 1 |
| `3300 (111 11 …)` | `SV_U`, cells | 1 |
| `3300 (112 11 …)` | `SV_V`, cells | 4 |
| `3300 (18 8 …)` | `SV_FLUX`, interior faces | 6 |

Kutta's reader refuses the file: "binary data body length disagrees with its header in section 3300; the file has no CR LF byte pair, consistent with a text-mode transfer that stripped them". It never repairs input.

`bench/fluent_restore_cr.hpp` is a separate, explicit forensic tool, used by `bench_restore_stripped_cr`, the legacy-data test and the slow gate:

- **Method.** Within each short body it puts back exactly as many CRs as are missing. Each goes before the unique LF whose restoration makes the following values plausible again (zero, or 1e-30 ≤ |v| ≤ 1e7). Misaligned doubles decode to absurd magnitudes, so the choice is sharp.
- **Uniqueness.** Each of the 15 insertions had exactly one candidate LF. The tool fails on zero or several.
- **Validation.** The restored file then reads with every body length matching its header.
- **Hashes and size.** Input SHA-256 `a00a56c0…f682bbf2d1`, 2,292,922 bytes. Output `cb7dce25827ca4aea45f6fc4bd9936d4ff55fa84ee14d361cec65c7d90ec57cb`, 2,292,937 bytes.

[`runs/restore.log`](runs/restore.log) lists every restored byte.

Independent checks on the restored file:

- Python `struct` decodes match Kutta's reader, including restored values: `SV_V` cell 4044 = 13.970779809010404.
- On the inlet, symmetry planes and wall, Fluent's face pressure equals the adjacent cell pressure exactly (difference 0.0), as standard pressure interpolation requires.
- The comparison errors above are small and smooth. A wrongly placed byte would scramble every following value of its field.

**Owner action (recommended):** add the restored file to the corpus manifest as a derived transform of the upstream file, with both hashes. The tests restore it on the fly, so nothing in the repository depends on that.

## Convergence and iteration count

| Metal iterations | Velocity L2 | Pressure drop | Wall force |
| ---: | ---: | ---: | ---: |
| 600 | 0.0675% | 0.37% | 0.073% |
| 1000 | 0.0645% | 0.509% | 0.0705% |
| 1500 | identical to 1000 at print precision | | |
| 3000 | 0.0645% | 0.506% | 0.0705% |

Kutta's float plateau sits at scaled continuity 4.6e-4 and x-momentum 2.6e-8. That is the float32 velocity-rounding floor, not a stall in the solution: the fields stop changing by 1000 iterations, and one random ulp on the velocities moves the continuity residual by as much as the plateau itself. Fluent's reference run was double precision. The 1e-5 criterion is not reachable in float32 under the documented scaling; see [the residual floor report](../residual_floor/REPORT.md). With the float-pair state, now the default for SIMPLE and SIMPLEC, the cylinder crosses 1e-5 at iteration 1700 and plateaus at 1.5e-6 to 1.9e-6 scaled ([float-pair state report](../paired_state/REPORT.md)). The tables above are main's float32 runs. With the pair, the 1000-iteration errors are 0.0645% (velocity), 0.517% (pressure), 0.506% (pressure drop) and 0.0706% (wall force), and \|ΔCl\|/Cd is 1.21e-4. The final 1000-iteration solves took 89 s on CPU and 21 s on Metal: about 0.089 and 0.021 s per iteration. The machine and GPU were shared, so these timings are noisy.

## Ablations

These were Metal runs, 3000 iterations. Each variant was selected with a temporary local switch that is not in the committed code.

| Variant | Velocity L2 | Pressure L2 | Pressure drop | Wall force | Shear force | Cd error | Cl |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **Green-Gauss, linear (committed)** | 0.0645% | 0.517% | 0.506% | 0.0705% | 1.07% | 0.069% | 0.339 |
| Least squares | 0.0701% | 0.535% | 0.012% | 0.129% | 1.37% | 0.129% | 0.445 |
| Green-Gauss with projected reconstruction | 0.0653% | 0.513% | 0.466% | 0.159% | 1.39% | 0.157% | 0.155 |

Fluent's Cl is 0.488. The committed configuration is the case's documented selection. It gives the best velocity, wall force, shear and drag. The one exception is the inlet pressure, where least squares is better: 50.955 Pa against Fluent's 50.949, where Green-Gauss gives 51.21.

The case also stores the undocumented `recon/cpld-use-lsq-deriv? #t`, as do `elbow_param` and VKI. One possible explanation is that Fluent's Coupled solver uses least-squares derivatives somewhere. That is unverified, and Kutta does not interpret the flag.

## 3D parity

Main `ccd9ed4` was built from `git archive` in a temporary directory, outside the main checkout. Both builds then ran `kutta solve … --iterations 40 --out DIR` on each case and backend:

- V6 `mixing_elbow.cas.h5 --approximate wfgc`;
- `elbow_param.cas.h5 --approximate compressibility,ke-options`;
- nozzle `nozzle_para_named_1.cas.h5 --approximate wfgc`.

| Case | CPU | Metal |
| --- | --- | --- |
| V6 | identical | identical |
| `elbow_param` | identical | identical |
| nozzle | identical | identical |

"Identical" means `fluxes.csv`, `forces.csv`, `loads.json` and `result.vtu` are byte-identical, and `residuals.csv` matches apart from its timing column. See [`runs/parity.log`](runs/parity.log).

## Tests

- **New quick tests:**
  - `tests/case/test_legacy_data.cpp`: synthetic ASCII and binary fields, residual rows, grid size, skipped sections and every rejection, including a one-byte loss in a field body and in a residual body. Corpus checks cover cylinder rejection, restoration and decoded values, plus reading Borealis and VKI.
  - `tests/fv/test_green_gauss.cpp`: matches the definition on three meshes, is exact for a linear field on a uniform hex box, and equals the planar gradient on an extrusion.
- **Extended quick test:**
  - "Fluent corpus pairs run only with their named approximations and map their numerics" now includes the cylinder. It checks the laminar mapping, the symmetry caps and the multiplier note. It also checks that the cylinder alone uses Green-Gauss and every other pair keeps least squares.
- **New slow gate:**
  - "Fluent laminar 2D cylinder matches converged Fluent fields" restores the data, verifies its SHA-256 and cold-starts 1000 iterations on CPU and Metal. Bounds: velocity < 9.7e-4, pressure < 7.8e-3, pressure drop < 7.6e-3, wall force < 1.1e-3, Cd < 1.05e-3, \|ΔCl\|/Cd < 1.85e-4, and no outlet inflow. It skips with a message when the corpus is absent. Result: **207/207 assertions pass on CPU and Metal**, 107 s ([`runs/slow-cylinder.xml.gz`](runs/slow-cylinder.xml.gz)).
- **Quick tier:** 246/246 test cases and 7,415,453 assertions pass, with fixtures present ([`runs/quick.log`](runs/quick.log)). The run reports 26 skipped, and those are the slow suite.
  - The first full run failed one assertion. "Case support is an allowlist for authoritative physics flags" pinned the old rule that a case with `rp-3d? #f` is rejected by the case reader.
  - That decision now belongs to setup, so the test now checks that `rp-3d? #f` is interpreted as dimension 2. The cylinder quick test checks that setup rejects a 2D mesh whose depth does not match `reference-depth`.
- **Existing slow gate:** "Fluent standard k-epsilon elbow and SST nozzle match converged Fluent fields" exercises the shared comparison code changed for laminar cases. It passes unchanged: 17,044/17,044 assertions on CPU and Metal ([`runs/slow-existing.xml.gz`](runs/slow-existing.xml.gz)).
- **App and lint:** the app builds with `-DKUTTA_BUILD_APP=ON`, and `ctest -L app` passes 10/10. `tools/lint_portability.py` passes.
- No existing gate, budget or tolerance changed.

## Borealis and VKI

**Borealis** (three CC0 tetrahedral fracture states):

- The mesh now reads: 511,466 cells; section 2040 was the blocker. All three data files read; the 26 MPa state has 8 cell fields, 13 face fields and 541 residual samples.
- Setup still rejects the cases, and these are physics, not reading:
  - gravity;
  - the Green-Gauss node-based gradient;
  - `swirl-model?`;
  - an absent `rough-bc`;
  - `bc-source?`;
  - a mesh zone of type 8;
  - `shear-bc`;
  - outlet zone 7's absent `prevent-reverse-flow?` and `p-profile-multiplier`, which have the same version caveat as the cylinder;
  - a non-boolean `hyb-init?`.

**VKI turbine:**

- The data read: 14 cell and 17 face fields, 7 residuals, last iteration 1493.
- The mesh uses nonconformal periodic interfaces: type codes 1012 and 1008, with zero-node parent faces in sections 3061/3067/3018. Kutta's reader now rejects these by name. Supporting them needs nonconformal interface intersection, which is out of scope.

## Open limits

- **Inlet pressure** is 0.5% high (51.21 Pa against 50.95), which accounts for most of the pressure L2 and the pressure-drop error. Least squares removes it but worsens forces. `cpld-use-lsq-deriv?` is a candidate explanation, not a verified one.
- **Lift** is near zero: Kutta gives 0.33 N against Fluent's 0.49 N, on a drag of 1271 N. A symmetric flow would have none, so the nonzero lift presumably comes from the unstructured mesh's asymmetry; that is not verified. Lift is gated relative to drag.
- **Green-Gauss with SST or SA** is rejected, not implemented. No corpus case needs it.
- **Unverified Fluent-default questions.** The default of an absent `prevent-reverse-flow?` and the meaning of `reversed-flow/pressure-bc?` are undocumented. Here they are harmless because the cylinder has no outlet inflow.
- **The corpus file** remains damaged until the owner replaces it; see "Owner action" above.
- **CUDA** compiles in CI only and was not run here.
