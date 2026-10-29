# Ghia, Ghia & Shin (1982): lid-driven cavity centerline velocities

- **Reference:** U. Ghia, K. N. Ghia, C. T. Shin, "High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method", J. Comput. Phys. 48(3):387-411, 1982. doi:10.1016/0021-9991(82)90058-4. Tables I (u on the vertical centerline) and II (v on the horizontal centerline).
- **Fetched:** 2026-10-29
- **Why there are several sources:** I found no machine-readable copy from the authors or publisher, and the paper is paywalled. Instead I downloaded six independent public transcriptions and built the CSVs by per-point majority vote for Re = 100, 400 and 1000.

| raw file | source URL (pinned) | license |
|---|---|---|
| `raw/gist_ivan-pi_ghiau.txt` (u only) | https://gist.githubusercontent.com/ivan-pi/3e9326d18a366ffe6a8e5bfda6353219/raw/1b746e02ecf16771b476ee82ebaf3a7e802ab9dc/ghiau.txt | gist, no license stated |
| `raw/pysph_ghia_cavity_data.py.txt` (u, v; a Python module saved as text, never executed) | https://raw.githubusercontent.com/pypr/pysph/1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5/pysph/examples/ghia_cavity_data.py | PySPH (BSD-style; GitHub reports NOASSERTION) |
| `raw/petibm_ghia_et_al_1982_lid_driven_cavity.dat` (u, v; Re 100/1000/3200/5000/10000) | https://raw.githubusercontent.com/barbagroup/PetIBM/be91f218702ada7fe569fe77c258fa748f075713/examples/data/ghia_et_al_1982_lid_driven_cavity.dat | BSD-3-Clause |
| `raw/shamsbasir_ghiau.txt`, `raw/shamsbasir_ghiav.txt` | https://raw.githubusercontent.com/shamsbasir/investigating_mitigating_failure_modes_in_pinns/266fcf5f613c52845d880edd4e70ca021825896b/Lid_driven_cavity/ghia{u,v}.txt | MIT |
| `raw/2dLBM_verticalGhia.dau`, `raw/2dLBM_horizGhia.dau` (Re=1000 only, with 129-grid index) | https://raw.githubusercontent.com/enstrophy/2dLBM/500979cd382d3e353ab03592f5360e8748635267/{verticalGhia,horizGhia}.dau | GPL-2.0 |
| `raw/psiboil_Ghia-u.dat`, `raw/psiboil_Ghia-v.dat` (grid index + all 7 Re) | https://raw.githubusercontent.com/Niceno/PSI-BOIL/ddc9e283d75429c43f496ca90d87dbc529b38289/Benchmarks/Lid/Ghia-{u,v}.dat | no license stated |

The tabulated *values* are factual data from a 1982 journal paper. The raw transcription files carry their repositories' licenses.

## Derived files

- `ghia1982_u_vertical_centerline.csv`: columns `grid_index_129, y_table, y_exact, u_Re100, u_Re400, u_Re1000`
- `ghia1982_v_horizontal_centerline.csv`: columns `grid_index_129, x_table, x_exact, v_Re100, v_Re400, v_Re1000`

`*_table` is the 4-decimal coordinate printed by Ghia. `*_exact` is `(grid_index-1)/128` on Ghia's uniform 129x129 grid; interpolate to this one. Velocities are scaled by the lid speed, and the lid moves in +x at y=1.

**Majority-vote outcome.** All sources agree at every u point. All v points agree except these three:

| point | majority (used) | dissent |
|---|---|---|
| v, Re=400, x=0.5 | 0.05186 (shamsbasir, psiboil + many other public copies) | PySPH 0.05188 |
| v, Re=1000, x=0.5 | 0.02526 (4/5; independently backed by Erturk 2005 601^2 value 0.0258) | PetIBM 0.02426 |
| v, Re=1000, x=0.9063 | -0.51550 (4/5) | shamsbasir -0.51500 |

**Suspect point kept as printed:** v, Re=400, x=0.9063 = -0.23827. It breaks the smooth profile between -0.22847 at x=0.9453 and -0.44993 at x=0.8594, and others have flagged it on CFD Online (thread 88848). It is probably a typo in the original table. Exclude it from pass/fail or give it a loose tolerance. The well-known u typos at Re=3200 (y=0.4531) and Re=10000 (y=0.5) fall outside the Re range extracted here.

## Files
| file | bytes | sha256 |
|---|---:|---|
| `ghia1982_u_vertical_centerline.csv` | 832 | `73698b6a0fbd59601421e858c2a6a6cbf6ca0d70ea8fd09549674c711472c7aa` |
| `ghia1982_v_horizontal_centerline.csv` | 831 | `4fc1d2726b447c25c0aa6a75cb2f837f5772d08f1fc404b30ffb9ccafcf8cd26` |
| `raw/2dLBM_horizGhia.dau` | 261 | `6780b68736ab360d637095080ced967c20701d3eb0e3ea3425fa4cfa86ecd4bd` |
| `raw/2dLBM_verticalGhia.dau` | 257 | `39835a7d874695ee82e9389c9f15dc1eeae0b8693913ab8503862a799784c77d` |
| `raw/gist_ivan-pi_ghiau.txt` | 1816 | `51b65a487eb5b424b4aecdd11ebfa612e3dc646a04217bdd8b14c20cf72f3b54` |
| `raw/petibm_ghia_et_al_1982_lid_driven_cavity.dat` | 2017 | `d575b8c4af5839dec1390b7037d4d2e7bbe42b9add8540aa128e1a31b3ef0298` |
| `raw/psiboil_Ghia-u.dat` | 1645 | `a216a0baf956befdc8ee108dca89a4a87a485989446c8400eb52c0cde896605b` |
| `raw/psiboil_Ghia-v.dat` | 1654 | `c100b4211f6bec3a8beb5cffc238c24bab70234984f74705a61599e18d382aee` |
| `raw/pysph_ghia_cavity_data.py.txt` | 3655 | `466d2a28b74304204481a7d8a709864c9c68c4c164907332a404daf98bf571a7` |
| `raw/shamsbasir_ghiau.txt` | 1189 | `56c76d7e0243ada526f9b1a2175137040f37d6ad2b93c1b671787cd438c61700` |
| `raw/shamsbasir_ghiav.txt` | 1190 | `33bc9e92605fbe03911fb36d13bd0c29e6aea19b7f36f2cb114a63c14993a1b0` |
