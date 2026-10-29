# Erturk, Corke & Gokcol (2005): fine-grid lid-driven cavity, Re=1000 subset

- **Reference:** E. Erturk, T. C. Corke, C. Gokcol, "Numerical solutions of 2-D steady incompressible driven cavity flow at high Reynolds numbers", Int. J. Numer. Meth. Fluids 48:747-774, 2005. Preprint arXiv:cs/0411047.
- **Source used:** https://arxiv.org/pdf/cs/0411047, fetched 2026-10-29. The PDF is **not stored** here (arXiv license "assumed-1991-2003", no redistribution grant). Read copy sha256: `293b42bd5149a01b0471e44c5f1aa44772a7c03a9cf970e7a97d03c71f83911d`.
- **Extraction:** text extracted with pypdf, then the Re=1000 column of Table 6 (u, vertical centerline, 23 points) and Table 7 (v, horizontal centerline, 23 points) was parsed. Both tables come from the 601x601 grid, which is second-order and very close to grid-converged. `erturk2005_primary_vortex_Re1000.csv` transcribes the Re=1000 rows of Table 1 (with the paper's citations): Erturk 601^2 and its Richardson-extrapolated value, Botella & Peyret 1998 (Chebyshev N=160, the most accurate published value: |psi|min = 0.1189366, omega = 2.067753, centre (0.5308, 0.5652)), Barragy & Carey 1997, and Ghia 1982.
- **Cross-check:** v(x=0.5) = 0.0258 agrees with the Ghia majority value 0.02526 to within Ghia's grid error. This is how PetIBM's 0.02426 was identified as a typo.

## Files
| file | bytes | sha256 |
|---|---:|---|
| `erturk2005_primary_vortex_Re1000.csv` | 366 | `9f97310ca574a87172b0c13269487790d30a3081e3c62754090106668d681b0a` |
| `erturk2005_u_vertical_centerline_Re1000.csv` | 321 | `b03833be10059b03bc15c65c886b5432080f08bb23336b8659d297b1b193a887` |
| `erturk2005_v_horizontal_centerline_Re1000.csv` | 320 | `b410221a7d441acf0ca3f87b74207bfb0769b676caeda9b26978328fb39ca802` |
