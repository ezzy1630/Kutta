# Generated tables (tools/fsae_aero/report.py)

### Zerihan / Tyrrell 026, sectional Cl vs h/c (level 2, Metal)

| h/c | Kutta Cl | CFL3D k-omega | diff | 5% gate | CFL3D SA | expt (p-integrated) | Kutta Cd | wing y+ max | cells | iterations | max scaled residual | force change over window | stop | solve s | peak GB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | 0.8612 | 0.872 | -1.2% | PASS | 0.885 | 0.766 | 0.0161 | 0.59 | 185,964 | 1800 | 2.1e-05 | 4.7e-04 | converged | 166 | 1.01 |
| 0.671 | 1.0273 | 1.030 | -0.3% | PASS | 1.052 | 0.902 | 0.0191 | 0.62 | 251,502 | 2850 | 1.3e-05 | 4.1e-04 | converged | 328 | 1.41 |
| 0.448 | 1.1403 | 1.136 | +0.4% | PASS | 1.174 | 1.009 | 0.0215 | 0.63 | 227,294 | 2850 | 2.0e-05 | 4.5e-04 | converged | 327 | 1.25 |
| 0.313 | 1.2754 | 1.253 | +1.8% | PASS | 1.320 | 1.145 | 0.0250 | 0.64 | 211,593 | 2750 | 1.4e-05 | 4.9e-04 | converged | 248 | 1.17 |
| 0.224 | 1.4015 | 1.352 | +3.7% | PASS | 1.475 | 1.286 | 0.0305 | 0.64 | 195,581 | 2400 | 1.4e-05 | 4.3e-04 | converged | 243 | 1.05 |
| 0.179 | 1.4555 | 1.393 | +4.5% | PASS | 1.565 | 1.362 | 0.0360 | 0.63 | 187,656 | 1800 | 8.1e-06 | 4.4e-04 | converged | 228 | 1.03 |
| 0.134 | 1.4675 | 1.386 | +5.9% | FAIL | 1.604 | 1.385 | 0.0452 | 0.61 | 179,126 | 2000 | 1.1e-05 | 4.3e-04 | converged | 317 | 0.98 |
| 0.090 | 1.4107 | 1.297 | +8.8% | FAIL | 1.539 | 1.371 | 0.0590 | 0.56 | 177,312 | 2600 | 2.4e-05 | 4.7e-04 | converged | 450 | 0.97 |

### Zerihan mesh study (Metal); * = unsettled, mean of the last 2000 iterations

| h/c | level | Kutta Cl | vs k-omega | Kutta Cd | O-block layers | first cell / c | cells | iterations | max scaled residual | force change over window | stop | solve s | peak GB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | L1 | 0.8647 | -0.8% | 0.0169 | 14 | 4.0e-05 | 60,376 | 1300 | 9.6e-06 | 3.0e-04 | converged | 228 | 0.42 |
| free | L2 | 0.8612 | -1.2% | 0.0161 | 29 | 2.0e-05 | 185,964 | 1800 | 2.1e-05 | 4.7e-04 | converged | 166 | 1.01 |
| free | L3 | 0.8526 | -2.2% | 0.0159 | 57 | 1.0e-05 | 655,996 | 3300 | 3.8e-05 | 3.5e-04 | converged | 994 | 3.38 |
| 0.224 | L1 | 1.3844 | +2.4% | 0.0324 | 14 | 4.0e-05 | 56,463 | 1500 | 4.2e-06 | 4.7e-04 | converged | 191 | 0.41 |
| 0.224 | L2 | 1.4015 | +3.7% | 0.0305 | 29 | 2.0e-05 | 195,581 | 2400 | 1.4e-05 | 4.3e-04 | converged | 243 | 1.05 |
| 0.224 | L3 | 1.3809 | +2.1% | 0.0302 | 57 | 1.0e-05 | 732,898 | 3200 | 1.4e-05 | 4.0e-04 | converged | 864 | 3.77 |
| 0.090 | L1 | 1.4618* | +12.7% | 0.0660* | 14 | 4.0e-05 | 51,709 | 6000 | 7.7e-03 | 7.6e-02 | not settled | 174 | 0.38 |
| 0.090 | L2 | 1.4107 | +8.8% | 0.0590 | 29 | 2.0e-05 | 177,312 | 2600 | 2.4e-05 | 4.7e-04 | converged | 450 | 0.97 |
| 0.090 | L3 | 1.4018 | +8.1% | 0.0583 | 57 | 1.0e-05 | 660,776 | 5100 | 1.2e-05 | 5.0e-04 | converged | 1437 | 3.40 |

### CPU vs Metal (same mesh)

| h/c | level | Cl CPU | Cl Metal | Cd CPU | Cd Metal | iterations CPU / Metal | solve s CPU / Metal | peak GB CPU / Metal |
|---|---|---|---|---|---|---|---|---|
| 0.224 | L2 | 1.40147 | 1.40147 | 0.03051 | 0.03051 | 2400 / 2400 | 2099 / 243 | 0.92 / 1.05 |

### Far-field distance (inlet 30c / outlet 40c / top 30c, then doubled)

| h/c | level | Cl | Cl, 2x domain | change | Cd | Cd, 2x domain |
|---|---|---|---|---|---|---|
| free | L2 | 0.8612 | 0.8626 | +0.15% | 0.01614 | 0.01592 |

### Selig S1223, Re 2e+05 (level 2, Metal); UIUC values interpolated at the run angle; UIUC CLmax 2.188 at 16.9 deg; * = unsettled, mean of the last 2000 iterations

| alpha | Kutta CL | UIUC CL | dCL | +-0.10 gate | Kutta CD | UIUC CD | dCD | +-25% gate | cells | iterations | max scaled residual | force change over window | stop | solve s | peak GB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1.129 | 1.083 | +0.046 | PASS | 0.0231 | 0.0191 | +21% | PASS | 184,806 | 2100 | 4.5e-05 | 4.6e-04 | converged | 490 | 1.01 |
| 4.76 | 1.585 | 1.510 | +0.075 | PASS | 0.0301 | 0.0245 | +23% | PASS | 184,470 | 2150 | 3.1e-05 | 2.5e-04 | converged | 504 | 1.00 |
| 7.93 | 1.822 | 1.753 | +0.069 | PASS | 0.0386 | 0.0314 | +23% | PASS | 184,566 | 2100 | 2.3e-05 | 3.7e-04 | converged | 641 | 1.01 |
| 10 | 1.934 | 1.902 | +0.031 | PASS | 0.0477 | 0.0357 | +34% | FAIL | 183,892 | 2100 | 3.1e-05 | 4.7e-04 | converged | 571 | 0.99 |
| 13.17 | 1.949 | 2.094 | -0.145 | trend | 0.0791 | 0.0413 | +91% | FAIL | 183,980 | 3800 | 1.8e-05 | 4.5e-04 | converged | 650 | 1.01 |
| 16 | 2.180* | 2.181 | -0.001 | trend | 0.1659* | n/a | n/a | n/a | 183,868 | 6000 | 4.7e-03 | 1.5e-01 | not settled | 648 | 1.01 |

### Selig S1223, Re 3e+05 (level 2, Metal); UIUC values interpolated at the run angle; UIUC CLmax 2.212 at 15.6 deg; * = unsettled, mean of the last 2000 iterations

| alpha | Kutta CL | UIUC CL | dCL | +-0.10 gate | Kutta CD | UIUC CD | dCD | +-25% gate | cells | iterations | max scaled residual | force change over window | stop | solve s | peak GB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.48 | 1.213 | 1.131 | +0.082 | PASS | 0.0221 | 0.0164 | +35% | FAIL | 184,792 | 2150 | 2.8e-05 | 3.5e-04 | converged | 323 | 1.01 |
| 4 | 1.551 | 1.461 | +0.090 | PASS | 0.0269 | 0.0203 | +32% | FAIL | 184,606 | 2150 | 3.8e-05 | 3.8e-04 | converged | 299 | 1.00 |
| 8.21 | 1.874 | 1.800 | +0.074 | PASS | 0.0371 | 0.0257 | +44% | FAIL | 184,216 | 2100 | 1.6e-05 | 2.7e-04 | converged | 282 | 1.00 |
| 12 | 2.036 | 2.064 | -0.028 | trend | 0.0572 | 0.0321 | +78% | FAIL | 183,860 | 2400 | 1.2e-05 | 3.8e-04 | converged | 325 | 1.01 |
| 14.74 | 1.839* | 2.175 | -0.336 | trend | 0.1139* | 0.0371 | +207% | FAIL | 183,454 | 6000 | 4.8e-04 | 7.5e-03 | not settled | 749 | 0.99 |
| 16 | 2.151* | 2.207 | -0.056 | trend | 0.1448* | 0.0435 | +233% | FAIL | 183,868 | 6000 | 4.5e-03 | 9.7e-02 | not settled | 606 | 1.00 |

### Ahmed body 35 deg, Cd mesh study (Metal); reference Cd 0.26, components front 0.020 / slant 0.090 / base 0.095

| level | Kutta Cd (final, or mean of last 2000 iterations if unsettled) | Cd range, last 2000 | vs 0.26 | 8% gate | friction Cd | front Cp drag | slant Cp drag | base Cp drag | Kutta CL (z) | slant area with reversed shear | body y+ mean / max | cells | iterations | max scaled residual | force change over window | stop | solve s | peak GB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L1 | 0.3531 | 0.3530-0.3532 | +35.8% | FAIL | 0.0469 | -0.0042 | 0.1816 | 0.1287 | 0.209 | 0.84 | 3.2 / 7 | 403,604 | 5000 | 1.1e-05 | 1.5e-03 | not settled | 715 | 2.16 |
| L2 | 0.3322 | 0.3269-0.3374 | +27.8% | FAIL | 0.0474 | -0.0059 | 0.1673 | 0.1180 | 0.169 | 0.88 | 3.0 / 7 | 690,104 | 5000 | 6.8e-03 | 3.2e-01 | not settled | 798 | 3.35 |
| L3 | 0.3137 | 0.3105-0.3170 | +20.7% | FAIL | 0.0481 | -0.0075 | 0.1505 | 0.1198 | 0.094 | 0.89 | 2.9 / 7 | 1,330,032 | 5000 | 5.7e-03 | 2.5e-01 | not settled | 1388 | 6.79 |

### Ahmed 35 deg rear Cp at the reference tap locations (Kutta L3, nearest wall face)

| reference | region | taps | reference mean Cp | Kutta mean Cp | RMS(Kutta - reference) | max tap-to-face distance mm |
|---|---|---|---|---|---|---|
| experiment | slant | 188 | -0.266 | -0.316 | 0.123 | 5.7 |
| experiment | base | 99 | -0.220 | -0.227 | 0.070 | 4.8 |
| CFX-5 SST | slant | 1110 | -0.233 | -0.342 | 0.139 | 6.4 |
| CFX-5 SST | base | 1041 | -0.173 | -0.253 | 0.133 | 6.3 |

### Ahmed slant topology on the centreline, stations x = -203, -163, -123, -83, -43, -3 mm (R reversed, A attached near-wall flow); reference rows in REPORT.md

| slant deg | source | stations | slant area with reversed near-wall flow |
|---|---|---|---|
| 25 | Kutta SST L2 | A R R R R R | 0.62 |
| 35 | Kutta SST L1 | - R R R R R | 0.84 |
| 35 | Kutta SST L2 | - R R R R A | 0.88 |
| 35 | Kutta SST L3 | - R R R R A | 0.89 |

### Verdicts

| case | point | measured vs reference | verdict |
|---|---|---|---|
| Zerihan 2D | h/c free | Cl 0.861 vs CFL3D k-omega 0.872 (-1.2%) | PASS |
| Zerihan 2D | h/c 0.671 | Cl 1.027 vs CFL3D k-omega 1.030 (-0.3%) | PASS |
| Zerihan 2D | h/c 0.448 | Cl 1.140 vs CFL3D k-omega 1.136 (+0.4%) | PASS |
| Zerihan 2D | h/c 0.313 | Cl 1.275 vs CFL3D k-omega 1.253 (+1.8%) | PASS |
| Zerihan 2D | h/c 0.224 | Cl 1.401 vs CFL3D k-omega 1.352 (+3.7%) | PASS |
| Zerihan 2D | h/c 0.179 | Cl 1.456 vs CFL3D k-omega 1.393 (+4.5%) | PASS |
| Zerihan 2D | h/c 0.134 | Cl 1.467 vs CFL3D k-omega 1.386 (+5.9%) | FAIL |
| Zerihan 2D | h/c 0.090 | Cl 1.411 vs CFL3D k-omega 1.297 (+8.8%) | FAIL |
| S1223 Re 2e+05 | alpha 0 | CL 1.129 vs 1.083 (+0.046) | PASS |
| S1223 Re 2e+05 | alpha 0 | CD 0.0231 vs 0.0191 (+21%) | PASS |
| S1223 Re 2e+05 | alpha 4.76 | CL 1.585 vs 1.510 (+0.075) | PASS |
| S1223 Re 2e+05 | alpha 4.76 | CD 0.0301 vs 0.0245 (+23%) | PASS |
| S1223 Re 2e+05 | alpha 7.93 | CL 1.822 vs 1.753 (+0.069) | PASS |
| S1223 Re 2e+05 | alpha 7.93 | CD 0.0386 vs 0.0314 (+23%) | PASS |
| S1223 Re 2e+05 | alpha 10 | CL 1.934 vs 1.902 (+0.031) | PASS |
| S1223 Re 2e+05 | alpha 10 | CD 0.0477 vs 0.0357 (+34%) | FAIL |
| S1223 Re 2e+05 | alpha 13.17 | CD 0.0791 vs 0.0413 (+91%) | FAIL |
| S1223 Re 3e+05 | alpha 0.48 | CL 1.213 vs 1.131 (+0.082) | PASS |
| S1223 Re 3e+05 | alpha 0.48 | CD 0.0221 vs 0.0164 (+35%) | FAIL |
| S1223 Re 3e+05 | alpha 4 | CL 1.551 vs 1.461 (+0.090) | PASS |
| S1223 Re 3e+05 | alpha 4 | CD 0.0269 vs 0.0203 (+32%) | FAIL |
| S1223 Re 3e+05 | alpha 8.21 | CL 1.874 vs 1.800 (+0.074) | PASS |
| S1223 Re 3e+05 | alpha 8.21 | CD 0.0371 vs 0.0257 (+44%) | FAIL |
| S1223 Re 3e+05 | alpha 12 | CD 0.0572 vs 0.0321 (+78%) | FAIL |
| S1223 Re 3e+05 | alpha 14.74 | CD 0.1139 vs 0.0371 (+207%) | FAIL |
| S1223 Re 3e+05 | alpha 16 | CD 0.1448 vs 0.0435 (+233%) | FAIL |
| Ahmed 35 deg | L3 (1,330,032 cells) | Cd 0.314 vs 0.26 (+20.7%) | FAIL |
