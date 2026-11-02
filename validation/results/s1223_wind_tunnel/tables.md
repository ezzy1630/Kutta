# S1223 coefficients at every reference angle

Lift gate: absolute error <= 0.10 for alpha <= 10 deg. Drag gate: relative error <= 25% wherever UIUC CD exists.

After: original meshes, molecular omega-production denominator ratio 0.001; transition also uses harmonic k. Solver defaults are unchanged. The refined-tail sweep is a separate mesh experiment.

`*` means forces or residuals did not meet the existing convergence criterion; the value is a last-2000-iteration mean.

| Re | alpha | UIUC CL | UIUC CD | SST CL before | SST CD before | SST CL after | SST CD after | gamma CL before | gamma CD before | Fluent-details CL after | Fluent-details CD after | dCL after | dCD after | after gates CL/CD |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 200000 | 0 | 1.0830 | 0.01910 | 1.1289 | 0.02305 | 1.1224 | 0.02362 | 1.1556* | 0.02144* | 1.1545* | 0.02143* | +0.0715 | +12.2% | PASS/PASS |
| 200000 | 4.76 | 1.5100 | 0.02450 | 1.5847 | 0.03007 | 1.5748 | 0.03095 | 1.6089* | 0.02848* | 1.6129 | 0.02812 | +0.1029 | +14.8% | FAIL/PASS |
| 200000 | 7.93 | 1.7530 | 0.03140 | 1.8224 | 0.03861 | 1.8081 | 0.03988 | 1.8544 | 0.03705 | 1.8600 | 0.03656 | +0.1070 | +16.4% | FAIL/PASS |
| 200000 | 10 | 1.9022 | 0.03568 | 1.9337 | 0.04774 | 1.9149 | 0.04955 | 1.9489 | 0.04828 | 1.9684 | 0.04623 | +0.0662 | +29.6% | PASS/FAIL |
| 200000 | 13.17 | 2.0940 | 0.04130 | 1.9493 | 0.07907 | 1.8965 | 0.08563 | 1.6848* | 0.16569* | 1.9011* | 0.08937* | -0.1929 | +116.4% | trend/FAIL |
| 200000 | 16 | 2.1814 | n/a | 2.1802* | 0.16587* | 2.1954* | 0.17264* | 2.0387* | 0.13634* | 2.0611* | 0.13290* | -0.1204 | n/a | trend/n/a |
| 300000 | 0.48 | 1.1310 | 0.01640 | 1.2133 | 0.02209 | 1.2084 | 0.02249 | 1.2514* | 0.01932* | 1.2509 | 0.01935 | +0.1199 | +18.0% | FAIL/PASS |
| 300000 | 4 | 1.4607 | 0.02032 | 1.5508 | 0.02692 | 1.5443 | 0.02746 | 1.5951 | 0.02292 | 1.5966 | 0.02288 | +0.1358 | +12.6% | FAIL/PASS |
| 300000 | 8.21 | 1.8000 | 0.02570 | 1.8739 | 0.03709 | 1.8636 | 0.03795 | 1.9525 | 0.03101 | 1.9553 | 0.03083 | +0.1553 | +20.0% | FAIL/PASS |
| 300000 | 12 | 2.0641 | 0.03210 | 2.0357 | 0.05715 | 2.0198 | 0.05881 | 2.1861 | 0.04367 | 2.1854 | 0.04386 | +0.1213 | +36.6% | trend/FAIL |
| 300000 | 14.74 | 2.1750 | 0.03710 | 1.8386* | 0.11393* | 1.7863* | 0.12485* | 2.2575 | 0.06098 | 2.2541 | 0.06145 | +0.0791 | +65.6% | trend/FAIL |
| 300000 | 16 | 2.2073 | 0.04348 | 2.1512* | 0.14482* | 2.1614* | 0.14991* | 2.2229 | 0.07655 | 2.2237 | 0.07662 | +0.0165 | +76.2% | trend/FAIL |
