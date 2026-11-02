# Separate eightfold trailing-edge mesh sweep

Both models use the molecular omega-production ratio 0.001; gamma also uses harmonic k. Original body spacing, wall-normal resolution and the polygon are retained; outer triangles are regenerated. This is not the source-only before/after comparison.

`*` means stopping criteria unmet and a final-2000 mean. A failed solve has no accuracy coefficient; its last finite sample remains in the receipt.

| Re | alpha | SST CL | SST CD | gamma CL | gamma CD | gamma dCL | gamma dCD | SST / gamma iterations |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| 200000 | 0 | failed | failed | failed | failed | failed | failed | 1050 failed / 1150 failed |
| 200000 | 4.76 | 1.5571 | 0.03065 | 1.5958 | 0.02788 | +0.0858 | +13.8% | 2650 / 2300 |
| 200000 | 7.93 | 1.7936 | 0.03963 | 1.8450 | 0.03632 | +0.0920 | +15.7% | 2800 / 2600 |
| 200000 | 10 | 1.8947 | 0.04969 | 1.9509 | 0.04621 | +0.0487 | +29.5% | 3250 / 2350 |
| 200000 | 13.17 | failed | failed | 1.8707* | 0.09167* | -0.2233 | +122.0% | 1050 failed / 6000* |
| 200000 | 16 | 2.1383* | 0.17166* | 2.1088* | 0.13262* | -0.0726 | n/a | 6000* / 6000* |
| 300000 | 0.48 | 1.1879 | 0.02223 | 1.2289* | 0.01909* | +0.0979 | +16.4% | 1850 / 6000* |
| 300000 | 4 | 1.5256 | 0.02717 | 1.5774 | 0.02260 | +0.1167 | +11.2% | 2300 / 2150 |
| 300000 | 8.21 | 1.8513 | 0.03766 | 1.9415 | 0.03055 | +0.1415 | +18.9% | 2850 / 2950 |
| 300000 | 12 | 2.0063 | 0.05849 | 2.1718 | 0.04354 | +0.1077 | +35.6% | 4250 / 4000 |
| 300000 | 14.74 | 1.8167* | 0.11995* | 2.2407 | 0.06070 | +0.0657 | +63.6% | 6000* / 3050 |
| 300000 | 16 | 2.0976* | 0.14628* | 2.2086* | 0.07711* | +0.0013 | +77.4% | 6000* / 6000* |
