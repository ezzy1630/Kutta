| case | backend | cells | SIMPLEC its / s | Coupled its / s | speedup | loads settled S → C | Coupled vs SIMPLEC (final) | Coupled B/cell | peak RSS MiB S / C |
|---|---|---|---|---|---|---|---|---|---|
| bump-sa | cpu | 56320 | not reached in 3000 (723.1 s) | not reached in 3000 (241.8 s) |  | no → 1300 / 115.6 s | n/a: not both settled | 340 | 266 / 297 |
| bump-sa | metal | 56320 | not reached in 3000 (90.1 s) | not reached in 3000 (24.2 s) |  | no → 1300 / 12.2 s | n/a: not both settled | 340 | 249 / 257 |
| bump-sst | cpu | 56320 | not reached in 3000 (810.2 s) | not reached in 3000 (417.9 s) |  | no → 1300 / 194.6 s | n/a: not both settled | 340 | 296 / 299 |
| bump-sst | metal | 56320 | not reached in 3000 (137.5 s) | not reached in 3000 (55.8 s) |  | no → 1300 / 26.6 s | n/a: not both settled | 340 | 249 / 266 |
| channel-983k | metal | 983040 | not reached in 200 (27.7 s) | not reached in 200 (17.9 s) |  | no → no | n/a: not both settled | 343 | 2224 / 2448 |
| channel | cpu | 4608 | diverging (raw 4.2e-01 at 1000) | diverging (raw 9.2e+23 at 1000) |  | - → - | n/a: not both settled | 344 | 22 / 26 |
| channel | metal | 4608 | diverging (raw 3.9e-01 at 1000) | diverged at 73 (0.6 s) |  | - → - | n/a: not both settled | 344 | 33 / 33 |
| elbow | cpu | 17822 | 550 / 55.36 s (raw 1.45e-06) | 350 / 25.21 s (raw 5.73e-08) | 2.20x |  | drag 6.2e-06 | 389 | 154 / 167 |
| elbow | metal | 17822 | 550 / 13.34 s (raw 3.26e-07) | 350 / 8.40 s (raw 5.73e-08) | 1.59x |  | drag 7.9e-06 | 389 | 156 / 158 |
| plate-sa | cpu | 13056 | 1700 / 95.42 s (raw 1.40e-07) | 300 / 15.25 s (raw 6.84e-08) | 6.26x | 1700 / 95.4 s → 300 / 15.2 s | drag 3.7e-05, cf0 8.4e-05 | 340 | 72 / 79 |
| plate-sa | metal | 13056 | 1700 / 39.90 s (raw 1.08e-07) | 300 / 6.09 s (raw 7.18e-08) | 6.55x | 1700 / 39.9 s → 300 / 6.1 s | drag 3.4e-05, cf0 7.7e-05 | 340 | 73 / 77 |
| plate-sst | cpu | 13056 | 1500 / 83.26 s (raw 1.02e-07) | 400 / 25.37 s (raw 6.98e-08) | 3.28x | 1500 / 83.3 s → 400 / 25.4 s | drag 5.4e-06, cf0 4.9e-06 | 340 | 77 / 86 |
| plate-sst | metal | 13056 | 1500 / 54.99 s (raw 1.46e-07) | 400 / 10.27 s (raw 7.09e-08) | 5.35x | 1500 / 55.0 s → 400 / 10.3 s | drag 6.1e-06, cf0 5.2e-06 | 340 | 75 / 82 |
