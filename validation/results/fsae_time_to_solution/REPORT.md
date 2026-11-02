# FSAE time to solution on Metal

Ahmed L1 SIMPLEC reaches the same converged answer faster in all four paired
runs: solve speedup 1.098 [1.003, 1.218], or about 8.9% less solve time at the
paired median. Whole-process speedup is 1.090 [1.006, 1.208]. These are shared-GPU
observations, with a smallest solve gain of only 0.3%; they do not establish an
uncontended or universal speedup. Wing SIMPLEC has a slower paired median,
wing Coupled is mixed, and Ahmed Coupled does not converge within the existing
budget. The million-cell fixed-workload results are also too variable to claim
a consistent gain.

The solver is frozen at `7b817f3`. All requested CPU/Metal parity checks and the
quick tier pass. The full slow tier preserves main's exact assertion set and
values, including its 62 known failures.


Integration note: the measured commits were built on `85c489b`. On main the code
commits are `c8cc008` (1a89f69), `254afa2` (cc7366c), `5cf7aed` (ea3c058) and
`4fabc2d` (7b817f3). The comparator change `cd798e9` was superseded by main's
equivalent `971fbb3`. Timings and hashes below were measured on the original
commits; integrated parity was rechecked on CPU and Metal before merge.

## Builds and inputs

The unchanged baseline is `85c489bbbcd685933f30188272315f9a61702474`.
The candidate numerical source is `7b817f3b0e463a582593f2cba6427c9d29adf940`.
The comparison-tool fix is `cd798e9` and changes no solver code.

Both builds use Release, C++20, HDF5 and Metal on an Apple M4 Pro with 24 GB
unified memory. [Metadata](raw/metadata.json) records the compiler, platform and
binary SHA-256 values. Main was built before any edits, in `build-baseline`.
Frozen copies remain in `<checkout>/build-baseline/frozen`; candidate copies
are in `<checkout>/build/final-frozen`.

Generated Metal source is embedded in each executable and compiled by the
system Metal runtime with safe math and precise floating-point functions.
There is no external Kutta library or metallib to drift behind a frozen binary;
its SHA-256 binds both host code and embedded kernel source.

| Case | Cells | Input |
| --- | ---: | --- |
| Wing in ground effect, h/c 0.224, L2 | 195,581 | Existing generated fixture |
| Ahmed body, 35 degrees, L1 | 403,560 | Existing generated fixture |
| Ahmed body, 35 degrees, L3 | 1,329,048 | Generated locally with `ahmed.py --slant 35 --level 3` |

[Mesh hashes](raw/meshes.json) pin the exact inputs. The current Ahmed fixture
has 403,560 cells; the older coupled-plateau report names a different generation
with 403,604. Every baseline/candidate pair uses the same pinned file. The L3
mesh and all restricted data stay out of git.

The FSAE settings are unchanged. All cases use SST and the standard limiter.
The wing uses Re 462,000, downforce, moving ground, turbulence intensity 0.001
and viscosity ratio 1. The Ahmed cases use Re 2,780,000, length 1.044, reference
area 0.112, body-wall forces, z lift and turbulence intensity 0.01. SIMPLEC uses
momentum relaxation 0.7 and pressure relaxation 1. Coupled uses 0.5 and 0.5 with
flow CFL 200.

## Fixed iteration measurements

| Case | Cells | Scheme | Build | s/iteration | Solve seconds, 200 iterations | Process seconds | Peak GB |
| --- | ---: | --- | --- | --- | --- | --- | ---: |
| wing-L2 | 195,581 | SIMPLEC | baseline | 0.18956 [0.17827, 0.21553] | 37.912 [35.655, 43.107] | 40.03 [37.49, 44.87] | 0.703 |
| wing-L2 | 195,581 | SIMPLEC | candidate | 0.18233 [0.17107, 0.20583] | 36.466 [34.214, 41.165] | 38.55 [36.09, 42.94] | 0.703 |
| wing-L2 | 195,581 | Coupled | baseline | 0.18280 [0.15918, 0.18980] | 36.559 [31.837, 37.959] | 38.38 [33.76, 40.12] | 0.759 |
| wing-L2 | 195,581 | Coupled | candidate | 0.17514 [0.15469, 0.20174] | 35.028 [30.939, 40.347] | 36.98 [32.94, 42.12] | 0.786 |
| ahmed-L1 | 403,560 | SIMPLEC | baseline | 0.24755 [0.20526, 0.27254] | 49.510 [41.053, 54.508] | 60.08 [51.94, 64.77] | 1.234 |
| ahmed-L1 | 403,560 | SIMPLEC | candidate | 0.24621 [0.22889, 0.26128] | 49.243 [45.778, 52.257] | 60.08 [56.14, 63.72] | 1.234 |
| ahmed-L1 | 403,560 | Coupled | baseline | 0.27457 [0.18365, 0.31031] | 54.914 [36.729, 62.062] | 65.53 [48.65, 72.73] | 1.372 |
| ahmed-L1 | 403,560 | Coupled | candidate | 0.25161 [0.17904, 0.27379] | 50.322 [35.809, 54.759] | 61.47 [46.99, 65.20] | 1.437 |
| ahmed-L3 | 1,329,048 | SIMPLEC | baseline | 0.66469 [0.42008, 0.68065] | 132.938 [84.017, 136.131] | 166.97 [125.11, 171.13] | 3.762 |
| ahmed-L3 | 1,329,048 | SIMPLEC | candidate | 0.60133 [0.40600, 0.84041] | 120.267 [81.201, 168.083] | 152.56 [121.63, 201.54] | 3.762 |
| ahmed-L3 | 1,329,048 | Coupled | baseline | 0.76762 [0.65469, 1.20450] | 153.524 [130.938, 240.900] | 188.70 [169.33, 271.25] | 4.242 |
| ahmed-L3 | 1,329,048 | Coupled | candidate | 0.90147 [0.48809, 1.11295] | 180.294 [97.618, 222.590] | 213.26 [137.88, 255.95] | 4.445 |

Paired ratios are baseline time divided by candidate time. Each entry is the
median [min, max] of four pairs. A ratio above 1 favors the candidate.

| Case | Scheme | Solve speedup | Process speedup | Faster solve pairs |
| --- | --- | --- | --- | ---: |
| wing-L2 | SIMPLEC | 1.045 [1.024, 1.056] | 1.042 [1.023, 1.055] | 4/4 |
| wing-L2 | Coupled | 1.034 [0.900, 1.097] | 1.030 [0.903, 1.097] | 3/4 |
| ahmed-L1 | SIMPLEC | 1.007 [0.897, 1.043] | 0.995 [0.925, 1.031] | 2/4 |
| ahmed-L1 | Coupled | 1.035 [1.026, 1.260] | 1.041 [1.015, 1.192] | 4/4 |
| ahmed-L3 | SIMPLEC | 1.051 [0.792, 1.172] | 1.041 [0.849, 1.138] | 3/4 |
| ahmed-L3 | Coupled | 1.067 [0.658, 1.406] | 1.050 [0.731, 1.267] | 3/4 |

![Paired fixed-workload speedups with full ranges](speedup.png)

All four pairs favor the candidate for wing SIMPLEC and Ahmed L1 Coupled.
The other combinations have slower pairs, so no consistent gain is claimed.
The largest Coupled case has a slower unpaired candidate median despite a
paired median ratio above 1. Its wide range prevents a reliable speedup claim.
All 24 pairs have identical final reported numerical metrics.
Peak memory is the maximum of four runs, in decimal GB (1e9 bytes).

Each comparison alternates baseline and candidate, with four pairs per case
and scheme. Pair order is baseline/candidate, candidate/baseline,
baseline/candidate, candidate/baseline. The fixed workload is 200 iterations with the usual checks every
50. Both processes inherit `MallocSpaceEfficient=1`. Solve time includes the
benchmark's monitoring and force sampling. Process time also includes mesh
reading, setup and final wall-report construction. Peak memory is BSD time's peak process footprint,
which includes touched Metal buffers on this unified-memory machine.

Tables report medians and the complete min-to-max range. Speedups use the
median of the four paired baseline/candidate ratios. Earlier development runs
shared the GPU with the full test tiers and were stopped before a complete
measurement series. They do not establish a speedup. Their completed numeric
receipts are retained separately.

## Time to the unchanged convergence criterion

The criterion and stopping iteration are identical between builds. The three
converged combinations each have four counterbalanced pairs. Times are medians
[min, max]; N is the actual stopping iteration, not a reduced iteration budget.

| Case | Scheme | N | Build | Solve seconds | Process seconds |
| --- | --- | ---: | --- | --- | --- |
| wing-L2 | SIMPLEC | 2,400 | baseline | 457.09 [304.66, 750.12] | 458.75 [306.42, 751.76] |
| wing-L2 | SIMPLEC | 2,400 | candidate | 474.71 [395.18, 531.80] | 476.28 [397.05, 534.12] |
| wing-L2 | Coupled | 700 | baseline | 94.44 [62.90, 128.61] | 96.11 [64.47, 131.05] |
| wing-L2 | Coupled | 700 | candidate | 91.63 [57.84, 128.95] | 93.47 [59.49, 130.62] |
| ahmed-L1 | SIMPLEC | 1,100 | baseline | 216.47 [140.11, 233.16] | 225.81 [149.85, 241.75] |
| ahmed-L1 | SIMPLEC | 1,100 | candidate | 186.46 [139.70, 213.75] | 195.19 [148.96, 225.25] |

| Case | Scheme | Paired solve speedup | Paired process speedup | Faster solve pairs |
| --- | --- | --- | --- | ---: |
| wing-L2 | SIMPLEC | 0.946 [0.771, 1.459] | 0.946 [0.772, 1.458] | 1/4 |
| wing-L2 | Coupled | 1.041 [0.896, 1.268] | 1.041 [0.897, 1.262] | 3/4 |
| ahmed-L1 | SIMPLEC | 1.098 [1.003, 1.218] | 1.090 [1.006, 1.208] | 4/4 |

![Paired time-to-convergence speedups with full ranges](converged-speedup.png)

Wing SIMPLEC loses three of four convergence pairs despite winning all four
200-iteration pairs. Its full-trajectory paired median is 0.946, with a broad
0.771–1.459 range. Wing Coupled wins three of four. Only Ahmed SIMPLEC has a
favorable result in every convergence pair, including process time. None of
these observations isolates the effects of individual commits.

Ahmed L1 Coupled fails both residual and force-settling gates at 6,000 iterations
in both builds. Its single budget probe takes 843.32 baseline solve seconds and
883.15 candidate seconds (854.45 and 894.18 process seconds). It is not a
converged-answer speedup measurement and has no repeated median. The protocol
stops repeating a combination if either build fails its first convergence pair.

Every one of the 13 convergence pairs has the same stopping iteration, final
reported values and exact reported force/residual checkpoint strings. The checkpoint
histories are hashed in [the analysis](summary.json). The timestamps
and scalar results stay in [the run receipts](raw/convergence.jsonl).

The benchmark requires every applicable Fluent-scaled residual to be at most
1e-3 and both force coefficients to change by at most 5e-4 over 500 iterations.
The run limit is the benchmark's existing 6,000 iterations. No convergence
threshold, iteration budget or numerical acceptance test was changed.

## Profile before and after

Values below are GPU elapsed milliseconds per outer iteration, median [min,
max] of three profiles. They are timestamp spans, not exclusive GPU busy cycles;
shared GPU work can affect them. The groups are kernel families, not exact
physical-phase partitions. CSR/reductions and replay serve both AMG and Krylov;
transport assembly serves both momentum and turbulence.

**Wing L2, SIMPLEC**

| Kernel family | Baseline ms/iteration | Candidate ms/iteration |
| --- | --- | --- |
| Krylov/CSR/reductions | 19.04 [19.04, 22.80] | 19.09 [18.95, 22.31] |
| AMG smoothing, transfers and shared residuals | 9.51 [9.22, 10.92] | 9.43 [8.53, 9.56] |
| Small-operation replay (mixed) | 17.09 [16.68, 19.14] | 16.46 [15.79, 18.98] |
| Operator refresh and validation | 0.31 [0.31, 0.40] | 0.32 [0.31, 0.38] |
| Gradients, limiter and face reconstruction | 9.43 [8.76, 9.43] | 9.23 [8.84, 9.44] |
| Transport assembly (momentum and turbulence) | 2.02 [1.94, 2.17] | 2.02 [1.96, 2.18] |
| Momentum, pressure and flux coupling | 2.48 [2.33, 3.15] | 2.45 [2.42, 3.07] |
| SST closure | 3.11 [2.93, 3.87] | 3.15 [3.13, 3.88] |
| Other FV/runtime | 1.64 [1.34, 1.68] | 1.65 [1.40, 1.75] |

**Wing L2, Coupled**

| Kernel family | Baseline ms/iteration | Candidate ms/iteration |
| --- | --- | --- |
| Krylov/CSR/reductions | 18.61 [18.19, 20.75] | 17.77 [17.37, 20.25] |
| AMG smoothing, transfers and shared residuals | 30.76 [29.85, 32.40] | 28.59 [28.04, 30.05] |
| Small-operation replay (mixed) | 7.28 [6.94, 8.17] | 7.45 [6.95, 7.89] |
| Operator refresh and validation | 3.37 [3.25, 3.42] | 3.75 [3.65, 3.84] |
| Gradients, limiter and face reconstruction | 9.30 [8.91, 9.67] | 9.37 [9.00, 9.70] |
| Transport assembly (momentum and turbulence) | 2.09 [2.00, 2.22] | 2.04 [2.02, 2.15] |
| Momentum, pressure and flux coupling | 2.31 [2.29, 2.93] | 2.36 [2.33, 2.93] |
| SST closure | 3.08 [3.02, 3.79] | 3.15 [2.94, 3.74] |
| Other FV/runtime | 0.52 [0.48, 0.55] | 0.59 [0.53, 0.61] |

**Ahmed L1, SIMPLEC**

| Kernel family | Baseline ms/iteration | Candidate ms/iteration |
| --- | --- | --- |
| Krylov/CSR/reductions | 33.23 [32.78, 35.73] | 33.05 [32.31, 35.65] |
| AMG smoothing, transfers and shared residuals | 14.91 [14.82, 15.30] | 12.91 [12.46, 13.09] |
| Small-operation replay (mixed) | 26.52 [25.95, 28.79] | 25.07 [24.89, 28.35] |
| Operator refresh and validation | 0.94 [0.92, 1.11] | 1.03 [0.91, 1.04] |
| Gradients, limiter and face reconstruction | 26.68 [26.53, 28.12] | 27.87 [26.66, 28.08] |
| Transport assembly (momentum and turbulence) | 4.91 [4.79, 5.12] | 5.01 [4.88, 5.24] |
| Momentum, pressure and flux coupling | 4.17 [3.92, 4.78] | 4.03 [3.99, 4.53] |
| SST closure | 4.65 [4.54, 5.73] | 4.53 [4.40, 5.93] |
| Other FV/runtime | 1.91 [1.82, 2.25] | 2.02 [1.97, 2.23] |

**Ahmed L1, Coupled**

| Kernel family | Baseline ms/iteration | Candidate ms/iteration |
| --- | --- | --- |
| Krylov/CSR/reductions | 23.40 [22.65, 27.08] | 23.36 [21.25, 25.37] |
| AMG smoothing, transfers and shared residuals | 43.57 [42.79, 45.26] | 41.02 [39.34, 42.72] |
| Small-operation replay (mixed) | 9.81 [9.15, 12.43] | 9.50 [8.00, 11.84] |
| Operator refresh and validation | 12.54 [12.16, 14.13] | 14.70 [14.24, 14.85] |
| Gradients, limiter and face reconstruction | 27.80 [26.89, 28.58] | 28.27 [27.76, 28.78] |
| Transport assembly (momentum and turbulence) | 5.20 [4.97, 5.37] | 5.19 [5.19, 5.61] |
| Momentum, pressure and flux coupling | 3.91 [3.41, 4.12] | 3.98 [3.23, 4.33] |
| SST closure | 6.02 [4.86, 6.48] | 5.97 [4.16, 6.27] |
| Other FV/runtime | 1.05 [0.93, 1.13] | 1.16 [0.93, 1.36] |

Coupled moves work into operator refresh: Ahmed's refresh median rises from
12.54 to 14.70 ms/iteration while the smoothing/transfer/shared-residual group
falls from 43.57 to 41.02. Wing's corresponding groups move from 3.37 to 3.75
and 30.76 to 28.59. Those profile changes do not by themselves establish net
native wall-time gains. SIMPLEC's Ahmed smoothing group falls from 14.91 to
12.91; shared CSR/reduction work remains about 33 ms and gradients about 27–28 ms.

The following counts come from the separate native-batching probes. Counts
are identical in all three repetitions; wait milliseconds show the full range.
A nonempty sync is an API synchronization that actually has queued GPU work.
A buffer wait is one wait call inside that sync; one sync can wait on several
buffers. Neither is equivalent to the public `synchronize()` call count.

| Case | Scheme | Build | Nonempty syncs/iteration | Buffer waits/iteration | Host wait ms/iteration | Uploads/iteration |
| --- | --- | --- | ---: | ---: | --- | ---: |
| wing-L2 | SIMPLEC | baseline | 30.34 | 38.98 | 75.03 [70.02, 95.88] | 0.02 |
| wing-L2 | SIMPLEC | candidate | 24.58 | 33.22 | 73.41 [68.20, 94.83] | 0.02 |
| wing-L2 | Coupled | baseline | 25.80 | 25.80 | 87.26 [82.21, 103.82] | 1.02 |
| wing-L2 | Coupled | candidate | 23.88 | 23.88 | 83.63 [79.51, 101.23] | 0.02 |
| ahmed-L1 | SIMPLEC | baseline | 27.40 | 31.68 | 128.88 [128.15, 136.57] | 0.02 |
| ahmed-L1 | SIMPLEC | candidate | 21.98 | 26.26 | 125.29 [123.10, 135.07] | 0.02 |
| ahmed-L1 | Coupled | baseline | 17.10 | 17.10 | 137.98 [130.48, 151.92] | 1.02 |
| ahmed-L1 | Coupled | candidate | 15.28 | 15.28 | 136.15 [128.27, 149.73] | 0.02 |

The predictor removes 5.76 and 5.42 nonempty synchronizations per SIMPLEC
iteration on wing and Ahmed, and 1.92 and 1.82 on Coupled. Public sync-call
counts can rise when a norm and status are read together: the second call is
empty. Coupled's former scalar-limit upload also caused an empty sync, so its
removal must not be counted as an additional real round trip. Native
backpressure waits are zero in every measured window.

All builds perform seven public Device allocations per 50 measured iterations:
one force-mask upload plus six pressure/viscous-force output buffers in
`SteadySolver::forces()`. Their total requested bytes are 10,985,884 for wing
and 755,328 for Ahmed. There is one such sample per measured window; no other
public Device allocation churn is observed inside these iterations. Downloads
per iteration are unchanged: wing SIMPLEC 10.42, wing Coupled 15.86, Ahmed
SIMPLEC 7.52 and Ahmed Coupled 7.22. These counters do not inspect host malloc,
Metal driver allocation or transfer byte counts.

The instrumented binaries use command-buffer GPU start/end timestamps, splitting
buffers at kernel-family boundaries. These probes change dispatch submission,
so their elapsed wall times are not performance measurements. Profiles run
50 warmup iterations, then measure iterations 51 through 100 with the same
50-iteration reporting cadence.

Changing kernel names can also split more timestamped buffers. Family totals
locate expensive work; they do not establish isolated speedups. The setup row
means operator refresh and validation inside the iteration loop. AMG hierarchy
construction happens once during solver setup; its CPU cost is included in
setup and process time, without a finer CPU partition.

A separate count-only mode keeps native command-buffer batching. It records
nonempty synchronizations, backpressure waits, time inside buffer waits,
Device allocations and transfers. Wait time overlaps GPU time. The profile
counts public Device allocations, not host malloc or driver allocations.
[Profiler metadata](raw/profilers.json) binds each instrumented binary to its
numerical source and hash. No probes remain in solver source.

## Changes

| Commit | Change | Numerical contract |
| --- | --- | --- |
| `1a89f69` | Cache the five directed-edge coefficients of the Coupled fine operator once per outer solve. Initialize its scalar limits with the existing kernel. | Same coefficient expressions and face order; same 24 scalar values |
| `cc7366c` | Prepare each fine and coarse smoother block's three reciprocals, three H-times-inverse products and Schur denominator once per outer solve. | Same products and parenthesization |
| `ea3c058` | Skip matrix products in AMG's first Jacobi sweep after its iterate has been cleared to zero. | Finite matrix entries give residual b plus zero. A nonfinite RHS takes the original CSR path. Signed zeros are tested. |
| `7b817f3` | On larger Metal systems, compute the reliable residual when the recorded prior exit lies in the current batch. Reuse dead Krylov work vectors. | Same residual, norm, projection and verification; CPU and CUDA keep their existing policy. |
| `cd798e9` | Normalize CLI timing cells, temporary directory IDs and file-label pointer addresses in slow-assertion contexts. | Residual values, assertion outcomes and acceptance criteria remain unchanged. |

The Coupled caches trade memory for repeated work. The edge cache stores five
floats per directed cell-face entry. The block cache adds seven floats per fine
and coarse row. SIMPLEC needs no new buffers for the zero-start sweep.

The zero-start oracle compares every output bit against the ordinary CSR sweep
on CPU and Metal, with 513 and 2,053 rows. Its 15,396 checks cover signed zeros,
subnormals, large finite values, nonfinite RHS values and frozen solver states.

## Parity and verification

The CLI comparison uses 40 iterations for each case and backend. It compares
`fluxes.csv`, `forces.csv`, `loads.json` and `result.vtu` byte for byte. For
`residuals.csv`, it removes only lines beginning with `#` and the
`seconds_per_iteration` column, retaining every other byte and line ending.

| Case | CPU | Metal |
| --- | --- | --- |
| Mixing elbow, `--approximate wfgc` | Identical | Identical |
| Parameter elbow, `--approximate compressibility,ke-options` | Identical | Identical |
| Mixing elbow converted with the SIMPLEC override, then `--approximate wfgc` | Identical | Identical |

[Full hashes and source identities](raw/parity.jsonl) record all six comparisons.
The converted case uses exactly
`{"v":1,"overrides":{"methods/pressure-velocity-coupling":"simplec"}}`.
Full fields and converted meshes remain in ignored local runs.

All four `bench_overhead parity` hashes match main on each backend: eight
matching hashes in total. [The receipt](raw/overhead-parity.json) records the
source and binary identities. These hashes compare each backend with its own
baseline; CPU and Metal are not asserted to match one another.

The final quick tier passes 278 cases and 7,989,882 assertions, with zero
failures. The predictor oracle contributes 524,316 checks of output bits,
residual-stat bits and actual readback counts on CPU and Metal. Its baseline
policy fails only the two expected Metal wait-count checks; prediction passes
them. An additional 40-iteration wing L2 check compares all cell and wall
outputs byte for byte with main and finds zero field change.
The headline Ahmed L1 SIMPLEC case is also exported at its actual convergence
point: both builds stop at 1,100 iterations, and `cells.csv` and `walls.csv`
are byte-identical. [The converged-field receipt](raw/converged-field-parity.json)
records zero exported field change, input and executable hashes, and the
unchanged 6,000-iteration limit. These extra exports are parity checks, not
additional timing observations.

The baseline slow tier completed all 37 cases with 201,571 passed assertions and
62 failed assertions. Its six failing cases are V6 mixing elbow, 25; V12 fan,
1; V2 tet refinement, 8; V5 SST plate, 6; V5 SST bump, 8; and V7 SA, 14.
The final `7b817f3` candidate has exactly the same assertion set, outcomes and
expanded values: 201,633 assertions across all 37 cases, zero regressions, recoveries,
changed values, missing assertions or new assertions. The 62 existing failures
remain. [The comparison](raw/slow-comparison.json) agrees between the original
XML and compact receipts; both original XML hashes are recorded in the status
receipts. The full test process exits 1 for those known failures; the comparison
exits 0. No Fluent-reference or slow-tier acceptance criterion changed.

The slow comparison tool initially treated volatile CLI diagnostics as assertion
identities. Four regression checks now verify that timing, temporary paths and
pointer addresses can vary, physical residual changes remain distinguishable,
and a failed CLI assertion remains a regression. The original failure and
passing checks are retained in `raw/comparator-*.log`.

Complete doctest XML embeds generated meshes inside passing assertions. The
original XML stays local. Published compact receipts retain assertion identities,
outcomes and SHA-256 hashes of values. The original XML is compared before any
compaction, then the receipts are checked with the same comparator.

CUDA compile-only validation passes in the CUDA 13.0.2 development container for
architectures 75, 86 and 89 with `--fmad=false`. The generated CUDA file matches
the Mac build byte for byte. CUDA execution is unverified; this Mac has no
NVIDIA GPU. Portability lint passes.

## Trials not retained

The coefficient-layout and replay-threshold trials below used `cc7366c` as
their starting point, with three alternating pairs each. Times are solve seconds, median [min, max].
A paired speedup below 1 means the trial was slower. These runs had other test
jobs sharing the GPU and are separate from the final measurement series.

| Trial, Coupled | Case | Existing | Trial | Paired speedup |
| --- | --- | --- | --- | --- |
| 32-row coefficient tiles, 200 iterations | Wing L2 | 16.498 [13.532, 17.469] | 15.915 [13.224, 18.544] | 1.023 [0.890, 1.098] |
| 32-row coefficient tiles, 200 iterations | Ahmed L1 | 31.259 [25.680, 32.129] | 30.446 [28.320, 34.445] | 0.933 [0.843, 1.104] |
| Replay threshold 2,048 to 8,192, 100 iterations | Wing L2 | 8.532 [7.757, 9.666] | 9.331 [9.172, 11.116] | 0.870 [0.846, 0.914] |
| Replay threshold 2,048 to 8,192, 100 iterations | Ahmed L1 | 12.655 [12.480, 14.666] | 19.301 [18.800, 20.419] | 0.673 [0.647, 0.718] |

The tiled layout had mixed results and lost two of three Ahmed pairs. The larger
replay threshold lost every pair. Both were reverted. Dispatch widths 64, 128
and 512 were screened in 18 individual observations; that incomplete screening
does not establish a winner. The original dispatch widths remain.
[Trial receipts](raw/trials.json) retain the pairs and source basis.

A combined matrix-refresh/solve API was also tried and rejected. It combined
coefficient validation with the solve's entry norms, saving one readback per
scalar solve. Its 205,924 bitwise, guard and readback checks passed, as did all
six CLI comparisons. Three pairs did not establish a reliable wall-time gain:

| SIMPLEC, 200 iterations | Existing cached candidate | Combined refresh | Paired speedup |
| --- | --- | --- | --- |
| Wing L2 | 49.799 [49.534, 52.752] | 50.099 [48.012, 51.876] | 1.017 [0.994, 1.032] |
| Ahmed L1 | 71.701 [64.980, 79.296] | 71.461 [68.649, 78.884] | 1.044 [0.824, 1.110] |

The extra public API and split validation interface were removed. The existing
`refresh()` and `solve()` interfaces remain. [The pilot receipts](raw/refresh-pilot.jsonl)
retain every pair; prototype metadata records the source basis, patch hash and
binary hash. The emitted commit label is the prototype's configured basis.

## Remaining work by measured cost

1. Sparse linear algebra remains the largest measured GPU opportunity. On Ahmed
   SIMPLEC, CSR/reductions account for 33.05 ms/iteration and mixed replay for
   25.07. On Ahmed Coupled, the smoothing/transfer/shared-residual group is
   41.02 ms; fine smoothing alone is 16.54 and coarse smoothing 9.96. These
   categories overlap in purpose, so their costs must not be double-counted as
   separate predicted savings. Preserve operation order when removing passes.
2. Gradient and face reconstruction is still 27.87 ms/iteration for Ahmed
   SIMPLEC and 28.27 for Coupled. `fv_jump_gradient` is about 8 ms and
   `fv_moment_hop` about 7.4–7.5. Reusing an unchanged reconstruction or fusing
   compatible face passes is a better measured target than the roughly 5 ms
   transport assembly or 4.5–6 ms SST-closure group.
3. Native host round trips remain at 15.28–24.58 per iteration after prediction.
   Wait spans are 73–136 ms at the medians, overlapping GPU execution. Further
   batching needs to save real round trips without adding costly sparse work;
   the rejected refresh API shows that fewer waits alone is insufficient.
4. Resolve the missing converged answer for Ahmed Coupled: it consumes all
   6,000 iterations and 843–883 solve seconds without satisfying either gate.
   Also profile late wing SIMPLEC iterations before claiming that case improved.
   The present timestamp window covers iterations 51–100 and cannot explain its
   slower full-trajectory median. Any numerical change needs a separate commit
   and the unchanged reference/slow gates.

There is no measured payoff estimate for these unimplemented changes. Initial
AMG setup is not rebuilt every iteration, and the measured allocation churn is
periodic force reporting. Those are weaker targets than the costs above.
The million-cell case has memory headroom but no consistent timing gain; CUDA
runtime behavior and performance remain unverified.

## Reproduce

Configure the baseline and candidate source revisions in separate owned build
trees with the same options below. Keep those trees as well as the frozen
binaries: the tests invoke the CLI path compiled into each test binary.

```sh
cmake -S <checkout> -B <build> -DCMAKE_BUILD_TYPE=Release \
  -DKUTTA_BUILD_APP=OFF -DKUTTA_BUILD_TESTS=ON \
  -DKUTTA_ENABLE_METAL=ON -DKUTTA_ENABLE_CUDA=OFF \
  -DKUTTA_ENABLE_HDF5=ON -DKUTTA_REQUIRE_HDF5=ON
cmake --build <build> --parallel 4
```

Run the comparisons from `<checkout>`. Point `KUTTA_FLUENT_EXAMPLES` and
`KUTTA_RUNS` at the read-only fixture checkout. Keep the two frozen binaries
separate. The commands below use the retained build names from this run.

```sh
python3 validation/results/fsae_time_to_solution/reproduce.py \
  --baseline build-baseline/frozen/bench_fsae_aero \
  --candidate build/final-frozen/bench_fsae_aero --out runs/fsae-recheck
python3 validation/results/fsae_time_to_solution/reproduce.py \
  --baseline build-baseline/frozen/bench_fsae_aero \
  --candidate build/final-frozen/bench_fsae_aero --convergence \
  --cases wing-L2 ahmed-L1 --out runs/fsae-convergence-recheck
python3 validation/results/fsae_time_to_solution/parity.py \
  --baseline build-baseline/frozen/kutta --candidate build/final-frozen/kutta \
  --out runs/fsae-parity-recheck
build/final-frozen/kutta-tests --test-suite-exclude=slow
build/final-frozen/kutta-tests -ts=slow -r=xml --success > runs/slow-recheck.xml
python3 tools/forensics/compare_slow_assertions.py BEFORE.xml AFTER.xml
python3 tools/forensics/test_compare_slow_assertions.py
python3 tools/lint_portability.py
python3 validation/results/fsae_time_to_solution/analyze.py
uv run --with matplotlib==3.11.2 python validation/results/fsae_time_to_solution/plot.py
```
