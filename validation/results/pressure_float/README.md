# Float pressure assembly evidence (2026-10-29)

Numerics: `be8b993`; capture/replay benchmark: `35042d5`. The full dated report,
frozen V5 gates and study outcomes are in
[`VERIFICATION.md`](../../../src/kutta/fv/VERIFICATION.md#float-pressure-m-matrix-repair-2026-10-29).
`identity.json` binds the native tests to their checkout and executable hash.
These artifacts establish a pressure definiteness repair; they do not establish
an accepted TMR bump solution.

## Original-state replay

`bump_iteration18.assembly.gz` stores the original pressure assembly inputs:
incidence/CSR maps, stored float face coefficients, provisional mass flux and
reference selection. A CPU-only build of archived source `152ac72` recaptured
the state after 18 completed outer steps. Its complete COO output is byte for
byte identical to the uncompressed, unmodified
[`tmr_v5/bump_pressure_iteration18.coo.gz`](../tmr_v5/bump_pressure_iteration18.coo.gz).
The archived numerics were unchanged; only the benchmark's argument positions
were adapted to the old kernel signature, which did not take CSR row offsets.
`capture_provenance.json` records that equality and its SHA-256.

Replay uses the actual new `fv_pressure_assemble` kernel with those original
inputs. CPU and Metal COO outputs are byte identical. The retained
`bump_iteration18_reassembled.coo.gz` preserves every original off-diagonal;
1,750 diagonals change. The largest change is 4.65661287308e-10. Small downward
changes are possible where the original independent face sum rounded upward
more than the ceiling of the final stored row.

| Independent double check | Original | Reassembled |
| --- | ---: | ---: |
| Symmetry error | 0 | 0 |
| Negative row sums (`math.fsum`) | 1,712 | 0 |
| Smallest eigenvalue | -1.10101827e-11 | +1.87656240e-17 |
| Cholesky | FAIL | PASS |

The raw reassembled spectral condition estimate is 4.75113e14. Mathematical
diagonal congruence reduces it to 3.85366e7; congruence of the original matrix
still has a negative eigenvalue. No float row/column scaling is applied in the
solver. Its existing L1/Jacobi smoothing and pressure AMG meet the unchanged
1e-3 relative linear tolerance in six iterations on both replay backends
(verified full stored-operator residual ratio .000955778245003). The original
state's RHS norm is already 1.33542359297e14, and its correction reaches
4.64774617491e25. Those are divergent-flow diagnostics, not evidence of a
physically usable solution. See `replay_{cpu,metal}.log` and both analysis logs.

`bump_current_iteration15.coo.gz` is an actual matrix captured from the repaired
CPU flow path after 15 completed steps. It too has exact symmetry, no negative
row sums and passing Cholesky (smallest eigenvalue 6.98154696e-14). The following
monitor shows enormous loads; restoring definiteness alone does not repair the
coarse-grid startup transient. Its capture and analysis logs are retained.

Compressed artifact hashes are in `artifact_sha256.json`; gzip mtimes are zero.
The assembly text begins with cell count, reference cell (-1 for no pin), and
completed steps. Nine following lines each contain a length then that vector's
values: incidence offsets, face IDs, signs, adjacent cells, matrix entries,
CSR row offsets, diagonal IDs, face coefficients, mass flux.

`pressure_spd_comparison.svg` and `.png` plot these independent checks. Regenerate
with `python3 tools/plot_pressure_float.py` (Matplotlib required).

```sh
mkdir -p build/pressure-replay
gzip -dc validation/results/pressure_float/bump_iteration18.assembly.gz \
  > build/pressure-replay/original.assembly
./build/bench_pressure replay build/pressure-replay/original.assembly cpu build/pressure-replay/cpu
./build/bench_pressure replay build/pressure-replay/original.assembly metal build/pressure-replay/metal
cmp build/pressure-replay/cpu.coo build/pressure-replay/metal.coo
python3 tools/analyze_tmr_pressure.py build/pressure-replay/cpu.coo \
  --compare validation/results/tmr_v5/bump_pressure_iteration18.coo.gz --require-spd
# Capture a new matrix from the actual coarse bump flow path:
./build/bench_pressure capture validation/data/nasa_tmr cpu 15 build/pressure-replay/current
```

## Verification logs

The repeated-neighbour regression fails with independent face-diagonal
summation on CPU and Metal (`regression_before.log`) and passes when restored
(`regression_after.log`). The quick suite also checks a valid 32-cell mesh with
cell-width ratio 1e7 and maximum aspect ratio 2.5e7, boundary/pinned/Neumann
pressure, independent double Cholesky, deterministic reassembly, momentum and
SST transport, actual pressure coarse matrices before/after nonuniform refresh,
and explicit constant-nullspace solves/incompatible forcing.

`quick.log`, `ctest.log` and `cuda.log` retain the native quick tier, CTest with
portability lint, and CUDA 13.0.2 Docker build/CPU tests. CUDA uses the CI
architecture list `75-real;86-real;89-real;120`; no NVIDIA runtime is available.
Native quick/slow runs explicitly use the external Fluent fixture directory
requested by the user. The container uses bundled fixtures; absent external
Fluent datasets are explicitly reported in its log.

`slow_v1_v4.log`, `slow_v6.log`, `slow_v5.log`, `slow_status.json` and
`measurements.csv` retain the frozen slow study outputs and timing. The V5
CSV is generated by `tools/report_tmr.py`, preserving failed states as NaN.
No V5 gate, study configuration, run budget or retained slow assertion changed.

The final slow outcomes are six retained V2 failures, thirty retained V6
failures, and twenty-six V5 failures. All six plate runs complete 3,000 steps
but fail normalized continuity. Coarse/middle bump runs fail on nonfinite
input; CPU fine reaches its 600-second check at 2,200 steps, Metal fine reaches
3,000 steps. Their last-100-step load changes are 2.0295% / 1.2115%, so both
remain unconverged. The extra V5 failure is the CPU fine front Cf at its earlier
budget endpoint; no accuracy or convergence gate changed.
