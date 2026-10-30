# Hybrid initialization and pressure startup evidence — 2026-10-30

Hybrid initialization is implemented and the uniform-start bump is stable.
The frozen accuracy/convergence studies still fail: V5 retains 16/83 failed
assertions, V6 retains 29 settings/accuracy failures, and V1–V4 retains exactly
six tet-force failures. No acceptance threshold or frozen case setup changed.

The full report is in
[`src/kutta/fv/VERIFICATION.md`](../../../src/kutta/fv/VERIFICATION.md#hybrid-initialization-and-bounded-pressure-skew-correction-2026-10-30).
The [hybrid contract](../../../docs/research/hybrid-initialization.md) records
public sources and all implementation choices where Fluent's description is
incomplete. Code commits: `8b7e954`, `f7e6ef3`, `4be291c`. Nothing was pushed.

## Artifact map

| Artifact | Meaning |
| --- | --- |
| `metadata.json` | Exact checkout, code commit, platform, source/build/input hashes, unchanged protocol hashes and check counts |
| `unbounded_startup_cpu.log` | Original `58aa804` pressure-jump kernel, current synchronizing benchmark; nonfinite after 16 completed steps |
| `bounded_startup_cpu.log`, `bounded_startup_metal.log` | Same uniform setup, 24 steps with bounded extrapolation; pressure RHS and gradient/skew terms |
| `startup_rhs.csv`, `startup_rhs.png` | Measured CPU pressure RHS comparison, indexed by completed steps before assembly |
| `startup_regression_before.log` | New 24-step regression run with original kernel: eight failed assertions across CPU/Metal |
| `verified_quick.log`, `verified_ctest.log` | Native complete quick tier: 165 cases, 1,539,822 assertions, 12.32 s; ctest 2/2, 12.22 s |
| `verified_cuda_build.log` | CUDA 13.0.2 image, exact CI configure/build, portability lint and CPU quick tier |
| `baseline_v1_v4.log`, `verified_v1_v4.log`, `laminar_comparison.json` | Complete before/after CPU/Metal V1–V4 + Couette; all 308 numerical output lines equal; six retained tet-force failures |
| `pure_link_v2_rejected.log` | Rejected experiment that removed skew correction everywhere; eight failures, including two added tet pressure-order failures |
| `verified_v5_plate.log` | Complete plate-case excerpt, six grids; force/trend gates pass, six convergence assertions fail |
| `incomplete_v5_before_expansion.log` | Original full scoped log: plate completed, then bump aborted at an initially unexpanded fine-grid fixture; not a complete bump result |
| `verified_v5_bump.log` | Complete rerun after checksum-verified fixture expansion: all six bump grids finite, ten failed assertions / 47 |
| `measurements.log`, `measurements.csv`, `measurements_report.txt` | Exactly twelve unique plate/bump endpoints and pinned reference values, checked by `tools/report_tmr.py` |
| `verified_v6_hybrid.log`, `verified_v6_proxy.log`, `v6_measurements.json` | Actual elbow, CPU/Metal, seeded/hybrid/explicit stored-standard proxy; unchanged V6 gates evaluated before optional continuation |
| `plateau_hybrid/*.csv.gz`, `plateau_proxy/*.csv.gz` | All 17,822 cold-start cell rows at 1,200 steps, corresponding owned RCM order, columns u/v/w/p/k/omega |
| `initializer_comparison.json` | Full-field differences without resampling or pressure-gauge alignment |
| `stored_hybrid_settings.txt` | Relevant observed stored hybrid controls from the exact external elbow case |
| `cli_hybrid.log`, `cli_summary.json` | Real elbow CLI initialization and one CPU step, no data file; only WFGC disabled in a local diagnostic case copy |

Text logs use LF line endings with trailing terminal whitespace removed.
Missing metrics in the measurement CSV are blank, not zero. `converged=0`
endpoints remain transient diagnostics even when individual force checks pass.
The original failed-matrix replay and definiteness proof remain unchanged in
[`pressure_float`](../pressure_float/README.md).

## Reproduction

From the repository root, prepare the bundled compressed grids before the slow
study. The preparation script verifies their pinned hashes and downloads nothing.
The external Fluent fixtures are supplied separately.

```sh
python3 tools/prepare_tmr.py
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
export KUTTA_FLUENT_EXAMPLES=/Users/ezzy/.Others/Kutta/validation/fluent_examples
./build/kutta-tests --test-suite-exclude=slow --duration
ctest --test-dir build -LE slow --output-on-failure

# These slow commands intentionally retain their real nonzero gate exits:
./build/kutta-tests --test-suite=slow --test-case-exclude='V5*,*V6*' --duration
./build/kutta-tests --test-suite=slow --test-case='V5*' --duration
KUTTA_V6_COMPARE_DIR=validation/results/hybrid_startup/plateau_hybrid \
  ./build/kutta-tests --test-suite=slow --test-case='*V6*' --duration
KUTTA_V6_STANDARD_PROXY=1 \
  KUTTA_V6_COMPARE_DIR=validation/results/hybrid_startup/plateau_proxy \
  ./build/kutta-tests --test-suite=slow --test-case='*V6*' --duration

./build/bench_startup validation/data/nasa_tmr cpu 0 24
./build/bench_startup validation/data/nasa_tmr metal 0 24
python3 validation/results/hybrid_startup/compare_initializers.py
python3 tools/report_tmr.py \
  --test-log validation/results/hybrid_startup/measurements.log \
  --output validation/results/hybrid_startup
# Plot reproduction additionally needs matplotlib (measured with 3.11.2):
python3 validation/results/hybrid_startup/plot_startup.py
```

The original-kernel trace used only `fv_flow.kern.h` from `58aa804` in place of
the bounded kernel when building `bench_startup`. The benchmark does not change
fields or coefficients; its reads deliberately synchronize, so its wall time
is not a solver-performance measurement. The ordinary test executable was
rebuilt with the committed bounded kernel afterward. CUDA GPU execution is
unverified; Docker compilation covers all four requested architectures.

The optional V6 comparison continues each solver to 1,200 steps **after** the
normal convergence/accuracy checks. Main `V6_RESULT` times exclude initialization
and the later continuation. Both starts meet V6 residual criteria at the final
comparison. Velocity/k/omega relative L2 differences are below `5e-7`; pressure
is `1.5224e-5` CPU / `4.7420e-6` Metal. Pressure is within the inner solver's
`1e-3` tolerance, but not a `1e-6` relative field-difference bound. Raw pressure
comparisons are retained, including their small mean differences.

For the actual CLI diagnostic, a local copy of the original elbow changed only
`(recon/cell-lsf? #t)` to `#f` through HDF5. Its SHA-256 is
`b3b550e8b9c2a2902c45319d7543cc5eada682367b3e53fb8427a2748e69a0ff`.
Stored `hyb-init? #t` and all hybrid controls were retained:

```sh
# Optional diagnostic-copy preparation requires h5py; it verifies both hashes:
python3 validation/results/hybrid_startup/prepare_cli_diagnostic.py
./build/kutta solve build/elbow-hybrid-cold.cas.h5 \
  --backend cpu --iterations 1 --out build/elbow-hybrid-cli
```

This proves actual CLI dispatch/outputs; it does not claim that one step is a
converged solution or that the original WFGC-enabled case is supported. The
input case copy and generated VTU are local build artifacts, not redistributed
fixtures. Original external files were unchanged and their hashes freshly
matched the prior V6 report.

The CUDA run uses `nvidia/cuda:13.0.2-devel-ubuntu24.04` under linux/amd64 Docker,
with the repository mounted read-only and an anonymous build volume removed
with the container. Dependencies match CI. Its exact configure line is:

```sh
cmake -S . -B build -G Ninja -DKUTTA_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES="75-real;86-real;89-real;120"
cmake --build build -j 2
./build/kutta-tests --test-suite-exclude=slow --test-case-exclude='compiled GPU backends*'
```

External Fluent fixtures were not mounted into that compile-check container;
its native CPU fixture checks are supplied by the explicit macOS run above.
