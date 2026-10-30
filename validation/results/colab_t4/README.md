# CUDA on a Colab T4 (2026-10-30)

Raw evidence for [`GPU-RUN-2026-10-30.md`](../../../src/kutta/runtime/cuda/GPU-RUN-2026-10-30.md).
Source `a0c5caa`, Tesla T4, driver 580.82.07, nvcc 13.0.88, GCC 13.3, Ubuntu 24.04.

- `notebook_cells.txt`: every notebook cell as executed, plus the commands typed into the Colab terminal.
- `logs/`: one log per command. `environment.log` holds `nvidia-smi`, `nvcc` and host details.
  `tests_fast.log` is the default build; `tests_fast_nofmad.log` and
  `tests_parity_nofmad.log` are the `--fmad=false` diagnostic build;
  `tests_fast_without_fixtures.log` ran before the Ansys files were fetched.
  `tests_slow_v3_re100_400.log` contains both V3 cases; `tests_slow_v3_re1000.log`
  is a duplicate run that was killed on purpose. `fma_probe.log` holds the
  compiler contraction probe and CMake flags.
- `mixing_elbow_{cuda,cpu}/`: `summary.json`, `forces.csv`, `residuals.csv` from
  `kutta solve --approximate wfgc --iterations 600` (cold start; `result.vtu` not kept).
- `m4pro_same_commit/`: `bench_bandwidth`, `bench_linalg`, `bench_fv` for the same
  source on the owner's M4 Pro (CPU + Metal), run while the machine was loaded
  by other jobs. The source came from `git archive`, so the commit field is blank.
