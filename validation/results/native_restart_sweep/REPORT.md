# Native restart sweep on the CPU

All 15 rotating-shaft runs met the same configured residual criteria. Warm
sweeps took 1,720 iterations each, compared with 1,580 for five independent
starts. The wall-time difference was 7–9% on this pass. Peak process footprint
stayed near 5.7 MB. The largest torque difference at a matched operating point
was 0.0059%.

| Sweep | Operating points (rad/s) | Total iterations | Elapsed, all 5 (s) | Median per point, min–max (s) | Peak process footprint (MB) | Worst torque difference from cold |
| --- | --- | ---: | ---: | --- | ---: | ---: |
| Cold | 8, 9, 10, 11, 12 | 1,580 | 5.02 | 1.01, 0.93–1.05 | 5.67 | — |
| Warm forward | 8, 9, 10, 11, 12 | 1,720 | 5.48 | 1.13, 0.96–1.18 | 5.62 | 0.0029% |
| Warm reverse | 12, 11, 10, 9, 8 | 1,720 | 5.38 | 1.06, 1.02–1.17 | 5.67 | 0.0059% |

The analytical Taylor–Couette setup uses the packaged rotating-shaft case on
one 16 × 72 annular grid. Shaft speed is the only varying input. Each warm
chain has one cold start and four restarts. All runs wrote a checkpoint, so
checkpoint output is included in each timing. The run cap was 4,000 iterations.
Residual checks used continuity 10⁻⁴ and both momentum equations 10⁻⁷. Each
condition met its threshold with no waivers. Warm starts each rebuild the
boundary fluxes and normalize continuity over their own first five
iterations. A matching scaled threshold can therefore correspond to different
absolute residuals. That affects the interpretation of iterations to converge.

The run used the Release CPU solver on an Apple M4 Pro. The measurement is one
interleaved pass. Elapsed time includes process startup, setup, solving and
output files. Peak is the largest sampled macOS process footprint reported by
the CLI's memory phases. [Raw measurements](sweep.json) include per-point
shaft and housing torque, the solver and mesh hashes and checkpoint hashes.
This experiment establishes same-mesh load repeatability for this laminar
case. Timing and memory savings on large or turbulent geometries remain open.

```sh
python3 validation/cases/native_restart_sweep.py build/kutta runs/native-restart-sweep
```
