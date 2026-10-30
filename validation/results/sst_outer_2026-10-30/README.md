# SST outer-iteration measurements, 2026-10-30

The dated Performance section in `src/kutta/fv/VERIFICATION.md` is the report.
`metadata.json` identifies the checkout, numerical source, benchmark harness,
artifacts, unchanged numerical protocols and measurement limits. Nothing was
pushed. No source under `src/kutta/app` or `src/kutta/setup` was changed.

`measurements.json` preserves throughput and synchronization/allocation counts;
`kernels.json` preserves the intrusive per-kernel timings and logical traffic.
`analyze.py` regenerates the summary CSVs. `plot.py` uses matplotlib 3.11.2 to
regenerate both figures. The Re_H=10,000 comparison is the primary benchmark;
the Re_H=1,000 pilot records incremental changes. `raw` contains the unedited measurement/test logs. Its scoped Git whitespace
attribute preserves terminal whitespace; it does not change source-code checks
or numerical acceptance gates. Generated CSVs use LF line endings.

The 1M-cell generated SST channel is a transient execution benchmark. Its smaller
4,608-tet counterpart fails to converge on baseline and final, at both tested
Reynolds numbers. The 122,880-tet half-resolution check also fails: baseline
triggers a nonfinite-input guard after checkpoint2840, final remains
unconverged at3000 with unphysical loads. That requested same-channel convergence result is unmet.
The separate 17,822-cell Fluent elbow converges at unchanged 350/550-step stops;
`v6_field_comparison.json` compares every corresponding cell after 1,200 steps,
without resampling or pressure-gauge alignment. The CSV field arrays remain in
the ignored local `build/performance/plateau-*` directories; their hashes are
recorded in the comparison.

The complete slow tier intentionally returns nonzero: six retained tet-force,
29 V6 and 15 V5 failed assertions. V5 previously retained 16; the faster CPU fine
bump reaches 2,300 rather than 2,200 steps within its unchanged 600-second budget,
and its later transient endpoint passes one old Cf gate. No steady-convergence
claim follows. All 308 V1-V4 numerical/progress lines match the retained main
baseline byte for byte. Gate sources and thresholds are unchanged.

Byte counts are the profiler's full-buffer-argument model, not measured DRAM
traffic. The available Metal counter set contains only GPUTimestamp. Allocation
counts cover public Device requests, not host heap or backend/driver allocations.
CUDA 13.0.2 compiles all four CI architectures and passes the CPU quick tier;
CUDA GPU execution and Windows/RTX performance are unverified.
