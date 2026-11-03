# 2026-11-03 evidence

Read the dated section in [REPORT.md](../REPORT.md) first. These are aggregate
solver outputs and reproduction sources; restricted meshes, references, field
exports and slow-tier XML are deliberately absent.

`metadata.json` binds main `e014e36` and the two source commits to frozen native
binaries and input meshes. `numerical-source.patch.gz` is the exact measured source
difference (gzip-compressed without changing the patch bytes). Use `gzip -cd`
to inspect or apply either archived patch. The binaries were built before the commits, so their embedded commit
string remains `e014e36`; use the recorded binary hashes and source patch.
The native executables remain under `<checkout>/build/speed/{baseline,candidate,candidate-final}`.
The intermediate `candidate` has only the face-jump change.

`convergence.jsonl` contains all 16 complete native runs. `convergence-*.log`
contains checkpoint history and BSD `time -l` memory/process measurements.
`convergence-summary.json` compares every printed checkpoint and reported final
numerical value within each pair. The original existing `../reproduce.py` started
the series; `converge_resume.py` resumes the same artifact hashes and settings,
retains complete receipts and streams subsequent logs. Its hard-coded output is
`build/speed/convergence`; do not overwrite retained measurements.

`profiles.json`, `kernels.csv` and `phases.csv` summarize two alternating early
profile pairs per case/scheme. `late-profiles.json` and `late-kernels.csv` summarize
one diagnostic wing SIMPLEC pair near convergence. Rates are buffer argument
footprints divided by intrusive timestamp spans, **not measured DRAM bandwidth**.
The phase families are defined by the existing `../analyze.py`, with native replay
kept as one measured mixed kernel. Both profile binaries contain the removed
`profiling.patch.gz` in addition to the corresponding numerical source.

To repeat profiling after building matching frozen artifacts and setting
`KUTTA_RUNS` to the read-only fixture checkout's `runs` directory:

```sh
python3 validation/results/fsae_time_to_solution/2026-11-03/profile_pairs.py --out build/speed/new-profiles --pairs 2
python3 validation/results/fsae_time_to_solution/2026-11-03/summarize_profiles.py --root build/speed/new-profiles --out build/speed/new-profiles.json
python3 validation/results/fsae_time_to_solution/2026-11-03/profile_pairs.py --out build/speed/new-late --pairs 1 --first 2350 --cases wing-L2 --schemes simplec
python3 validation/results/fsae_time_to_solution/2026-11-03/setup_pairs.py --out build/speed/new-setup --pairs 3
```

Profiling uses `build/speed/profile-baseline` and `profile-candidate`, not the
native timing binaries. `setup_pairs.py` compares the face-jump-only candidate
with the final candidate to isolate the pivot scan. A one-iteration setup probe
is diagnostic and does not change convergence acceptance.

`micro_jump.cpp` and `micro_inverse.cpp` are the exact standalone component
experiment sources, compiled with the project's Release flags and library.
`micro-summary.json` reduces all five alternating pairs; the original logs retain
every sample. The inverse benchmark consumes all 100 outputs and checks every
matrix entry; an earlier single-output pilot is excluded.

`parity-cli.jsonl` and `parity-wing.jsonl` record file comparison receipts for the
40-iteration CPU/Metal runs. `overhead-parity.json` adds four solver-state hashes
per backend. `verification.json` records quick-tier and direct-oracle counts,
artifact hashes and unverified scope. It does not claim a full slow-tier run.
`fixture-cleanup.json` records removal of all 268 temporary fixture symlinks.

Build configuration: Release, C++20, Metal/HDF5/tests enabled, CUDA and Qt app
disabled. Use the repository's CMake targets `kutta`, `kutta-tests`,
`bench_fsae_aero`, `bench_fv` and `bench_overhead`. CUDA generation is checked;
CUDA compilation and execution require another machine.
