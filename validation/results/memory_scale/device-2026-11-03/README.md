# Device-storage evidence, 2026-11-03

The parent REPORT.md explains the measurements and their limits.

- `baseline-*` is e014e36 with allocation diagnostics; `candidate-*` is the
  solver source committed as e63e157. Exact commits and binary hashes are in
  `manifest.json`. The frozen binaries retain their configure-time e014e36
  stamp; the stamp in benchmark output is not the candidate source identity.
- Each `*-buffers.csv.gz` has one row per live, nonempty solver-owned device
  allocation after the inventory iterations. Bytes are exact; the ratios are
  rounded to 12 significant digits. Source/column plus instance distinguishes
  repeated allocations. Backend reduction/command storage is excluded.
- The matching `*-bench.log.gz` records setup and iteration snapshots, mesh
  sizes, totals, host heap and process footprints. `inventory-summary.json`
  contains exact payload sums. Hex is Ahmed L1, poly is V6, native is Ahmed
  tet/prism. Both SIMPLEC and Coupled use plain SST and default paired state.
- `phases.jsonl` contains the eight full CLI scale runs, including process
  peaks and live device payload at each phase. It is not a device high-water
  counter. All child exits are zero; L3 ran once per binary, serially.
- `input-hashes.json` identifies the actual input bytes, including the native
  40 m/s air/SST case. Meshes, case wrappers and computed flow outputs are not
  included here.
- `parity.jsonl` covers ten 40-iteration case/backend pairs. Four full files
  pass cmp; residuals pass grep/cut/cmp. Hashes use that exact normalization.
- `checks.json` records the quick tier, focused storage test, build, lint and
  the unexecuted slow/CUDA scope. No slow XML or restricted reference values
  are included.

CSV generation uses tools/memory_scale/inventory.py. The allocating source is
read from e014e36 for baseline logs and from e63e157 for candidate logs.
