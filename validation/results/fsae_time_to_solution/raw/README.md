# Evidence index

The report and `../analyze.py` use only the final `measurements.jsonl`,
`convergence.jsonl`, `profile-quiet-*.log` and `counts-quiet-*.log` series.
Their baseline is 85c489b and numerical candidate is 7b817f3. JSON receipts
bind runs to source identities, executable hashes and mesh hashes.

`parity.jsonl`, `overhead-parity-*`, `zero-oracle.log`, `quick-candidate.log`,
the compact slow assertion receipts and `slow-comparison.json` are the final
numerical gates. `cuda-7b817f3-build.log` and `cuda-generation.json` establish
CUDA compilation, not execution. `profilers.json` describes the temporary
GPU timestamp and native-count probes; no instrumentation is retained in
solver source.

Development receipts are retained to distinguish experiments from final
measurements. `development/` contains the interrupted, contended timing logs
and earlier cold profiles. `edge-only-measurements.jsonl` and
`contended-measurements-*.jsonl` are incomplete series, not speedup evidence.
`cache-pilot-*`, `block-pilot-*`, `tiled-pilot-*`, `sequence-pilot-*`,
`tuning.jsonl`, `tune-*` and `initial-*` are development probes. The numerical
`RESULT` source tag identifies each version. The tiles and threshold probes
are based on cc7366c; their source basis and completed pairs are in
`trials.json`. None replaces a missing observation in the final series.

`development/ea3c058-verification/` preserves the complete passing gates for
the earlier cached candidate. `development/ea3c058-profiles/` holds its
interrupted profile series, which exposed the remaining readback costs.
`refresh-pilot*` records the rejected combined refresh/solve API. Prototype
`RESULT` tags show the configured source basis; metadata's patch hash and
binary hash bind the actual uncommitted experiment. Those tags do not claim
the prototypes are pristine builds of their basis commit.

Full doctest XML is ignored because passing assertions embed generated
meshes. Compact XML receipts preserve identities and hash every expanded
value and INFO context. Original XML hashes bind the receipts to the local
full run. No generated mesh or restricted reference dataset is published.
Home directory prefixes in text logs are replaced by `<checkout>` or
`<fixtures>`.

`converged-field-parity.json` verifies complete exported Ahmed L1 SIMPLEC cell
and wall fields at the actual 1,100-iteration convergence point. It is a parity
check, not an extra timing observation. Full fields remain in ignored local runs.
