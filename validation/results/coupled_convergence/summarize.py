"""Summarize measure.sh output as Markdown tables.

Usage: python3 summarize.py <measure outdir>

fluent: the `done` event of `kutta solve --progress jsonl` with the case's
Fluent residual criteria (iterations, solve seconds). strict: the first
bench_convergence STEP with converged=1 (V6 criterion). Iterations are
deterministic per build and backend; seconds are given as the range over
repetitions, since the CPU and GPU were shared.
"""
import collections
import json
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1])
runs = collections.defaultdict(lambda: {"iterations": set(), "seconds": [], "bytes": None, "converged": set()})


def key(path, prefix):
    name, backend, side, _ = path.stem[len(prefix):].split("_")
    return name, backend, side


for path in sorted(root.glob("fluent_*.jsonl")):
    done = [json.loads(line) for line in path.read_text().splitlines() if '"event":"done"' in line]
    run = runs[("fluent",) + key(path, "fluent_")]
    if not done:
        run["converged"].add("error")
        continue
    run["converged"].add(done[-1]["converged"])
    run["iterations"].add(done[-1]["iterations"])
    run["seconds"].append(done[-1]["elapsed"])

for path in sorted(root.glob("strict_*.log")):
    run = runs[("strict",) + key(path, "strict_")]
    hit = None
    for line in path.read_text(errors="replace").splitlines():
        fields = dict(re.findall(r"(\w+)=(\S+)", line))
        if line.startswith("CASE "):
            run["bytes"] = float(fields["coupled_bytes_per_cell"])
        if line.startswith("STEP ") and fields.get("converged") == "1" and hit is None:
            hit = (int(fields["it"]), float(fields["seconds"]))
    run["converged"].add(hit is not None)
    if hit:
        run["iterations"].add(hit[0])
        run["seconds"].append(hit[1])


def cell(run):
    if run["converged"] != {True}:
        return "not reached"
    its = "/".join(str(i) for i in sorted(run["iterations"]))
    lo, hi = min(run["seconds"]), max(run["seconds"])
    return f"{its} / {lo:.1f}–{hi:.1f} s"


for criterion, title in (("fluent", "Fluent residual criteria (kutta solve)"), ("strict", "V6 criterion (bench_convergence)")):
    print(f"\n### {title}\n")
    print("| case | backend | Coupled before | Coupled after | SIMPLEC |")
    print("|---|---|---|---|---|")
    for name in ("elbow", "nozzle"):
        for backend in ("cpu", "metal"):
            row = [cell(runs[(criterion, name, backend, side)]) for side in ("base", "cand", "simplec")]
            print(f"| {name} | {backend} | " + " | ".join(row) + " |")

print("\n### Coupled solver bytes per cell (CoupledSolver::bytes)\n")
print("| case | before | after |")
print("|---|---|---|")
for name in ("elbow", "nozzle"):
    b, a = runs[("strict", name, "cpu", "base")]["bytes"], runs[("strict", name, "cpu", "cand")]["bytes"]
    print(f"| {name} | {b:.1f} | {a:.1f} |")
