"""Summarize bench_convergence logs from run_comparison.sh as a Markdown table.

Usage: python3 summarize.py <log directory>

Elbow: first STEP with bench's V6 criterion (converged=1). TMR: first STEP
whose transport defect max(u, v, k, omega, nu~) <= 1e-5 with the 100-iteration
load change <= 0.1% ("transport+load"); the full bench gate also needs scaled
continuity <= 1e-5, reported as gate=1 when it is reached. "Loads settled" is
the first STEP after which that load change stays <= 0.1% to the end of the
run. Accuracy is the relative difference of Coupled's final drag and cf0 from
SIMPLEC's; it is a fixed-point comparison only when both runs settled. A run
whose final raw continuity is over 10x its lowest is reported as diverging even
when it did not throw (a finite blown-up state can be stationary).
"""
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1])


def fields(line):
    return dict(m.groups() for m in re.finditer(r"(\w+)=(\S+)", line))


def parse(path):
    run = {"steps": [], "failure": None, "rss": None, "case": {}}
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("CASE "):
            run["case"] = fields(line)
        elif line.startswith("STEP "):
            run["steps"].append({k: v for k, v in fields(line).items()})
        elif line.startswith("FAILURE "):
            run["failure"] = fields(line)
        elif "maximum resident set size" in line:
            run["rss"] = int(line.split()[0])
    return run


def first(run, predicate):
    for step in run["steps"]:
        if predicate(step):
            return step
    return None


def transport(step):
    keys = ("u", "v", "k", "omega", "sa")
    return max(float(step[k]) for k in keys) <= 1e-5 and float(step["change"]) <= 1e-3


def bench_gate(step):
    return step["converged"] == "1"


def blown_up(run):
    raws = [float(step["raw"]) for step in run["steps"]]
    return bool(raws) and raws[-1] > 10 * min(raws)


def settled(run):
    if run["failure"] or not run["steps"] or blown_up(run):
        return "-"
    hit = None
    for step in run["steps"]:
        if float(step["change"]) <= 1e-3:
            hit = hit or step
        else:
            hit = None
    return f"{hit['it']} / {float(hit['seconds']):.1f} s" if hit else "no"


def describe(run, elbow):
    if run["failure"]:
        f = run["failure"]
        return f"diverged at {f['completed_iterations']} ({float(f['seconds']):.1f} s)", None
    if blown_up(run):
        last = run["steps"][-1]
        return f"diverging (raw {float(last['raw']):.1e} at {last['it']})", None
    hit = first(run, bench_gate if elbow else transport)
    last = run["steps"][-1] if run["steps"] else None
    if hit is None:
        return f"not reached in {last['it']} ({float(last['seconds']):.1f} s)" if last else "no output", None
    gate = "" if elbow else (" gate=1" if first(run, bench_gate) else "")
    return f"{hit['it']} / {float(hit['seconds']):.2f} s (raw {float(hit['raw']):.2e}){gate}", hit


rows = []
for simplec in sorted(root.glob("*_simplec.log")):
    stem = simplec.name[: -len("_simplec.log")]
    coupled = root / f"{stem}_coupled.log"
    if not coupled.exists():
        continue
    name, backend = stem.rsplit("_", 1)
    a, b = parse(simplec), parse(coupled)
    elbow = name == "elbow"
    text_a, hit_a = describe(a, elbow)
    text_b, hit_b = describe(b, elbow)
    speedup = ""
    if hit_a and hit_b:
        speedup = f"{float(hit_a['seconds']) / float(hit_b['seconds']):.2f}x"
    comparable = (hit_a and hit_b) if elbow else all(settled(r) not in ("-", "no") for r in (a, b))
    if not comparable:
        accuracy = "n/a: not both settled"
    else:
        fa, fb = a["steps"][-1], b["steps"][-1]
        diffs = []
        for key in ("drag", "cf0"):
            ref = float(fa[key])
            if ref != 0:
                diffs.append(f"{key} {abs(float(fb[key]) - ref) / abs(ref):.1e}")
        accuracy = ", ".join(diffs)
    cells = b["case"].get("cells", "")
    per_cell = b["case"].get("coupled_bytes_per_cell", "")
    rss = ""
    if a["rss"] and b["rss"]:
        rss = f"{a['rss'] / 2**20:.0f} / {b['rss'] / 2**20:.0f}"
    load = "" if elbow else f"{settled(a)} \u2192 {settled(b)}"
    rows.append((name, backend, cells, text_a, text_b, speedup, load, accuracy,
                 f"{float(per_cell):.0f}" if per_cell else "", rss))

print("| case | backend | cells | SIMPLEC its / s | Coupled its / s | speedup | loads settled S \u2192 C | Coupled vs SIMPLEC (final) | Coupled B/cell | peak RSS MiB S / C |")
print("|---|---|---|---|---|---|---|---|---|---|")
for row in rows:
    print("| " + " | ".join(row) + " |")
