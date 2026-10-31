#!/usr/bin/env python3
"""Markdown tables from tools/memory_scale/phases.py records.

For each case: bytes per cell by phase, before and after, as the phase peak
(highest footprint since the previous phase) and the device bytes allocated
at the phase end. Repeated runs report their maximum.

  python3 tools/memory_scale/tables.py validation/results/memory_scale/raw/phases.jsonl.gz
"""
import collections
import gzip
import json
import sys


def main():
    path = sys.argv[1]
    lines = (gzip.open(path, "rt") if path.endswith(".gz") else open(path)).read().splitlines()
    runs = [json.loads(line) for line in lines if line.strip()]
    by_case = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in runs:
        if r["exit"] == 0:
            by_case[(r["case"], r["cells"])][r["label"]].append(r)
    for (case, cells), labels in sorted(by_case.items(), key=lambda kv: kv[0][1]):
        names = [n for n in ("before", "before-returned", "after") if n in labels]
        print(f"\n### {case}, {cells:,} cells\n")
        header = " | ".join(f"{n} peak B/cell" for n in names)
        print(f"| phase | {header} | device B/cell (after) |")
        print("| --- |" + " ---: |" * (len(names) + 1))
        phases = [p["phase"] for p in labels[names[-1]][0]["phases"]]
        for phase in phases:
            row = []
            for n in names:
                values = [p["peak_per_cell"] for r in labels[n] for p in r["phases"] if p["phase"] == phase]
                row.append(f"{max(values):,.0f}" if values else "")
            device = [p["device_per_cell"] for r in labels[names[-1]] for p in r["phases"] if p["phase"] == phase]
            print(f"| {phase} | " + " | ".join(row) + f" | {max(device):,.0f} |")
        totals = " | ".join(
            f"**{max(r['peak_footprint_per_cell'] for r in labels[n]):,.0f}** "
            f"({max(r['peak_footprint'] for r in labels[n]) / 1e9:.2f} GB, {len(labels[n])} runs)" for n in names)
        print(f"| **process peak** | {totals} | |")


if __name__ == "__main__":
    main()
