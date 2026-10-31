#!/usr/bin/env python3
"""Per-phase memory of one `kutta solve` run (KUTTA_MEMORY_PHASES probe).

Runs the solve under /usr/bin/time -l (macOS), parses the probe's
`kutta-memory` lines and the process peak footprint, and appends one JSON
record to --jsonl with bytes and bytes per cell for every phase.

  python3 tools/memory_scale/phases.py --kutta build/kutta --label candidate \\
      --case runs/mesh/ahmed35_L1.cas --cells 403604 --iterations 5 --jsonl out.jsonl
"""
import argparse
import json
import os
import pathlib
import re
import subprocess
import tempfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kutta", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--cells", type=int, required=True)
    ap.add_argument("--iterations", type=int, default=5)
    ap.add_argument("--backend", default="metal")
    ap.add_argument("--approximate", default="wfgc")
    ap.add_argument("--jsonl", required=True)
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as out:
        cmd = ["/usr/bin/time", "-l", a.kutta, "solve", a.case, "--backend", a.backend, "--iterations",
               str(a.iterations), "--approximate", a.approximate, "--out", out]
        env = dict(os.environ, KUTTA_MEMORY_PHASES="1")
        proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, env=env)
    text = proc.stderr
    phases = []
    for line in text.splitlines():
        if line.startswith("kutta-memory "):
            fields = line.split()
            values = dict(f.split("=", 1) for f in fields[2:])
            record = {"phase": fields[1], "seconds": float(values["seconds"])}
            for key in ("current", "peak", "heap", "device"):
                record[key] = int(values[key])
                record[key + "_per_cell"] = round(int(values[key]) / a.cells, 1)
            phases.append(record)
    peak = re.search(r"(\d+)\s+peak memory footprint", text)
    rss = re.search(r"(\d+)\s+maximum resident set size", text)
    real = re.search(r"([\d.]+) real", text)
    result = {"label": a.label, "case": pathlib.Path(a.case).name, "cells": a.cells, "backend": a.backend,
              "iterations": a.iterations, "exit": proc.returncode,
              "peak_footprint": int(peak.group(1)) if peak else None,
              "peak_footprint_per_cell": round(int(peak.group(1)) / a.cells, 1) if peak else None,
              "max_rss": int(rss.group(1)) if rss else None,
              "wall_seconds": float(real.group(1)) if real else None, "phases": phases}
    if proc.returncode:
        result["error"] = text[-2000:]
    with open(a.jsonl, "a") as f:
        f.write(json.dumps(result) + "\n")
    print(f"{a.label} {result['case']} exit={proc.returncode} peak={result['peak_footprint']} "
          f"({result['peak_footprint_per_cell']} B/cell) wall={result['wall_seconds']}s")


if __name__ == "__main__":
    main()
