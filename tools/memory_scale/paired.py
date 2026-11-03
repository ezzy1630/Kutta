#!/usr/bin/env python3
"""Alternating paired CLI timing, process peak and exact-output checks (macOS).

The timed sample is the final residual CSV batch of a fixed-iteration solve.
Setup, export and process peak are recorded separately. Runs are sequential;
other machine users can still contend. No quiet-machine speed claim is made.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile

from parity import FILES, digest


def measure(binary, case, backend, iterations, approximate):
    with tempfile.TemporaryDirectory() as output:
        command = ["/usr/bin/time", "-l", str(binary), "solve", str(case), "--backend", backend,
                   "--iterations", str(iterations), "--out", output]
        if approximate:
            command += ["--approximate", approximate]
        result = subprocess.run(command, env=dict(os.environ, KUTTA_MEMORY_PHASES="1"),
                                capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stderr[-4000:])
        directory = Path(output)
        actual = json.loads((directory / "loads.json").read_text())["iterations"]
        metadata = json.loads((directory / "summary.json").read_text())
        if actual != iterations:
            raise RuntimeError(f"completed {actual} iterations, expected {iterations}")
        with (directory / "residuals.csv").open() as stream:
            residuals = list(csv.DictReader(line for line in stream if not line.startswith("#")))
        phases = []
        for line in result.stderr.splitlines():
            if not line.startswith("kutta-memory "):
                continue
            words = line.split()
            fields = dict(word.split("=", 1) for word in words[2:])
            phases.append(dict(phase=words[1], **{key: float(value) if key == "seconds" else int(value)
                                                for key, value in fields.items()}))
        match = re.search(r"(\d+)\s+peak memory footprint", result.stderr)
        if not match or not phases or not residuals:
            raise RuntimeError("missing process/phase/timing measurements")
        return {"iterations": actual, "cells": metadata["cells"], "embedded_commit": metadata["kutta_commit"],
                "seconds_per_iteration": float(residuals[-1]["seconds_per_iteration"]),
                "timed_batch_start": int(residuals[-2]["iteration"]) + 1 if len(residuals) > 1 else 1,
                "timed_batch_end": int(residuals[-1]["iteration"]), "peak_bytes": int(match.group(1)),
                "phases": phases, "sha256": {name: digest(directory / name) for name in FILES}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--case", required=True, type=Path)
    parser.add_argument("--cells", required=True, type=int)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--pairs", type=int, default=3)
    parser.add_argument("--approximate", default="wfgc")
    parser.add_argument("--backend", default="metal")
    parser.add_argument("--peak-limit-gb", type=float, default=3, help="stop after an oversized run; 0 disables")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if (args.pairs < 1 or args.iterations < 1 or args.cells < 1 or
            not math.isfinite(args.peak_limit_gb) or args.peak_limit_gb < 0):
        parser.error("pairs, iterations and cells must be positive; peak limit must be nonnegative")
    args.out.mkdir(parents=True, exist_ok=True)
    if (args.out / "runs.jsonl").exists():
        parser.error("use a fresh output directory; runs.jsonl already exists")
    runs, ratios = [], []
    reference = None
    for pair in range(args.pairs):
        measured = {}
        order = ("baseline", "candidate") if pair % 2 == 0 else ("candidate", "baseline")
        for label in order:
            record = measure(getattr(args, label), args.case, args.backend, args.iterations, args.approximate)
            if record["cells"] != args.cells:
                raise RuntimeError(f"mesh has {record['cells']} cells, expected {args.cells}")
            if (record["timed_batch_end"] != args.iterations or
                    not math.isfinite(record["seconds_per_iteration"]) or record["seconds_per_iteration"] <= 0):
                raise RuntimeError("missing final timing batch")
            record.update(label=label, pair=pair, case=args.case.name, cells=args.cells, backend=args.backend)
            record["binary_sha256"] = hashlib.sha256(getattr(args, label).read_bytes()).hexdigest()
            if reference is None:
                reference = record["sha256"]
            if record["sha256"] != reference:
                raise RuntimeError(f"output mismatch in {label} pair {pair}")
            record["case_sha256"] = hashlib.sha256(args.case.read_bytes()).hexdigest()
            measured[label] = record["seconds_per_iteration"]
            runs.append(record)
            with (args.out / "runs.jsonl").open("a") as stream:
                stream.write(json.dumps(record) + "\n")
            print(f"pair={pair} {label} {record['seconds_per_iteration']:.6f} s/it "
                  f"peak={record['peak_bytes']/1e9:.3f} GB", flush=True)
            if args.peak_limit_gb and record["peak_bytes"] > args.peak_limit_gb * 1e9:
                raise RuntimeError("process peak exceeded --peak-limit-gb; stopped before the next solve")
        ratios.append(measured["candidate"] / measured["baseline"])
    summary = {"case": args.case.name, "cells": args.cells, "backend": args.backend, "pairs": args.pairs,
               "iterations": args.iterations, "paired_ratios": ratios, "median_paired_ratio": statistics.median(ratios)}
    for label in ("baseline", "candidate"):
        samples = [record for record in runs if record["label"] == label]
        times = [record["seconds_per_iteration"] for record in samples]
        summary[label] = {"seconds_per_iteration_median": statistics.median(times),
                          "seconds_per_iteration_min": min(times), "seconds_per_iteration_max": max(times),
                          "peak_bytes_max": max(record["peak_bytes"] for record in samples),
                          "device_bytes_max": max(phase["device"] for record in samples for phase in record["phases"])}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
