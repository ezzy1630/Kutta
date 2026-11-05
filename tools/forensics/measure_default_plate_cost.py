#!/usr/bin/env python3
"""Time default-option effects where SIMPLEC controls are active.

Replay an existing checked bench_convergence cost receipt on the coarse NASA
plate, with its SST and SA models and ordinary solver controls. Every run
keeps the original 3,000-step budget and convergence gates. Baseline and
candidate alternate AB, BA; candidate early convergence is retained, not
replaced with extra iterations. All files belong in an ignored build folder.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys

from check_solver_controls import check_controls


CONTROL_FIELDS = (
    "compensated_momentum", "paired_state_slack", "momentum_absolute",
    "pressure_absolute", "scalar_compensation", "scalar_pseudo_time",
    "scalar_steady_tolerance", "scalar_absolute", "dt",
    "momentum_relaxation", "pressure_relaxation",
)
BUDGET = 3000
ORDINARY_CONTROLS = {
    "compensated_momentum": 0.0, "paired_state_slack": 1.0,
    "momentum_absolute": 1e-10, "pressure_absolute": 1e-12,
    "scalar_compensation": 0.0, "scalar_pseudo_time": 1.0,
    "scalar_steady_tolerance": 0.0, "scalar_absolute": 1e-12, "dt": 0.0,
}
VARIANTS = {
    "momentum-compensation": {"compensated_momentum": 1.0},
    "state-slack-16": {"paired_state_slack": 16.0},
    "momentum-absolute": {"momentum_absolute": 1e-20},
    "pressure-absolute": {"pressure_absolute": 1e-20},
    "scalar-compensation": {"scalar_compensation": 1.0},
    "scalar-no-pseudotime": {"scalar_pseudo_time": 0.0},
    "scalar-steady": {"scalar_steady_tolerance": 1e-6},
    "scalar-absolute": {"scalar_absolute": 1e-20},
    "fixed-flow-time": {"dt": 0.1},
}
VARIANTS["inner-absolute-group"] = {
    **VARIANTS["momentum-absolute"], **VARIANTS["pressure-absolute"],
    **VARIANTS["scalar-absolute"],
}
VARIANTS["scalar-group"] = {
    **VARIANTS["scalar-compensation"], **VARIANTS["scalar-no-pseudotime"],
    **VARIANTS["scalar-steady"], **VARIANTS["scalar-absolute"],
}
VARIANTS["precision-group"] = {
    **VARIANTS["momentum-compensation"], **VARIANTS["state-slack-16"],
    **VARIANTS["inner-absolute-group"], **VARIANTS["scalar-group"],
}
# Qualify the six individually eligible defaults together.
VARIANTS["eligible-defaults"] = {
    **VARIANTS["inner-absolute-group"], **VARIANTS["scalar-group"],
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(lines, prefix):
    return [dict(word.split("=", 1) for word in line.split()[1:])
            for line in lines if line.startswith(prefix + " ")]


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--data", type=Path, default=Path("validation/data/nasa_tmr"))
    parser.add_argument("--pairs", type=int, default=2)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    if sys.platform != "darwin":
        parser.error("physical-footprint measurement currently requires macOS")
    if not 1 <= args.pairs <= 4:
        parser.error("--pairs must be 1..4")
    source = json.loads(args.receipt.read_text())
    if source.get("settings_checked") is not True:
        parser.error("receipt must have checked every emitted solver control")
    if source.get("variant") not in VARIANTS:
        parser.error("receipt must name one of the frozen default-study variants")
    binaries = {"baseline": args.baseline.resolve(), "candidate": args.candidate.resolve()}
    hashes = {name: sha(path) for name, path in binaries.items()}
    for name in binaries:
        if hashes[name] != source[name + "_sha256"]:
            parser.error(f"{name} binary does not match the checked receipt")
    expected = {}
    for row in source["late_iteration_cost"]:
        key = (row["model"], row["backend"], row["variant"])
        if key in expected or row["failure"] is not None or row["exit"] not in (0, 2):
            parser.error("receipt has duplicate or failed coarse-plate observations")
        fields = records([row["settings"]], "CASE")
        if len(fields) != 1 or fields[0]["name"] != "plate" or fields[0]["scheme"] != "simplec":
            parser.error("receipt must describe the ordinary SIMPLEC plate protocol")
        case = fields[0]
        if case["model"] != key[0] or case["backend"] != key[1] or int(case["cells"]) != 816:
            parser.error("receipt case identity does not match its row")
        controls = {field: float(case[field]) for field in CONTROL_FIELDS}
        intended = {**ORDINARY_CONTROLS, "momentum_relaxation": 0.7,
                    "pressure_relaxation": 1.0}
        if row["variant"] == "candidate":
            intended.update(VARIANTS[source["variant"]])
        if controls != intended:
            parser.error(f"{key} does not match the frozen baseline and named option")
        expected[key] = {"controls": controls, "monitor": row["last_monitor"],
                         "work": row["cumulative_inner_iterations"], "exit": row["exit"]}
    required = {(model, backend, name) for model in ("sst", "sa")
                for backend in ("cpu", "metal") for name in binaries}
    if expected.keys() != required:
        parser.error("receipt must cover both models and both backends for each binary")
    mesh = args.data / "2DZP_flatplate/grids/flatplate_clust2_4levelsdown_35x25.p2dfmt"
    frozen = {str(path.resolve()): sha(path) for path in (*binaries.values(), args.receipt, mesh)}
    if args.validate_only:
        print("Verified binary hashes, all eight checked control profiles and mesh identity")
        return 0
    args.directory.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, MallocSpaceEfficient="1")
    env.pop("KUTTA_MEMORY_PHASES", None)
    result = {"variant": source["variant"], "complete": False, "pairs": args.pairs,
              "order": "AB, BA", "outer_budget": BUDGET, "frozen_inputs_sha256": frozen,
              "protocol": "Coarse NASA plate, SST/SA, ordinary SIMPLEC controls; solver-loop time. "
                          "Early convergence is retained; both status and iteration count are reported. "
                          "Late per-iteration timing is available only when both runs reach step3000.",
              "contended": True, "rows": [], "summary": []}
    progress = args.directory / "cost.json"
    save(progress, result)
    for model in ("sst", "sa"):
        for backend in ("cpu", "metal"):
            for pair in range(args.pairs):
                order = ("baseline", "candidate") if pair % 2 == 0 else ("candidate", "baseline")
                for name in order:
                    for path, digest in frozen.items():
                        if sha(Path(path)) != digest:
                            raise ValueError(f"frozen input changed: {path}")
                    reference = expected[(model, backend, name)]
                    command = [str(binaries[name]), "--case", "plate", "--model", model,
                               "--backend", backend, "--data", str(args.data), "--level", "0",
                               "--iterations", str(BUDGET), "--interval", "100",
                               "--max-seconds", "1000000000"]
                    override = source.get("candidate_compensation_override")
                    if name == "candidate" and override is not None:
                        command += ["--compensated-momentum", str(override)]
                    log = args.directory / f"{model}-{backend}-{pair}-{name}.log"
                    with log.open("w") as output:
                        run = subprocess.run(["/usr/bin/time", "-l", *command], env=env,
                                             stdout=output, stderr=subprocess.STDOUT)
                    lines = log.read_text().splitlines()
                    check_controls(lines, reference["controls"])
                    cases = records(lines, "CASE")
                    if any(cases[0][key] != value for key, value in
                           {"name": "plate", "model": model, "backend": backend,
                            "scheme": "simplec", "cells": "816"}.items()):
                        raise ValueError("replayed case identity changed")
                    steps, work, memory = (records(lines, prefix) for prefix in ("STEP", "WORK", "MEMORY"))
                    final = dict(steps[-1])
                    seconds = float(final.pop("seconds"))
                    if run.returncode != reference["exit"] or final != reference["monitor"]:
                        raise ValueError(f"numerical endpoint changed in {log}")
                    if work[-1] != reference["work"]:
                        raise ValueError(f"cumulative inner work changed in {log}")
                    times = {int(row["it"]): float(row["seconds"]) for row in steps}
                    late = (times[BUDGET] - times[1500]) / 1500 if BUDGET in times and 1500 in times else None
                    peak = int(next(line.split()[0] for line in lines if "peak memory footprint" in line))
                    if peak > 3_000_000_000:
                        raise ValueError(f"resource guidance exceeded in {log}: {peak}")
                    row = {"model": model, "backend": backend, "pair": pair, "variant": name,
                           "solver_seconds_at_stop": seconds, "outer_iterations": int(final["it"]),
                           "converged": final["converged"] == "1", "maximum_defect": float(final["defect"]),
                           "late_seconds_per_iteration_1501_3000": late, "inner_work": work[-1],
                           "device_bytes_at_stop": int(memory[-1]["device"]), "process_peak_bytes": peak,
                           "controls_verified": True, "numerical_endpoint_and_work_reproduced": True}
                    result["rows"].append(row)
                    save(progress, result)
                    print(model, backend, pair, name, "steps=" + final["it"], "seconds=" + str(seconds), flush=True)
            rows = [row for row in result["rows"] if row["model"] == model and row["backend"] == backend]
            ratios, late_ratios = [], []
            for pair in range(args.pairs):
                selected = {row["variant"]: row for row in rows if row["pair"] == pair}
                before, after = selected["baseline"], selected["candidate"]
                ratios.append(after["solver_seconds_at_stop"] / before["solver_seconds_at_stop"])
                if before["late_seconds_per_iteration_1501_3000"] is not None and after["late_seconds_per_iteration_1501_3000"] is not None:
                    late_ratios.append(after["late_seconds_per_iteration_1501_3000"] / before["late_seconds_per_iteration_1501_3000"])
            result["summary"].append({"model": model, "backend": backend,
                "paired_solver_time_to_stop_ratio_median": statistics.median(ratios),
                "paired_solver_time_to_stop_ratio_range": [min(ratios), max(ratios)],
                "paired_late_time_per_iteration_ratio_median": statistics.median(late_ratios) if late_ratios else None,
                "paired_late_time_per_iteration_ratio_range": [min(late_ratios), max(late_ratios)] if late_ratios else None,
                "baseline_outer_iterations": rows[0]["outer_iterations"],
                "candidate_outer_iterations": next(row["outer_iterations"] for row in rows if row["variant"] == "candidate"),
                "baseline_converged": rows[0]["converged"],
                "candidate_converged": next(row["converged"] for row in rows if row["variant"] == "candidate")})
            save(progress, result)
    for path, digest in frozen.items():
        if sha(Path(path)) != digest:
            raise ValueError(f"frozen input changed: {path}")
    result["complete"] = True
    save(progress, result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
