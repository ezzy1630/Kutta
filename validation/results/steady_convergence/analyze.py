#!/usr/bin/env python3
"""Read the retained logs, keeping failed/censored solves separate from convergence."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import xml.etree.ElementTree as ET


def read(path):
    if path.exists():
        return path.read_text()
    return gzip.decompress(Path(str(path) + ".gz").read_bytes()).decode()


def fields(line):
    result = {}
    for key, value in re.findall(r"(\w+)=([^\s]*)", line):
        try:
            value = float(value)
            if not math.isfinite(value):
                value = None
        except ValueError:
            pass
        if key in ("cf", "cp", "Cf", "Cp"):
            result.setdefault(key, []).append(value)
        else:
            result[key] = value
    return result


def steps(path):
    return [fields(line) for line in read(path).splitlines() if line.startswith("STEP ")]


def assertions(path):
    tree = ET.fromstring(read(path))
    result = {}
    totals = {}
    for case in tree.iter("TestCase"):
        if case.get("skipped") == "true":
            continue
        counts = collections.Counter()
        case_name = case.get("name")
        for expression in case.iter("Expression"):
            file = expression.get("filename").split("/tests/")[-1]
            original = " ".join(expression.findtext("Original", "").split())
            info = [" ".join(info.itertext()).strip() for info in expression.findall("Info")]
            # Optional files change the inventory loop's order. Match its
            # assertions by source-file identity, not a shifted loop position.
            scope = " | ".join(info) if case_name == "Fluent research downloads retain published byte identities" else ""
            base = (case_name, file, expression.get("line"), original, scope)
            index = counts[base]
            counts[base] += 1
            key = (*base, index)
            result[key] = {
                "success": expression.get("success") == "true",
                "expanded": " ".join(expression.findtext("Expanded", "").split()),
                "info": info,
            }
        totals[case_name] = {
            "assertions": sum(counts.values()),
            "failed": sum(not value["success"] for key, value in result.items() if key[0] == case_name),
            "exceptions": [" ".join(error.itertext()).strip() for error in case.iter("Exception")],
        }
    return result, totals


def is_v1_v6(name):
    return name.startswith(("FV V1 ", "FV V2 ", "FV V3 ", "FV V4 ", "FV V6 ", "V5 ",
                            "FV SIMPLEC develops planar Couette "))


def compare(before, after, destination):
    old, old_totals = assertions(before)
    new, new_totals = assertions(after)
    changes, rows = [], []
    for key in sorted(old.keys() | new.keys()):
        a, b = old.get(key), new.get(key)
        state = "unchanged" if a and b and a["success"] == b["success"] else "changed"
        row = {"test_case": key[0], "file": key[1], "line": key[2], "expression": key[3],
               "inventory_file": key[4], "occurrence": key[5], "before_pass": a["success"] if a else None,
               "after_pass": b["success"] if b else None,
               "before_value": a["expanded"] if a else None,
               "after_value": b["expanded"] if b else None,
               "before_context": " | ".join(a["info"]) if a else None,
               "after_context": " | ".join(b["info"]) if b else None, "state": state}
        rows.append(row)
        if state == "changed":
            changes.append(row)
    with (destination / "assertions.csv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    v1_v6 = lambda row: is_v1_v6(row["test_case"])
    scoped = {}
    for phase, totals in (("before", old_totals), ("after", new_totals)):
        selected = [case for name, case in totals.items() if is_v1_v6(name)]
        scoped[phase] = {"cases": len(selected), "assertions": sum(case["assertions"] for case in selected),
                         "failed": sum(case["failed"] for case in selected)}
    return {"before": old_totals, "after": new_totals, "changes": changes,
            "new_failures": [row for row in changes if row["before_pass"] and row["after_pass"] is False],
            "missing_assertions": [row for row in changes if row["after_pass"] is None],
            "added_assertions": [row for row in changes if row["before_pass"] is None],
            "v1_v6_totals": scoped,
            "v1_v6_changes": [row for row in changes if v1_v6(row)],
            "v1_v6_missing_assertions": [row for row in changes if v1_v6(row) and row["after_pass"] is None],
            "v1_v6_new_failures": [row for row in changes if v1_v6(row) and row["before_pass"] and row["after_pass"] is False]}


def spatial(root, stem):
    faces = list(csv.DictReader(read(root / (stem + "-faces.csv")).splitlines()))
    cells = list(csv.DictReader(read(root / (stem + "-cells.csv")).splitlines()))
    residual = [0.0] * len(cells)
    patches, regions = collections.defaultdict(float), collections.defaultdict(float)
    for face in faces:
        flux = float(face["flux"])
        residual[int(face["owner"])] += flux
        neighbour = int(face["neighbour"])
        if neighbour >= 0:
            residual[neighbour] -= flux
        else:
            patches[face["patch"]] += flux
    for cell, mass in zip(cells, residual):
        region = ("upstream" if float(cell["x"]) < 0 else "plate")
        region += "_low" if float(cell["y"]) < .01 else "_outer"
        regions[region] += abs(mass)
    worst = sorted(range(len(cells)), key=lambda index: abs(residual[index]), reverse=True)[:10]
    return {"sum_abs_double_from_stored_faces": sum(map(abs, residual)), "region_l1": regions,
            "patch_signed_mass": patches,
            "worst_cells": [{**cells[index], "mass_imbalance": residual[index]} for index in worst]}


def trace_summary(path):
    rows = list(csv.DictReader(read(path).splitlines()))
    result = {"samples": len(rows), "first_iteration": int(rows[0]["iteration"]),
              "last_iteration": int(rows[-1]["iteration"]),
              "last_raw": float(rows[-1]["continuity_raw_if_monitored"])}
    for key in ("cell_u", "cell_v", "cell_p", "delta_u", "delta_v", "delta_p",
                "max_delta_u", "max_delta_v", "max_delta_p", "pressure_rhs_l2"):
        values = [float(row[key]) for row in rows]
        result[key] = {"min": min(values), "max": max(values), "range": max(values) - min(values),
                       "mean": statistics.mean(values), "max_abs": max(map(abs, values))}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, default=Path(__file__).parent / "raw")
    parser.add_argument("--out", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    summary = {"full_task_accepted": False, "gate_protocol_changed": False, "runs": [], "cycles": {}}
    for phase in ("baseline", "final"):
        for line in read(args.raw / (phase + "-slow.log")).splitlines():
            if line.startswith("V5 plate ") or line.startswith("V5 bump measurement ") or line.startswith("V7 SA ") or line.startswith("V6_RESULT "):
                values = fields(line)
                values.update(phase=phase, record=line.split(" backend=")[0])
                summary["runs"].append(values)
    for phase in ("before", "after"):
        for backend in ("cpu", "metal"):
            name = f"channel-million-{phase}-{backend}"
            data = steps(args.raw / (name + ".log"))
            summary["runs"].append({"record": "million_cell_channel", "phase": phase, "backend": backend,
                                    "last_checkpoint": data[-1], "converged": bool(data[-1]["converged"]),
                                    "convergence_seconds": data[-1]["seconds"] if data[-1]["converged"] else None})
            name = f"sa-cycle-{phase}-{backend}"
            data = steps(args.raw / (name + ".log"))
            tail = [row for row in data if row["it"] >= 900]
            transport = [row for row in data if max(row["u"], row["v"], row["sa"]) <= 1e-5 and row["change"] is not None and row["change"] <= .001]
            summary["cycles"][name] = {"last_checkpoint": data[-1],
                "tail_max_step_change": max(row["step_change"] for row in tail),
                "first_transport_and_load_gate": transport[0] if transport else None,
                "all_gates_pass": bool(data[-1]["converged"])}
    summary["channel_checkpoint_comparison"] = {}
    for backend in ("cpu", "metal"):
        before = {row["it"]: row for row in steps(args.raw / f"channel-million-before-{backend}.log")}
        after = {row["it"]: row for row in steps(args.raw / f"channel-million-after-{backend}.log")}
        common = sorted(before.keys() & after.keys())
        differences = [{"iteration": iteration, "field": key, "before": before[iteration].get(key),
                        "after": after[iteration].get(key)}
            for iteration in common for key in before[iteration].keys() | after[iteration].keys()
            if key != "seconds" and before[iteration].get(key) != after[iteration].get(key)]
        summary["channel_checkpoint_comparison"][backend] = {"common_checkpoints": len(common),
                                                            "differences_excluding_seconds": differences}
    summary["spatial"] = {stem: spatial(args.raw, stem) for stem in ("plate0", "plate2-metal")}
    summary["late_plate_trace"] = {backend: trace_summary(args.raw / f"plate-trace-{backend}.csv")
                                   for backend in ("cpu", "metal")}
    summary["regression"] = compare(args.raw / "baseline-slow.xml", args.raw / "final-slow.xml", args.out)
    summary["slow_totals"] = {}
    for phase in ("baseline", "final"):
        tree = ET.fromstring(read(args.raw / (phase + "-slow.xml")))
        summary["slow_totals"][phase] = {name: {key: int(value) for key, value in tree.find(name).attrib.items()}
            for name in ("OverallResultsAsserts", "OverallResultsTestCases")}
    timing_rows = []
    for run in summary["runs"]:
        end = run.get("last_checkpoint", run)
        converged = bool(run.get("converged", end.get("converged")))
        timing_rows.append({"case": run["record"], "phase": run["phase"], "backend": run["backend"],
            "mode": run.get("mode", ""), "level": run.get("level", ""),
            "cells": run.get("cells", 983040 if run["record"] == "million_cell_channel" else ""),
            "iterations_observed": end.get("iterations", end.get("steps", end.get("it"))),
            "seconds_observed": end["seconds"], "converged": converged,
            "seconds_to_convergence": end["seconds"] if converged else "",
            "continuity_at_stop": end.get("continuity", "")})
    with (args.out / "timings.csv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=timing_rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(timing_rows)
    (args.out / "measurements.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    hashes = {str(path.relative_to(args.out)): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in sorted(args.out.rglob("*")) if path.is_file() and path.name != "sha256.json"
              and "__pycache__" not in path.parts and path.suffix != ".pyc"}
    (args.out / "sha256.json").write_text(json.dumps(hashes, indent=2) + "\n")


if __name__ == "__main__":
    main()
