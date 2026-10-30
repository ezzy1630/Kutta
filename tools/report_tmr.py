#!/usr/bin/env python3
"""Summarize bench_tmr logs without converting failed solves into valid data."""
import argparse
import csv
import re
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--runs", type=Path, default=Path("build/v5"))
parser.add_argument("--test-log", type=Path, help="Use final slow-suite measurement lines instead of benchmark logs")
parser.add_argument("--output", type=Path, default=Path("validation/results/tmr_v5"))
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
fields = ["case", "backend", "level", "cells", "iterations", "setup_s", "solve_s", "converged", "budget", "defect", "continuity_scaled", "change", "drag", "cf0", "cf1", "cf2", "cp0", "cp1", "cp2", "failure"]
def measurements(line):
    row = {}
    counts = {"cf": 0, "cp": 0}
    for key, value in re.findall(r"(\w+)=([^\s]+)", line):
        if key in counts:
            index = counts[key]
            counts[key] += 1
            key = f"{key}{index}"
        key = {"seconds": "solve_s", "continuity": "continuity_scaled"}.get(key, key)
        if key in fields:
            row[key] = value
    return row

rows = []
if args.test_log:
    for line in args.test_log.read_text().splitlines():
        case = "plate" if line.startswith("V5 plate backend=") else ("bump" if line.startswith("V5 bump measurement backend=") else None)
        if case is None:
            continue
        row = measurements(line)
        row["case"] = case
        row["level"] = str(({816: 0, 3264: 1, 13056: 2} if case == "plate" else {3520: 0, 14080: 1, 56320: 2})[int(row["cells"])])
        row["failure"] = line.partition(" failure=")[2]
        rows.append(row)
    if len(rows) != 12:
        raise SystemExit(f"Expected twelve measured grids; found {len(rows)}")
else:
    for case in ("plate", "bump"):
        for backend in ("cpu", "metal"):
            for level in range(3):
                path = args.runs / f"{case}-{backend}-{level}.log"
                if not path.exists():
                    raise SystemExit(f"Missing run: {path}")
                lines = path.read_text().splitlines()
                results = [line for line in lines if line.startswith("RESULT ")]
                if len(results) != 1:
                    raise SystemExit(f"Expected one completed result in {path}")
                row = measurements(results[0])
                row.update(case=case, backend=backend, level=level)
                row["failure"] = next((line for line in lines if line.startswith("FAILURE ")), "")
                rows.append(row)
expected = {(case, backend, str(level)) for case in ("plate", "bump") for backend in ("cpu", "metal") for level in range(3)}
observed = {(row["case"], row["backend"], str(row["level"])) for row in rows}
if observed != expected or len(rows) != len(expected):
    raise SystemExit("Measurements must cover each case/backend/level exactly once")
with (args.output / "measurements.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print("| Case | Backend | Cells | Iterations | Solve s | Converged | Cf (station order) | CD |")
print("| --- | --- | ---: | ---: | ---: | --- | --- | ---: |")
for row in rows:
    values = [row.get(f"cf{i}", "nan") for i in range(1 if row["case"] == "plate" else 3)]
    print(f'| {row["case"]} | {row["backend"]} | {row["cells"]} | {row["iterations"]} | {float(row["solve_s"]):.2f} | {row["converged"]} | {", ".join(values)} | {row["drag"]} |')
# Coefficients in the fixtures are already normalized by freestream quantities.
root = Path(__file__).resolve().parents[1] / "validation/data/nasa_tmr"
def zones(path):
    result = []
    for line in path.read_text().splitlines():
        if line.lower().startswith("zone"):
            result.append([])
        else:
            try:
                values = [float(word) for word in line.split()]
            except ValueError:
                continue
            if values and result:
                result[-1].append(values)
    return result

def interp(points, x):
    points = sorted(points)
    for a, b in zip(points, points[1:]):
        if a[0] <= x <= b[0]:
            return a[1] + (b[1] - a[1]) * (x - a[0]) / (b[0] - a[0])
    raise ValueError("station outside profile")

print("\nFinest CFL3D references and relative model spreads:")
for variant in ("sst2003", "sst", "sstv"):
    cf = zones(root / f"2DBC_bump/{variant}/cf_convergence_{variant}.dat")
    cp = zones(root / f"2DBC_bump/{variant}/cp_bump_{variant}.dat")
    print(variant, "cf", [cf[i][0][-1] for i in (2, 0, 4)], "cp", [interp(cp[0], x) for x in (.6321975, .75, .8678025)])
# Interpolate in x within each original structured cell row, then in y.
profiles = []
velocity = zones(root / "2DZP_flatplate/sst/flatplate_u_sstv.dat")[0]
viscosity = zones(root / "2DZP_flatplate/sst/mut_0.97_sstv.dat")[0]
for backend in ("cpu", "metal"):
    path = args.runs / f"plate-{backend}-2/cells.csv"
    if not path.exists():
        continue
    by_row = {}
    for cell in csv.DictReader(path.open()):
        row = int(cell["original_cell"]) // 136
        by_row.setdefault(row, []).append({key: float(value) for key, value in cell.items()})
    u, mut = [], []
    for cells in by_row.values():
        y = cells[0]["y"]
        u.append([y, interp([[c["x"], c["u"]] for c in cells], .970084071)])
        mut.append([y, interp([[c["x"], c["mut_over_mu"]] for c in cells], .970084071)])
    for y in (.0001, .001, .003, .01, .02):
        profiles.append({"backend": backend, "y": y, "u": interp(u, y), "tmr_u": interp([[r[1], r[0]] for r in velocity], y), "mut_over_mu": interp(mut, y), "tmr_mut_over_mu": interp([[r[1], r[2]] for r in viscosity], y)})
if profiles:
    with (args.output / "flatplate_profiles.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(profiles[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(profiles)
    print("\nFinest plate profiles:")
    for profile in profiles:
        print(profile)
