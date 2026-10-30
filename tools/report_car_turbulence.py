#!/usr/bin/env python3
"""Extract complete V7/V8 logs and plot accuracy diagnostics (requires matplotlib)."""
import argparse
import bisect
import csv
import json
import re
from pathlib import Path


def zones(path):
    result = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.lower().startswith("zone"):
            result.append([])
        elif line and line[0] not in "#" and not line.lower().startswith("variables"):
            result[-1].append([float(x) for x in line.split()])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--ke-log", type=Path, help="Use a separately pinned V8 rerun instead of V8 entries in log")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1] / "validation/data/nasa_tmr"
    sa, ke = [], []
    lines = args.log.read_text().splitlines()
    if args.ke_log:
        lines = [line for line in lines if not line.startswith("V8 KE channel ")]
        lines += [line for line in args.ke_log.read_text().splitlines() if line.startswith("V8 KE channel ")]
    for line in lines:
        if not line.startswith(("V7 SA ", "V8 KE channel ")):
            continue
        pairs = re.findall(r"(\w+)=([^ ]+)", line)
        row = {k: (v if k == "backend" else float(v)) for k, v in pairs if k not in {"Cf", "Cp", "failure"}}
        if line.startswith("V7 SA "):
            row["case"] = line.split()[2]
            row["cf"] = [float(v) for k, v in pairs if k == "Cf"]
            row["cp"] = [float(v) for k, v in pairs if k == "Cp"]
            folder = root / ("2DZP_flatplate/sa" if row["case"] == "plate" else "2DBC_bump/sa")
            cf = zones(folder / ("cf_convergence.dat" if row["case"] == "plate" else "cf_convergence_sa.dat"))
            cd = zones(folder / ("drag_convergence.dat" if row["case"] == "plate" else "force_convergence.dat"))
            indices = [0] if row["case"] == "plate" else [2, 0, 4]
            row["cf_reference"] = [cf[z][0][3] for z in indices]
            row["cd_reference"] = cd[0][0][3 if row["case"] == "plate" else 4]
            row["cf_error_percent"] = [100 * (row["cf"][j] / ref - 1) for j, ref in enumerate(row["cf_reference"])]
            row["cd_error_percent"] = 100 * (row["Cd"] / row["cd_reference"] - 1)
            row["cp_reference"], row["cp_error"] = [], []
            if row["case"] == "bump":
                distribution = sorted(zones(folder / "cp_bump.dat")[0])
                xvalues = [point[0] for point in distribution]
                for j, station in enumerate([.6321975, .75, .8678025]):
                    high = bisect.bisect_right(xvalues, station)
                    lo, hi = distribution[high - 1], distribution[high]
                    ref = lo[1] + (hi[1] - lo[1]) * (station - lo[0]) / (hi[0] - lo[0])
                    row["cp_reference"].append(ref)
                    row["cp_error"].append(abs(row["cp"][j] - ref))
            sa.append(row)
        else:
            row["u_error_percent"] = 100 * row["max_U_error"]
            ke.append(row)
    # Missing backend/grid results are a failed study, not a shortened report.
    assert len(sa) == 12 and {(r["case"], r["backend"], int(r["level"])) for r in sa} == {
        (case, backend, level) for case in ("plate", "bump") for backend in ("cpu", "metal") for level in range(3)
    }, "incomplete V7 study"
    assert len(ke) == 6 and {(r["backend"], int(r["cells"])) for r in ke} == {
        (backend, cells) for backend in ("cpu", "metal") for cells in (20, 40, 80)
    }, "incomplete V8 study"
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "measurements.json").write_text(json.dumps({"sa": sa, "ke": ke}, indent=2, allow_nan=False) + "\n")
    with (args.out / "measurements.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["study", "backend", "cells", "steps", "seconds", "converged", "continuity_or_defect", "load_or_profile_change", "cf_or_u_error_percent", "cd_error_percent"])
        for r in sa:
            writer.writerow(["SA " + r["case"], r["backend"], r["cells"], r["steps"], r["seconds"], r["converged"], r["continuity"], r["load_change"], max(map(abs, r["cf_error_percent"])), r["cd_error_percent"]])
        for r in ke:
            writer.writerow(["KE channel", r["backend"], r["cells"], r["steps"], r["seconds"], r["converged"], r["defect"], r["change"], r["u_error_percent"], ""])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullLocator
    figure, axes = plt.subplots(1, 3, figsize=(11, 3.7), layout="constrained")
    for panel, (title, rows, limit) in enumerate([
        ("SA plate Cf", [r for r in sa if r["case"] == "plate"], 5),
        ("SA bump max station Cf", [r for r in sa if r["case"] == "bump"], 10),
        ("Realizable k-epsilon channel U", ke, 10),
    ]):
        axis = axes[panel]
        for backend, color, style in [("cpu", "#2d5f91", "-"), ("metal", "#a35520", "--")]:
            values = sorted([r for r in rows if r["backend"] == backend], key=lambda r: r["cells"])
            x = [r["cells"] for r in values]
            y = [r["u_error_percent"] if panel == 2 else max(map(abs, r["cf_error_percent"])) for r in values]
            axis.plot(x, y, style, color=color, label=backend.upper())
            for xi, yi, r in zip(x, y, values):
                axis.plot(xi, yi, "o", color=color, markerfacecolor=color if r["converged"] else "white", markersize=6)
        axis.axhline(limit, color="#777777", linestyle=":", label="Accuracy limit")
        axis.set(xscale="log", xlabel="Cells", ylabel="Absolute relative error (%)", title=title)
        cells = sorted({int(r["cells"]) for r in rows})
        axis.set_xticks(cells, [f"{n:,}" for n in cells], fontsize=9)
        axis.xaxis.set_minor_locator(NullLocator())
        axis.set_ylim(bottom=0)
        axis.grid(alpha=.18)
    axes[0].legend(frameon=False, fontsize=8)
    figure.suptitle("V7/V8 accuracy diagnostics — hollow markers fail iterative convergence", fontsize=11)
    figure.savefig(args.out / "accuracy.png", dpi=180)
    figure.savefig(args.out / "accuracy.svg")
    print(f"Extracted {len(sa)} SA and {len(ke)} channel measurements")


if __name__ == "__main__":
    main()
