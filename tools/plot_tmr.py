#!/usr/bin/env python3
"""Plot the recorded V5 plate trend (requires matplotlib; no solver runs)."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

root = Path(__file__).resolve().parents[1]
results = root / "validation/results/tmr_v5"
rows = list(csv.DictReader((results / "measurements.csv").open()))

def reference(name, zone=0):
    current = -1
    points = []
    path = root / "validation/data/nasa_tmr/2DZP_flatplate/sst" / name
    for line in path.read_text().splitlines():
        if line.lower().startswith("zone"):
            current += 1
        elif current == zone:
            try:
                values = [float(word) for word in line.split()]
            except ValueError:
                continue
            if len(values) == 4:
                points.append((values[0], values[3] * 1000))
    return sorted(points)

plt.rcParams.update({"font.size": 10, "svg.fonttype": "none", "svg.hashsalt": "kutta-tmr-v5"})
fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.2), layout="constrained")
for ax, field, title, data in zip(axes, ("cf0", "drag"), ("Skin friction at x=0.970084071", "Drag; reference area = 2 × span"), ("cf_convergence_sstv.dat", "drag_convergence_sstv.dat")):
    for backend, marker, color in (("cpu", "o", "#2563eb"), ("metal", "+", "#f97316")):
        subset = sorted((int(r["cells"]), float(r[field]) * 1000) for r in rows if r["case"] == "plate" and r["backend"] == backend)
        ax.plot(*zip(*subset), marker=marker, markersize=7, linewidth=1.5, color=color, label=f"Kutta {backend.upper()} (unconverged)")
    for zone, label, style in ((0, "TMR CFL3D M=0.2", "--"), (1, "TMR FUN3D M=0.2", ":")):
        ax.plot(*zip(*reference(data, zone)), style, color="#64748b", label=label)
    if field == "cf0":
        ax.plot(*zip(*reference("cf_incomp_results_sstv.dat")), "-.", color="#15803d", label="TMR FUN3D incompressible")
    ax.set_xscale("log", base=4)
    ax.set_xticks([816, 3264, 13056, 52224, 208896])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.tick_params(axis="x", labelsize=8)
    ax.set_xlabel("Cells (finer →)")
    ax.set_ylabel("Coefficient × 1000")
    ax.set_title(title)
    ax.grid(alpha=.2)
axes[0].legend(fontsize=8, loc="lower right")
axes[1].legend(fontsize=8, loc="lower right")
fig.suptitle("V5 flat plate: force trend agrees; normalized continuity gate remains red", fontsize=12)
svg = results / "flatplate_convergence.svg"
fig.savefig(svg, metadata={"Date": None})
svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
fig.savefig(results / "flatplate_convergence.png", dpi=180)
