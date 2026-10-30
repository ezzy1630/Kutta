#!/usr/bin/env python3
"""Plot the independent original/reassembled pressure checks (matplotlib)."""
from pathlib import Path
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[1] / "validation/results/pressure_float"
logs = [(root / name).read_text() for name in ("original_analysis.log", "reassembly_analysis.log")]
negative = [int(re.search(r"negative row sums=(\d+)", log)[1]) for log in logs]
minimum = [float(re.search(r"smallest eigenvalues=\[([^\s\]]+)", log)[1]) for log in logs]
names = ["Original", "Reassembled"]
colors = ["#b34440", "#24775a"]
fig, axes = plt.subplots(1, 2, figsize=(8, 2.8), constrained_layout=True)
for ax in axes:
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xticks(range(2), names)
axes[0].bar(range(2), negative, color=colors, width=.55)
axes[0].set_title("Negative stored row sums")
axes[0].set_ylim(0, 2100)
for i, value in enumerate(negative):
    axes[0].text(i, value + 65, f"{value:,}", ha="center")
axes[1].set_yscale("symlog", linthresh=1e-18)
axes[1].axhline(0, color=".65", linewidth=.7)
axes[1].scatter(range(2), minimum, c=colors, s=55, zorder=3)
axes[1].set_title("Smallest eigenvalue (symlog)")
for i, value in enumerate(minimum):
    axes[1].annotate(f"{value:+.3e}", (i, value), xytext=(0, 12),
                     textcoords="offset points", ha="center", fontsize=9)
axes[1].set_xlim(-.5, 1.5)
axes[1].margins(y=.35)
fig.suptitle("Identical TMR pressure assembly inputs: 3,520 cells", fontsize=12)
fig.savefig(root / "pressure_spd_comparison.png", dpi=180)
fig.savefig(root / "pressure_spd_comparison.svg")
