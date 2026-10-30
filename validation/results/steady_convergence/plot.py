#!/usr/bin/env python3
"""Publication/export figures; logs and spatial captures are the data source."""
import argparse
import csv
import math
import statistics
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from analyze import read, steps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, default=Path(__file__).parent / "raw")
    parser.add_argument("--out", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.6), sharey=True, layout="constrained")
    for axis, backend in zip(axes, ("cpu", "metal")):
        for phase, label, color in (("before", "fb827e0", "#ba5728"), ("after", "Full local SA derivative", "#24618b")):
            data = [row for row in steps(args.raw / f"sa-cycle-{phase}-{backend}.log") if 980 <= row["it"] <= 1000]
            axis.plot([row["it"] for row in data], [1000 * row["cf0"] for row in data],
                      "o-", ms=3, lw=1.3, label=label, color=color)
        axis.set(title=backend.upper(), xlabel="Outer iteration", xticks=[980, 985, 990, 995, 1000])
        axis.grid(alpha=.15)
    axes[0].set_ylabel("Station Cf × 1000")
    axes[0].legend(loc="upper left", frameon=False, fontsize=8)
    fig.suptitle("Coarse SA plate: the two-step load cycle is removed; continuity still fails", fontsize=12)
    for suffix in ("png", "svg"):
        fig.savefig(args.out / ("sa-cycle." + suffix), dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 5.3), sharex=True, sharey="row", layout="constrained")
    for column, backend in enumerate(("cpu", "metal")):
        data = list(csv.DictReader(read(args.raw / f"plate-trace-{backend}.csv").splitlines()))
        iteration = [int(row["iteration"]) for row in data]
        pressure = [float(row["cell_p"]) for row in data]
        axes[0, column].plot(iteration, [1e6 * float(row["cell_v"]) for row in data], lw=1, color="#24618b")
        axes[1, column].plot(iteration, [1e6 * (p - statistics.mean(pressure)) for p in pressure], lw=1, color="#ba5728")
        axes[0, column].set_title(backend.upper())
        axes[1, column].set_xlabel("Outer iteration")
        for axis in axes[:, column]:
            axis.grid(alpha=.15)
    axes[0, 0].set_ylabel("Cell v × 10⁶")
    axes[1, 0].set_ylabel("Cell pressure deviation × 10⁶")
    fig.suptitle("Baseline SST plate: local fields still fluctuate at the global residual floor\nRCM cell 795, x=-0.266386, y=0.000123023", fontsize=12)
    for suffix in ("png", "svg"):
        fig.savefig(args.out / ("pressure-velocity-trace." + suffix), dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True, layout="constrained")
    for axis, stem, label in zip(axes, ("plate0", "plate2-metal"), ("816 cells, CPU", "13,056 cells, Metal")):
        faces = list(csv.DictReader(read(args.raw / (stem + "-faces.csv")).splitlines()))
        cells = list(csv.DictReader(read(args.raw / (stem + "-cells.csv")).splitlines()))
        residual = [0.] * len(cells)
        for face in faces:
            flux = float(face["flux"])
            residual[int(face["owner"])] += flux
            if int(face["neighbour"]) >= 0:
                residual[int(face["neighbour"])] -= flux
        image = axis.scatter([float(row["x"]) for row in cells], [float(row["y"]) for row in cells],
                             c=[max(abs(value), 1e-16) for value in residual], s=5,
                             norm=LogNorm(vmin=1e-15, vmax=1e-7), cmap="magma")
        axis.set(yscale="log", title=label, xlabel="x", ylim=(1e-6, 1.1))
        axis.axvline(0, color="grey", ls=":", lw=.7)
        axis.axhline(.01, color="grey", ls=":", lw=.7)
    axes[0].set_ylabel("y, logarithmic spacing")
    fig.colorbar(image, ax=axes, label="Absolute cell mass imbalance, kg/s")
    fig.suptitle("SST plate at 3,000 steps: divergence of stored provisional face flux", fontsize=12)
    for suffix in ("png", "svg"):
        fig.savefig(args.out / ("continuity-spatial." + suffix), dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
