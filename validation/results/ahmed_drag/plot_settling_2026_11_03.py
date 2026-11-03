#!/usr/bin/env python3
"""Plot the archived solver monitors; no restricted reference data is needed."""
import gzip
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(evidence, output):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True,
                             layout="constrained")
    controls = [
        ("85c489b_continuation", "Before shielding, pressure relaxation 1", "#427699"),
        ("df59f85", "First bad commit, pressure relaxation 1", "#b86d45"),
        ("e014e36_pressure_urf05_continuation", "Shielding retained, pressure relaxation 0.5", "#388271"),
    ]
    for name, label, color in controls:
        with gzip.open(evidence / (name + ".history.csv.gz"), "rt") as f:
            history = np.genfromtxt(f, delimiter=",", names=True)
        monitors = json.loads((evidence / (name + ".monitors.json")).read_text())
        iteration = np.array([r[0] for r in monitors])
        change = np.array([r[4] for r in monitors], dtype=float)
        residual = np.array([r[3] for r in monitors])
        valid = np.isfinite(change)
        axes[0].plot(history["iteration"], history["lift"], color=color,
                     lw=1.1, label=label)
        axes[1].plot(iteration[valid], change[valid], color=color, lw=1.1)
        passed = valid & (change <= 5e-4) & (residual <= 1e-3)
        axes[1].scatter(iteration[passed], change[passed], s=13, color=color,
                        zorder=3)
    axes[0].set_ylim(.2088, .20985)
    axes[0].set_ylabel("Lift coefficient, Cl")
    axes[0].legend(frameon=False, loc="upper right", fontsize=9)
    axes[0].set_title("Old solver's 1,100-iteration pass was a quiet window")
    axes[1].axhline(5e-4, color="#444444", ls="--", lw=1,
                    label="Unchanged both-force gate")
    axes[1].set_yscale("log")
    axes[1].set_ylim(5e-5, .01)
    axes[1].set_ylabel("Maximum relative force change\nover 500 iterations")
    axes[1].set_xlabel("Steady solver iteration (50-iteration monitor interval)")
    axes[1].set_xlim(800, 6000)
    axes[1].legend(frameon=False, loc="upper right")
    for ax in axes:
        ax.axvline(1100, color="#777777", ls=":", lw=1)
        ax.grid(alpha=.18)
    fig.suptitle("Ahmed 35°: persistent lift variation predates df59f85", fontsize=14)
    fig.supxlabel("403,560 cells; same mesh, settings and 6,000-iteration cap. "
                  "Dots meet residual and both-force gates.\n"
                  "Continuations disable early stopping only; they do not loosen either gate.",
                  fontsize=9)
    fig.savefig(output, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    here = pathlib.Path(__file__).resolve().parent
    plot(here / "evidence/2026-11-03/settling", here / "settling_2026-11-03.png")
