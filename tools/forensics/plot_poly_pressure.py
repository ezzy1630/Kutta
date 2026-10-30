#!/usr/bin/env python3
"""Plot the recorded polyhedral pressure diagnostics (requires matplotlib)."""
import argparse
import math
import pathlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def rows(path, prefix):
    return [dict(item.split("=", 1) for item in line.split()[1:] if "=" in item)
            for line in path.read_text().splitlines() if line.startswith(prefix)]


def cpu(path, prefix):
    return [row for row in rows(path, prefix) if row.get("backend") == "cpu"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    args = parser.parse_args()
    root = args.directory
    coarse = float(cpu(root / "refine_baseline8.log", "POLY_SOLVE")[0]["pressure_L2"])
    middle = float(cpu(root / "refine_baseline16.log", "POLY_SOLVE")[0]["pressure_L2"])
    fine = float([r for r in cpu(root / "baseline_v72.log", "FSAE V7.2")
                  if r["poly"] == "1"][0]["pressure_relative_L2"])
    # The two rejected trial meshes were run independently; no middle level
    # was measured for this trial, so only its two observed points are shown.
    side = [float(cpu(root / name, "POLY_SOLVE")[0]["pressure_L2"])
            for name in ("side8.log", "side32.log")]
    mu_errors = [float(cpu(root / name, "POLY_SOLVE")[0]["pressure_L2"])
                 for name in ("mu01_8.log", "refine_baseline8.log", "mu2_8.log")]
    viscous = [r for r in cpu(root / "quadratic_viscous.log", "POLY_VISCOUS")
               if r["component"] == "0" and r["exact_grad"] == "0"]
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    fig, ax = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    ax[0, 0].loglog([8, 16, 32], [100*x for x in (coarse, middle, fine)], "o-",
                    color="#b83a38", label="Retained solver")
    ax[0, 0].scatter([8, 32], [100*x for x in side], color="#285d91", marker="s",
                     label="Rejected side-velocity trial")
    ax[0, 0].axhline(1, color="#555555", linestyle="--", label="Unchanged 1% gate")
    ax[0, 0].set(title="Taylor–Couette pressure", xlabel="Radial rings (8 sectors/ring)",
                 ylabel="Pressure relative L2 (%)", xticks=[8, 16, 32], xticklabels=[8, 16, 32])
    ax[0, 0].legend(fontsize=8)
    ax[0, 1].plot([.1, 1, 2], [100*x for x in mu_errors], "o-", color="#b83a38")
    ax[0, 1].set(title="Coarse pressure error vs viscosity", xlabel="Dynamic viscosity (ρ = 1)",
                 ylabel="Pressure relative L2 (%)", ylim=(0, 20))
    for key, label in (("normal_derivative0", "Azimuthal faces"),
                       ("normal_derivative1", "Fine/coarse radial faces"),
                       ("normal_derivative3", "Outer wall")):
        ax[1, 0].loglog([int(r["nr"]) for r in viscous], [float(r[key]) for r in viscous],
                        "o-", label=label)
    ax[1, 0].set(title="Quadratic viscous flux, u = x²",
                 xlabel="Radial rings", ylabel="Area-weighted RMS", xticks=[8, 16, 32], xticklabels=[8, 16, 32])
    ax[1, 0].legend(fontsize=8)
    momentum = [r for r in cpu(root / "final_momentum_decomposition.log", "POLY_MOMENTUM") if r["nr"] == "32"]
    force = [r for r in cpu(root / "final_momentum_decomposition.log", "POLY_PRESSURE_FORCE")
             if r["nr"] == "32" and r["mode"] == "0"][0]
    def tangential(mode, convection):
        return float(next(r for r in momentum if r["gradient_mode"] == mode and r["convection"] == convection)["tangential_defect"])
    ax[1, 1].barh(["Viscous, LSQ", "Viscous, exact ∇u", "Pressure force", "Convection"],
                   [tangential("0", "0"), tangential("1", "0"), float(force["tangential_defect"]), tangential("0", "1")],
                   color=["#b83a38", "#b83a38", "#777777", "#777777"])
    ax[1, 1].invert_yaxis()
    ax[1, 1].set(title="Tangential defect, exact Taylor-Couette", xlabel="Volume-weighted RMS (μ = 1)")
    for axis in (ax[0, 0], ax[1, 0]):
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.set_xlim(7, 36)
    fig.savefig(root / "diagnostics.png", dpi=160)
    fig.savefig(root / "diagnostics.svg")
    print("Pressure refinement orders:", math.log2(coarse/middle), math.log2(middle/fine))
    for key in ("normal_derivative0", "normal_derivative1", "normal_derivative2", "normal_derivative3"):
        print(key, "finest-pair order", math.log2(float(viscous[-2][key])/float(viscous[-1][key])))


if __name__ == "__main__":
    main()
