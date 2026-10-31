"""Figures and tables for validation/results/transition_verification/REPORT.md.

Inputs: Kutta runs under runs/verify (tools/transition/flatplate.py run <case>
--out runs/verify/flatplate; bench_tmr ... plate metal 2 30000 runs/verify/tmr/<mode> <mode>)
and the tracked independent_<case>.csv from boundary_layer.py.

  build-py/bin/python tools/transition/verify/figures.py
"""
import importlib.util
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tools/transition"))
import boundary_layer as bl  # noqa: E402
import channel  # noqa: E402
import flatplate  # noqa: E402

# tools/fsae_aero/report.py (palette and rcParams), loaded by path as tools/transition/report.py does.
_spec = importlib.util.spec_from_file_location("fsae_aero_report", ROOT / "tools/fsae_aero/report.py")
aero = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(aero)

RUNS = ROOT / "runs/verify"
OUT = ROOT / "validation/results/transition_verification"
BLUE, ORANGE, AQUA, VIOLET, INK, MUTED, GRID = aero.BLUE, aero.ORANGE, aero.AQUA, aero.VIOLET, aero.INK, aero.MUTED, aero.GRID
CASES = {"t3a": "T3A", "t3b": "T3B", "t3am": "T3A-"}


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def independent(case):
    return np.genfromtxt(OUT / f"independent_{case}.csv", delimiter=",", names=True)


def t3(lines):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    lines.append("## T3 flat plates: Kutta vs the independent solution of the documented equations\n")
    for ax, (case, title) in zip(axes, CASES.items()):
        run = RUNS / "flatplate" / f"{case}_L1_gamma"
        comp, exp = flatplate.cf_table(case, run)
        ind = independent(case)
        rex = np.geomspace(3e3, ind["re_x"][-1], 300)
        ax.plot(rex, 0.664 / np.sqrt(rex), color=GRID, linewidth=1.5, label="Blasius")
        ax.plot(ind["re_x"], ind["cf_sst"], color=MUTED, linewidth=1, linestyle=":", label="independent, SST")
        ax.plot(comp[:, 0], comp[:, 1], color=BLUE, linewidth=2, label="Kutta, SST + gamma (L1)")
        ax.plot(ind["re_x"], ind["cf_gamma"], color=ORANGE, linewidth=1.5, linestyle="--", label="independent, SST + gamma")
        ax.plot(exp[:, 1], exp[:, 3], "o", color=INK, markersize=4, label="experiment")
        ax.set_xscale("log")
        ax.set_xlim(1e4 if case != "t3am" else 1e5, ind["re_x"][-1])
        ax.set_ylim(0, 0.009 if case != "t3am" else 0.006)
        ax.set_title(title)
        ax.set_xlabel("Re_x")
        if case == "t3a":
            ax.set_ylabel("skin friction Cf")
            ax.legend(fontsize=7.5)
        on_k, end_k = bl.transition(comp[:, 0], comp[:, 1])
        on_i, end_i = bl.transition(ind["re_x"], ind["cf_gamma"])
        on_e, end_e = bl.transition(exp[:, 1], exp[:, 3])
        lines.append(f"### {title}\n")
        lines.append(f"Onset (Cf minimum) Kutta {on_k:.3e}, independent {on_i:.3e} ({100*(on_k/on_i-1):+.1f}%), "
                     f"experiment {on_e:.3e}. End (first Cf maximum) Kutta {end_k:.3e}, independent {end_i:.3e}"
                     f" ({100*(end_k/end_i-1):+.1f}%), experiment {end_e:.3e}"
                     f"{' (last station; still rising)' if end_e >= exp[-1, 1] else ''}.\n")
        # Measured Cf values are restricted (validation/data/ercoftac_t3_flat_plate/SOURCE.md): only
        # the comparison's summary is quoted; per-station rows compare Kutta with the independent solution.
        pre = [float(np.interp(r[1], ind["re_x"], ind["cf_gamma"]) / r[3] - 1) for r in exp if r[1] < on_e]
        lines.append(f"Largest pre-transition Cf difference, independent solution vs experiment: {100*max(pre, key=abs):+.1f}%.\n")
        lines.append("| Re_x (measured stations) | Cf Kutta | Cf independent | Kutta vs independent | Blasius |")
        lines.append("|---:|---:|---:|---:|---:|")
        for row in exp:
            k = float(np.interp(row[1], comp[:, 0], comp[:, 1]))
            i = float(np.interp(row[1], ind["re_x"], ind["cf_gamma"]))
            lines.append(f"| {row[1]:.3e} | {k:.5f} | {i:.5f} | {100*(k/i-1):+.1f}% | {0.664/np.sqrt(row[1]):.5f} |")
        lines.append("")
    fig.tight_layout()
    save(fig, "t3_independent")


def tmr_cf(mode):
    d = np.genfromtxt(RUNS / "tmr" / mode / "wall_cf.csv", delimiter=",", names=True, dtype=None, encoding="ascii")
    return np.sort(d[d["x"] > 0], order="x")


def tmr(lines):
    s, g = tmr_cf("high-tu"), tmr_cf("high-tu-transition")
    ind = independent("tmr")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.8))
    ax1.plot(s["x"], s["cf"], color=MUTED, linewidth=1.5, label="Kutta, SST")
    ax1.plot(g["x"], g["cf"], color=BLUE, linewidth=2, label="Kutta, SST + gamma")
    ax1.plot(ind["x"], ind["cf_sst"], color=MUTED, linewidth=1, linestyle=":", label="independent, SST")
    ax1.plot(ind["x"], ind["cf_gamma"], color=ORANGE, linewidth=1.5, linestyle="--", label="independent, SST + gamma")
    ax1.set_xscale("log")
    ax1.set_xlim(1e-3, 2)
    ax1.set_ylim(0, 0.008)
    ax1.set_xlabel("x (Re 5e6 per unit length)")
    ax1.set_ylabel("skin friction Cf")
    ax1.set_title("TMR plate, inlet Tu 5%, mu_t/mu 100")
    ax1.legend(fontsize=7.5)
    xk = s["x"]
    ax2.plot(xk, 100 * (np.interp(xk, g["x"], g["cf"]) / s["cf"] - 1), color=BLUE, linewidth=2, label="Kutta")
    ax2.plot(ind["x"], 100 * (ind["cf_gamma"] / ind["cf_sst"] - 1), color=ORANGE, linewidth=1.5, linestyle="--",
             label="independent")
    ax2.axhline(0, color=MUTED, linewidth=1)
    ax2.set_xscale("log")
    ax2.set_xlim(0.1, 2)
    ax2.set_xticks([0.1, 0.2, 0.5, 1, 2], ["0.1", "0.2", "0.5", "1", "2"])
    ax2.minorticks_off()
    ax2.set_ylim(-2, 12)
    ax2.set_xlabel("x")
    ax2.set_ylabel("Cf, SST + gamma vs SST (%)")
    ax2.set_title("Turbulent offset")
    ax2.legend(fontsize=7.5)
    fig.tight_layout()
    save(fig, "tmr_independent")
    lines.append("## TMR plate: turbulent offset\n")
    lines.append("| x | Cf SST Kutta | Cf gamma Kutta | offset Kutta | Cf SST independent | Cf gamma independent | offset independent |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|")
    for x in (0.05, 0.1, 0.2, 0.5, 0.97, 1.5, 1.9):
        ks, kg = np.interp(x, s["x"], s["cf"]), np.interp(x, g["x"], g["cf"])
        i_s, i_g = np.interp(x, ind["x"], ind["cf_sst"]), np.interp(x, ind["x"], ind["cf_gamma"])
        lines.append(f"| {x} | {ks:.6f} | {kg:.6f} | {100*(kg/ks-1):+.2f}% | {i_s:.6f} | {i_g:.6f} | {100*(i_g/i_s-1):+.2f}% |")
    lines.append("")


def column(run, x0):
    d = np.genfromtxt(run / "cells.csv", delimiter=",", names=True)
    xs = np.unique(np.round(d["x"], 9))
    c = d[abs(d["x"] - xs[np.argmin(abs(xs - x0))]) < 1e-8]
    return c[np.argsort(c["y"])]


def profiles():
    nu = 1 / 5e6
    cg, cs = column(RUNS / "tmr/high-tu-transition", 0.95), column(RUNS / "tmr/high-tu", 0.95)
    utau = np.sqrt(nu * cg["u"][0] / cg["y"][0])
    ypk = cg["y"] * utau / nu
    c = bl.CASES["tmr"]
    _, _, pg = bl.run(*c, True, profile_x=(0.95,))
    _, _, ps = bl.run(*c, False, profile_x=(0.95,))
    pg, ps = pg[0.95], ps[0.95]
    ypi = pg["y"] * np.sqrt(pg["cf"] / 2) / nu
    ch_s, ch_g = channel.solve(2000), channel.solve(2000, gamma_model=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.8))
    ax1.plot(ypk, cg["intermittency"], color=BLUE, linewidth=2, label="Kutta, TMR plate x = 0.95")
    ax1.plot(ypi[1:], pg["g"][1:], color=ORANGE, linewidth=1.5, linestyle="--", label="independent BL, x = 0.95")
    ax1.plot(ch_g["y"], ch_g["g"], color=VIOLET, linewidth=1.5, linestyle=":", label="independent channel, Re_tau 2000")
    ax1.axhline(1 / 50, color=MUTED, linewidth=1)
    ax1.text(200, 0.045, "1/c_e2", color=MUTED, fontsize=8)
    ax1.set_xscale("log")
    ax1.set_xlim(0.1, 1000)
    ax1.set_ylim(0, 1.05)
    ax1.set_xlabel("y+")
    ax1.set_ylabel("intermittency gamma")
    ax1.set_title("Gamma across the turbulent boundary layer")
    ax1.legend(fontsize=7.5, loc="upper left")
    ax2.plot(ypk, cg["mut_over_mu"] / cs["mut_over_mu"], color=BLUE, linewidth=2, label="Kutta")
    ax2.plot(ypi[1:], pg["nut"][1:] / np.maximum(ps["nut"][1:], 1e-300), color=ORANGE, linewidth=1.5, linestyle="--",
             label="independent BL")
    ax2.plot(ch_g["y"], ch_g["nut"] / ch_s["nut"], color=VIOLET, linewidth=1.5, linestyle=":", label="independent channel")
    ax2.axhline(1, color=MUTED, linewidth=1)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(0.1, 1000)
    ax2.set_ylim(0.8, 1e3)
    ax2.set_xlabel("y+")
    ax2.set_ylabel("nu_t, SST + gamma / SST")
    ax2.set_title("Eddy viscosity ratio to plain SST")
    ax2.legend(fontsize=7.5)
    fig.tight_layout()
    save(fig, "gamma_profile")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["# Tables (generated by tools/transition/verify/figures.py)\n"]
    t3(lines)
    tmr(lines)
    profiles()
    (OUT / "tables.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
