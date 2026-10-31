#!/usr/bin/env python3
"""Plots, tables and gzipped evidence for validation/results/fsae_aero.

Reads the run directories written by run.py (runs/fsae_aero/results*/<tag>).
Reference values are read at run time from the restricted files under
validation/data. Plots and tables carry only the handful of reference values
quoted in REPORT.md (Table 15 points, UIUC values at the run angles, Ahmed Cd
and region-mean Cp), never the reference curves themselves.

  build-py/bin/python tools/fsae_aero/report.py
"""
import csv
import gzip
import json
import math
import pathlib
import shutil

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs/fsae_aero"
OUT = ROOT / "validation/results/fsae_aero"
DATA = ROOT / "validation/data"
# Categorical slots in fixed order (validated: CVD and normal-vision separation
# pass; the aqua slot is low-contrast, so it always has its own marker shape,
# a legend entry and a table row).
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                     "lines.linewidth": 2, "lines.markersize": 7, "svg.hashsalt": "kutta-fsae-aero",
                     "axes.titlesize": 9.5, "axes.titlelocation": "left"})

AHMED_CD = 0.26                     # Ahmed et al. 1984, 35 deg (as cited by Liu & Moser)
AHMED_COMPONENTS = {"front": 0.020, "base": 0.095, "slant": 0.090}
S1223_LINEAR = 10.0                 # deg; CL gated up to here, trend only above
L_AHMED, SLANT_LEN, ROOF = 1.044, 0.222, 0.338


def results():
    rows = []
    for path in sorted(RUNS.glob("results*/*/result.json")):
        r = json.loads(path.read_text())
        r["dir"] = path.parent
        if "result" in r:
            rows.append(r)
    return rows


def pick(rows, tag="", **want):
    """Runs matching `want`; tagged runs (sensitivity checks) only on request."""
    return [r for r in rows if (r.get("tag") or "") == tag and all(r.get(k) == v for k, v in want.items())]


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def settle(res):
    vals = [res.get("settle_lift"), res.get("settle_drag")]
    return max(np.inf if v is None else v for v in vals)


def run_cols(r):
    res = r["result"]
    status = "converged" if res["residuals_met"] and res["settled"] else \
        "settled forces" if res["settled"] else "budget" if res["budget"] else "FAILED" if res["failure"] else "not settled"
    return (f"{res['cells']:,} | {res['iterations']} | {res['max_residual']:.1e} | {settle(res):.1e} | {status} | "
            f"{res['solve_s']:.0f} | {r.get('peak_footprint_bytes', 0) / 1e9:.2f}")


RUN_HEAD = "cells | iterations | max scaled residual | force change over window | stop | solve s | peak GB"


def table(lines, title, head, rows):
    lines += ["", f"### {title}", "", "| " + head + " |", "|" + "---|" * (head.count("|") + 1)]
    lines += ["| " + row + " |" for row in rows]


# ----------------------------------------------------------------------------- Zerihan
def zerihan_refs():
    rows = csv.DictReader(open(DATA / "zerihan_tyrrell026_ground_effect/table15_2d_sectional_Cl_vs_hc.csv"))
    return {r["h_over_c"]: (float(r["Cl_CFL3D_kw"]), float(r["Cl_CFL3D_SA"]), float(r["Cl_expt_centre_pressure_integrated"]))
            for r in rows}


def hc(point):
    return math.inf if point == "free" else float(point)


def zref(refs, point):
    if point == "free":
        return refs["freestream"]
    return next(v for k, v in refs.items() if k != "freestream" and abs(float(k) - float(point)) < 1e-6)


def zerihan(rows, lines, verdicts):
    refs = zerihan_refs()
    l2 = sorted(pick(rows, case="zerihan", level=2, backend="metal"), key=lambda r: -hc(r["point"]))
    if not l2:
        return
    body = []
    xs, ys, kw, sa, ex = [], [], [], [], []
    for r in l2:
        res, ref = r["result"], zref(refs, r["point"])
        d = res["lift"] / ref[0] - 1
        ok = abs(d) <= 0.05 and not res["failure"]
        verdicts.append(("Zerihan 2D", f"h/c {r['point']}", f"Cl {res['lift']:.3f} vs CFL3D k-omega {ref[0]:.3f} ({100 * d:+.1f}%)",
                         "PASS" if ok else "FAIL"))
        body.append(f"{r['point']} | {res['lift']:.4f} | {ref[0]:.3f} | {100 * d:+.1f}% | {'PASS' if ok else 'FAIL'} | "
                    f"{ref[1]:.3f} | {ref[2]:.3f} | {res['drag']:.4f} | {res['yplus_max']:.2f} | " + run_cols(r))
        x = 0.8 if r["point"] == "free" else float(r["point"])
        xs.append(x), ys.append(res["lift"]), kw.append(ref[0]), sa.append(ref[1]), ex.append(ref[2])
    table(lines, "Zerihan / Tyrrell 026, sectional Cl vs h/c (level 2, Metal)",
          "h/c | Kutta Cl | CFL3D k-omega | diff | 5% gate | CFL3D SA | expt (p-integrated) | Kutta Cd | wing y+ max | " + RUN_HEAD,
          body)

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    ax.plot(xs, kw, "s", color=ORANGE, label="Thesis 2D CFL3D k-omega (Table 15)", zorder=3)
    ax.plot(xs, sa, "D", color=VIOLET, markersize=5.5, label="Thesis 2D CFL3D SA (Table 15)", zorder=3)
    ax.plot(xs, ex, "^", color=AQUA, label="Experiment, pressure-integrated (Table 15)", zorder=3)
    # Freestream sits at a nominal x position, so it is not joined to the ground curve.
    g = [i for i, r in enumerate(l2) if r["point"] != "free"]
    ax.plot([xs[i] for i in g], [ys[i] for i in g], "o-", color=BLUE, label="Kutta SST, level 2", zorder=4)
    ax.plot([x for i, x in enumerate(xs) if i not in g], [y for i, y in enumerate(ys) if i not in g], "o", color=BLUE, zorder=4)
    for i in range(len(xs)):
        ax.plot([xs[i]] * 2, [kw[i] * 0.95, kw[i] * 1.05], color=ORANGE, alpha=0.35, linewidth=6, solid_capstyle="butt",
                label="k-omega +-5% (acceptance)" if i == 0 else None, zorder=2)
    ax.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    ax.set_xticklabels(["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "free"])
    ax.set_xlabel("ride height h/c")
    ax.set_ylabel("sectional downforce coefficient Cl")
    ax.set_title("Tyrrell 026 main element, alpha_ref 1 deg, Re 4.62e5, moving ground")
    ax.legend(loc="upper right", fontsize=8.5)
    save(fig, "zerihan_cl_vs_hc")

    # Mesh study and backend parity.
    body = []
    for point in ("free", "0.224", "0.090"):
        for level in (1, 2, 3):
            for r in pick(rows, case="zerihan", level=level, backend="metal", point=point):
                res, ref = r["result"], zref(refs, point)
                lift, drag, mark = value(r, "lift"), value(r, "drag"), "" if res["settled"] else "*"
                body.append(f"{point} | L{level} | {lift:.4f}{mark} | {100 * (lift / ref[0] - 1):+.1f}% | {drag:.4f}{mark} | "
                            f"{r['mesh']['oblock_layers']} | {r['mesh']['ds_wall']:.1e} | " + run_cols(r))
    if body:
        table(lines, "Zerihan mesh study (Metal); * = unsettled, mean of the last 2000 iterations",
              "h/c | level | Kutta Cl | vs k-omega | Kutta Cd | O-block layers | first cell / c | " + RUN_HEAD, body)
    cpu = []
    for c in pick(rows, case="zerihan", backend="cpu"):
        for g in pick(rows, case="zerihan", backend="metal", level=c["level"], point=c["point"]):
            cpu.append(f"{c['point']} | L{c['level']} | {c['result']['lift']:.5f} | {g['result']['lift']:.5f} | "
                       f"{c['result']['drag']:.5f} | {g['result']['drag']:.5f} | {c['result']['iterations']} / {g['result']['iterations']} | "
                       f"{c['result']['solve_s']:.0f} / {g['result']['solve_s']:.0f} | "
                       f"{c.get('peak_footprint_bytes', 0) / 1e9:.2f} / {g.get('peak_footprint_bytes', 0) / 1e9:.2f}")
    if cpu:
        table(lines, "CPU vs Metal (same mesh)",
              "h/c | level | Cl CPU | Cl Metal | Cd CPU | Cd Metal | iterations CPU / Metal | solve s CPU / Metal | peak GB CPU / Metal",
              cpu)
    far = []
    for f in pick(rows, tag="_far2", case="zerihan"):
        for g in pick(rows, case="zerihan", backend=f["backend"], level=f["level"], point=f["point"]):
            far.append(f"{f['point']} | L{f['level']} | {g['result']['lift']:.4f} | {f['result']['lift']:.4f} | "
                       f"{100 * (f['result']['lift'] / g['result']['lift'] - 1):+.2f}% | {g['result']['drag']:.5f} | {f['result']['drag']:.5f}")
    if far:
        table(lines, "Far-field distance (inlet 30c / outlet 40c / top 30c, then doubled)",
              "h/c | level | Cl | Cl, 2x domain | change | Cd | Cd, 2x domain", far)


# ----------------------------------------------------------------------------- S1223
def s1223_reference(re):
    rows = [r for r in csv.DictReader(open(DATA / "uiuc_lsat_s1223/s1223_vol1_drag_polars.csv")) if r["config"] == "clean"]
    groups = {2e5: ("200400", "198500"), 3e5: ("299100",)}[re]
    pts = sorted((float(r["alpha_deg"]), float(r["cl"]), float(r["cd"])) for r in rows if r["re_avg"] in groups)
    lift = [r for r in csv.DictReader(open(DATA / "uiuc_lsat_s1223/s1223_vol1_lift_curves.csv"))
            if r["config"] == "clean" and r["re_avg"] == {2e5: "198100", 3e5: "302200"}[re]]
    lift = [(float(r["alpha_deg"]), float(r["cl"])) for r in lift]
    a, cl, cd = (np.array(v) for v in zip(*pts))

    def at(alpha):
        if alpha > a.max():  # beyond the polar: CL from the lift-curve sweep, no CD
            up = sorted(lift[: len(lift) // 2])
            return float(np.interp(alpha, [p[0] for p in up], [p[1] for p in up])), None
        return float(np.interp(alpha, a, cl)), float(np.interp(alpha, a, cd))
    clmax = max(lift, key=lambda p: p[1])
    return at, clmax


def s1223(rows, lines, verdicts):
    runs = [r for r in pick(rows, case="s1223", level=2, backend="metal")]
    if not runs:
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0))
    for re, color, marker in ((2e5, BLUE, "o"), (3e5, VIOLET, "s")):
        at, clmax = s1223_reference(re)
        pts = sorted((r for r in runs if abs(r["reynolds"] - re) < 1), key=lambda r: float(r["point"]))
        if not pts:
            continue
        body, al, kl, kd, rl, rd = [], [], [], [], [], []
        for r in pts:
            res, alpha = r["result"], float(r["point"])
            lift, drag = value(r, "lift"), value(r, "drag")
            cl_ref, cd_ref = at(alpha)
            dl = lift - cl_ref
            dd = None if cd_ref is None else drag / cd_ref - 1
            linear = alpha <= S1223_LINEAR
            ok_l = abs(dl) <= 0.10 and not res["failure"]
            ok_d = dd is not None and abs(dd) <= 0.25 and not res["failure"]
            gate = ("PASS" if ok_l else "FAIL") if linear else "trend"
            dgate = "n/a" if dd is None else ("PASS" if ok_d else "FAIL")
            if linear:
                verdicts.append((f"S1223 Re {re:.0e}", f"alpha {alpha:g}", f"CL {lift:.3f} vs {cl_ref:.3f} ({dl:+.3f})", gate))
            if dd is not None:
                verdicts.append((f"S1223 Re {re:.0e}", f"alpha {alpha:g}", f"CD {drag:.4f} vs {cd_ref:.4f} ({100 * dd:+.0f}%)", dgate))
            mark = "" if res["settled"] else "*"
            body.append(f"{alpha:g} | {lift:.3f}{mark} | {cl_ref:.3f} | {dl:+.3f} | {gate} | {drag:.4f}{mark} | "
                        + ("n/a | n/a" if dd is None else f"{cd_ref:.4f} | {100 * dd:+.0f}%") + f" | {dgate} | " + run_cols(r))
            al.append(alpha), kl.append(lift), kd.append(drag), rl.append(cl_ref), rd.append(cd_ref)
        table(lines, f"Selig S1223, Re {re:.0e} (level 2, Metal); UIUC values interpolated at the run angle; "
                     f"UIUC CLmax {clmax[1]:.3f} at {clmax[0]:.1f} deg; * = unsettled, mean of the last 2000 iterations",
              "alpha | Kutta CL | UIUC CL | dCL | +-0.10 gate | Kutta CD | UIUC CD | dCD | +-25% gate | " + RUN_HEAD, body)
        label = f"Re {re / 1e5:.0f}e5"
        ax1.plot(al, rl, marker, color=AQUA if re == 2e5 else ORANGE, fillstyle="none", markeredgewidth=1.6,
                 label=f"UIUC experiment, {label} (at run angles)")
        ax1.plot(al, kl, marker + "-", color=color, label=f"Kutta SST, {label}")
        sel = [i for i, d in enumerate(rd) if d is not None]
        ax2.plot([rd[i] for i in sel], [rl[i] for i in sel], marker, color=AQUA if re == 2e5 else ORANGE,
                 fillstyle="none", markeredgewidth=1.6, label=f"UIUC experiment, {label}")
        ax2.plot(kd, kl, marker + "-", color=color, label=f"Kutta SST, {label}")
    ax1.axvspan(S1223_LINEAR, 17, color=GRID, alpha=0.5, linewidth=0)
    ax1.text(S1223_LINEAR + 0.3, 0.97, "trend only", color=MUTED, fontsize=8.5, va="top",
             transform=ax1.get_xaxis_transform())
    ax1.set_xlabel("angle of attack, deg")
    ax1.set_ylabel("lift coefficient CL")
    ax1.set_title("Lift curve (fully turbulent SST)")
    ax1.legend(fontsize=7.5, loc="lower right")
    ax2.set_xlabel("drag coefficient CD")
    ax2.set_ylabel("lift coefficient CL")
    ax2.set_xlim(0.01, 0.09)
    ax2.text(0.98, 0.5, "post-stall points beyond\nCD 0.09: see table", color=MUTED, fontsize=8, ha="right",
             transform=ax2.transAxes)
    ax2.set_title("Drag polar")
    ax2.legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    save(fig, "s1223_lift_polar")


# ----------------------------------------------------------------------------- Ahmed
def walls(run_dir, patch="body"):
    """walls.csv rows for `patch`. Its normals are face area vectors, which point
    out of the fluid (into the body), and its wall shear is the stress on the
    fluid, opposite to the near-wall flow (both checked on the Ahmed roof)."""
    data = np.genfromtxt(run_dir / "walls.csv", delimiter=",", names=True, dtype=None, encoding="ascii")
    return data[data["patch"] == patch]


def regions(w, slant_deg=35):
    """Front (nose + fillets), slant and base faces, from outward body normals."""
    x, nx, nz = w["x"], -w["nx"], -w["nz"]
    base = (x > -1e-4) & (nx > 0.99)
    slant_x = -SLANT_LEN * np.cos(np.radians(slant_deg))
    slant = (x > slant_x - 1e-4) & (nz > 0.3) & (nx > 0.3) & ~base
    front = x < -L_AHMED + 0.1 + 1e-4
    return {"front": front, "slant": slant, "base": base}


def ahmed_breakdown(r):
    w = walls(r["dir"])
    area = 0.389 * 0.288
    # Area vectors point into the body, the direction pressure pushes it.
    fx = w["cp"] * 0.5 * w["nx"] * w["area"]
    return {k: float(fx[m].sum() / (0.5 * area)) for k, m in regions(w).items()}, w


def lienhart(path):
    pts = []
    for line in open(path, encoding="latin-1"):
        s = line.split()
        if len(s) == 4 and not line.startswith(("#", "%")):
            try:
                pts.append([float(v) for v in s])
            except ValueError:
                pass
    p = np.array(pts)
    return p[:, :3] / 1000.0, p[:, 3]


def cp_compare(w, slant):
    from scipy.spatial import cKDTree
    tree = cKDTree(np.column_stack([w["x"], w["y"], w["z"]]))
    out = []
    for name, path in (("experiment", f"exp_lienhart_becker/Ahmed_{slant}_Cp.dat"),
                       ("CFX-5 SST", f"cfd11_cfx5_sst/CFD11_SST_Ahmed_{slant}_Cp.dat")):
        xyz, cp = lienhart(DATA / "ercoftac_ahmed_body" / path)
        dist, idx = tree.query(xyz)
        kutta = w["cp"][idx]
        groups = {"slant": (xyz[:, 0] < -1e-3) & (xyz[:, 2] < ROOF - 1e-3) & (xyz[:, 2] > slant_bottom(slant)),
                  "base": xyz[:, 0] > -1e-3}
        out.append((name, xyz, cp, kutta, dist, groups))
    return out


def slant_bottom(slant):
    return ROOF - SLANT_LEN * np.sin(np.radians(slant)) - 1e-3


def separated_fraction(w, slant):
    """Share of slant area whose near-wall flow points back up the slant."""
    m = regions(w, slant)["slant"]
    down = np.array([np.cos(np.radians(slant)), 0, -np.sin(np.radians(slant))])
    flow = -np.column_stack([w["tau_x"], w["tau_y"], w["tau_z"]])[m] @ down
    return float(w["area"][m][flow < 0].sum() / w["area"][m].sum())


def window(r, key, last=2000):
    """Mean, min and max of `key` over the last `last` iterations; steady runs
    that do not settle oscillate about a stationary mean (bluff-body wake)."""
    h = np.genfromtxt(r["dir"] / "history.csv", delimiter=",", names=True)
    v = h[key][h["iteration"] > h["iteration"][-1] - last]
    return float(v.mean()), float(v.min()), float(v.max())


def value(r, key):
    """Final value of a settled run, else the mean over the last 2000 iterations."""
    return r["result"][key] if r["result"]["settled"] else window(r, key)[0]


def ahmed_cd(r):
    return value(r, "drag")


STATIONS = (-203, -163, -123, -83, -43, -3)  # mm, the reference slant profile stations


def centreline_topology(w, slant):
    """'R' (reversed) or 'A' (attached) near-wall flow on the slant centreline
    at each reference station, from the wall-shear direction of the nearest faces."""
    m = regions(w, slant)["slant"] & (np.abs(w["y"]) < 0.03)
    down = np.array([np.cos(np.radians(slant)), 0, -np.sin(np.radians(slant))])
    flow = -np.column_stack([w["tau_x"], w["tau_y"], w["tau_z"]])[m] @ down
    out = []
    for x in STATIONS:
        near = np.abs(w["x"][m] - x / 1000) < 0.012
        out.append("-" if not near.any() else "R" if np.average(flow[near], weights=w["area"][m][near]) < 0 else "A")
    return " ".join(out)


def ahmed_topology(rows, lines):
    body = []
    for slant in ("25", "35"):
        for r in sorted(pick(rows, case="ahmed", point=slant, backend="metal"), key=lambda r: r["level"]):
            if (r["dir"] / "walls.csv").exists() and not r["result"]["failure"]:
                w = walls(r["dir"])
                body.append(f"{slant} | Kutta SST L{r['level']} | {centreline_topology(w, int(slant))} | "
                            f"{separated_fraction(w, int(slant)):.2f}")
    if body:
        table(lines, "Ahmed slant topology on the centreline, stations x = " + ", ".join(map(str, STATIONS)) +
              " mm (R reversed, A attached near-wall flow); reference rows in REPORT.md",
              "slant deg | source | stations | slant area with reversed near-wall flow", body)


def ahmed(rows, lines, verdicts):
    runs = sorted(pick(rows, case="ahmed", point="35", backend="metal"), key=lambda r: r["level"])
    if not runs:
        return
    body, cells, cds = [], [], []
    for r in runs:
        res = r["result"]
        cd = ahmed_cd(r)
        mean, lo, hi = window(r, "drag")
        d = cd / AHMED_CD - 1
        comp, w = ahmed_breakdown(r) if (r["dir"] / "walls.csv").exists() else ({}, None)
        sep = separated_fraction(w, 35) if w is not None else float("nan")
        body.append(f"L{r['level']} | {cd:.4f} | {lo:.4f}-{hi:.4f} | {100 * d:+.1f}% | {'PASS' if abs(d) <= 0.08 and not res['failure'] else 'FAIL'} | "
                    f"{res['drag'] - res['drag_pressure']:.4f} | {comp.get('front', math.nan):.4f} | {comp.get('slant', math.nan):.4f} | "
                    f"{comp.get('base', math.nan):.4f} | {res['lift']:.3f} | {sep:.2f} | {res['yplus_mean']:.1f} / {res['yplus_max']:.0f} | "
                    + run_cols(r))
        if not res["failure"]:
            cells.append(res["cells"]), cds.append(cd)
    finest = runs[-1]["result"]
    cd = ahmed_cd(runs[-1])
    d = cd / AHMED_CD - 1
    verdicts.append(("Ahmed 35 deg", f"L{runs[-1]['level']} ({finest['cells']:,} cells)",
                     f"Cd {cd:.3f} vs 0.26 ({100 * d:+.1f}%)", "PASS" if abs(d) <= 0.08 and not finest["failure"] else "FAIL"))
    table(lines, "Ahmed body 35 deg, Cd mesh study (Metal); reference Cd 0.26, components front 0.020 / slant 0.090 / base 0.095",
          "level | Kutta Cd (final, or mean of last 2000 iterations if unsettled) | Cd range, last 2000 | vs 0.26 | 8% gate | friction Cd | front Cp drag | slant Cp drag | base Cp drag | Kutta CL (z) | "
          "slant area with reversed shear | body y+ mean / max | " + RUN_HEAD, body)

    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.axhspan(AHMED_CD * 0.92, AHMED_CD * 1.08, color=AQUA, alpha=0.12, linewidth=0, label="0.26 +-8% (acceptance)")
    ax.axhline(AHMED_CD, color=AQUA, linewidth=1.5, linestyle="--", label="Ahmed 1984, Cd 0.26")
    ax.plot(cells, cds, "o-", color=BLUE, label="Kutta SST, levels 1-3")
    ax.set_xscale("log")
    ax.set_xlabel("cells (full model, hexahedra)")
    ax.set_ylabel("drag coefficient Cd")
    ax.set_title("Ahmed body, 35 deg slant, Re_L 2.78e6")
    ax.legend(fontsize=8.5)
    save(fig, "ahmed_cd_mesh_study")

    # Rear-surface Cp aggregates on the finest finished run.
    done = [r for r in runs if not r["result"]["failure"] and (r["dir"] / "walls.csv").exists()]
    if done:
        w = walls(done[-1]["dir"])
        body = []
        for name, xyz, cp, kutta, dist, groups in cp_compare(w, 35):
            for g, m in groups.items():
                body.append(f"{name} | {g} | {m.sum()} | {cp[m].mean():.3f} | {kutta[m].mean():.3f} | "
                            f"{np.sqrt(np.mean((kutta[m] - cp[m]) ** 2)):.3f} | {dist[m].max() * 1000:.1f}")
        table(lines, f"Ahmed 35 deg rear Cp at the reference tap locations (Kutta L{done[-1]['level']}, nearest wall face)",
              "reference | region | taps | reference mean Cp | Kutta mean Cp | RMS(Kutta - reference) | max tap-to-face distance mm", body)


# ----------------------------------------------------------------------------- convergence
def convergence(rows):
    chosen = [("Zerihan h/c 0.224, L2", pick(rows, case="zerihan", level=2, backend="metal", point="0.224")),
              ("S1223 Re 2e5, alpha 4.76, L2", [r for r in pick(rows, case="s1223", point="4.76")]),
              ("Ahmed 35 deg, finest level", sorted(pick(rows, case="ahmed", point="35"), key=lambda r: r["level"])[-1:])]
    chosen = [(t, r[0]) for t, r in chosen if r and (r[0]["dir"] / "history.csv").exists()]
    if not chosen:
        return
    fig, axes = plt.subplots(2, len(chosen), figsize=(3.6 * len(chosen), 5.4), squeeze=False)
    for i, (title, r) in enumerate(chosen):
        h = np.genfromtxt(r["dir"] / "history.csv", delimiter=",", names=True)
        h = h[np.isfinite(h["drag"])]
        for name, color in (("continuity", BLUE), ("u", ORANGE), ("k", AQUA), ("omega", VIOLET)):
            axes[0, i].semilogy(h["iteration"], h[name], color=color, linewidth=1.5, label=name)
        axes[0, i].axhline(1e-3, color=MUTED, linewidth=1, linestyle="--")
        axes[0, i].set_title(title)
        lift_name = "CL" if r["case"] != "ahmed" else "CL (z)"
        axes[1, i].plot(h["iteration"], h["drag"] / value(r, "drag"), color=BLUE, linewidth=1.5, label="Cd")
        axes[1, i].plot(h["iteration"], h["lift"] / value(r, "lift"), color=ORANGE, linewidth=1.5, label=lift_name)
        axes[1, i].set_ylim(0.9, 1.1)
        axes[1, i].set_xlabel("iteration")
    axes[0, 0].set_ylabel("Fluent-scaled residual")
    axes[1, 0].set_ylabel("coefficient / reported value")
    axes[0, 0].legend(fontsize=8)
    axes[1, 0].legend(fontsize=8)
    fig.tight_layout()
    save(fig, "convergence")


def evidence():
    dest = OUT / "runs"
    if dest.exists():
        shutil.rmtree(dest)
    for path in sorted(RUNS.glob("results*/*/result.json")):
        target = dest / path.parent.name
        target.mkdir(parents=True, exist_ok=True)
        text = path.read_text().replace(str(ROOT), "<checkout>")
        (target / "result.json").write_text(text)
        for name in ("history.csv", "solve.log"):
            src = path.parent / name
            if src.exists():
                data = src.read_bytes().replace(str(ROOT).encode(), b"<checkout>")
                with gzip.GzipFile(target / (name + ".gz"), "wb", mtime=0) as gz:
                    gz.write(data)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = results()
    lines, verdicts = ["# Generated tables (tools/fsae_aero/report.py)"], []
    zerihan(rows, lines, verdicts)
    s1223(rows, lines, verdicts)
    ahmed(rows, lines, verdicts)
    ahmed_topology(rows, lines)
    convergence(rows)
    table(lines, "Verdicts", "case | point | measured vs reference | verdict", [" | ".join(v) for v in verdicts])
    (OUT / "tables.md").write_text("\n".join(lines) + "\n")
    evidence()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
