#!/usr/bin/env python3
"""Plots, tables and gzipped evidence for validation/results/transition.

Reads the runs written by flatplate.py (runs/transition/flatplate), run.py with
--transition (runs/transition/s1223) and bench_tmr's high-tu modes
(runs/transition/tmr). The fully turbulent S1223 baseline is the tracked
evidence in validation/results/fsae_aero/runs. Reference values are read at
run time from the restricted files under validation/data; plots and tables
carry only the values quoted in REPORT.md (T3A summary Cf, UIUC values at the
run angles, measured T3B/T3A- onset).

  build-py/bin/python tools/transition/report.py
"""
import gzip
import importlib.util
import json
import pathlib
import shutil
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/transition"))
import flatplate  # noqa: E402

# tools/fsae_aero/report.py (S1223 reference, styling); loaded by path because
# this module has the same name.
_spec = importlib.util.spec_from_file_location("fsae_aero_report", ROOT / "tools/fsae_aero/report.py")
aero = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(aero)

RUNS = ROOT / "runs/transition"
OUT = ROOT / "validation/results/transition"
BASELINE = ROOT / "validation/results/fsae_aero/runs"
BLUE, ORANGE, AQUA, VIOLET, INK, MUTED, GRID = aero.BLUE, aero.ORANGE, aero.AQUA, aero.VIOLET, aero.INK, aero.MUTED, aero.GRID
CASES = {"t3a": "T3A", "t3b": "T3B", "t3am": "T3A-"}


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def result(run_dir):
    """The RESULT line of a flatplate.py run log."""
    for line in reversed((pathlib.Path(str(run_dir) + ".log")).read_text().splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[7:])["result"]
    return None


# Re_x windows: the onset search ends before full turbulence, the end search
# covers the transition region. T3A- transitions late (measured Cf is still
# rising at its last station, so its measured end is a lower bound).
WINDOW = {"t3a": (4e5, 8e5), "t3b": (4e5, 8e5), "t3am": (1.9e6, 2.1e6)}


def onset(rex, cf, case):
    """Re_x of the Cf minimum between the leading-edge region and full turbulence."""
    sel = (rex > 2e4) & (rex < WINDOW[case][0])
    return float(rex[sel][np.argmin(cf[sel])])


def end(rex, cf, start, case):
    """Re_x of the Cf maximum after onset (end of transition)."""
    sel = (rex > start) & (rex < WINDOW[case][1])
    return float(rex[sel][np.argmax(cf[sel])])


# ----------------------------------------------------------------------------- flat plates
def plates(lines):
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9), sharey=False)
    rows = []
    for ax, (case, name) in zip(axes, CASES.items()):
        base = RUNS / "flatplate"
        exp = flatplate.summary(case)
        m_on = onset(exp[:, 1], exp[:, 3], case)
        m_end = end(exp[:, 1], exp[:, 3], m_on, case)
        fit = json.loads((base / f"{case}_fit.json").read_text()) if (base / f"{case}_fit.json").exists() else None
        rex = np.geomspace(5e3, 2.5e6, 200)
        ax.plot(rex, 0.664 / np.sqrt(rex), "--", color=MUTED, linewidth=1.2, label="Blasius laminar")
        ax.plot(rex, 0.0592 * rex ** -0.2, ":", color=MUTED, linewidth=1.2, label="turbulent 0.0592 Re_x^-0.2")
        for model, color, label in (("gamma", BLUE, "Kutta SST + intermittency"), ("sst", ORANGE, "Kutta SST, fully turbulent")):
            d = base / f"{case}_L1_{model}"
            if not (d / "walls.csv").exists():
                continue
            comp, _ = flatplate.cf_table(case, d)
            comp = comp[comp[:, 0] > 3e3]
            ax.plot(comp[:, 0], comp[:, 1], color=color, label=label)
            res = result(d)
            row = {"case": name, "model": model, "result": res}
            if model == "gamma":
                on = onset(comp[:, 0], comp[:, 1], case)
                row.update(onset=on, end=end(comp[:, 0], comp[:, 1], on, case), measured_onset=m_on, measured_end=m_end,
                           end_bound=m_end >= exp[-1, 1])
                errs = [(r, np.interp(r, comp[:, 0], comp[:, 1]) / c - 1) for r, c in zip(exp[:, 1], exp[:, 3])]
                row["laminar"] = max(abs(e) for r, e in errs if r < m_on)
                after = [abs(e) for r, e in errs if r > (2.6e5 if case == "t3a" else m_end * 1.05)]
                row["turbulent"] = max(after) if after else None
                row["errors"] = errs
            rows.append(row)
        if case == "t3a":
            ax.plot(exp[:, 1], exp[:, 3], "o", color=INK, fillstyle="none", markeredgewidth=1.4, label="Coupland, measured")
        else:
            ax.axvline(m_on, color=INK, linewidth=1)
            ax.text(m_on * 1.06, 0.0075, "measured onset\n(Cf minimum)", color=INK, fontsize=8, va="top")
        ax.set_xscale("log")
        ax.set_xlim(1e4, 2.5e6 if case == "t3am" else 1.1e6)
        ax.set_ylim(0, 0.008)
        ax.set_xlabel("Re_x (local freestream)")
        tu = f", fitted Tu at LE {100 * fit['tu_leading_edge']:.2f}%" if fit else ""
        ax.set_title(f"{name}{tu}")
    axes[0].set_ylabel("skin friction Cf")
    axes[0].legend(fontsize=7.2, loc="upper right")
    fig.tight_layout()
    save(fig, "t3_cf")
    lines.append("### Flat plates (L1, Metal): transition onset and Cf error vs measured\n")
    lines.append("| case | onset Re_x Kutta | onset measured | error | end Re_x Kutta | end measured | max laminar Cf error | max turbulent Cf error | iterations | max scaled residual | drag change |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        if r["model"] != "gamma":
            continue
        res = r["result"]
        turbulent = "n/a" if r["turbulent"] is None else f"{100 * r['turbulent']:.1f}%"
        bound = ">= " if r["end_bound"] else ""
        lines.append(f"| {r['case']} | {r['onset']:.3g} | {r['measured_onset']:.3g} | {100 * (r['onset'] / r['measured_onset'] - 1):+.0f}% | "
                     f"{r['end']:.3g} | {bound}{r['measured_end']:.3g} | {100 * r['laminar']:.1f}% | {turbulent} | "
                     f"{res['iterations']} | {res['max_residual']:.1e} | {res['settle_drag']:.1e} |")
    t3a = next((r for r in rows if r["case"] == "T3A" and r["model"] == "gamma"), None)
    if t3a:
        sst = next((r for r in rows if r["case"] == "T3A" and r["model"] == "sst"), None)
        sst_comp = flatplate.cf_table("t3a", RUNS / "flatplate/t3a_L1_sst")[0] if sst else None
        lines.append("\n### T3A Cf at the measured stations\n")
        lines.append("| Re_x | measured Cf | Kutta gamma | error | Kutta SST | error |")
        lines.append("|---|---|---|---|---|---|")
        exp = flatplate.summary("t3a")
        for (rex, e), c in zip(t3a["errors"], exp[:, 3]):
            g = c * (1 + e)
            s = "| n/a | n/a |" if sst_comp is None else \
                f"| {np.interp(rex, sst_comp[:, 0], sst_comp[:, 1]):.5f} | {100 * (np.interp(rex, sst_comp[:, 0], sst_comp[:, 1]) / c - 1):+.1f}% |"
            lines.append(f"| {rex:.4g} | {c:.5f} | {g:.5f} | {100 * e:+.1f}% {s}")
    return rows


# ----------------------------------------------------------------------------- S1223
def s1223_runs(directory, tag):
    out = {}
    for path in sorted(directory.glob(f"s1223_re*_L2_metal{tag}/result.json")):
        r = json.loads(path.read_text())
        r["dir"] = path.parent
        if "result" in r:
            out[(float(r["reynolds"]), float(r["point"]))] = r
    return out


def s1223(lines):
    gamma = s1223_runs(RUNS / "s1223", "_gamma")
    sst = s1223_runs(BASELINE, "")
    if not gamma:
        return
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.6))
    for row, re in enumerate((2e5, 3e5)):
        at, clmax = aero.s1223_reference(re)
        ax1, ax2 = axes[row]
        body = []
        series = {}
        for name, runs in (("sst", sst), ("gamma", gamma)):
            pts = sorted((a, r) for (rr, a), r in runs.items() if abs(rr - re) < 1)
            series[name] = [(a, aero.value(r, "lift"), aero.value(r, "drag"), r) for a, r in pts]
        for a, cl, cd, r in series["gamma"]:
            cl_ref, cd_ref = at(a)
            s = next((x for x in series["sst"] if abs(x[0] - a) < 1e-6), None)
            res = r["result"]
            mark = "" if res["settled"] else "*"
            linear = a <= aero.S1223_LINEAR
            gl = ("PASS" if abs(cl - cl_ref) <= .10 else "FAIL") if linear else "trend"
            gd = "n/a" if cd_ref is None else ("PASS" if abs(cd / cd_ref - 1) <= .25 else "FAIL")
            lsb = bubble(r["dir"])
            body.append(f"{a:g} | {cl:.3f}{mark} | {'' if s is None else f'{s[1]:.3f}'} | {cl_ref:.3f} | {cl - cl_ref:+.3f} | {gl} | "
                        f"{cd:.4f}{mark} | {'' if s is None else f'{s[2]:.4f}'} | "
                        + ("n/a | n/a" if cd_ref is None else f"{cd_ref:.4f} | {100 * (cd / cd_ref - 1):+.0f}%")
                        + f" | {gd} | {lsb} | " + aero.run_cols(r))
        lines.append(f"\n### Selig S1223, Re {re:.0e} (level 2, Metal), SST + intermittency vs fully turbulent SST and UIUC; "
                     f"UIUC CLmax {clmax[1]:.3f} at {clmax[0]:.1f} deg; * = unsettled, mean of the last 2000 iterations\n")
        lines.append("| alpha | gamma CL | SST CL | UIUC CL | dCL (gamma) | +-0.10 gate | gamma CD | SST CD | UIUC CD | dCD (gamma) | +-25% gate | "
                     "suction-side reversed flow (x/c) | " + aero.RUN_HEAD + " |")
        lines.append("|" + "---|" * 18)
        lines.extend("| " + b + " |" for b in body)
        label = f"Re {re / 1e5:.0f}e5"
        ref = [(a, *at(a)) for a, *_ in series["gamma"]]
        ax1.plot([p[0] for p in ref], [p[1] for p in ref], "o", color=INK, fillstyle="none", markeredgewidth=1.4,
                 label="UIUC experiment (at run angles)")
        sel = [p for p in ref if p[2] is not None]
        ax2.plot([p[2] for p in sel], [p[1] for p in sel], "o", color=INK, fillstyle="none", markeredgewidth=1.4,
                 label="UIUC experiment")
        for name, color, text in (("sst", ORANGE, "Kutta SST, fully turbulent"), ("gamma", BLUE, "Kutta SST + intermittency")):
            p = series[name]
            ax1.plot([x[0] for x in p], [x[1] for x in p], "-", marker="s" if name == "sst" else "D", color=color, label=text)
            ax2.plot([x[2] for x in p], [x[1] for x in p], "-", marker="s" if name == "sst" else "D", color=color, label=text)
        ax1.axvspan(aero.S1223_LINEAR, 17, color=GRID, alpha=0.5, linewidth=0)
        ax1.text(aero.S1223_LINEAR + 0.3, 0.97, "trend only", color=MUTED, fontsize=8.5, va="top", transform=ax1.get_xaxis_transform())
        ax1.set_xlabel("angle of attack, deg")
        ax1.set_ylabel("lift coefficient CL")
        ax1.set_title(f"Lift curve, {label}")
        ax2.set_xlabel("drag coefficient CD")
        ax2.set_ylabel("lift coefficient CL")
        ax2.set_xlim(0.0, 0.09)
        ax2.set_title(f"Drag polar, {label}")
        ax1.legend(fontsize=7.5, loc="lower right")
        ax2.legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    save(fig, "s1223_transition")


def bubble(run_dir):
    """Chord-fraction spans of reversed near-wall flow on the suction (upper)
    surface aft of 2% chord: a "bubble" reattaches before 98% chord, "TE
    separation" reaches it; "none" if attached. The chord runs from the
    trailing edge (largest x) to the wall point farthest from it; walls.csv
    shear is the stress on the fluid, so the near-wall flow is along -tau."""
    path = run_dir / "walls.csv"
    if not path.exists():
        return "n/a"
    w = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="ascii")
    w = w[w["patch"] == "wing"]
    p = np.stack([w["x"], w["y"]], axis=1)
    te = p[np.argmax(p[:, 0])]
    le = p[np.argmax(np.hypot(*(p - te).T))]
    c = te - le
    xc = (p - le) @ c / (c @ c)
    # Suction side: the into-body normal points to the chord's right (the
    # cambered pressure side lies partly above the chord line, so position
    # alone cannot separate the surfaces).
    upper = -c[1] * w["nx"] + c[0] * w["ny"] < 0
    # Surface tangent oriented from leading to trailing edge.
    t = np.stack([-w["ny"], w["nx"]], axis=1)
    t *= np.sign(t @ c)[:, None]
    flow = -(w["tau_x"] * t[:, 0] + w["tau_y"] * t[:, 1])
    sel = upper & (xc > 0.02)
    order = np.argsort(xc[sel])
    xs, rev = xc[sel][order], (flow[sel] < 0)[order]
    spans, start, prev = [], None, None
    for x, r in zip(xs, rev):
        if r and start is None:
            start = x
        if not r and start is not None:
            spans.append((start, prev))
            start = None
        prev = x
    if start is not None:
        spans.append((start, prev))
    out = [f"{a:.2f}-{b:.2f} " + ("TE separation" if b >= 0.98 else "bubble") for a, b in spans if b - a > 0.005]
    return ", ".join(out) or "none"


# ----------------------------------------------------------------------------- TMR plate
def tmr(lines):
    d = RUNS / "tmr"
    runs = {}
    for mode in ("high-tu", "high-tu-transition"):
        log = d / f"{mode}.log"
        if not log.exists() or not (d / mode / "wall_cf.csv").exists():
            return
        res = next(line for line in reversed(log.read_text().splitlines()) if line.startswith("RESULT"))
        values = {}
        for kv in res.split()[1:]:  # repeated keys (cf, cp per station): keep the first, the plate station
            if "=" in kv:
                values.setdefault(*kv.split("=", 1))
        cf = np.genfromtxt(d / mode / "wall_cf.csv", delimiter=",", names=True, dtype=None, encoding="ascii")
        cf = np.sort(cf[cf["x"] > 0], order="x")
        runs[mode] = (values, cf)
    (s, cs), (g, cg) = runs["high-tu"], runs["high-tu-transition"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.8))
    ax1.plot(cs["x"], cs["cf"], color=ORANGE, label="SST, fully turbulent")
    ax1.plot(cg["x"], cg["cf"], color=BLUE, label="SST + intermittency")
    ax1.set_xscale("log")
    ax1.set_ylim(0, 0.008)
    ax1.set_xlabel("x (plate length 2, Re 5e6 per unit length)")
    ax1.set_ylabel("skin friction Cf")
    ax1.set_title("TMR flat plate, inlet Tu 5%, mu_t/mu 100")
    ax1.legend(fontsize=8)
    ratio = np.interp(cs["x"], cg["x"], cg["cf"]) / cs["cf"] - 1
    ax2.plot(cs["x"], 100 * ratio, color=BLUE)
    ax2.axhline(0, color=MUTED, linewidth=1)
    ax2.set_xscale("log")
    ax2.set_ylim(-20, 10)
    ax2.set_xlabel("x")
    ax2.set_ylabel("Cf difference, intermittency vs SST (%)")
    ax2.set_title("Difference")
    fig.tight_layout()
    save(fig, "tmr_high_tu")
    lines.append("\n### TMR flat plate at inlet Tu 5%, mu_t/mu 100 (137x97, Metal): intermittency vs fully turbulent SST\n")
    lines.append("| quantity | SST | SST + intermittency | difference |")
    lines.append("|---|---|---|---|")
    for key, name in (("cf", "Cf at x = 0.97"), ("drag", "plate drag coefficient")):
        a, b = float(s[key]), float(g[key])
        lines.append(f"| {name} | {a:.6f} | {b:.6f} | {100 * (b / a - 1):+.2f}% |")
    for x in (0.01, 0.05, 0.1, 0.5, 1.0, 1.9):
        a, b = np.interp(x, cs["x"], cs["cf"]), np.interp(x, cg["x"], cg["cf"])
        lines.append(f"| Cf at x = {x:g} | {a:.6f} | {b:.6f} | {100 * (b / a - 1):+.2f}% |")
    for name, v in (("SST", s), ("SST + intermittency", g)):
        lines.append(f"\n{name}: {v['iterations']} iterations, defect {float(v['defect']):.2e}, force change {float(v['change']):.1e}, "
                     f"solve {float(v['solve_s']):.0f} s.")


# ----------------------------------------------------------------------------- evidence
def scrub(data):
    """Run logs and results name absolute paths; keep them relative to the checkout."""
    return data.replace(str(ROOT).encode() + b"/", b"./").replace(str(ROOT).encode(), b".")


def store(src, dst, compress):
    data = scrub(pathlib.Path(src).read_bytes())
    if compress:
        with gzip.GzipFile(dst, "wb", mtime=0) as g:
            g.write(data)
    else:
        pathlib.Path(dst).write_bytes(data)


def evidence():
    """Histories, logs and results next to REPORT.md (no meshes, no fields)."""
    dest = OUT / "runs"
    shutil.rmtree(dest, ignore_errors=True)
    dest.mkdir()
    for d in sorted((RUNS / "s1223").glob("s1223_*")):
        out = dest / d.name
        out.mkdir()
        for name, compress in (("result.json", False), ("history.csv", True), ("solve.log", True)):
            if (d / name).exists():
                store(d / name, out / (name + (".gz" if compress else "")), compress)
    for d in sorted((RUNS / "flatplate").glob("t3*_L*_*")):
        if not d.is_dir():
            continue
        out = dest / d.name
        out.mkdir()
        for src, name in ((d / "history.csv", "history.csv.gz"), (pathlib.Path(str(d) + ".log"), "solve.log.gz")):
            if src.exists():
                store(src, out / name, True)
    for f in sorted((RUNS / "flatplate").glob("*_fit.json")):
        store(f, dest / f.name, False)
    for mode in ("high-tu", "high-tu-transition"):
        out = dest / f"tmr_{mode}"
        out.mkdir()
        for src, name in ((RUNS / "tmr" / f"{mode}.log", "solve.log.gz"), (RUNS / "tmr" / mode / "history.csv", "history.csv.gz"),
                          (RUNS / "tmr" / mode / "wall_cf.csv", "wall_cf.csv.gz")):
            if src.exists():
                store(src, out / name, True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    plates(lines)
    s1223(lines)
    tmr(lines)
    evidence()
    (OUT / "tables.md").write_text("<!-- generated by tools/transition/report.py -->\n" + "\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
