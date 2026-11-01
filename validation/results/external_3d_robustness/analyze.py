#!/usr/bin/env python3
"""Tables for the Ahmed y+ ~ 1 mesh study (REPORT.md), from run.py directories.

Uses the FSAE aero report's definitions (tools/fsae_aero/report.py): force
components by body region, slant topology at the reference stations, the
reversed share of slant area and rear Cp at the reference tap locations.
Reference curves are read at run time from validation/data and never written.

  build-py/bin/python validation/results/external_3d_robustness/analyze.py runs/fsae_aero/study
"""
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/fsae_aero"))
import report  # noqa: E402

# Experiment Cd by slant: 35 deg as in report.py, 25 deg as quoted in the FSAE aero report.
EXPERIMENT_CD = {35: report.AHMED_CD, 25: 0.285}


def run_rows(study):
    rows = []
    for path in sorted(pathlib.Path(study).glob("*/result.json")):
        r = json.loads(path.read_text())
        r["dir"] = path.parent
        if "result" in r:
            rows.append(r)
    return rows


def history(r):
    return np.genfromtxt(r["dir"] / "history.csv", delimiter=",", names=True)


def summary(r):
    res, h = r["result"], history(r)
    settled = res["settled"]
    last = h["iteration"] > h["iteration"][-1] - 2000
    cd = res["drag"] if settled else float(h["drag"][last].mean())
    cl = res["lift"] if settled else float(h["lift"][last].mean())
    friction = res["drag"] - res["drag_pressure"]
    row = dict(scheme="Coupled" if "coupled" in r.get("tag", "") else "SIMPLEC", slant=r["point"], level=r["level"],
               cells=res["cells"], poor=res.get("poor_cells", 0), iterations=res["iterations"], seconds=res["solve_s"],
               gb=r.get("peak_footprint_bytes", 0) / 1e9, cd=cd, cd_lo=float(h["drag"][last].min()),
               cd_hi=float(h["drag"][last].max()), cl=cl, friction=friction, residual=res["max_residual"],
               raw=float(h["continuity_raw"][-1]), settled=settled, met=res["residuals_met"],
               yplus=(res["yplus_mean"], res["yplus_max"]), failure=res["failure"])
    if (r["dir"] / "walls.csv").exists() and not res["failure"]:
        comp, w = report.ahmed_breakdown(r)
        slant = int(float(r["point"]))
        row.update(components=comp, topology=report.centreline_topology(w, slant),
                   reversed=report.separated_fraction(w, slant), cp=cp_rms(w, slant))
    return row


def cp_rms(w, slant):
    out = {}
    for name, xyz, cp, kutta, dist, groups in report.cp_compare(w, slant):
        for g, m in groups.items():
            m = m & (dist < 0.006)
            out[f"{name} {g}"] = (int(m.sum()), float(cp[m].mean()), float(kutta[m].mean()),
                                  float(np.sqrt(np.mean((kutta[m] - cp[m]) ** 2))))
    return out


def main():
    rows = [summary(r) for r in run_rows(sys.argv[1])]
    rows.sort(key=lambda s: (s["slant"], s["scheme"], s["level"]))
    print("| slant | scheme | level | cells | flagged cells | iterations | solve s | peak GB | Cd | Cd range, last 2000 | "
          "vs experiment | friction Cd | front / slant / base pressure Cd | CL | max residual | raw continuity | stop | y+ mean / max |")
    print("|" + " --- |" * 18)
    for s in rows:
        c = s.get("components", {})
        stop = "failed: " + s["failure"][:60] if s["failure"] else (
            "converged" if s["met"] and s["settled"] else "residuals met, forces drifting" if s["met"] else "not converged")
        print(f"| {s['slant']} | {s['scheme']} | L{s['level']:g} | {s['cells']:,} | {s['poor']:,} | {s['iterations']} | "
              f"{s['seconds']:.0f} | {s['gb']:.2f} | {s['cd']:.4f} | {s['cd_lo']:.4f}-{s['cd_hi']:.4f} | "
              f"{100 * (s['cd'] / EXPERIMENT_CD[int(float(s['slant']))] - 1):+.1f}% | {s['friction']:.3f} | "
              f"{c.get('front', float('nan')):.3f} / {c.get('slant', float('nan')):.3f} / {c.get('base', float('nan')):.3f} | "
              f"{s['cl']:.3f} | {s['residual']:.1e} | {s['raw']:.1e} | {stop} | {s['yplus'][0]:.2f} / {s['yplus'][1]:.2f} |")
    print()
    print("| slant | scheme | level | centreline stations " + " ".join(map(str, report.STATIONS)) + " mm | reversed slant area |")
    print("| --- | --- | --- | --- | --- |")
    for s in rows:
        if "topology" in s:
            print(f"| {s['slant']} | {s['scheme']} | L{s['level']:g} | {s['topology']} | {s['reversed']:.2f} |")
    print()
    print("| slant | scheme | level | reference | region | taps | reference mean Cp | Kutta mean Cp | RMS difference |")
    print("| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: |")
    for s in rows:
        for key, (n, ref, kutta, rms) in s.get("cp", {}).items():
            name, region = key.rsplit(" ", 1)
            print(f"| {s['slant']} | {s['scheme']} | L{s['level']:g} | {name} | {region} | {n} | {ref:.3f} | {kutta:.3f} | {rms:.3f} |")


if __name__ == "__main__":
    main()
