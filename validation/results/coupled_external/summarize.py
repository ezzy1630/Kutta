#!/usr/bin/env python3
"""Coupled vs SIMPLEC table from run_fsae.sh output (stdlib only).

SIMPLEC runs stop on bench_fsae_aero's criterion: every Fluent-scaled residual
<= 1e-3 and lift and drag settled to 5e-4 (relative) over the window (500
iterations, 1000 for Ahmed). Coupled's Fluent-scaled continuity is normalized
by a start-up imbalance 25-230 times smaller than SIMPLEC's, so that criterion
asks Coupled for a lower raw continuity than SIMPLEC finishes at (REPORT.md).
Coupled's convergence iteration is therefore the first sample where:
  - lift and drag are settled by the same rule and window;
  - u, v (w in 3D), k and omega Fluent-scaled residuals are <= 1e-3;
  - raw continuity is at or below SIMPLEC's final raw continuity on that mesh.

  python3 validation/results/coupled_external/summarize.py <run_fsae.sh out dir>
"""
import csv
import json
import pathlib
import sys


def history(path):
    with open(path) as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def result(log):
    for line in open(log):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT "):])["result"]
    return None


def settled(rows, i, window, tol=5e-4):
    last = rows[i]
    if rows[0]["iteration"] > last["iteration"] - window:
        return False
    for row in rows[: i + 1]:
        if last["iteration"] - row["iteration"] > window:
            continue
        for key in ("lift", "drag"):
            if abs(row[key] - last[key]) > tol * max(abs(last[key]), 1e-12):
                return False
    return True


def converged_at(rows, window, continuity, three_d):
    keys = ["u", "v", "k", "omega"] + (["w"] if three_d else [])
    for i, row in enumerate(rows):
        if all(row[k] <= 1e-3 for k in keys) and row["continuity_raw"] <= continuity and settled(rows, i, window):
            return row
    return None


def main(out):
    out = pathlib.Path(out)
    tags = sorted({p.name.rsplit("_", 1)[0] for p in out.glob("*_simplec.log")})
    print("| case | cells | SIMPLEC its / s | Coupled its / s | speedup (time) | SIMPLEC lift / drag | Coupled lift / drag | lift diff | drag diff | raw continuity S / C |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for tag in tags:
        s_log, c_log = out / f"{tag}_simplec.log", out / f"{tag}_coupled.log"
        s, c = result(s_log), result(c_log) if c_log.exists() else None
        if s is None or c is None:
            print(f"| {tag} | | {'running' if s is None else 'done'} | {'running' if c is None else 'done'} | | | | | | |")
            continue
        sh, ch = history(out / f"{tag}_simplec/history.csv"), history(out / f"{tag}_coupled/history.csv")
        if c["failure"] or not ch:
            last = int(ch[-1]["iteration"]) if ch else 0
            print(f"| {tag} | {c['cells']} | {s['iterations']} / {s['solve_s']:.0f} | failed before iteration {last + 50}: "
                  f"{c['failure'].rsplit(': ', 1)[-1]} | | {s['lift']:.6f} / {s['drag']:.6f} | | | | |")
            continue
        window = 1000 if tag.startswith("ahmed") else 500
        three_d = tag.startswith("ahmed")
        s_ok = s["residuals_met"] and s["settled"]
        s_text = f"{s['iterations']} / {s['solve_s']:.0f}" + ("" if s_ok else " (not converged)")
        hit = converged_at(ch, window, sh[-1]["continuity_raw"], three_d)
        c_text = f"{int(hit['iteration'])} / {hit['seconds']:.0f}" if hit else f"not in {c['iterations']}"
        speed = f"{s['solve_s'] / hit['seconds']:.2f}x" if hit and s_ok else ""
        dl = c["lift"] / s["lift"] - 1
        dd = c["drag"] / s["drag"] - 1
        print(f"| {tag} | {c['cells']} | {s_text} | {c_text} | {speed} | {s['lift']:.6f} / {s['drag']:.6f} | "
              f"{c['lift']:.6f} / {c['drag']:.6f} | {dl:+.1e} | {dd:+.1e} | "
              f"{sh[-1]['continuity_raw']:.2e} / {ch[-1]['continuity_raw']:.2e} |")


if __name__ == "__main__":
    main(sys.argv[1])
