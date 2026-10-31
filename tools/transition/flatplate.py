"""ERCOFTAC T3A / T3B / T3A- zero-pressure-gradient flat-plate transition.

Subcommands:
  fit   fit inlet Tu and viscosity ratio to the measured freestream Tu decay
  mesh  write a structured 2D Fluent .msh of the plate (tools/fsae_aero/fluent_msh.py)
  run   mesh (if missing), fit and solve through build/bench_fsae_aero
  cf    skin friction vs Re_x from a run's walls.csv, with the experiment

Nondimensional as in validation/cases/fsae_aero.hpp: U = rho = 1, lengths in
metres, nu = 1/Re_unit. Re_unit is the median of REX/X in the experiment's
summary table, so the computed Re_x = x*Re_unit uses the same velocity scale
as the measured Re_x. The plate leading edge is x = 0; the inlet is X_IN
upstream, over a symmetry plane. The top boundary is a symmetry plane.

Inlet turbulence is not given by the experiment, only its decay. `fit` fits
(Tu_in, mu_t/mu) to the measured TUFS(x) with SST's freestream decay
(F1 = 0: beta = .0828, beta* = .09), the usual practice for these cases.

The reference data are restricted (CC BY-NC-SA, see
validation/data/ercoftac_t3_flat_plate/SOURCE.md) and are read at run time.
"""
import argparse
import json
import math
import pathlib
import subprocess
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/fsae_aero"))
import fluent_msh  # noqa: E402

DATA = ROOT / "validation/data/ercoftac_t3_flat_plate/raw"
FILES = {"t3a": "t3ay.dat", "t3b": "t3by.dat", "t3am": "t3amy.dat"}
X_IN = 0.1      # inlet distance upstream of the leading edge (m)
LENGTH = 1.7    # plate length (m)
HEIGHT = 0.3    # domain height (m)
BETA, BETA_STAR = 0.0828, 0.09


def summary(case):
    """Rows of the experiment's summary table: x (m), Re_x, U0, Cf, H, Re_theta, Tu (fraction)."""
    rows = []
    for line in (DATA / FILES[case]).read_text().splitlines():
        parts = line.split()
        if len(parts) == 9 and parts[0].startswith("1W"):
            x, rex, u0, cf, h, rt, _, tu = map(float, parts[1:])
            rows.append((x / 1000, rex, u0, cf, h, rt, tu / 100))
    if not rows:
        raise SystemExit(f"no summary rows in {DATA / FILES[case]}")
    return np.array(rows)


def reynolds_unit(case):
    s = summary(case)
    return float(np.median(s[:, 1] / s[:, 0]))


def decay(x, tu_in, ratio, nu):
    """Freestream SST decay from the inlet at x = -X_IN with U = 1."""
    k0 = 1.5 * tu_in**2
    w0 = k0 / (nu * ratio)
    t = np.asarray(x) + X_IN
    return np.sqrt(2.0 / 3.0 * k0 * (1 + BETA * w0 * t) ** (-BETA_STAR / BETA))


def fit(case):
    s = summary(case)
    nu = 1 / reynolds_unit(case)
    x, tu = s[:, 0], s[:, 6]
    best = None
    # Coarse log grid, then a local refinement; the error is smooth in both.
    for tu_in in np.geomspace(0.5 * tu[0], 4 * tu[0], 120):
        for ratio in np.geomspace(0.5, 2000, 160):
            err = float(np.sum(np.log(decay(x, tu_in, ratio, nu) / tu) ** 2))
            if best is None or err < best[0]:
                best = (err, tu_in, ratio)
    _, tu_in, ratio = best
    for _ in range(3):
        for a in np.linspace(0.95, 1.05, 41) * tu_in:
            for b in np.geomspace(0.9, 1.1, 41) * ratio:
                err = float(np.sum(np.log(decay(x, a, b, nu) / tu) ** 2))
                if err < best[0]:
                    best = (err, a, b)
        _, tu_in, ratio = best
    rms = math.sqrt(best[0] / len(x))
    return {"case": case, "re_unit": 1 / nu, "x_inlet": -X_IN, "tu_inlet": float(tu_in), "viscosity_ratio": float(ratio),
            "tu_leading_edge": float(decay(0.0, tu_in, ratio, nu)), "rms_log_error": rms}


def stretched(n, first, length):
    """n cells over `length`, first cell `first`, geometric growth (solved)."""
    lo, hi = 1.0 + 1e-9, 2.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        total = first * (r**n - 1) / (r - 1)
        lo, hi = (r, hi) if total < length else (lo, r)
    r = 0.5 * (lo + hi)
    d = first * r ** np.arange(n)
    return np.concatenate([[0.0], np.cumsum(d * length / d.sum())])


def mesh(case, level, out):
    re = reynolds_unit(case)
    # First cell y+ about 0.3 for a turbulent Cf of .006, at every Re_unit.
    first_y = 0.3 / (re * math.sqrt(0.003)) / level
    scale = level
    ny, nx_plate, nx_up = 120 * scale, 400 * scale, 40 * scale
    y = stretched(ny, first_y, HEIGHT)
    ds_le = 5e-4 / level
    xp = stretched(nx_plate, ds_le, LENGTH)
    xu = -stretched(nx_up, ds_le, X_IN)[::-1]
    x = np.concatenate([xu[:-1], xp])
    nxn, nyn = len(x), len(y)
    X, Y = np.meshgrid(x, y, indexing="ij")
    points = np.stack([X.ravel(), Y.ravel(), np.zeros(X.size)], axis=1)
    node = lambda i, j: i * nyn + j  # noqa: E731
    I, J = np.meshgrid(np.arange(nxn - 1), np.arange(nyn - 1), indexing="ij")
    I, J = I.ravel(), J.ravel()
    quads = np.stack([node(I, J), node(I + 1, J), node(I + 1, J + 1), node(I, J + 1)], axis=1)
    le = nx_up  # node index of the leading edge
    i = np.arange(nxn - 1)
    j = np.arange(nyn - 1)
    bottom = np.stack([node(i, 0), node(i + 1, 0)], axis=1)
    boundaries = [
        ("plate", "wall", bottom[le:]),
        ("upstream", "symmetry", bottom[:le]),
        ("top", "symmetry", np.stack([node(i, nyn - 1), node(i + 1, nyn - 1)], axis=1)),
        ("inlet", "velocity-inlet", np.stack([node(0, j), node(0, j + 1)], axis=1)),
        ("outlet", "pressure-outlet", np.stack([node(nxn - 1, j), node(nxn - 1, j + 1)], axis=1)),
    ]
    if abs(x[le]) > 1e-15:
        raise SystemExit("leading edge node is not at x = 0")
    stats = fluent_msh.write(out, points, [(3, quads)], boundaries, comment=f"transition flatplate {case}")
    meta = dict(stats, case=case, level=level, re_unit=re, first_y=first_y, ds_le=ds_le,
                cells=int(len(quads)), nx=nxn - 1, ny=nyn - 1)
    pathlib.Path(str(out) + ".json").write_text(json.dumps(meta, indent=1))
    return meta


def edge_velocity(run_dir):
    """Local freestream speed U_e(x): the largest u in each cell column below
    y = 0.1 m (outside every boundary layer here, inside the domain's
    slightly accelerating core)."""
    d = np.genfromtxt(pathlib.Path(run_dir) / "cells.csv", delimiter=",", names=True, usecols=("x", "y", "u"))
    d = d[d["y"] < 0.1]
    xs, inverse = np.unique(d["x"], return_inverse=True)
    ue = np.full(len(xs), -np.inf)
    np.maximum.at(ue, inverse, d["u"])
    return xs, ue


def cf_table(case, run_dir):
    """Computed (Re_x, Cf) on plate faces with the experiment's definitions:
    both use the local freestream velocity U_e(x), as the summary's U0 and REX
    do. Returns the computed rows and the experiment's summary."""
    import csv
    re = reynolds_unit(case)
    xs, ue = edge_velocity(run_dir)
    rows = []
    with open(pathlib.Path(run_dir) / "walls.csv") as f:
        for r in csv.DictReader(f):
            if r["patch"] == "plate":
                x = float(r["x"])
                u = float(np.interp(x, xs, ue))
                # tau is the stress on the fluid (antiparallel to the flow).
                rows.append((x * re * u, -2.0 * float(r["tau_x"]) / u**2))
    rows.sort()
    return np.array(rows), summary(case)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=["fit", "mesh", "run", "cf"])
    ap.add_argument("case", choices=list(FILES))
    ap.add_argument("--level", type=int, default=1, help="mesh refinement factor in each direction")
    ap.add_argument("--out", default=str(ROOT / "runs/transition/flatplate"))
    ap.add_argument("--backend", default="metal")
    ap.add_argument("--turbulent", action="store_true", help="fully turbulent SST (no transition model)")
    ap.add_argument("--iterations", type=int, default=8000)
    ap.add_argument("--max-seconds", type=float, default=3600)
    ap.add_argument("--extra", default="")
    args = ap.parse_args()
    out = pathlib.Path(args.out)
    if args.command == "fit":
        print(json.dumps(fit(args.case), indent=1))
        return
    msh = out / f"{args.case}_L{args.level}.msh"
    if args.command == "mesh" or (args.command == "run" and not msh.exists()):
        out.mkdir(parents=True, exist_ok=True)
        print(json.dumps(mesh(args.case, args.level, msh)))
        if args.command == "mesh":
            return
    tag = f"{args.case}_L{args.level}_{'sst' if args.turbulent else 'gamma'}"
    if args.command == "cf":
        comp, exp = cf_table(args.case, out / tag)
        for rex, cf in comp[:: max(1, len(comp) // 40)]:
            print(f"{rex:12.4e} {cf:10.6f}")
        return
    f = fit(args.case)
    (out / f"{args.case}_fit.json").write_text(json.dumps(f, indent=1))
    cmd = [str(ROOT / "build/bench_fsae_aero"), "--mesh", str(msh), "--backend", args.backend, "--out", str(out / tag),
           "--re", repr(f["re_unit"]), "--length", "1", "--area", "1", "--walls", "plate",
           "--tu", repr(f["tu_inlet"]), "--ratio", repr(f["viscosity_ratio"]),
           "--iterations", str(args.iterations), "--max-seconds", str(args.max_seconds), "--cells"]
    if not args.turbulent:
        cmd.append("--transition")
    cmd += args.extra.split()
    print(" ".join(cmd), flush=True)
    log = out / f"{tag}.log"
    with open(log, "w") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    print(log.read_text().splitlines()[-1])
    sys.exit(p.returncode)


if __name__ == "__main__":
    main()
