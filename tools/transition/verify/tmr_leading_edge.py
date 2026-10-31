"""Leading-edge mesh check for the TMR high-Tu transition case.

The TMR plate grids put 2-8 mm cells at the leading edge (Re_dx 1e4-4e4 at
Re 5e6 per unit length). This writes flatplate.py's T3 layout at Re 5e6 with a
chosen leading-edge spacing. The inlet sits 0.1 m upstream instead of TMR's
0.33333, so it computes the inlet Tu and mu_t/mu that give the same leading-edge
k and omega as bench_tmr's high-tu mode after SST freestream decay. It then
prints the bench_fsae_aero command.

  python3 tools/transition/verify/tmr_leading_edge.py runs/verify/tmrle/plate_le50um.msh --ds-le 5e-5
"""
import argparse
import json
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/transition"))
sys.path.insert(0, str(ROOT / "tools/fsae_aero"))
import flatplate as fp  # noqa: E402
import fluent_msh  # noqa: E402

RE = 5e6
BETA, BETA_STAR = 0.0828, 0.09


def inlet():
    """bench_tmr high-tu: Tu 5%, mu_t/mu 100 at x = -0.33333; matched at x = 0 from x = -X_IN."""
    nu, k0 = 1 / RE, 1.5 * 0.05**2
    w0 = k0 / (nu * 100)
    t = 0.33333
    w_le = w0 / (1 + BETA * w0 * t)
    k_le = k0 * (1 + BETA * w0 * t) ** (-BETA_STAR / BETA)
    w_in = w_le / (1 - BETA * w_le * fp.X_IN)
    k_in = k_le * (1 + BETA * w_in * fp.X_IN) ** (BETA_STAR / BETA)
    return math.sqrt(2 * k_in / 3), k_in / (nu * w_in)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--ds-le", type=float, default=5e-5, help="streamwise spacing at the leading edge (m)")
    args = ap.parse_args()
    y = fp.stretched(120, 0.3 / (RE * math.sqrt(0.003)), fp.HEIGHT)
    xp = fp.stretched(400, args.ds_le, fp.LENGTH)
    xu = -fp.stretched(40, args.ds_le, fp.X_IN)[::-1]
    x = np.concatenate([xu[:-1], xp])
    nxn, nyn = len(x), len(y)
    X, Y = np.meshgrid(x, y, indexing="ij")
    points = np.stack([X.ravel(), Y.ravel(), np.zeros(X.size)], axis=1)
    node = lambda i, j: i * nyn + j  # noqa: E731
    I, J = np.meshgrid(np.arange(nxn - 1), np.arange(nyn - 1), indexing="ij")
    I, J = I.ravel(), J.ravel()
    quads = np.stack([node(I, J), node(I + 1, J), node(I + 1, J + 1), node(I, J + 1)], axis=1)
    i, j, le = np.arange(nxn - 1), np.arange(nyn - 1), 40
    bottom = np.stack([node(i, 0), node(i + 1, 0)], axis=1)
    boundaries = [
        ("plate", "wall", bottom[le:]),
        ("upstream", "symmetry", bottom[:le]),
        ("top", "symmetry", np.stack([node(i, nyn - 1), node(i + 1, nyn - 1)], axis=1)),
        ("inlet", "velocity-inlet", np.stack([node(0, j), node(0, j + 1)], axis=1)),
        ("outlet", "pressure-outlet", np.stack([node(nxn - 1, j), node(nxn - 1, j + 1)], axis=1)),
    ]
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    print(json.dumps(fluent_msh.write(out, points, [(3, quads)], boundaries, comment="TMR-conditions plate")))
    tu, ratio = inlet()
    print(" ".join(["build/bench_fsae_aero", "--mesh", str(out), "--backend", "metal", "--out", str(out.parent / "gamma"),
                    "--re", repr(RE), "--length", "1", "--area", "1", "--walls", "plate", "--tu", repr(tu),
                    "--ratio", repr(ratio), "--iterations", "6000", "--cells", "--transition"]))


if __name__ == "__main__":
    main()
