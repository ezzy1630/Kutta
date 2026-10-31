"""Double-precision evaluation of the documented gamma terms on one column of a Kutta TMR plate run.

Reads bench_tmr's cells.csv (u, v, k, omega, mu_t/mu, intermittency) for the
2DZP plate (Re 5e6 per unit length), takes S = |du/dy| and dV/dy = dv/dy from
second-order differences along the structured column, prints the TG 4.8.2 terms
at a few cells, and solves the wall-normal gamma equation on Kutta's frozen k,
omega, nu_t and S (zero flux at the wall, Kutta's gamma at the column top).
Streamwise convection is neglected in that solve; boundary_layer.py includes it.

  build-py/bin/python tools/transition/verify/column.py RUN_DIR [--x 0.95] [--sst SST_RUN_DIR] [--json OUT]
"""
import argparse
import json

import numpy as np
from scipy.linalg import solve_banded

import model as m

NU = 1 / 5e6


def column(run, x0):
    d = np.genfromtxt(f"{run}/cells.csv", delimiter=",", names=True)
    xs = np.unique(np.round(d["x"], 9))
    xx = xs[np.argmin(abs(xs - x0))]
    c = d[abs(d["x"] - xx) < 1e-8]
    return xx, c[np.argsort(c["y"])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--x", type=float, default=0.95)
    ap.add_argument("--sst", help="matching fully turbulent run, for k and nu_t ratios")
    ap.add_argument("--json")
    args = ap.parse_args()
    xx, c = column(args.run, args.x)
    y, u, v, k, w, g = c["y"], c["u"], c["v"], c["k"], c["omega"], c["intermittency"]
    nut = c["mut_over_mu"] * NU
    n = len(y)

    def ddy(q):
        out = np.gradient(q, y, edge_order=2)
        out[0] = q[0] / y[0]  # one-sided to the wall value 0 (u and v vanish there)
        return out

    S, dvdy = np.abs(ddy(u)), ddy(v)
    utau = np.sqrt(NU * u[0] / y[0])
    yp = y * utau / NU
    o = m.onset(k, w, y, NU, S, dvdy)
    a, b = m.gamma_rates(S, S, o)
    # Wall-normal equilibrium of the gamma equation on Kutta's frozen fields.
    yf = np.concatenate([[0.0], 0.5 * (y[1:] + y[:-1]), [1.5 * y[-1] - 0.5 * y[-2]]])
    vol = np.diff(yf)
    coef = 0.5 * ((NU + nut / m.SIGMA_GAMMA)[1:] + (NU + nut / m.SIGMA_GAMMA)[:-1]) / np.diff(y)
    gs = np.ones(n)
    for _ in range(2000):
        lo, di, up = np.zeros(n), np.zeros(n), np.zeros(n)
        di[1:] += coef; lo[1:] -= coef; di[:-1] += coef; up[:-1] -= coef
        di += (a + m.CE2 * b) * gs * vol
        rhs = (a + b) * gs * vol
        di[-1], lo[-1], rhs[-1] = 1.0, 0.0, g[-1]
        ab = np.zeros((3, n)); ab[0, 1:] = up[:-1]; ab[1] = di; ab[2, :-1] = lo[1:]
        new = np.clip(solve_banded((1, 1), ab, rhs), 1e-7, 1)
        done = np.max(abs(new - gs)) < 1e-13
        gs = 0.5 * new + 0.5 * gs
        if done:
            break
    ratios = None
    if args.sst:
        _, s = column(args.sst, args.x)
        ratios = dict(k=(k / s["k"]).tolist(), nut=(c["mut_over_mu"] / s["mut_over_mu"]).tolist(), u=(u / s["u"]).tolist())
    print(f"column x = {xx:.4f}, u_tau = {utau:.5f} (from the wall cell)")
    print("    y+     R_T  F_turb   Tu_L lambda Re_th_c    Re_V F_onset  P/(rho S)  E/(rho S)  gamma Kutta  gamma eq.")
    for t in (0.2, 1, 3, 5, 8, 10, 13, 16, 20, 30, 100):
        i = int(np.argmin(abs(yp - t)))
        print(f"{yp[i]:7.2f} {o['RT'][i]:7.3f} {o['Fturb'][i]:7.4f} {o['TuL'][i]:6.2f} {o['lam'][i]:6.3f} {o['Rec'][i]:7.1f} "
              f"{o['ReV'][i]:7.1f} {o['Fon'][i]:7.4f} {a[i]/S[i]*g[i]*(1-g[i]):10.3e} {b[i]/S[i]*g[i]*(m.CE2*g[i]-1):10.3e} "
              f"{g[i]:12.4f} {gs[i]:9.4f}")
    near = yp < 40
    print(f"max |gamma Kutta - gamma equation| for y+ < 40: {np.max(abs(g[near] - gs[near])):.4f}")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(dict(x=float(xx), utau=float(utau), yplus=yp.tolist(), gamma_kutta=g.tolist(),
                           gamma_equation=gs.tolist(), ratios=ratios), f)


if __name__ == "__main__":
    main()
