"""Fully developed channel, SST vs SST + intermittency, from the documented equations alone.

Wall units (nu = u_tau = 1, half-height h = Re_tau, dp/dx = -1/h), cell-centred
finite volumes on a stretched grid, wall at y = 0, symmetry at y = h, double
precision. In fully developed 1D shear S = Omega and dV/dy = 0. Omega is fixed in
the wall cell at 6 nu/(beta1 y^2); k has zero flux (Fluent) or k = 0 (Menter) at
the wall. Prints Cf for both models on three grids and the term ablation that
locates the offset.

  build-py/bin/python tools/transition/verify/channel.py [--re-tau 2000]
"""
import argparse

import numpy as np
from scipy.linalg import solve_banded

import model as m


def stretched_faces(h, first, n):
    lo, hi = 1.0 + 1e-12, 1.5
    for _ in range(200):
        r = 0.5 * (lo + hi)
        lo, hi = (r, hi) if first * (r**n - 1) / (r - 1) < h else (lo, r)
    d = first * r ** np.arange(n)
    return np.concatenate([[0.0], np.cumsum(d)]) * h / d.sum()


def solve(h=2000.0, first=0.2, n=160, gamma_model=False, k_wall="flux0", floor=0.1,
          gamma_fixed=None, iters=60000, tol=1e-11):
    yf = stretched_faces(h, first, n)
    y = 0.5 * (yf[1:] + yf[:-1])
    dy, dc = np.diff(yf), np.diff(y)
    U = np.minimum(y, np.log(1 + 0.41 * y) / 0.41 + 5.0)
    k = 3.3 * np.minimum(1.0, (y / 30) ** 2) + 1e-8
    w = np.maximum(1.0 / (0.3 * m.KAPPA * y), 6 / (m.B1 * y**2))
    g = np.ones(n) if gamma_fixed is None else gamma_fixed.copy()

    def face(q):
        qf = np.empty(n + 1)
        qf[1:-1] = q[:-1] + (q[1:] - q[:-1]) * (yf[1:-1] - y[:-1]) / dc
        qf[0], qf[-1] = q[0], q[-1]
        return qf

    def ddy(q, wall=None):
        qf = face(q)
        if wall is not None:
            qf[0] = wall
        return np.diff(qf) / dy

    def solve_tri(G, wall_coef, sink, src, fix_first=None):
        c = G[1:-1] / dc
        lo, di, up = np.zeros(n), np.zeros(n), np.zeros(n)
        di[:-1] += c; up[:-1] -= c; di[1:] += c; lo[1:] -= c
        di[0] += wall_coef
        di = di + sink * dy
        rhs = src * dy
        if fix_first is not None:
            di[0], up[0], rhs[0] = 1.0, 0.0, fix_first
        ab = np.zeros((3, n)); ab[0, 1:] = up[:-1]; ab[1] = di; ab[2, :-1] = lo[1:]
        return solve_banded((1, 1), ab, rhs)

    prev = 0.0
    for it in range(iters):
        S = np.abs(ddy(U, 0.0))
        dk, dw = ddy(k, 0.0 if k_wall == "dirichlet" else None), ddy(w)
        b = m.sst_blend(k, w, y, 1.0, dk, dw, S, gamma_model)
        G = face(1 + b["nut"])
        U = solve_tri(G, G[0] / y[0], np.zeros(n), np.full(n, 1.0 / h))
        S = np.abs(ddy(U, 0.0))
        Pk = np.minimum(b["nut"] * S * S, 10 * m.BSTAR * k * w)
        prod, dmul = Pk, 1.0
        if gamma_model:
            prod, dmul = m.k_coupling(g, Pk, b["nut"], S, S, y, 1.0, floor)
        Gk = face(1 + b["sk"] * b["nut"])
        k_new = solve_tri(Gk, Gk[0] / y[0] if k_wall == "dirichlet" else 0.0, dmul * m.BSTAR * w, prod)
        k = np.maximum(0.7 * k_new + 0.3 * k, 1e-20)
        cross = 2 * (1 - b["F1"]) * m.SW2 / w * dk * dw
        src = b["alpha"] * S * S + np.maximum(cross, 0) + b["beta"] * w * w
        sink = 2 * b["beta"] * w + np.maximum(-cross, 0) / w
        w_new = solve_tri(face(1 + b["sw"] * b["nut"]), 0.0, sink, src, fix_first=6 / (m.B1 * y[0] ** 2))
        w = np.maximum(0.7 * w_new + 0.3 * w, 1e-20)
        if gamma_model and gamma_fixed is None:
            nut = m.A1 * k / np.maximum(m.A1 * w, S * b["F2"])
            o = m.onset(k, w, y, 1.0, S, np.zeros(n))
            a, bb = m.gamma_rates(S, S, o)
            g_new = solve_tri(face(1 + nut / m.SIGMA_GAMMA), 0.0, (a + m.CE2 * bb) * g, (a + bb) * g)
            g = np.clip(0.7 * g_new + 0.3 * g, 1e-7, 1.0)
        if it % 2000 == 0:
            ub = np.sum(U * dy) / h
            if it and abs(ub - prev) < tol * ub:
                break
            prev = ub
    nut = m.A1 * k / np.maximum(m.A1 * w, S * b["F2"])
    ub = np.sum(U * dy) / h
    return dict(y=y, U=U, k=k, w=w, g=g, nut=nut, Ub=ub, Cf=2 / ub**2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--re-tau", type=float, default=2000)
    h = ap.parse_args().re_tau
    print(f"Fully developed channel, Re_tau {h:g}: Cf of SST + intermittency relative to SST")
    for kw in ("flux0", "dirichlet"):
        for first, n in ((0.4, 120), (0.2, 160), (0.1, 220)):
            s, g = solve(h, first, n, False, kw), solve(h, first, n, True, kw)
            print(f"  k wall {kw:9s} y1+ {first:.1f}: Cf SST {s['Cf']:.6e}, gamma {g['Cf']:.6e}, {100*(g['Cf']/s['Cf']-1):+.2f}%")
    s, g = solve(h), solve(h, gamma_model=True)
    one = solve(h, gamma_model=True, gamma_fixed=np.ones_like(g["g"]))
    prod_only = solve(h, gamma_model=True, gamma_fixed=g["g"], floor=1.0)
    print("Ablation (y1+ 0.2, zero k flux), Cf relative to SST:")
    print(f"  full model                                {100*(g['Cf']/s['Cf']-1):+.2f}%")
    print(f"  gamma = 1 everywhere (only F3 differs)    {100*(one['Cf']/s['Cf']-1):+.2f}%")
    print(f"  model gamma, destruction not scaled       {100*(prod_only['Cf']/s['Cf']-1):+.2f}%")
    print("Profiles (model / SST):")
    print("    y+    gamma   k ratio  nu_t ratio")
    for t in (0.2, 1, 2, 3, 5, 7, 9, 11, 13, 16, 20, 30, 100, 500):
        i = int(np.argmin(abs(g["y"] - t)))
        print(f"  {g['y'][i]:6.2f} {g['g'][i]:7.4f} {g['k'][i]/s['k'][i]:9.3f} {g['nut'][i]/s['nut'][i]:9.4f}")


if __name__ == "__main__":
    main()
