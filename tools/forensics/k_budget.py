"""Kutta's k equation at Fluent's saved V6 state, cell by cell, in one wall zone's pipe (default wall-inlet).

The question is which k-equation term disagrees with Fluent where Kutta's k is high. Fluent's own k residual
at its saved state is ~0, so evaluating Kutta's terms there shows a disagreement directly:
  R = C + D - P + Y,  C convection (saved SV_FLUX, SOU, Barth-Jespersen), D diffusion (net outflow),
  P = min(mu_t S^2, 10 rho beta* k omega) with the TG 4-419..422 blend in wall cells, Y = rho beta* k omega.
R < 0 means Kutta's equation has a net source that Fluent's does not. Every sum is normalised by Y.

Kutta's strain S at Fluent's saved velocity comes from the V6 slow test with KUTTA_V6_WALL_DIR set
(<backend>-saved-sst.csv). In mu_t-limited cells Fluent's own S is exact: S_F = rho a1 k / (mu_t F2).
Elsewhere S is not saved, so section 3 infers it from the balance, S_need = sqrt((C + D + Y) / (mu_t V)),
and checks that inference on the limited cells first.

Sections: 1 S against S_F; 2 R by wall layer and pipe band; 3 implied S; 4 standard k-omega's f_beta*
(TG 4.4.2) as a test of TG 4.4.3's "f_beta* = 1 for SST"; 5 optional Kutta k / Fluent k (volume-weighted) from solved runs,
also by distance from the symmetry plane.
Usage: python -I k_budget.py <case_dir> <stem> <cache_dir> <wall_dir> <backend> [zone] [run_dir/run ...]
  run_dir/run names a <run>-cells.csv written by the V6 test (Fluent cell order).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import fvops
import momentum as mo
import mut_limiter as ml

BS, A1, KAPPA = 0.09, 0.31, 0.4187
BANDS, LAYERS = 5, (1, 2, 3, 4, 5, 6)


def circle(p):
    """Least-squares circle through 2D points (Kasa): centre, radius."""
    A = np.column_stack([2 * p, np.ones(len(p))])
    sol = np.linalg.lstsq(A, (p**2).sum(1), rcond=None)[0]
    return sol[:2], np.sqrt(sol[2] + sol[:2] @ sol[:2])


def ratio_sum(R, Y, s):
    return R[s].sum() / Y[s].sum() if s.any() else np.nan


def main():
    case_dir, stem, cache, wall_dir, backend = sys.argv[1:6]
    zone = sys.argv[6] if len(sys.argv) > 6 else "wall-inlet"
    runs = sys.argv[7:]
    st = mo.State(case_dir, stem, cache)
    ctx = st.ctx
    rho, mu, mut, k, w, y = st.rho, st.mu, st.mut, st.k, st.w, st.y
    kb, U, Ub = st.kb, st.U, st.Ub
    K = np.genfromtxt(os.path.join(wall_dir, backend + "-saved-sst.csv"), delimiter=",", names=True)
    SK = K["strain"]
    layer = ml.wall_layers(st.m, 8)

    # The zone's pipe: axis along its longest extent, circle fitted to the wall-face centroids.
    z = next(z for z in st.m["zones"] if z["name"] == zone)
    fc = ctx.fc[z["lo"]:z["hi"]]
    ax = int(np.argmax(fc.max(0) - fc.min(0)))
    other = [i for i in range(3) if i != ax]
    centre, radius = circle(fc[:, other])
    x = ctx.cen[:, ax]
    pipe = (np.linalg.norm(ctx.cen[:, other] - centre, axis=1) < radius) & (x >= fc[:, ax].min()) & (x <= fc[:, ax].max())
    edges = np.quantile(fc[:, ax], np.linspace(0, 1, BANDS + 1))
    band = np.clip(np.digitize(x, edges) - 1, 0, BANDS - 1)
    print(f"{zone}: axis {'xyz'[ax]} [{edges[0]:+.4f}, {edges[-1]:+.4f}], radius {radius:.5f}; {pipe.sum()} cells; "
          f"bands {' '.join(f'{e:+.4f}' for e in edges)}")

    F2 = ml.f2(k, w, y, mu, rho)
    limited = mut * w / (rho * k) < 1 - 1e-9
    SF = rho * A1 * k / (mut * F2)
    print("\n1. Kutta's S / Fluent's S_F in mu_t-limited cells (median; S_F is exact there)")
    print("  whole domain by layer: " + " ".join(f"L{L} {np.median((SK / SF)[limited & (layer == L)]):.3f}" for L in (2, 3, 4, 5)) +
          f" deep {np.median((SK / SF)[limited & (layer == 99)]):.3f}")
    print(f"  {zone} pipe, layers 2-6: {np.median((SK / SF)[limited & pipe & (layer >= 2)]):.3f} (n={(limited & pipe & (layer >= 2)).sum()})")

    # k equation at Fluent's state with Kutta's operators: LSQ 1/d^2 without wall rows (zero-flux k), BJ-limited SOU.
    G = fvops.lsq(ctx, 2, wall_rows=False)
    gk, gw = ctx.apply(G, k, kb), ctx.apply(fvops.lsq(ctx, 2), w, st.wb)
    cd = np.maximum(2 * rho * np.einsum("ij,ij->i", gk, gw) / (1.168 * w), 1e-10)
    arg1 = np.minimum(np.maximum(np.sqrt(k) / (BS * w * y), 500 * mu / (rho * y * y * w)), 4 * rho * k / (1.168 * cd * y * y))
    F1 = np.tanh(arg1**4)
    Gk = mu + mut * (F1 / 1.176 + (1 - F1))
    C = mo.scalar_convection(st, k, kb, mo.bj_limit(ctx, k, kb, gk)[0])
    D = mo.scalar_diffusion(st, k, kb, gk, Gk)
    Y = rho * BS * k * w * ctx.vol
    walls = st.zone(3)
    cw, nh = ctx.c0[walls], ctx.nhat[walls]
    ut = np.linalg.norm(U[cw] - np.einsum("ij,ij->i", U[cw], nh)[:, None] * nh, axis=1)
    tau = np.linalg.norm(st.wall_shear[walls], axis=1) / ctx.area[walls]
    dn = np.abs(st.dn[walls])
    yp = rho[cw] * np.sqrt((mu[cw] * ut / dn + 0.3 * rho[cw] * k[cw]) / rho[cw]) * dn / mu[cw]
    plam, plog = mut[cw] * (tau / mu[cw]) ** 2, tau**2 / (mu[cw] * KAPPA * yp)
    wall_pk = plam * plog / (plam + plog)  # TG 4-422; one wall face per wall cell here

    def production(mu_t, S):
        P = np.minimum(mu_t * S**2, 10 * rho * BS * k * w)
        P[cw] = wall_pk
        return P * ctx.vol

    mut_kutta = rho * k / np.maximum(w, SK * F2 / A1)
    P_kutta = production(mut_kutta, SK)
    P_exact = production(np.where(limited, mut, mut_kutta), np.where(limited, SF, SK))
    print(f"  production limiter: max P/(10 rho beta* k omega) in the pipe {(mut * SK**2 / (10 * rho * BS * k * w))[pipe].max():.3f}")

    print("\n2. R / Y at Fluent's state (R < 0: Kutta's k equation has a net source Fluent's does not), layer x band")
    for lab, P in (("Kutta's S and mu_t", P_kutta), ("S_F and Fluent's mu_t where limited", P_exact)):
        R = C + D - P + Y
        print(f"  {lab}:")
        for L in LAYERS:
            print(f"    L{L}: " + " ".join(f"{ratio_sum(R, Y, pipe & (layer == L) & (band == b)):+.3f}" for b in range(BANDS)))
        print(f"    bands 2-5, all layers: {ratio_sum(R, Y, pipe & (band >= 1)):+.4f}; wall cells of every zone: {ratio_sum(R, Y, np.isin(np.arange(ctx.n), cw)):+.4f}")
    s = pipe & (band >= 1)
    print(f"  bands 2-5 totals / Y: P {P_kutta[s].sum() / Y[s].sum():.3f}, C {C[s].sum() / Y[s].sum():+.3f}, D {D[s].sum() / Y[s].sum():+.3f}")

    print("\n3. implied S: S_need = sqrt((C + D + Y) / (mu_t V)), medians")
    S_need = np.sqrt(np.maximum((C + D + Y) / ctx.vol, 0) / mut)
    ok = S_need > 0
    print("  control, limited cells (S_need / S_F): " + " ".join(f"L{L} {np.median((S_need / SF)[ok & limited & pipe & (layer == L)]):.3f}" for L in (2, 3, 4, 5)) + " (pipe); "
          f"layers 2-4 elsewhere {np.median((S_need / SF)[ok & limited & ~pipe & (layer >= 2) & (layer <= 4)]):.3f}")
    production_matters = SK**2 / (BS * w**2) >= 0.3
    for L in (2, 3):
        print(f"  unlimited L{L} (S_need / S_Kutta) by band: " + " ".join(
            f"{np.median((S_need / SK)[q]):.3f} (n={q.sum()})" if q.sum() > 3 else "-"
            for q in (ok & ~limited & pipe & (layer == L) & (band == b) for b in range(BANDS))))
    q = ok & ~limited & (layer >= 2) & (layer <= 3) & production_matters
    print(f"  unlimited L2-3 with P/Y >= 0.3: pipe {np.median((S_need / SK)[q & pipe]):.3f} (n={(q & pipe).sum()}), elsewhere {np.median((S_need / SK)[q & ~pipe]):.3f} (n={(q & ~pipe).sum()})")

    print("\n4. standard k-omega f_beta* = (1 + 680 chi^2)/(1 + 400 chi^2), chi = grad k . grad omega / omega^3 > 0, on the limited control")
    chi = np.einsum("ij,ij->i", gk, gw) / w**3
    f = np.where(chi > 0, (1 + 680 * chi**2) / (1 + 400 * chi**2), 1.0)
    S_need_f = np.sqrt(np.maximum((C + D + Y * f) / ctx.vol, 0) / mut)
    for region, sel in (("layers 2-4", (layer >= 2) & (layer <= 4)), ("layers >= 5", (layer >= 5) | (layer == 99))):
        for lo, hi in ((1, 1.0001), (1.01, 1.05), (1.05, 1.15), (1.15, 1.3), (1.3, 1.8)):
            q = limited & sel & (f >= lo) & (f < hi)
            if q.sum() >= 5:
                print(f"  {region:11s} f_beta* [{lo}, {hi}) n={q.sum():5d}: S_need/S_F without {np.median((S_need / SF)[q]):.3f}, "
                      f"with {np.median((S_need_f / SF)[q]):.3f} (a real f_beta* would give 1/sqrt(f) = {np.median(1 / np.sqrt(f[q])):.3f} without)")

    symmetry_hop = np.full(ctx.n, 99)
    symmetry_hop[ctx.c0[st.zone(7)]] = 1
    I, c0, c1 = ctx.I, ctx.c0, ctx.c1
    for h in (2, 3):
        front = symmetry_hop == h - 1
        reached = np.zeros(ctx.n, bool)
        reached[c1[I][front[c0[I]]]] = True
        reached[c0[I][front[c1[I]]]] = True
        symmetry_hop[reached & (symmetry_hop == 99)] = h
    for run in runs:
        c = np.genfromtxt(run + "-cells.csv", delimiter=",", names=True)
        kv, fv = c["k"] * ctx.vol, k * ctx.vol
        print(f"\n5. {run}: Kutta k / Fluent k, volume-weighted by layer x band")
        for L in LAYERS:
            print(f"    L{L}: " + " ".join(f"{ratio_sum(kv, fv, pipe & (layer == L) & (band == b)):.3f}" for b in range(BANDS)))
        wv, Fw = c["omega"] * ctx.vol, w * ctx.vol
        print("  bands 3-5 by cell hops from the symmetry plane (1 = has a symmetry face), k and omega ratios")
        for L in (1, 2, 3, 4):
            print(f"    L{L}: " + " | ".join(f"hop {'>3' if h == 99 else h} n={q.sum():3d} k {ratio_sum(kv, fv, q):.3f} omega {ratio_sum(wv, Fw, q):.3f}"
                                         for h in (1, 2, 3, 99) for q in [pipe & (band >= 2) & (layer == L) & (symmetry_hop == h)]))


if __name__ == "__main__":
    main()
