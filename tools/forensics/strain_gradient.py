"""Which cell velocity gradient gives Fluent's strain rate S?

Where the SST mu_t limiter is active, S_F = rho a1 k / (SV_MU_T F2) is Fluent's own S
per cell (mut_limiter.py establishes the limiter form). Each candidate gradient operator
is scored by |S/S_F - 1| on all limited cells and on "clean" ones (>= 3 face hops from
any boundary), by wall layer, and by the eddy viscosity it implies everywhere.
The per-cell root t* of S(G_LSQ + t C) = S_F tests the structure of a correction C:
the right structure gives one t for every cell.
Usage: python -I strain_gradient.py [case_dir] [stem] [cache_dir]   (see cff.py)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import cff
import fvops
import mut_limiter as ml

A1 = ml.A1


def quantiles(e):
    return f"median {np.median(e):.2e}  p90 {np.quantile(e, .9):.2e}  p99 {np.quantile(e, .99):.2e}"


def root_t(g_base, C, target, ts=np.linspace(-1.0, 3.0, 401)):
    """Per cell, the t nearest 0.5 with S(g_base + t C) = target (NaN if none in range)."""
    vals = np.stack([fvops.strain_rate(g_base + t * C) - target for t in ts])
    t_star = np.full(target.size, np.nan)
    for it, ic in zip(*np.nonzero(np.diff(np.sign(vals), axis=0) != 0)):
        tc = ts[it] - vals[it, ic] * (ts[1] - ts[0]) / (vals[it + 1, ic] - vals[it, ic])
        if np.isnan(t_star[ic]) or abs(tc - 0.5) < abs(t_star[ic] - 0.5):
            t_star[ic] = tc
    return t_star


def main():
    case_dir, stem, cache = cff.args()
    m = cff.load(case_dir, stem, cache)
    ctx = fvops.Ctx(m)
    n = ctx.n
    g = lambda name: cff.cell_field(m["cells"], name, n)  # noqa: E731
    k, w, mut, rho, mu, y = (g(s) for s in ("SV_K", "SV_O", "SV_MU_T", "SV_DENSITY", "SV_MU_LAM", "SV_WALL_DIST"))
    U = np.stack([g("SV_U"), g("SV_V"), g("SV_W")], 1)
    Ub = np.stack([cff.cell_field(m["faces"], s, ctx.nbf) for s in ("SV_U", "SV_V", "SV_W")], 1)
    assert np.isfinite(Ub).all()

    F2 = ml.f2(k, w, y, mu, rho)
    limited = mut * w / (rho * k) < 1 - 1e-9
    S_F = rho * A1 * k / (mut * F2)
    layer = ml.wall_layers(m, 8)
    clean = limited & (fvops.hops_from_boundary(ctx) >= 3)
    print(f"limited cells {limited.sum()}, clean subset (>= 3 hops from boundaries) {clean.sum()}")

    G2 = fvops.lsq(ctx, 2)
    G2_nowall = fvops.lsq(ctx, 2, wall_rows=False)
    cands = {
        "LSQ 1/d^0": fvops.lsq(ctx, 0),
        "LSQ 1/d^1": fvops.lsq(ctx, 1),
        "LSQ 1/d^2 (Kutta today)": G2,
        "LSQ 1/d^3": fvops.lsq(ctx, 3),
        "Green-Gauss cell (arith. faces)": fvops.green_gauss(ctx, fvops.arithmetic_faces(ctx)),
        "Green-Gauss of 2nd-order faces (LSQ2)": fvops.green_gauss(ctx, fvops.second_order_faces(ctx, G2)),
        "LSQ2 + 1/2 jump fit": (G2 + 0.5 * fvops.jump_correction(ctx, G2)).tocsr(),
        "LSQ2 + 1/2 jump fit, no wall row in base": (G2_nowall + 0.5 * fvops.jump_correction(ctx, G2_nowall)).tocsr(),
    }
    print("\n|S/S_F - 1| per candidate; mu_t columns use every cell")
    grads = {}
    for name, G in cands.items():
        gu = ctx.apply(G, U, Ub)
        grads[name] = gu
        S = fvops.strain_rate(gu)
        e = np.abs(S / S_F - 1)
        mut_pred = rho * A1 * k / np.maximum(A1 * w, S * F2)
        l2 = np.sqrt(np.sum((mut_pred - mut) ** 2) / np.sum(mut**2))
        agree = np.mean((S * F2 > A1 * w) == limited)
        print(f"{name:42s} all: {quantiles(e[limited])} | clean: median {np.median(e[clean]):.2e} p90 {np.quantile(e[clean], .9):.2e}"
              f" | mu_t rel L2 {l2:.2%} on/off agree {agree:.4f}")

    print("\nBy wall layer (median |S/S_F - 1|; layer 1 = wall-adjacent, never limited):")
    for name in ("LSQ 1/d^2 (Kutta today)", "LSQ2 + 1/2 jump fit", "LSQ2 + 1/2 jump fit, no wall row in base"):
        r = fvops.strain_rate(grads[name]) / S_F - 1
        print(f"  {name:42s} " + " ".join(f"L{L}:{np.median(np.abs(r[limited & (layer == L)])):.2e}" for L in (2, 3, 4, 5, 6, 7, 8, 99)))

    print("\nPer-cell root t* of S(G_LSQ2 + t C) = S_F on the clean subset (structure test):")
    g2 = grads["LSQ 1/d^2 (Kutta today)"]
    corrections = {
        "C = jump fit, face rows 1/|r_f|^2": ctx.apply(fvops.jump_correction(ctx, G2, 2), U, Ub),
        "C = jump fit, face rows 1/|r_f|^0": ctx.apply(fvops.jump_correction(ctx, G2, 0), U, Ub),
        "C = jump fit, face rows 1/|r_f|^4": ctx.apply(fvops.jump_correction(ctx, G2, 4), U, Ub),
        "C = GG(2nd-order faces) - LSQ2": grads["Green-Gauss of 2nd-order faces (LSQ2)"] - g2,
        "C = LSQ 1/d^3 - LSQ2": grads["LSQ 1/d^3"] - g2,
    }
    for name, C in corrections.items():
        t = root_t(g2[clean], C[clean], S_F[clean])
        ok = ~np.isnan(t)
        q = np.quantile(t[ok], [.1, .25, .5, .75, .9])
        print(f"  {name:36s} root found {ok.mean():.3f}  t* p10/p25/p50/p75/p90 {np.round(q, 3)}  IQR {q[3] - q[1]:.3f}")


if __name__ == "__main__":
    main()
