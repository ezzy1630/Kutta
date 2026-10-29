"""Which quantity does Fluent's SST eddy-viscosity limiter use: strain S or vorticity Omega?

mu_t = rho a1 k / max(a1 omega, X F2).  Where Fluent's limiter is active,
rho a1 k / SV_MU_T recovers X*F2 exactly, so each hypothesis X in {S, Omega}
is tested by predicting mu_t from own velocity gradients and comparing with SV_MU_T.
Also reports how often the k-production limiter (C_lim = 10) could be active.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import cff
import grad
import walldist

A1, BSTAR = 0.31, 0.09


def stats(x):
    a = np.abs(x)
    return f"median {np.median(a):.3e}  p95 {np.quantile(a, .95):.3e}  p99 {np.quantile(a, .99):.3e}"


def f2(k, w, y, mu, rho):
    phi2 = np.maximum(2 * np.sqrt(k) / (BSTAR * w * y), 500 * mu / (rho * y * y * w))
    return np.tanh(phi2 * phi2)


def wall_layers(m, depth=5):
    """Cell layer index counted from walls (1 = wall-adjacent); 99 beyond `depth`."""
    n, c0, c1 = m["n_cells"], m["c0"], m["c1"]
    inn = c1 >= 0
    a, b = c0[inn], c1[inn]
    layer = np.full(n, 99)
    layer[c0[walldist.wall_faces(m)]] = 1
    for L in range(2, depth + 1):
        prev = layer == L - 1
        nb = np.zeros(n, bool)
        nb[b[prev[a]]] = True
        nb[a[prev[b]]] = True
        layer[nb & (layer == 99)] = L
    return layer


def main():
    case_dir, stem, cache = cff.args()
    m = cff.load(case_dir, stem, cache)
    n = m["n_cells"]
    C, F = m["cells"], m["faces"]
    g = lambda k: cff.cell_field(C, k, n)
    k, w, mut, rho, mu = g("SV_K"), g("SV_O"), g("SV_MU_T"), g("SV_DENSITY"), g("SV_MU_LAM")
    U = np.stack([g("SV_U"), g("SV_V"), g("SV_W")], axis=1)
    nbf = max(z["hi"] for z in m["zones"] if z["type"] != 2)
    Ub = np.stack([cff.cell_field(F, c, nbf) for c in ("SV_U", "SV_V", "SV_W")], axis=1)
    assert np.isfinite(Ub).all()
    y_saved = g("SV_WALL_DIST")
    wd_cache = os.path.join(cache, stem + "_walldist.npz")
    y_own = np.load(wd_cache)["d_exact"] if os.path.exists(wd_cache) else walldist.wall_distances(m)["d_exact"]

    R = mut * w / (rho * k)
    limited = R < 1 - 1e-9
    free = ~limited
    print(f"cells {n}: mu_t*omega/(rho k): min {R.min():.4f} max {R.max():.17g}")
    print(f"limiter inactive (|R-1| < 1e-9) {free.sum()} ({free.mean():.1%}); max |R-1| there {np.abs(R[free] - 1).max():.1e}; active {limited.sum()}")
    print("H0 no limiter, mu_t = rho k/omega:")
    e0 = rho * k / w / mut - 1
    print("   inactive cells:", stats(e0[free]))
    print("   active cells:  ", stats(e0[limited]))

    layer = wall_layers(m)
    print("limited cells per wall layer (layer: limited/total):",
          {L: f"{int((limited & (layer == L)).sum())}/{int((layer == L).sum())}" for L in (1, 2, 3, 4, 5, 99)})
    XF2_fluent = rho * A1 * k / mut  # exact where limited

    pairs = grad.neighbour_pairs(m, nbf)
    methods = {
        "LSQ unweighted": lambda: grad.lsq_gradient(m, U, Ub, 0, pairs),
        "LSQ 1/d": lambda: grad.lsq_gradient(m, U, Ub, 1, pairs),
        "LSQ 1/d^2": lambda: grad.lsq_gradient(m, U, Ub, 2, pairs),
        "Green-Gauss cell": lambda: grad.green_gauss(m, U, Ub),
    }
    print("\nHypothesis tables: |mu_t,pred/SV_MU_T - 1|")
    for gname, gfun in methods.items():
        S, W, _ = grad.strain_vorticity(gfun())
        for yname, y in (("y=SV_WALL_DIST", y_saved), ("y=own exact", y_own)):
            if yname == "y=own exact" and gname != "LSQ 1/d^2":
                continue
            F2 = f2(k, w, y, mu, rho)
            for xname, X in (("S", S), ("Omega", W)):
                pred = rho * A1 * k / np.maximum(A1 * w, X * F2)
                e = pred / mut - 1
                agree = np.mean((X * F2 > A1 * w) == limited)
                print(f"[{gname:16s} {yname:15s}] X={xname:5s} limiter on/off agreement {agree:.4f}")
                print(f"      active   (n={limited.sum():5d}):", stats(e[limited]))
                print(f"      inactive (n={free.sum():5d}):", stats(e[free]))

    print("\nLSQ distance-weight exponent scan, limited cells, ratio r = (X F2)_own / (X F2)_Fluent:")
    for p in (0, 1, 1.5, 2, 2.5, 3, 4):
        S, W, _ = grad.strain_vorticity(grad.lsq_gradient(m, U, Ub, p, pairs))
        F2 = f2(k, w, y_saved, mu, rho)
        rs, rw = (S * F2 / XF2_fluent)[limited], (W * F2 / XF2_fluent)[limited]
        print(f"   w=1/d^{p:<3}: S: |r-1| median {np.median(abs(rs - 1)):.4f} p95 {np.quantile(abs(rs - 1), .95):.4f} signed median {np.median(rs - 1):+.4f}"
              f" | Omega: |r-1| median {np.median(abs(rw - 1)):.4f} p95 {np.quantile(abs(rw - 1), .95):.4f} signed median {np.median(rw - 1):+.4f}")

    S, W, _ = grad.strain_vorticity(grad.lsq_gradient(m, U, Ub, 2, pairs))
    F2 = f2(k, w, y_saved, mu, rho)
    rS, rW = S * F2 / XF2_fluent, W * F2 / XF2_fluent
    sub = limited & (np.abs(S - W) / np.maximum(S, W) > 0.2)
    print(f"\nDiscriminating subset (LSQ 1/d^2): limiter active and |S-Omega|/max(S,Omega) > 20 %: {sub.sum()} cells")
    for xname, r in (("S", rS), ("Omega", rW)):
        print(f"   X={xname:5s}: r median {np.median(r[sub]):.4f} IQR [{np.quantile(r[sub], .25):.4f}, {np.quantile(r[sub], .75):.4f}] "
              f"p5-p95 [{np.quantile(r[sub], .05):.4f}, {np.quantile(r[sub], .95):.4f}]")
    print("By wall layer (X=S, LSQ 1/d^2):")
    for L in (2, 3, 4, 5, 99):
        s = limited & (layer == L)
        if s.sum():
            print(f"   layer {L:2d} n={s.sum():5d} r median {np.median(rS[s]):.3f} p5-p95 [{np.quantile(rS[s], .05):.3f}, {np.quantile(rS[s], .95):.3f}]")
    print("By F2 (X=S, LSQ 1/d^2):")
    for lo, hi in ((0, .5), (.5, .8), (.8, .95), (.95, .999), (.999, 1.01)):
        s = limited & (F2 >= lo) & (F2 < hi)
        if s.sum():
            print(f"   F2 in [{lo}, {hi}) n={s.sum():5d} r median {np.median(rS[s]):.3f} IQR [{np.quantile(rS[s], .25):.3f}, {np.quantile(rS[s], .75):.3f}]")

    Pk = mut * S**2
    ratio = Pk / (10 * rho * BSTAR * k * w)
    print(f"\nProduction limiter: cells with mu_t S^2 > 10 rho beta* k omega: {(ratio > 1).sum()} of {n}; "
          f"max ratio {ratio.max():.3f}, p99 {np.quantile(ratio, .99):.3f}")
    print(f"viscosity ratio max mu_t/mu = {np.max(mut / mu):.2f} (turb-viscosity/max-limit is 1e5)")


if __name__ == "__main__":
    main()
