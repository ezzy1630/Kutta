"""Energy-equation residual of Fluent's saved state: the cleanest test of the scalar SOU operator.

Constant-property air (cp 1006.43, k 0.0242, Pr_t 0.85), adiabatic walls, no viscous heating or
pressure work, so h = cp (T - 298.15) and the equation is sum_f F h_f = sum_f k_eff,f grad T . A.
No pressure, transpose, Rhie-Chow or wall-function term enters; only reconstruction, limiter and diffusion.
Sections: 1 normalisation (a_P), 2 operator candidates, 3 limiter variants, 4 reconstruction point,
5 diffusion scale, 6 gradient-operator basis fit, 7 where / spatial character, 8 bounded limiter fit.
Usage: python -I energy_residual.py [case_dir] [stem] [cache_dir]   (see cff.py)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5py
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize
from scipy.sparse.linalg import spsolve

import cff
import fvops
import momentum as mo

CP, KLAM, PRT, TREF = 1006.43, 0.0242, 0.85, 298.15


def main():
    case_dir, stem, cache = cff.args()
    st = mo.State(case_dir, stem, cache)
    ctx, m = st.ctx, st.m
    n, I, B, c0, c1 = ctx.n, ctx.I, ctx.B, ctx.c0, ctx.c1
    F = st.F
    T = cff.cell_field(m["cells"], "SV_T", n)
    Tb = cff.cell_field(m["faces"], "SV_T", ctx.nbf)
    H = cff.cell_field(m["cells"], "SV_H", n)
    keff = KLAM + CP * st.mut / PRT
    dtau = mo.pseudo_dt(case_dir, stem)
    with h5py.File(os.path.join(case_dir, stem + ".dat.h5"), "r") as f:
        e0, e1 = f["results/residuals/phase-1/energy/data"][-1, :2]
    print(f"Fluent energy: raw {e0:.4e} W, normalisation {e1:.6g}, scaled {e0 / e1:.3e};  max|SV_H - cp(T - {TREF})| = {np.abs(H - CP * (T - TREF)).max():.1e}")
    w3 = st.zone(3)
    print(f"wall faces: max|T_b - T_c| {np.abs(Tb[w3] - T[c0[w3]]).max():.1e} (adiabatic); inlet T_b {np.unique(np.round(Tb[st.zone(10)], 4))}")

    print("\n1. normalisation = sum_P a_P T_P with absolute T (W):")
    out = mo.outflow_coefficient(st)
    dif = mo.face_sum(st, 0.5 * (keff[c0[I]] + keff[c1[I]]) * st.alpha)
    pt = st.rho * ctx.vol * CP / dtau
    So, Sd, Sp = np.sum(out * CP * T), np.sum(dif * T), np.sum(pt * T)
    print(f"  cp*outflow {So:.6g} + k_eff alpha {Sd:.6g} + rho V cp/dtau {Sp:.6g} = {So + Sd + Sp:.6g}; Fluent {e1:.6g} (ratio {e1 / (So + Sd + Sp):.5f}); "
          f"without pseudo-time ratio {e1 / (So + Sd):.4f}")

    G2 = fvops.lsq(ctx, 2)
    gT = ctx.apply(G2, T, Tb)
    glT, psiT = mo.bj_limit(ctx, T, Tb, gT)
    h, hb = CP * (T - TREF), CP * (Tb - TREF)
    D = mo.scalar_diffusion(st, T, Tb, gT, keff)
    hops = fvops.hops_from_boundary(ctx)
    deep = hops >= 3

    def conv(gT_rec):
        return mo.scalar_convection(st, h, hb, CP * gT_rec)

    def show(lab, R):
        print(f"  {lab:50s} sum|R| {np.abs(R).sum():.3e} W  scaled {np.abs(R).sum() / e1:.2e}  deep {np.abs(R[deep]).sum():.3e}")

    print(f"\n2. operator candidates (Fluent: {e0:.2e} W)")
    show("first-order upwind", conv(0 * gT) + D)
    show("SOU, LSQ2 unlimited", conv(gT) + D)
    R = conv(glT) + D
    show("SOU, LSQ2 + Barth-Jespersen (baseline)", R)
    for pw in (1, 1.5, 2.5, 3):
        show(f"SOU, LSQ 1/d^{pw} + BJ", conv(mo.bj_limit(ctx, T, Tb, ctx.apply(fvops.lsq(ctx, pw), T, Tb))[0]) + D)
    J = ctx.apply(fvops.jump_correction(ctx, G2), T, Tb)
    for t in (0.25, 0.5):
        show(f"SOU, LSQ2 + {t} jump correction + BJ", conv(mo.bj_limit(ctx, T, Tb, gT + t * J)[0]) + D)
    show("baseline, harmonic k_eff", conv(glT) + mo.scalar_diffusion(st, T, Tb, gT, keff, gamma_face="harm"))
    show("baseline, no non-orthogonal correction", conv(glT) + mo.scalar_diffusion(st, T, Tb, gT, keff, nonorth=False))
    s5 = st.zone(5)
    print(f"  global balance: signed sum R {R.sum():+.3e} W against enthalpy throughput {np.sum(F[s5] * h[c0[s5]]):.3f} W")

    print("\n3. limiter variants (energy)")
    show("BJ, boundary values excluded from min/max", conv(mo.bj_limit(ctx, T, Tb, gT, include_boundary=False)[0]) + D)
    show("smooth Venkatakrishnan function, no epsilon", conv(mo.bj_limit(ctx, T, Tb, gT, smooth=True)[0]) + D)
    allT = np.concatenate([T, Tb])
    cellf = np.concatenate([c0[I], c1[I], c0[B]])
    othf = np.concatenate([c1[I], c0[I], n + B])
    facesf = np.concatenate([I, I, B])
    lo, hi = T.copy(), T.copy()
    np.minimum.at(lo, cellf, allT[othf])
    np.maximum.at(hi, cellf, allT[othf])

    def limited(dhi, dlo, rvec=None):
        rv = ctx.fc[facesf] - ctx.cen[cellf] if rvec is None else rvec
        dd = np.einsum("ri,ri->r", rv, gT[cellf])
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(dd > 0, dhi / dd, np.where(dd < 0, dlo / dd, np.inf))
        psi = np.ones(n)
        np.minimum.at(psi, cellf, np.minimum(ratio, 1))
        return gT * np.clip(psi, 0, 1)[:, None]

    Dn = allT[othf] - T[cellf]
    show("two-cell bounds per face (allneighbors? #f reading)", conv(limited(np.maximum(Dn, 0), np.minimum(Dn, 0))) + D)
    show("BJ, cell-to-cell distance (r = d)", conv(limited(hi[cellf] - T[cellf], lo[cellf] - T[cellf],
                                                              np.concatenate([ctx.d, -ctx.d, ctx.fc[B] - ctx.cen[c0[B]]]))) + D)
    for e in (0.01, 0.05):
        rng = (hi - lo)[cellf]
        show(f"BJ, bounds widened by {e} x local range", conv(limited(hi[cellf] - T[cellf] + e * rng, lo[cellf] - T[cellf] - e * rng)) + D)

    print("\n4. reconstruction point (BJ-limited LSQ2 gradient)")
    n_ = ctx.nhat[I]
    x0 = ctx.cen[c0[I]]
    tt = np.einsum("ij,ij->i", ctx.fc[I] - x0, n_) / np.einsum("ij,ij->i", ctx.d, n_)
    xip = x0 + tt[:, None] * ctx.d

    def conv_r(r0, r1):
        p0 = h[c0[I]] + CP * np.einsum("ij,ij->i", glT[c0[I]], r0)
        p1 = h[c1[I]] + CP * np.einsum("ij,ij->i", glT[c1[I]], r1)
        pf = np.where(F[I] >= 0, p0, p1)
        Rr = np.zeros(n)
        np.add.at(Rr, c0[I], F[I] * pf)
        np.add.at(Rr, c1[I], -F[I] * pf)
        np.add.at(Rr, c0[B], F[B] * np.where(F[B] >= 0, h[c0[B]], hb))
        return Rr

    for lab, (r0, r1) in (("face centroid x_f", (ctx.r0, ctx.r1)), ("centroid line meets face plane", (xip - x0, xip - ctx.cen[c1[I]])),
                          ("centroid midpoint", (0.5 * ctx.d, -0.5 * ctx.d))):
        show(lab, conv_r(r0, r1) + D)

    print("\n5. diffusion scale (R = conv + c * diffusion)")
    C = conv(glT)
    print("  " + " ".join(f"c={c}: {np.abs(C + c * D).sum():.3e}" for c in (0, 0.8, 0.9, 1.0, 1.1, 1.25)))

    print("\n6. gradient-operator basis fit on cells whose own and upwind reconstructions are all unlimited (linear there)")
    limc = psiT < 1 - 1e-12
    up = np.where(F[I] >= 0, c0[I], c1[I])
    inv = limc.copy()
    np.logical_or.at(inv, c0[I], limc[up])
    np.logical_or.at(inv, c1[I], limc[up])
    mk = ~inv
    AF = fvops.arithmetic_faces(ctx)
    ops = {"jump correction J": fvops.jump_correction(ctx, G2), "LSQ1 - LSQ2": fvops.lsq(ctx, 1) - G2, "LSQ3 - LSQ2": fvops.lsq(ctx, 3) - G2,
           "GG(2nd-order faces) - LSQ2": fvops.green_gauss(ctx, fvops.second_order_faces(ctx, G2)) - G2,
           "GG cell - LSQ2": fvops.green_gauss(ctx, AF) - G2}
    Tk = {k: conv(ctx.apply(Dop, T, Tb)) for k, Dop in ops.items()}
    y = -R[mk]
    for k, v in Tk.items():
        c = np.dot(v[mk], y) / np.dot(v[mk], v[mk])
        print(f"  single {k:28s} coefficient {c:+.3f}; |R| {np.abs(y).sum():.3e} -> {np.abs(y - c * v[mk]).sum():.3e}")
    X = np.stack([v[mk] for v in Tk.values()], 1)
    cj, *_ = np.linalg.lstsq(X, y, rcond=None)
    print(f"  joint: {', '.join(f'{k} {v:+.3f}' for k, v in zip(Tk, cj))}; |R| -> {np.abs(y - X @ cj).sum():.3e}")
    Tsou = C - conv(0 * glT)
    ok = mk & (np.abs(Tsou) > 1e-3 * np.abs(Tsou).max())
    print(f"  per-cell R / SOU correction on these cells: p10/25/50/75/90 {np.round(np.quantile(R[ok] / Tsou[ok], [.1, .25, .5, .75, .9]), 4)}")

    print("\n7. where, and spatial character")
    gm = np.linalg.norm(gT, axis=1)
    for a_, b_ in ((0, 1), (1, 10), (10, 100), (100, 1e9)):
        s = (gm >= a_) & (gm < b_)
        print(f"  |grad T| in [{a_:g}, {b_:g}) K/m: cells {s.sum():5d} sum|R| {np.abs(R[s]).sum():.3e}  limited {(psiT[s] < 1 - 1e-12).mean():.2f}")
    print(f"  unlimited-stencil cells: {mk.sum()} carry {np.abs(R[mk]).sum():.3e}; neighbour correlation corr(R_c0, R_c1) {np.corrcoef(R[c0[I]], R[c1[I]])[0, 1]:+.3f}")
    Fi = F[I]
    Dk = 0.5 * (keff[c0[I]] + keff[c1[I]]) * st.alpha
    an0, an1 = Dk + np.maximum(-Fi, 0) * CP, Dk + np.maximum(Fi, 0) * CP
    ap = mo.face_sum(st, 0 * Dk)
    np.add.at(ap, c0[I], an0)
    np.add.at(ap, c1[I], an1)
    np.add.at(ap, c0[B], np.abs(F[B]) * CP)
    A1 = sp.csr_matrix((-np.concatenate([an0, an1]), (np.concatenate([c0[I], c1[I]]), np.concatenate([c1[I], c0[I]]))), shape=(n, n))
    e = spsolve((A1 + sp.diags(ap + pt)).tocsc(), R)
    print(f"  implied state error (first-order + pseudo-time operator)^-1 R: max {np.abs(e).max():.2f} K, median {np.median(np.abs(e)):.1e} K, "
          f"neighbour correlation {np.corrcoef(e[c0[I]], e[c1[I]])[0, 1]:+.2f}")

    print("\n8. bounded fit of a per-cell limiter factor psi in [0, 1] on the cells BJ limits (does any limiter explain it?)")
    val = np.abs(Fi) * CP * np.einsum("ij,ij->i", gT[up], np.where((Fi >= 0)[:, None], ctx.r0, ctx.r1))
    dnw = np.where(Fi >= 0, c1[I], c0[I])
    M = sp.csr_matrix((np.concatenate([val, -val]), (np.concatenate([up, dnw]), np.concatenate([up, up]))), shape=(n, n)).tocsc()
    base = conv(0 * gT) + D
    lim_idx = np.nonzero(limc & (np.abs(M).sum(0).A1 > 0))[0]
    fixed = base + M @ np.where(limc, 0.0, 1.0)
    Ml = M[:, lim_idx]
    scale = 1.0 / np.abs(fixed).sum()

    def fg(x):
        r = fixed + Ml @ x
        return scale * r @ r * 1e3, scale * 2e3 * (Ml.T @ r)

    sol = minimize(fg, psiT[lim_idx], jac=True, method="L-BFGS-B", bounds=[(0, 1)] * lim_idx.size, options={"maxiter": 2000})
    psi_fit = np.where(limc, 0.0, 1.0)
    psi_fit[lim_idx] = sol.x
    Rfit = base + M @ psi_fit
    print(f"  unknowns {lim_idx.size}; sum|R| {np.abs(R).sum():.3e} -> {np.abs(Rfit).sum():.3e}; corr(psi_fit, psi_BJ) {np.corrcoef(sol.x, psiT[lim_idx])[0, 1]:.3f}; "
          f"psi_fit p10/50/90 {np.round(np.quantile(sol.x, [.1, .5, .9]), 3)} vs BJ {np.round(np.quantile(psiT[lim_idx], [.1, .5, .9]), 3)}")


if __name__ == "__main__":
    main()
