"""k and omega residuals of Fluent's saved state with the SST sources (docs/research/fluent-momentum-forensics.md section 6).

SST as in src/kutta/kernels/sst.kern.h and fluent-sst-forensics.md: sigma as divisors, F1 from CD_kw,
production min(mu_t S^2, 10 rho beta* k omega), G_omega = alpha rho P_k / mu_t, cross diffusion
2 (1 - F1) rho grad k . grad omega / (1.168 omega); omega fixed in wall-adjacent cells (excluded).
Sections: 1 normalisations (a_P incl. linearised sinks), 2 wall-cell k production, 3 strain-rate gradient,
4 Fluent's own S in limiter-active cells, 5 omega near-wall ingredients.
Usage: python -I turbulence_residual.py [case_dir] [stem] [cache_dir]   (see cff.py)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import cff
import fvops
import momentum as mo
import mut_limiter as ml

BS, A1, KAPPA = 0.09, 0.31, 0.4187


def main():
    case_dir, stem, cache = cff.args()
    st = mo.State(case_dir, stem, cache)
    ctx = st.ctx
    n, I, c0 = ctx.n, ctx.I, ctx.c0
    rho, mu, mut, k, w, y = st.rho, st.mu, st.mut, st.k, st.w, st.y
    U, Ub, kb, wb = st.U, st.Ub, st.kb, st.wb
    dtau = mo.pseudo_dt(case_dir, stem)
    layer = ml.wall_layers(st.m, 8)
    wallcell = layer == 1
    res = st.res
    G2 = fvops.lsq(ctx, 2)
    gk, gw = ctx.apply(G2, k, kb), ctx.apply(G2, w, wb)
    g2 = ctx.apply(G2, U, Ub)
    G2nw = fvops.lsq(ctx, 2, wall_rows=False)
    gS = g2 + 0.5 * ctx.apply(fvops.jump_correction(ctx, G2), U, Ub)
    gSnw = ctx.apply(G2nw, U, Ub) + 0.5 * ctx.apply(fvops.jump_correction(ctx, G2nw), U, Ub)

    def F1of(gk_, gw_):
        cd = np.maximum(2 * rho * np.einsum("ij,ij->i", gk_, gw_) / (1.168 * w), 1e-10)
        a = np.sqrt(k) / (BS * w * y)
        b = 500 * mu / (rho * y * y * w)
        c = 4 * rho * k / (1.168 * cd * y * y)
        return np.tanh(np.minimum(np.maximum(a, b), c) ** 4)

    F1 = F1of(gk, gw)
    Gk = mu + mut * (F1 / 1.176 + (1 - F1))
    Gw = mu + mut * (F1 / 2.0 + (1 - F1) / 1.168)
    beta = 0.075 * F1 + 0.0828 * (1 - F1)
    alpha = F1 * (0.075 / 0.09 - 0.41**2 / (2.0 * 0.3)) + (1 - F1) * (0.0828 / 0.09 - 0.41**2 / (1.168 * 0.3))

    print("1. normalisations (col1 = sum a_P |phi_P|)")
    out = mo.outflow_coefficient(st)
    pt = rho * ctx.vol / dtau
    difk = mo.face_sum(st, 0.5 * (Gk[c0[I]] + Gk[ctx.c1[I]]) * st.alpha)
    difw = mo.face_sum(st, 0.5 * (Gw[c0[I]] + Gw[ctx.c1[I]]) * st.alpha)
    for lab, ap in (("out + dif + pt", out + difk + pt), ("out + dif + pt + rho beta* omega V", out + difk + pt + rho * BS * w * ctx.vol),
                    ("out + dif + rho beta* omega V (no pt)", out + difk + rho * BS * w * ctx.vol)):
        print(f"  k:     {lab:46s} ratio col1/sum {res['k'][1] / np.sum(ap * k):.4f}")
    nw = ~wallcell
    for lab, ap in (("out + dif + pt", out + difw + pt), ("out + dif + pt + rho beta omega V", out + difw + pt + rho * beta * w * ctx.vol),
                    ("out + dif + pt + 2 rho beta omega V (Newton)", out + difw + pt + 2 * rho * beta * w * ctx.vol)):
        print(f"  omega: {lab:46s} ratio col1/sum, wall cells excluded {res['omega'][1] / np.sum((ap * w)[nw]):.4f}")

    w3 = st.zone(3)
    cw, nh, area = c0[w3], ctx.nhat[w3], ctx.area[w3]
    Ut = U[cw] - np.einsum("ij,ij->i", U[cw], nh)[:, None] * nh
    ut = np.linalg.norm(Ut, axis=1)
    tau = np.linalg.norm(st.wall_shear[w3], axis=1) / area  # SV_WALL_SHEAR is the face force
    dn = np.abs(st.dn[w3])
    us = np.sqrt((mu[cw] * ut / dn + 0.3 * rho[cw] * k[cw]) / rho[cw])
    yp = rho[cw] * us * dn / mu[cw]
    utau = np.sqrt(tau / rho[cw])
    plam = mut[cw] * (tau / mu[cw]) ** 2
    plog = tau**2 / (mu[cw] * KAPPA * yp)

    def harmonic(a, b):
        return a * b / (a + b)

    def k_res(S, wall_pk=None, gk_=gk):
        glk = mo.bj_limit(ctx, k, kb, gk_)[0]
        Pk = np.minimum(mut * S**2, 10 * rho * BS * k * w)
        if wall_pk is not None:
            Pk = Pk.copy()
            Pk[cw] = wall_pk  # one wall face per wall cell in this mesh
        return mo.scalar_convection(st, k, kb, glk) + mo.scalar_diffusion(st, k, kb, gk_, Gk) - (Pk - rho * BS * k * w) * ctx.vol

    def w_res(S, gw_=gw, wb_=wb, cdk=1.0, gam="arith", gform="nut"):
        glw = mo.bj_limit(ctx, w, wb_, gw_)[0]
        Pk = np.minimum(mut * S**2, 10 * rho * BS * k * w)
        cd = cdk * 2 * (1 - F1) * rho * np.einsum("ij,ij->i", gk, gw_) / (1.168 * w)
        R = mo.scalar_convection(st, w, wb_, glw) + mo.scalar_diffusion(st, w, wb_, gw_, Gw, gamma_face=gam) \
            - ((alpha * rho * Pk / mut if gform == "nut" else alpha * w / k * Pk) - rho * beta * w * w + cd) * ctx.vol
        R[wallcell] = 0
        return R

    S2 = fvops.strain_rate(g2)
    print(f"\n2. wall-cell k production (k sum|R| over the {wallcell.sum()} wall cells; Fluent's whole-domain k residual {res['k'][0]:.2e})")
    variants = {"mu_t S^2 from LSQ2 (no wall treatment)": None, "0": 0 * tau,
                "tau_w tau_w/mu": tau * tau / mu[cw], "tau_w u*/(kappa dy)": tau * us / (KAPPA * dn), "tau_w u_tau/(kappa dy)": tau * utau / (KAPPA * dn),
                "mu_t (tau_w/mu)^2": plam, "tau_w^2/(mu kappa y+)": plog,
                "harmonic(mu_t (tau_w/mu)^2, tau_w^2/(mu kappa y+))  [Kutta sst_wall]": harmonic(plam, plog),
                "  same, kappa 0.41": harmonic(plam, tau**2 / (mu[cw] * 0.41 * yp)),
                "  same, y+ from u_tau": harmonic(plam, tau**2 / (mu[cw] * KAPPA * rho[cw] * utau * dn / mu[cw])),
                "  power-2 blend": (plam**-2 + plog**-2) ** -0.5, "  power-4 blend": (plam**-4 + plog**-4) ** -0.25,
                "  min": np.minimum(plam, plog)}
    for lab, pv in variants.items():
        print(f"  {lab:66s} {np.abs(k_res(S2, pv)[wallcell]).sum():.3e}")

    pk_wall = harmonic(plam, plog)
    print("\n3. strain-rate gradient in P_k = mu_t S^2 (wall cells use the harmonic wall production)")
    for lab, gU in (("LSQ2", g2), ("LSQ2 + 1/2 J (WFGC candidate)", gS), ("LSQ2 without wall rows + 1/2 J", gSnw)):
        S = fvops.strain_rate(gU)
        Rk, Rw = k_res(S, pk_wall), w_res(S)
        print(f"  {lab:34s} k sum|R| {np.abs(Rk).sum():.3e} (scaled {np.abs(Rk).sum() / res['k'][1]:.2e}; Fluent {res['k'][0] / res['k'][1]:.2e}) | "
              f"omega sum|R| {np.abs(Rw).sum():.3e} (scaled {np.abs(Rw).sum() / res['omega'][1]:.2e}; Fluent {res['omega'][0] / res['omega'][1]:.2e})")
    Rk, Rw = k_res(S2, pk_wall), w_res(S2)
    print("  by wall layer (LSQ2):", " ".join(f"L{L}: k {np.abs(Rk[layer == L]).sum():.1e} omega {np.abs(Rw[layer == L]).sum():.1e}" for L in (1, 2, 3, 4, 99)))
    glk, glw = mo.bj_limit(ctx, k, kb, gk)[0], mo.bj_limit(ctx, w, wb, gw)[0]
    souk = mo.scalar_convection(st, k, kb, glk) - mo.scalar_convection(st, k, kb, 0 * glk)
    souw = mo.scalar_convection(st, w, wb, glw) - mo.scalar_convection(st, w, wb, 0 * glw)
    souw[wallcell] = 0
    print(f"  residual / SOU correction: k {np.abs(Rk).sum() / np.abs(souk).sum():.2f}, omega {np.abs(Rw).sum() / np.abs(souw).sum():.2f} "
          "(momentum and energy: ~0.07)")

    print("\n4. Fluent's own S in limiter-active cells: S_F = rho a1 k / (mu_t F2)")
    F2 = ml.f2(k, w, y, mu, rho)
    limited = mut * w / (rho * k) < 1 - 1e-9
    SF = rho * A1 * k / (mut * F2)
    for lab, gU in (("LSQ2", g2), ("LSQ2 without wall rows + 1/2 J", gSnw)):
        for useF in (False, True):
            S = fvops.strain_rate(gU).copy()
            if useF:
                S[limited] = SF[limited]
            Rk, Rw = k_res(S, pk_wall), w_res(S)
            print(f"  {lab:32s} {'+ S_F where limited' if useF else '                   '}: k |R| limited cells {np.abs(Rk[limited]).sum():.3e}, all {np.abs(Rk).sum():.3e}"
                  f" | omega |R| limited cells {np.abs(Rw[limited]).sum():.3e}, all {np.abs(Rw).sum():.3e}")

    print("\n5. omega near-wall ingredients (S from LSQ2; wall cells fixed)")
    wb_cell = wb.copy()
    wb_cell[w3] = w[cw]

    def show(lab, R):
        print(f"  {lab:46s} omega sum|R| {np.abs(R).sum():.3e}  L2 {np.abs(R[layer == 2]).sum():.2e} L3 {np.abs(R[layer == 3]).sum():.2e} "
              f"L4 {np.abs(R[layer == 4]).sum():.2e} deep {np.abs(R[layer == 99]).sum():.2e}")

    show("saved wall-face omega in LSQ rows (baseline)", w_res(S2))
    show("omega LSQ without wall rows", w_res(S2, gw_=ctx.apply(G2nw, w, wb)))
    show("wall-face omega replaced by cell omega", w_res(S2, gw_=ctx.apply(G2, w, wb_cell), wb_=wb_cell))
    show("harmonic Gamma_omega", w_res(S2, gam="harm"))
    show("no cross diffusion", w_res(S2, cdk=0.0))
    show("S from LSQ2 without wall rows + 1/2 J", w_res(fvops.strain_rate(gSnw)))
    show("G_omega = alpha (omega/k) P_k", w_res(S2, gform="wk"))
    show("G_omega = alpha (omega/k) P_k, S from LSQ2nw+1/2J", w_res(fvops.strain_rate(gSnw), gform="wk"))
    show("  + no cross diffusion", w_res(fvops.strain_rate(gSnw), gform="wk", cdk=0.0))


if __name__ == "__main__":
    main()
