"""Reconstruct Fluent's saved face mass flux SV_FLUX from the saved cell fields.

SV_FLUX is positive from c0 to c1 (out of c0); sum|div| = 7.7e-8 kg/s, so it is the
solver's post-solve flux. Fluent's own final continuity residual (2.48e-6 kg/s) is the
imbalance of its flux formula recomputed at nearly these fields, so the true formula
should miss SV_FLUX by ~1e-5 of |F| per face. The tests:
  1. boundary faces;
  2. regressions of F - (mean of unlimited LSQ side reconstructions) on candidate terms,
     on interior faces >= 3 hops from boundaries whose two cells the limiter leaves alone
     (plus one on all deep faces for the limiter itself);
  3. the assembled formula against Kutta's current kernel formula on every face (error, divergence);
  4. the pressure-outlet formula.
a_P is a proxy built from the saved fluxes (upwind outflow + mu_eff diffusion), not Fluent's.
Usage: python -I face_flux.py [case_dir] [stem] [cache_dir]   (see cff.py)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import cff
import fvops
import mut_limiter as ml

FLUENT_CONTINUITY_L1 = 2.475884062e-6  # kg/s, last residual sample x normalisation (residuals.py)


def momentum_ap_proxy(ctx, flux, mu_eff, dtau):
    """(convective outflow, interior diffusion, boundary diffusion, pseudo-time rho V/dtau) per cell."""
    c0, c1, I, B, n = ctx.c0, ctx.c1, ctx.I, ctx.B, ctx.n
    out = np.zeros(n)
    np.add.at(out, c0, np.maximum(flux, 0))
    np.add.at(out, c1[I], np.maximum(-flux[I], 0))
    dif = np.zeros(n)
    Df = 0.5 * (mu_eff[c0[I]] + mu_eff[c1[I]]) * ctx.area[I] ** 2 / np.abs(np.einsum("ij,ij->i", ctx.A[I], ctx.d))
    np.add.at(dif, c0[I], Df)
    np.add.at(dif, c1[I], Df)
    dn = np.abs(np.einsum("ij,ij->i", ctx.fc[B] - ctx.cen[c0[B]], ctx.nhat[B]))
    np.add.at(dif, c0[B], np.where(ctx.zone_type[B] != 7, mu_eff[c0[B]] * ctx.area[B] / dn, 0))
    return out, dif, ctx.vol / dtau


def fit(y, cols, names, label):
    M = np.stack(cols, 1)
    c, *_ = np.linalg.lstsq(M, y, rcond=None)
    res = y - M @ c
    coef = ", ".join(f"{nm} {ci:+.3f}" for nm, ci in zip(names, c))
    print(f"  {label:58s} {coef:60s} R2 {1 - np.sum(res**2) / np.sum(y**2):.3f}")


def divergence(ctx, flux_all):
    dv = np.zeros(ctx.n)
    np.add.at(dv, ctx.c0, flux_all)
    np.add.at(dv, ctx.c1[ctx.I], -flux_all[ctx.I])
    return np.abs(dv)


def main():
    case_dir, stem, cache = cff.args()
    m = cff.load(case_dir, stem, cache)
    ctx = fvops.Ctx(m)
    n, c0, c1, I, B = ctx.n, ctx.c0, ctx.c1, ctx.I, ctx.B
    g = lambda name: cff.cell_field(m["cells"], name, n)  # noqa: E731
    rho, p, mu, mut = g("SV_DENSITY"), g("SV_P"), g("SV_MU_LAM"), g("SV_MU_T")
    U = np.stack([g("SV_U"), g("SV_V"), g("SV_W")], 1)
    Ub = np.stack([cff.cell_field(m["faces"], s, ctx.nbf) for s in ("SV_U", "SV_V", "SV_W")], 1)
    flux = cff.cell_field(m["faces"], "SV_FLUX", c0.size)
    pface = cff.cell_field(m["faces"], "SV_P", c0.size)
    pb = pface[: ctx.nbf]
    settings = cff_settings(case_dir, stem)
    dtau = settings["pseudo-auto-time-step"]
    S = -ctx.A  # out of c0
    print(f"faces {c0.size}, interior {I.size}; saved sum|div| {divergence(ctx, flux).sum():.2e} kg/s; "
          f"Fluent final continuity L1 {FLUENT_CONTINUITY_L1:.3e} kg/s; stored global pseudo dt {dtau}")

    print("\n1. Boundary faces: SV_FLUX vs rho U_face . S")
    for z in m["zones"]:
        if z["type"] == 2:
            continue
        s = np.arange(z["lo"], z["hi"])
        pred = rho[c0[s]] * np.einsum("ij,ij->i", Ub[s], S[s])
        scale = max(np.abs(flux[s]).max(), 1e-300)
        print(f"  {z['name']:18s} n={s.size:5d} max|F - rho Ub.S| / max|F| = {np.abs(pred - flux[s]).max() / scale:.1e} (max|F| {np.abs(flux[s]).max():.2e})")

    G2 = fvops.lsq(ctx, 2)
    g2 = ctx.apply(G2, U, Ub)
    gl = fvops.barth_jespersen(ctx, U, Ub, g2)
    unlimited_cell = np.all(np.abs(g2 - gl) < 1e-14 * np.abs(g2).max(), axis=(1, 2))
    hops = fvops.hops_from_boundary(ctx)
    deep = np.minimum(hops[c0[I]], hops[c1[I]]) >= 3
    quiet = deep & unlimited_cell[c0[I]] & unlimited_cell[c1[I]]
    rr = rho[c0[I]]
    Fi = flux[I]

    def side_fluxes(gu):
        u0 = U[c0[I]] + np.einsum("ni,nij->nj", ctx.r0, gu[c0[I]])
        u1 = U[c1[I]] + np.einsum("ni,nij->nj", ctx.r1, gu[c1[I]])
        return rr * np.einsum("ij,ij->i", u0, S[I]), rr * np.einsum("ij,ij->i", u1, S[I])

    out, dif, pt = momentum_ap_proxy(ctx, flux, mu + mut, dtau)
    ap = out + dif
    a0, a1 = ap[c0[I]], ap[c1[I]]
    w = a0 / (a0 + a1)
    geo = ctx.area[I] ** 2 / np.abs(np.einsum("ij,ij->i", ctx.A[I], ctx.d))
    d_f = rr * 0.5 * (ctx.vol[c0[I]] + ctx.vol[c1[I]]) / (0.5 * (a0 + a1)) * geo
    gp_lsq = ctx.apply(G2, p, pb)
    gp_mom = ctx.apply(fvops.green_gauss(ctx, fvops.second_order_faces(ctx, G2)), p, pb)
    gp_ggcb = ctx.apply(fvops.green_gauss(ctx, fvops.arithmetic_faces(ctx)), p, pb)

    def jump_rec(gp):
        return (p[c0[I]] + np.einsum("ij,ij->i", gp[c0[I]], ctx.r0)) - (p[c1[I]] + np.einsum("ij,ij->i", gp[c1[I]], ctx.r1))

    def jump_avg(gp):
        return (p[c0[I]] - p[c1[I]]) + 0.5 * np.einsum("ij,ij->i", gp[c0[I]] + gp[c1[I]], ctx.d)

    print(f"\n2. Interior regressions. deep faces {deep.sum()}, deep faces with both cells unlimited {quiet.sum()}")
    v0, v1 = side_fluxes(g2)
    D = Fi - 0.5 * (v0 + v1)
    print(f"  |F - mean(unlimited LSQ2 side recon)| / |F| median on quiet faces {np.median(np.abs(D[quiet]) / np.abs(Fi[quiet])):.2e}")
    print(" 2a. velocity interpolant (no RC), quiet faces, median |err|/|F|:")
    G0, G1 = g2[c0[I]], g2[c1[I]]
    xmid = 0.5 * (ctx.cen[c0[I]] + ctx.cen[c1[I]])
    wd = np.linalg.norm(ctx.r1, axis=1) / (np.linalg.norm(ctx.r0, axis=1) + np.linalg.norm(ctx.r1, axis=1))
    xip = wd[:, None] * ctx.cen[c0[I]] + (1 - wd)[:, None] * ctx.cen[c1[I]]
    forms = {
        "mean of side reconstructions u_c + G_c r_c": None,
        "mean values + mean gradient . (x_f - x_mid)": 0.5 * (U[c0[I]] + U[c1[I]]) + np.einsum("ni,nij->nj", ctx.fc[I] - xmid, 0.5 * (G0 + G1)),
        "dist-weighted values + weighted grad . (x_f - x_ip)": wd[:, None] * U[c0[I]] + (1 - wd)[:, None] * U[c1[I]]
        + np.einsum("ni,nij->nj", ctx.fc[I] - xip, wd[:, None, None] * G0 + (1 - wd)[:, None, None] * G1),
        "mean cell values only": 0.5 * (U[c0[I]] + U[c1[I]]),
    }
    for name, uf in forms.items():
        Fm = 0.5 * (v0 + v1) if uf is None else rr * np.einsum("ij,ij->i", uf, S[I])
        print(f"    {name:52s} {np.median(np.abs(Fi - Fm)[quiet] / np.abs(Fi[quiet])):.2e}")
    print(" 2b. reconstruction-gradient weight, quiet faces, after the RC term below (median |err|/|F|):")
    for pw in (1, 1.5, 2, 2.5, 3):
        gu = ctx.apply(fvops.lsq(ctx, pw), U, Ub)
        s0, s1 = side_fluxes(gu)
        e = np.abs(Fi - (w * s0 + (1 - w) * s1) - d_f * jump_rec(gp_mom))
        print(f"    LSQ 1/d^{pw:<3} {np.median(e[quiet] / np.abs(Fi[quiet])):.2e}")
    print(" 2c. single-term coefficients (1.0 = the term is present at full strength):")
    y = D[quiet]
    lv0, lv1 = side_fluxes(gl)
    lim_term = 0.5 * (lv0 + lv1) - 0.5 * (v0 + v1)
    yd = D[deep]
    fit(yd, [lim_term[deep], (d_f * jump_rec(gp_lsq))[deep]], ["limited-recon", "RC"], "all deep faces: limited minus unlimited recon + RC")
    asym = 0.5 * (a0 - a1) / (a0 + a1) * (v0 - v1)
    for label, J in (("RC jump, momentum gradient GG(2nd-order p_f)", jump_rec(gp_mom)), ("RC jump, LSQ 1/d^2 gradient", jump_rec(gp_lsq)),
                     ("RC jump, Green-Gauss cell gradient", jump_rec(gp_ggcb)), ("RC averaged-gradient form, momentum grad", jump_avg(gp_mom)),
                     ("RC averaged-gradient form, LSQ grad", jump_avg(gp_lsq))):
        fit(y, [(d_f * J)[quiet], asym[quiet]], ["RC", "a_P-weighted recon"], label)
    jfit = ctx.apply(fvops.jump_correction(ctx, G2), U, Ub)
    j0, j1 = side_fluxes(jfit)
    fit(y, [(d_f * jump_rec(gp_mom))[quiet], asym[quiet], (0.5 * (j0 + j1) - 0.5 * rr * np.einsum("ij,ij->i", U[c0[I]] + U[c1[I]], S[I]))[quiet]],
        ["RC", "a_P", "jump-fit (WFGC cand.)"], "does the S-derived jump correction enter the flux?")
    for label, apv in (("convection only", out), ("convection + diffusion", out + dif), ("+ pseudo-time rho V/dtau", out + dif + pt)):
        b0, b1 = apv[c0[I]], apv[c1[I]]
        dd = rr * 0.5 * (ctx.vol[c0[I]] + ctx.vol[c1[I]]) / (0.5 * (b0 + b1)) * geo
        fit(y, [(dd * jump_rec(gp_mom))[quiet], (0.5 * (b0 - b1) / (b0 + b1) * (v0 - v1))[quiet]], ["RC", "a_P-weighted recon"], "a_P proxy: " + label)

    print("\n3. Assembled formulas on every interior face (saved boundary fluxes kept for the divergence):")
    layer = ml.wall_layers(m, 8)
    s0, s1 = side_fluxes(gl)
    rx = w[:, None] * ctx.r0 + (1 - w)[:, None] * ctx.r1
    uk = w[:, None] * U[c0[I]] + (1 - w)[:, None] * U[c1[I]] + 0.5 * np.einsum("ni,nij->nj", rx, G0 + G1)
    d_kutta = rr * 0.5 * (ctx.vol[c0[I]] / a0 + ctx.vol[c1[I]] / a1) * geo
    g_nowall = fvops.barth_jespersen(ctx, U, Ub, ctx.apply(fvops.lsq(ctx, 2, wall_rows=False), U, Ub))
    n0, n1 = side_fluxes(g_nowall)
    formulas = {
        "Fluent candidate: a_P-wtd limited recon + RC(momentum grad)": w * s0 + (1 - w) * s1 + d_f * jump_rec(gp_mom),
        "  same, wall faces dropped from wall-cell LSQ": w * n0 + (1 - w) * n1 + d_f * jump_rec(gp_mom),
        "  same, unlimited velocity gradient": w * v0 + (1 - w) * v1 + d_f * jump_rec(gp_mom),
        "  same, RC jump with LSQ pressure gradient": w * s0 + (1 - w) * s1 + d_f * jump_rec(gp_lsq),
        "Kutta kernels today (fv_momentum_faces + fv_rhie_chow)": rr * np.einsum("ij,ij->i", uk, S[I]) + d_kutta * jump_rec(gp_lsq),
    }
    for name, Fint in formulas.items():
        Fa = flux.copy()
        Fa[I] = Fint
        dv = divergence(ctx, Fa)
        e = np.abs(Fint - Fi) / np.abs(Fi)
        by = " ".join(f"L{L}:{dv[layer == L].sum():.1e}" for L in (1, 2, 3, 99))
        print(f"  {name:62s} |err|/|F| median {np.median(e):.2e} p90 {np.quantile(e, .9):.2e} | sum|div| {dv.sum():.3e} kg/s ({by})")

    print("\n4. Pressure outlet: F = rho S.u_b + d_b [(p_c + G r_b) - p_b], d_b = rho V/a_P |A|^2/(A.r_b)")
    z = next(z for z in m["zones"] if z["type"] == 5)
    s = np.arange(z["lo"], z["hi"])
    cc = c0[s]
    rb = ctx.fc[s] - ctx.cen[cc]
    db = rho[cc] * ctx.vol[cc] / ap[cc] * ctx.area[s] ** 2 / np.abs(np.einsum("ij,ij->i", ctx.A[s], rb))
    for vname, ub in (("cell velocity", U[cc]), ("reconstructed velocity (Kutta)", U[cc] + np.einsum("ni,nij->nj", rb, g2[cc]))):
        for jname, gp in (("momentum grad", gp_mom), ("LSQ grad", gp_lsq)):
            Jb = p[cc] + np.einsum("ij,ij->i", gp[cc], rb) - pface[s]
            base_ = flux[s] - rho[cc] * np.einsum("ij,ij->i", ub, S[s])
            c = np.dot(db * Jb, base_) / np.dot(db * Jb, db * Jb)
            print(f"  u_b = {vname:31s} J with {jname:14s}: |err|/|F| median no RC {np.median(np.abs(base_) / np.abs(flux[s])):.2e}, "
                  f"RC x1 {np.median(np.abs(base_ - db * Jb) / np.abs(flux[s])):.2e}, fitted RC coef {c:.3f}")


def cff_settings(case_dir, stem):
    """Numeric rpvars needed here, parsed (never evaluated) from the case file's Rampant Variables."""
    import re

    import h5py

    with h5py.File(os.path.join(case_dir, stem + ".cas.h5"), "r") as f:
        raw = f["settings/Rampant Variables"][()]
    text = (raw[0] if not isinstance(raw, bytes) else raw).decode(errors="replace")
    out = {}
    for key in ("pseudo-auto-time-step",):
        hit = re.search(r"\(" + re.escape(key) + r" ([-0-9.e+]+)\)", text)
        out[key] = float(hit.group(1))
    return out


if __name__ == "__main__":
    main()
