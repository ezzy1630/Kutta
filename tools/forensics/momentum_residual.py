"""Momentum residual of Fluent's saved state under candidate discretisations (docs/research/fluent-momentum-forensics.md).

Sections:
  1. a_P from Fluent's own residual normalisation (col1 = sum a_P |U| per component).
  2. Baseline residual, term sizes and where it sits.
  3. Convective (SOU) reconstruction gradient and limiter.
  4. Pressure face reconstruction and boundary pressure.
  5. Viscous terms, including the jump-corrected (WFGC) gradient's role.
  6. Boundary treatments, scored on the cells next to each boundary type.
  7. Assembled best candidate vs a Kutta-like assembly.
  8. The a_P inside the Rhie-Chow coefficient d_f.
Usage: python -I momentum_residual.py [case_dir] [stem] [cache_dir]   (see cff.py)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import cff
import fvops
import momentum as mo
import mut_limiter as ml


def fmt(v):
    return " ".join(f"{x:.3e}" for x in v)


def main():
    case_dir, stem, cache = cff.args()
    st = mo.State(case_dir, stem, cache)
    ctx = st.ctx
    n, I, c0, c1 = ctx.n, ctx.I, ctx.c0, ctx.c1
    U, Ub, p, pb = st.U, st.Ub, st.p, st.pb
    me = st.mu + st.mut
    Um = np.linalg.norm(U, axis=1)
    dtau = mo.pseudo_dt(case_dir, stem)
    names = ("x-velocity", "y-velocity", "z-velocity")
    col0 = np.array([st.res[k][0] for k in names])
    col1 = np.array([st.res[k][1] for k in names])
    print(f"Fluent last iteration: raw sum|R| x/y/z {fmt(col0)} N; normalisation {fmt(col1)}; scaled {fmt(col0 / col1)}")

    import h5py

    with h5py.File(os.path.join(case_dir, stem + ".dat.h5"), "r") as f:
        hist = {e: f["results/residuals/phase-1/" + e + "/data"][:] for e in ("continuity", "x-velocity", "energy", "k", "omega")}
    print("history (scaled = col0/col1); the run has three phases, 1-100, 101-200, 201-294:")
    for it in (90, 91, 92, 93, 98, 99, 100, 190, 199, 292, 293, 294):
        print(f"  it {it:3d} " + " ".join(f"{e[:5]} {hist[e][it - 1, 0] / hist[e][it - 1, 1]:.3e}" for e in hist))

    print("\n0. saved boundary-face values")
    for t, lab in ((3, "wall"), (7, "symmetry"), (10, "inlet"), (5, "outlet")):
        s = st.zone(t)
        print(f"  {lab:8s} max|p_b - p_c| {np.abs(pb[s] - p[c0[s]]).max():.2e}")
    s7, s5 = st.zone(7), st.zone(5)
    nh7, nh5 = ctx.nhat[s7], ctx.nhat[s5]
    proj = U[c0[s7]] - np.einsum("ij,ij->i", U[c0[s7]], nh7)[:, None] * nh7
    print(f"  symmetry U_b vs tangential projection of U_P: max|diff| {np.abs(Ub[s7] - proj).max():.1e}")
    dU = Ub[s5] - U[c0[s5]]
    tang = dU - np.einsum("ij,ij->i", dU, nh5)[:, None] * nh5
    print(f"  outlet U_b - U_P: tangential part max {np.abs(tang).max():.1e}; normal part median |.|/|U_P| "
          f"{np.median(np.abs(np.einsum('ij,ij->i', dU, nh5)) / Um[c0[s5]]):.1e} (flux-consistent normal velocity)")

    print("\n1. a_P from Fluent's normalisation col1 = sum_P a_P,i |U_P|")
    out = mo.outflow_coefficient(st)
    dif = mo.face_sum(st, mo.face_mu(st, "arith") * st.alpha)
    pt = st.rho * ctx.vol / dtau
    w3, s7, s10 = st.zone(3), st.zone(7), st.zone(10)
    cw, nhw = c0[w3], ctx.nhat[w3]
    Ut = U[cw] - np.einsum("ij,ij->i", U[cw], nhw)[:, None] * nhw
    cwA = np.linalg.norm(st.wall_shear[w3], axis=1) / np.linalg.norm(Ut, axis=1)  # tau_w A / |U_t| (SV_WALL_SHEAR is a force)
    wall_i = np.array([np.sum(cwA * (1 - nhw[:, i] ** 2) * Um[cw]) for i in range(3)])
    sym_i = np.array([np.sum(me[c0[s7]] * ctx.area[s7] / np.abs(st.dn[s7]) * ctx.nhat[s7][:, i] ** 2 * Um[c0[s7]]) for i in range(3)])
    inlet = np.sum(me[c0[s10]] * ctx.area[s10] / np.abs(st.dn[s10]) * Um[c0[s10]])
    print(f"  pieces: outflow {np.sum(out * Um):.6f}, interior mu_eff alpha {np.sum(dif * Um):.6f}, rho V/dtau {np.sum(pt * Um):.6f} (dtau {dtau:.6g})")
    print(f"  wall c_w A (1 - n_i^2) x/y/z {fmt(wall_i)}; symmetry mu_eff A/dn n_i^2 x/y/z {fmt(sym_i)}; inlet mu_eff A/dn {inlet:.3e}")
    print(f"  (col1_x - col1_y) / (wall_x - wall_y) = {(col1[0] - col1[1]) / (wall_i[0] - wall_i[1]):.5f}   (walls are the only x/y-asymmetric term)")
    sym_from_data = (col1[2] - col1[0]) / 2 - (wall_i[2] - wall_i[0])
    print(f"  symmetry term implied by (col1_z - col1_x)/2 - wall: {sym_from_data:.6e} vs mu_eff A/dn n_z^2: {sym_i[2]:.6e} (ratio {sym_from_data / sym_i[2]:.5f})")
    for lab, tot in (("2 (out + dif + pt + wall_i + inlet)", 2 * (np.sum((out + dif + pt) * Um) + wall_i + inlet)),
                     ("2 (out + dif + wall_i + inlet), no pseudo-time", 2 * (np.sum((out + dif) * Um) + wall_i + inlet)),
                     ("1 (out + dif + pt + wall_i + inlet)", np.sum((out + dif + pt) * Um) + wall_i + inlet)):
        print(f"  {lab:48s} x/y/z {fmt(tot)}  ratio to col1 {' '.join(f'{r:.4f}' for r in col1 / tot)}")

    G2 = fvops.lsq(ctx, 2)
    g2 = ctx.apply(G2, U, Ub)
    gl, psi = mo.bj_limit(ctx, U, Ub, g2)
    gp = ctx.apply(G2, p, pb)
    layer = ml.wall_layers(st.m, 8)
    hops = fvops.hops_from_boundary(ctx)
    deep = hops >= 3
    lim = (psi < 1 - 1e-12).any(1)
    tb = mo.boundary_cells(st)
    nfc = st.m["nfaces_per_cell"]

    def show(lab, R):
        print(f"  {lab:52s} scaled x/y/z {fmt(mo.scaled(st, R))} | deep sum|R| {np.abs(R[deep]).sum():.3e}")

    print("\n2. Baseline: SOU with BJ-limited LSQ 1/d^2, Second Order pressure (LSQ 1/d^2), arithmetic mu_eff, "
          "over-relaxed non-orthogonal correction, transpose term, saved wall shear, symmetry mu_eff A/dn on u_n")
    Rc = mo.convection(st, gl)
    Rp = mo.pressure(st, gp)[0]
    Rv = mo.interior_viscous(st, g2, g2)
    Rb = mo.boundary_viscous(st, g2)
    for lab, R in (("term: convection", Rc), ("term: pressure", Rp), ("term: interior viscous", Rv), ("term: boundary viscous", Rb)):
        show(lab, R)
    R = Rc + Rp + Rv + Rb
    show("TOTAL", R)
    show("TOTAL, wall shear with the opposite sign", Rc + Rp + Rv + mo.boundary_viscous(st, g2, wall="shear_neg"))
    a = np.abs(R).sum(1)
    tot = a.sum()
    sc = out * Um
    print(f"  Fluent's own: sum|R| {col0.sum():.3e} N (x+y+z); here {tot:.3e} N  -> {tot / col0.sum():.0f}x")
    for lab, mk in [(f"wall layer {L}", layer == L) for L in (1, 2, 3, 4, 5)] + [
            ("deeper than layer 8", layer == 99), (">= 3 hops from any boundary", deep), (">= 3 hops, limiter inactive", deep & ~lim),
            ("limiter active (any component)", lim), ("cells with 6 faces", nfc == 6), ("cells with > 6 faces", nfc > 6)]:
        print(f"    {lab:34s} cells {mk.sum():6d}  share {a[mk].sum() / tot:6.1%}  mean |R|/cell {a[mk].mean():.2e}  median |R|/(out|U|) {np.median(a[mk] / sc[mk]):.2e}")
    T1 = mo.convection(st, 0 * gl)
    for lab, T in (("SOU correction", Rc - T1), ("pressure SO correction", Rp - mo.pressure(st, None)[0])):
        print(f"  scale scan of the {lab} (deep sum|R| for R + (c-1) T):",
              " ".join(f"{c}:{np.abs((R + (c - 1) * T)[deep]).sum():.3e}" for c in (0.9, 0.95, 1.0, 1.05, 1.1)))

    rest_c = Rp + Rv + Rb
    print("\n3. Convective reconstruction gradient / limiter (other terms as baseline)")
    J = ctx.apply(fvops.jump_correction(ctx, G2), U, Ub)
    for pw in (1, 1.5, 2, 2.5, 3):
        show(f"LSQ 1/d^{pw}, BJ", mo.convection(st, mo.bj_limit(ctx, U, Ub, ctx.apply(fvops.lsq(ctx, pw), U, Ub))[0]) + rest_c)
    for t in (-0.25, 0.25, 0.5):
        show(f"LSQ2 + {t} x jump correction (WFGC candidate), BJ", mo.convection(st, mo.bj_limit(ctx, U, Ub, g2 + t * J)[0]) + rest_c)
    show("LSQ2 unlimited", mo.convection(st, g2) + rest_c)
    show("LSQ2 BJ, boundary values excluded from min/max", mo.convection(st, mo.bj_limit(ctx, U, Ub, g2, include_boundary=False)[0]) + rest_c)
    show("LSQ2 BJ, one factor for all components", mo.convection(st, mo.bj_limit(ctx, U, Ub, g2, vector=True)[0]) + rest_c)
    show("LSQ2 smooth (Venkatakrishnan-function) limiter", mo.convection(st, mo.bj_limit(ctx, U, Ub, g2, smooth=True)[0]) + rest_c)
    show("Green-Gauss of second-order faces, BJ", mo.convection(st, mo.bj_limit(ctx, U, Ub, ctx.apply(fvops.green_gauss(ctx, fvops.second_order_faces(ctx, G2)), U, Ub))[0]) + rest_c)
    show("Green-Gauss cell-based, BJ", mo.convection(st, mo.bj_limit(ctx, U, Ub, ctx.apply(fvops.green_gauss(ctx, fvops.arithmetic_faces(ctx)), U, Ub))[0]) + rest_c)
    show("velocity LSQ without wall rows, BJ", mo.convection(st, mo.bj_limit(ctx, U, Ub, ctx.apply(fvops.lsq(ctx, 2, wall_rows=False), U, Ub))[0]) + rest_c)

    rest_p = Rc + Rv + Rb
    print("\n4. Pressure face reconstruction (Second Order) and boundary pressure")
    for pw in (1, 2, 3):
        show(f"p LSQ 1/d^{pw}", mo.pressure(st, ctx.apply(fvops.lsq(ctx, pw), p, pb))[0] + rest_p)
    show("p LSQ2, BJ-limited", mo.pressure(st, mo.bj_limit(ctx, p, pb, gp)[0])[0] + rest_p)
    show("p Green-Gauss of second-order faces", mo.pressure(st, ctx.apply(fvops.green_gauss(ctx, fvops.second_order_faces(ctx, G2)), p, pb))[0] + rest_p)
    show("p Green-Gauss cell-based", mo.pressure(st, ctx.apply(fvops.green_gauss(ctx, fvops.arithmetic_faces(ctx)), p, pb))[0] + rest_p)
    show("p LSQ2 + 1/2 jump correction", mo.pressure(st, gp + 0.5 * ctx.apply(fvops.jump_correction(ctx, G2), p, pb))[0] + rest_p)
    show("p arithmetic mean (no reconstruction)", mo.pressure(st, None)[0] + rest_p)
    Rpx = mo.pressure(st, gp, boundary="extrap")[0]
    show("p LSQ2, boundary faces p_c + G.r_b (not outlet)", Rpx + rest_p)
    s7 = st.zone(7)
    rb7, nh7 = st.rb[s7], ctx.nhat[s7]
    t7 = rb7 - np.einsum("ij,ij->i", rb7, nh7)[:, None] * nh7
    Rpt = Rpx.copy()  # Kutta's symmetry face pressure uses only the tangential offset: p_c + G.t
    np.add.at(Rpt, c0[s7], -(np.einsum("ij,ij->i", gp[c0[s7]], t7 - rb7))[:, None] * ctx.A[s7])
    symc = (mo.boundary_cells(st) & 2) > 0
    print(f"  symmetry cells sum|R|: p_c + G.r_b {np.abs((Rpx + rest_p)[symc]).sum():.3e}; p_c + G.t (Kutta) {np.abs((Rpt + rest_p)[symc]).sum():.3e}; "
          f"saved p_c {np.abs((Rp + rest_p)[symc]).sum():.3e}")

    rest_v = Rc + Rp
    print("\n5. Viscous terms (boundary viscous held as baseline)")
    gS = (g2 + 0.5 * J)
    G2nw = fvops.lsq(ctx, 2, wall_rows=False)
    gSnw = ctx.apply(G2nw, U, Ub) + 0.5 * ctx.apply(fvops.jump_correction(ctx, G2nw), U, Ub)
    for lab, kw in (("baseline", {}), ("no transpose", dict(transpose=False)), ("+ -2/3 mu div(u)", dict(div23=True)),
                    ("harmonic mu_eff", dict(mu_mode="harm")), ("normal-distance-weighted mu_eff (Kutta)", dict(mu_mode="dist")),
                    ("distance-weighted face gradient (Kutta)", dict(grad_mode="dist")), ("no non-orthogonal correction", dict(nonorth="none")),
                    ("minimum-correction non-orthogonal split", dict(nonorth="minimum")),
                    ("face gradient with compact normal part in transpose", dict(face_grad_correct=True))):
        show(lab, rest_v + mo.interior_viscous(st, g2, g2, **kw) + Rb)
    for lab, gg in (("WFGC gradient LSQ2 + 1/2 J in Laplacian + transpose", gS), ("WFGC, no wall rows in base", gSnw)):
        show(lab, rest_v + mo.interior_viscous(st, gg, gg) + mo.boundary_viscous(st, gg))
        show(lab + " (transpose only)", rest_v + mo.interior_viscous(st, g2, gg) + Rb)

    print("\n6. Boundary treatments: sum|R| (x+y+z) on the cells next to each boundary type")
    regions = {"wall cells": (tb & 1) > 0, "layer 2": layer == 2, "sym cells": (tb & 2) > 0, "inlet cells": (tb & 4) > 0,
               "outlet cells": (tb & 8) > 0, "all": np.ones(n, bool)}

    def build(vel_drop=(), p_drop=None, pbnd="saved", outlet_conv="cell", wall="shear", wall_normal=False, sym="normal", inlet="full", outlet_visc="none"):
        gu = ctx.apply(mo.lsq_operator(ctx, 2, drop_types=vel_drop), U, Ub)
        gg = mo.bj_limit(ctx, U, Ub, gu)[0]
        gpp = gp if p_drop is None else ctx.apply(mo.lsq_operator(ctx, 2, drop_types=p_drop), p, pb)
        return (mo.convection(st, gg, outlet=outlet_conv) + mo.pressure(st, gpp, boundary=pbnd)[0] + mo.interior_viscous(st, gu, gu)
                + mo.boundary_viscous(st, gu, wall=wall, wall_normal=wall_normal, sym=sym, inlet=inlet, outlet=outlet_visc))

    cands = {"baseline": {}, "boundary p = p_c + G.r_b": dict(pbnd="extrap"),
             "  + pressure LSQ without wall rows": dict(pbnd="extrap", p_drop=(3,)),
             "  + pressure LSQ without wall/sym rows": dict(pbnd="extrap", p_drop=(3, 7)),
             "  + pressure LSQ outlet rows only (Kutta)": dict(pbnd="extrap", p_drop=(3, 7, 10)),
             "velocity LSQ without wall rows": dict(vel_drop=(3,)),
             "velocity LSQ without symmetry rows": dict(vel_drop=(7,)), "wall normal viscous stress added": dict(wall_normal=True),
             "no symmetry normal force": dict(sym="none"), "symmetry 2 mu_eff A/dn": dict(sym="normal2"),
             "outlet convects saved face velocity": dict(outlet_conv="face"), "outlet convects u_P + G.r_b": dict(outlet_conv="recon"),
             "outlet carries mu (grad u)^T.A": dict(outlet_visc="transpose"), "outlet carries full mu (grad u + grad u^T).A": dict(outlet_visc="full"),
             "inlet: orthogonal diffusion only": dict(inlet="orth"), "inlet: no diffusion": dict(inlet="none")}
    print(f"  {'candidate':40s}" + "".join(f"{k:>13s}" for k in regions))
    for name, kw in cands.items():
        aa = np.abs(build(**kw)).sum(1)
        print(f"  {name:40s}" + "".join(f"{aa[mk].sum():13.3e}" for mk in regions.values()))

    print("\n7. Assembled")
    best = build(pbnd="extrap", outlet_visc="transpose", sym="none")
    show("best: + extrapolated boundary p, outlet transpose, no symmetry force", best)
    # Kutta-like: velocity LSQ rows only at walls/inlets, pressure LSQ rows only at outlets, distance-weighted mu and
    # face gradient, extrapolated boundary pressure, symmetry face value = projected reconstruction with full diffusion.
    gk_u = ctx.apply(mo.lsq_operator(ctx, 2, drop_types=(5, 7)), U, Ub)
    gk_p = ctx.apply(mo.lsq_operator(ctx, 2, drop_types=(3, 7, 10)), p, pb)
    Rk = (mo.convection(st, mo.bj_limit(ctx, U, Ub, gk_u)[0], outlet="recon") + mo.pressure(st, gk_p, boundary="extrap")[0]
          + mo.interior_viscous(st, gk_u, gk_u, mu_mode="dist", grad_mode="dist") + mo.boundary_viscous(st, gk_u, sym="normal", outlet="full"))
    show("Kutta-like assembly (see doc section 7)", Rk)

    print("\n8. a_P in the Rhie-Chow coefficient (quiet deep faces; coefficient 1 = this a_P is the one)")
    S = -ctx.A
    rr = st.rho[c0[I]]
    geo = ctx.area[I] ** 2 / np.abs(np.einsum("ij,ij->i", ctx.A[I], ctx.d))
    gpm = ctx.apply(fvops.green_gauss(ctx, fvops.second_order_faces(ctx, G2)), p, pb)
    Jp = (p[c0[I]] + np.einsum("ij,ij->i", gpm[c0[I]], ctx.r0)) - (p[c1[I]] + np.einsum("ij,ij->i", gpm[c1[I]], ctx.r1))
    u0 = U[c0[I]] + np.einsum("ni,nij->nj", ctx.r0, gl[c0[I]])
    u1 = U[c1[I]] + np.einsum("ni,nij->nj", ctx.r1, gl[c1[I]])
    s0, s1 = rr * np.einsum("ij,ij->i", u0, S[I]), rr * np.einsum("ij,ij->i", u1, S[I])
    quiet = (np.minimum(hops[c0[I]], hops[c1[I]]) >= 3) & ~lim[c0[I]] & ~lim[c1[I]]
    for lab, ap in (("outflow + mu_eff alpha", out + dif), ("outflow + mu_eff alpha + rho V/dtau", out + dif + pt),
                    ("2 x (outflow + mu_eff alpha + rho V/dtau)", 2 * (out + dif + pt))):
        a0, a1 = ap[c0[I]], ap[c1[I]]
        w = a0 / (a0 + a1)
        X = rr * 0.5 * (ctx.vol[c0[I]] + ctx.vol[c1[I]]) / (0.5 * (a0 + a1)) * geo * Jp
        base = st.F[I] - (w * s0 + (1 - w) * s1)
        c = np.dot(X[quiet], base[quiet]) / np.dot(X[quiet], X[quiet])
        print(f"  {lab:44s} RC coefficient {c:.3f}; median |err|/|F| at coefficient 1: {np.median(np.abs(base - X)[quiet] / np.abs(st.F[I][quiet])):.3e}")


if __name__ == "__main__":
    main()
