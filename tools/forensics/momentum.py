"""Steady momentum residual at Fluent's saved state, for candidate discretisations.

The saved SV_FLUX is the convecting face flux (it is the solver's own flux), so the
residual of the momentum equations isolates the convection, diffusion, pressure and
boundary operators. Per cell P and component i (outward area A_f from P):

  R_P = sum_f F_f u_f                                   SOU upwind face value
      + sum_f p_f A_f                                   pressure force
      - sum_f mu_f [alpha_f (u_N - u_P) + g_f.(A_f - alpha_f d_f)]   Laplacian, alpha = |A|^2/(A.d)
      - sum_f mu_f (g_f^T A_f)                          transpose part of the stress
      - wall / symmetry / inlet / outlet boundary forces

Fluent's scaled residual is sum_P |R_P| / sum_P a_P |U_P| (TG 37-28). The data file
stores both, per iteration, in /results/residuals (col0, col1).
Building blocks are here; momentum_residual.py runs the candidates.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5py
import numpy as np
import scipy.sparse as sp

import cff
import fvops


def fluent_residuals(case_dir, stem, eqs=("x-velocity", "y-velocity", "z-velocity", "k", "omega", "continuity")):
    """{eq: (col0 raw, col1 normalisation)} at the last stored iteration."""
    with h5py.File(os.path.join(case_dir, stem + ".dat.h5"), "r") as f:
        g = f["results/residuals/phase-1"]
        return {eq: tuple(g[eq]["data"][-1, :2]) for eq in eqs}


class State:
    """Saved Fluent fields in fvops conventions (A_f points into c0; SV_FLUX positive c0 -> c1)."""

    def __init__(self, case_dir, stem, cache):
        m = cff.load(case_dir, stem, cache)
        self.m = m
        self.ctx = ctx = fvops.Ctx(m)
        n, nbf, nf = ctx.n, ctx.nbf, ctx.c0.size
        g = lambda name: cff.cell_field(m["cells"], name, n)  # noqa: E731
        fb = lambda name: cff.cell_field(m["faces"], name, nbf)  # noqa: E731
        self.rho, self.mu, self.mut = g("SV_DENSITY"), g("SV_MU_LAM"), g("SV_MU_T")
        self.k, self.w, self.y = g("SV_K"), g("SV_O"), g("SV_WALL_DIST")
        self.p = g("SV_P")
        self.U = np.stack([g("SV_U"), g("SV_V"), g("SV_W")], 1)
        self.Ub = np.stack([fb("SV_U"), fb("SV_V"), fb("SV_W")], 1)
        self.kb, self.wb = fb("SV_K"), fb("SV_O")
        self.F = cff.cell_field(m["faces"], "SV_FLUX", nf)
        self.pf = cff.cell_field(m["faces"], "SV_P", nf)
        self.pb = self.pf[:nbf]
        sh = cff.cell_field(m["faces"], "SV_WALL_SHEAR", nbf)
        self.wall_shear = np.where(np.isfinite(sh), sh, 0.0)
        self.zt = ctx.zone_type[:nbf]
        self.rb = ctx.fc[: nbf] - ctx.cen[ctx.c0[:nbf]]
        self.dn = np.einsum("ij,ij->i", self.rb, ctx.nhat[:nbf])  # > 0: RH normal points into the fluid
        self.alpha = ctx.area[ctx.I] ** 2 / np.einsum("ij,ij->i", ctx.A[ctx.I], -ctx.d)  # A into c0, d c0->c1
        self.res = fluent_residuals(case_dir, stem)
        self.case_dir, self.stem = case_dir, stem

    def zone(self, ztype):
        return np.nonzero(self.zt == ztype)[0]


# ---------------------------------------------------------------- gradients

def lsq_operator(ctx, power=2.0, drop_types=(), wall_rows=True):
    """LSQ 1/d^power over face neighbours; boundary faces of zone types in drop_types are not rows."""
    keep = np.ones(ctx.nbf, bool)
    if not wall_rows:
        keep &= ~ctx.is_wall_face
    for t in drop_types:
        keep &= ctx.zone_type[: ctx.nbf] != t
    Bk = ctx.B[keep]
    I, c0, c1 = ctx.I, ctx.c0, ctx.c1
    cell = np.concatenate([c0[I], c1[I], c0[Bk]])
    dx = np.concatenate([ctx.d, -ctx.d, ctx.fc[Bk] - ctx.cen[c0[Bk]]])
    V = sp.vstack([ctx.unit_rows(c1[I]), ctx.unit_rows(c0[I]), ctx.unit_rows(ctx.n + Bk)]).tocsr()
    return fvops.rows_to_operator(ctx, cell, fvops.lsq_coef(ctx, cell, dx, np.linalg.norm(dx, axis=1) ** -float(power)), V)


def bj_limit(ctx, phi, phib, g, include_boundary=True, eps=0.0, smooth=False, vector=False):
    """Barth-Jespersen cell-to-face limiter variants on g[n,3(,m)].

    include_boundary: boundary-face values enter the min/max stencil.
    eps: Fluent's recon/limit/epsilon style guard, added as |dmax| < eps*|phi range| -> no limiting (scaled).
    smooth: Venkatakrishnan-like smooth min(1, y) replaced by (y^2+2y)/(y^2+y+2).
    vector: one factor per cell for all components (min over components).
    """
    phi = np.asarray(phi, float)
    scalar = phi.ndim == 1
    if scalar:
        phi, phib, g = phi[:, None], np.asarray(phib)[:, None], g[:, :, None]
    c0, c1, I, B, n = ctx.c0, ctx.c1, ctx.I, ctx.B, ctx.n
    allv = np.concatenate([phi, phib])
    cell = np.concatenate([c0[I], c1[I]] + ([c0[B]] if include_boundary else []))
    oth = np.concatenate([c1[I], c0[I]] + ([n + B] if include_boundary else []))
    lo, hi = phi.copy(), phi.copy()
    np.minimum.at(lo, cell, allv[oth])
    np.maximum.at(hi, cell, allv[oth])
    # The face projections always cover every face of the cell (boundary faces included).
    cellf = np.concatenate([c0[I], c1[I], c0[B]])
    faces = np.concatenate([I, I, B])
    dd = np.einsum("ri,rim->rm", ctx.fc[faces] - ctx.cen[cellf], g[cellf])
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(dd > 0, (hi[cellf] - phi[cellf]) / dd, np.where(dd < 0, (lo[cellf] - phi[cellf]) / dd, np.inf))
    if smooth:
        y = np.where(np.isfinite(ratio), ratio, 1e30)
        lim = np.where(y < 1e29, (y * y + 2 * y) / (y * y + y + 2), 1.0)
    else:
        lim = np.minimum(ratio, 1.0)
    psi = np.ones_like(phi)
    np.minimum.at(psi, cellf, lim)
    psi = np.clip(psi, 0, 1)
    if eps > 0:
        rng = np.abs(allv).max(axis=0)
        psi = np.where((hi - lo) < eps * rng, 1.0, psi)
    if vector:
        psi = np.repeat(psi.min(axis=1, keepdims=True), psi.shape[1], axis=1)
    out = g * psi[:, None, :]
    return (out[:, :, 0] if scalar else out), (psi[:, 0] if scalar else psi)


# ---------------------------------------------------------------- residual terms
# Each returns [n, 3] momentum contributions to R (net outflow + pressure - viscous forces).

def scatter(ctx, faces, owner_val, neigh_val=None):
    """Add owner_val[k] to c0[faces[k]] and neigh_val[k] to c1[faces[k]] (interior faces)."""
    out = np.zeros((ctx.n,) + owner_val.shape[1:])
    np.add.at(out, ctx.c0[faces], owner_val)
    if neigh_val is not None:
        np.add.at(out, ctx.c1[faces], neigh_val)
    return out


def convection(st, g_up, outlet="cell", inlet="face", g_outlet=None):
    """sum_f F_f u_f with u_f the upwind side's u + g.r (g_up: [n,3,3] already limited if wanted)."""
    ctx, U, F = st.ctx, st.U, st.F
    I, c0, c1 = ctx.I, ctx.c0, ctx.c1
    Fi = F[I]
    u0 = U[c0[I]] + np.einsum("ni,nij->nj", ctx.r0, g_up[c0[I]])
    u1 = U[c1[I]] + np.einsum("ni,nij->nj", ctx.r1, g_up[c1[I]])
    uf = np.where((Fi >= 0)[:, None], u0, u1)
    R = scatter(ctx, I, Fi[:, None] * uf, -Fi[:, None] * uf)
    B = ctx.B
    Fb = F[B]
    cb = c0[B]
    if outlet == "cell":
        ub_out = U[cb]
    elif outlet == "face":
        ub_out = st.Ub
    else:  # reconstructed
        go = g_up if g_outlet is None else g_outlet
        ub_out = U[cb] + np.einsum("ni,nij->nj", st.rb, go[cb])
    ub = np.where((Fb >= 0)[:, None], ub_out, st.Ub)
    np.add.at(R, cb, Fb[:, None] * ub)
    return R


def mass_imbalance_term(st):
    """u_P * sum_f F_f: subtract to use the non-conservative a_P = sum a_nb convention."""
    ctx = st.ctx
    dv = np.zeros(ctx.n)
    np.add.at(dv, ctx.c0, st.F)
    np.add.at(dv, ctx.c1[ctx.I], -st.F[ctx.I])
    return st.U * dv[:, None]


def pressure(st, gp=None, boundary="saved", gp_boundary=None):
    """sum_f p_f A_out. Interior: mean of side reconstructions with gp (None: arithmetic mean).
    boundary: 'saved' face SV_P (= cell p at walls/symmetry/inlets), or 'extrap' p_c + gp.r_b except outlets."""
    ctx, p = st.ctx, st.p
    I, B, c0, c1 = ctx.I, ctx.B, ctx.c0, ctx.c1
    if gp is None:
        pfi = 0.5 * (p[c0[I]] + p[c1[I]])
    else:
        pfi = 0.5 * (p[c0[I]] + np.einsum("ij,ij->i", gp[c0[I]], ctx.r0) + p[c1[I]] + np.einsum("ij,ij->i", gp[c1[I]], ctx.r1))
    A = ctx.A
    R = scatter(ctx, I, -pfi[:, None] * A[I], pfi[:, None] * A[I])
    pb = st.pb.copy()
    if boundary == "extrap":
        gb = gp if gp_boundary is None else gp_boundary
        ext = p[c0[B]] + np.einsum("ij,ij->i", gb[c0[B]], st.rb)
        pb = np.where(st.zt == 5, st.pb, ext)
    np.add.at(R, c0[B], -pb[:, None] * A[B])
    return R, pfi


def face_mu(st, mode="arith"):
    ctx = st.ctx
    me = st.mu + st.mut
    a, b = me[ctx.c0[ctx.I]], me[ctx.c1[ctx.I]]
    if mode == "arith":
        return 0.5 * (a + b)
    if mode == "harm":
        return 2 * a * b / (a + b)
    if mode == "dist":  # normal-distance weights (Kutta's mesh.weight)
        d0 = np.abs(np.einsum("ij,ij->i", ctx.r0, ctx.nhat[ctx.I]))
        d1 = np.abs(np.einsum("ij,ij->i", ctx.r1, ctx.nhat[ctx.I]))
        w = d1 / (d0 + d1)
        return w * a + (1 - w) * b
    if mode == "mut_arith":  # mu + arithmetic mean of mu_t (same as arith for constant mu)
        return st.mu[ctx.c0[ctx.I]] + 0.5 * (st.mut[ctx.c0[ctx.I]] + st.mut[ctx.c1[ctx.I]])
    raise ValueError(mode)


def face_gradient(st, g, mode="arith"):
    ctx = st.ctx
    a, b = g[ctx.c0[ctx.I]], g[ctx.c1[ctx.I]]
    if mode == "arith":
        return 0.5 * (a + b)
    if mode == "dist":
        d0 = np.abs(np.einsum("ij,ij->i", ctx.r0, ctx.nhat[ctx.I]))
        d1 = np.abs(np.einsum("ij,ij->i", ctx.r1, ctx.nhat[ctx.I]))
        w = (d1 / (d0 + d1))[:, None, None]
        return w * a + (1 - w) * b
    raise ValueError(mode)


def interior_viscous(st, g_lap, g_tr, mu_mode="arith", grad_mode="arith", nonorth="overrelaxed",
                     transpose=True, div23=False, g_div=None, face_grad_correct=False):
    """-sum_f [Laplacian flux + transpose flux] on interior faces (outward area from each side).

    nonorth: 'overrelaxed' alpha = |A|^2/(A.d) with correction g.(A - alpha d) (Fluent TG 23-?, Kutta);
             'none' (orthogonal part only); 'minimum' alpha = (A.d)/|d|^2.
    face_grad_correct: replace the face-normal part of the averaged gradient by the compact difference
             before using it in the transpose term (g_f + [(u_N-u_P)/|d| - g_f.e] e).
    """
    ctx, U = st.ctx, st.U
    I, c0, c1 = ctx.I, ctx.c0, ctx.c1
    mf = face_mu(st, mu_mode)
    A0 = -ctx.A[I]  # outward from c0
    d = ctx.d
    du = U[c1[I]] - U[c0[I]]
    gf = face_gradient(st, g_lap, grad_mode)
    if nonorth == "overrelaxed":
        al = st.alpha
        flux = al[:, None] * du + np.einsum("ni,nij->nj", A0 - al[:, None] * d, gf)
    elif nonorth == "minimum":
        al = np.einsum("ij,ij->i", A0, d) / np.einsum("ij,ij->i", d, d)
        flux = al[:, None] * du + np.einsum("ni,nij->nj", A0 - al[:, None] * d, gf)
    elif nonorth == "none":
        flux = st.alpha[:, None] * du
    else:
        raise ValueError(nonorth)
    flux = flux * mf[:, None]
    if transpose:
        gt = face_gradient(st, g_tr, grad_mode)
        if face_grad_correct:
            e = d / np.linalg.norm(d, axis=1)[:, None]
            corr = du / np.linalg.norm(d, axis=1)[:, None] - np.einsum("ni,nij->nj", e, gt)
            gt = gt + e[:, :, None] * corr[:, None, :]
        flux = flux + mf[:, None] * np.einsum("nij,ni->nj", np.transpose(gt, (0, 2, 1)), A0)  # sum_j du_j/dx_i A_j
    if div23:
        gd = face_gradient(st, g_tr if g_div is None else g_div, grad_mode)
        div = np.einsum("nii->n", gd)
        flux = flux - (2.0 / 3.0) * (mf * div)[:, None] * A0
    return scatter(ctx, I, -flux, flux)


def boundary_viscous(st, g, wall="shear", sym="normal", inlet="full", outlet="none", transpose=True,
                     sym_mu="cell_eff", wall_normal=False):
    """Viscous boundary forces (subtracted: R -= force on fluid)."""
    ctx, U = st.ctx, st.U
    B, c0 = ctx.B, ctx.c0
    R = np.zeros((ctx.n, 3))
    me = st.mu + st.mut
    nh = ctx.nhat[: ctx.nbf]  # into the fluid
    Aout = -ctx.A[: ctx.nbf]
    area = ctx.area[: ctx.nbf]
    # walls
    w = st.zone(3)
    cw = c0[w]
    if wall == "shear":
        np.add.at(R, cw, -st.wall_shear[w])
    elif wall == "shear_neg":
        np.add.at(R, cw, st.wall_shear[w])
    if wall_normal:  # normal viscous stress 2 mu du_n/dn on the wall face: force on fluid -2 mu A u_n/dn n
        un = np.einsum("ij,ij->i", U[cw], nh[w])
        np.add.at(R, cw, (2 * me[cw] * area[w] * un / st.dn[w])[:, None] * nh[w])
    # symmetry: normal velocity driven to zero through mu A/dn on the normal component
    s = st.zone(7)
    cs = c0[s]
    if sym != "none":
        mus = me[cs] if sym_mu == "cell_eff" else st.mu[cs]
        un = np.einsum("ij,ij->i", U[cs], nh[s])
        fac = {"normal": 1.0, "normal2": 2.0}[sym]
        force = -fac * (mus * area[s] / st.dn[s] * un)[:, None] * nh[s]  # force on fluid, opposes u_n
        np.add.at(R, cs, -force)
    # inlets: Dirichlet velocity
    for t in (10,):
        b = st.zone(t)
        cb = c0[b]
        if inlet == "none":
            continue
        rb = st.rb[b]
        alpha_b = area[b] ** 2 / np.einsum("ij,ij->i", Aout[b], rb)
        du = st.Ub[b] - U[cb]
        flux = alpha_b[:, None] * du
        if inlet == "full":
            flux = flux + np.einsum("ni,nij->nj", Aout[b] - alpha_b[:, None] * rb, g[cb])
        if transpose:  # component i: sum_j du_j/dx_i A_j
            flux = flux + np.einsum("nij,nj->ni", g[cb], Aout[b])
        np.add.at(R, cb, -(me[cb][:, None] * flux))
    # pressure outlet
    o = st.zone(5)
    co = c0[o]
    if outlet in ("transpose", "full"):
        flux = np.einsum("nij,nj->ni", g[co], Aout[o])  # transpose part
        if outlet == "full":  # plus the normal-derivative (Laplacian) part
            flux = flux + np.einsum("nji,nj->ni", g[co], Aout[o])
        np.add.at(R, co, -(me[co][:, None] * flux))
    return R


def scaled(st, R, comps=(0, 1, 2)):
    names = ("x-velocity", "y-velocity", "z-velocity")
    return [np.abs(R[:, i]).sum() / st.res[names[i]][1] for i in comps]


# ---------------------------------------------------------------- scalar transport (energy, k, omega)

def scalar_convection(st, phi, phib, g):
    """sum_f F_f phi_f, SOU upwind face value phi_up + g_up.r_up; boundary inflow takes the face value,
    outflow the cell value (Fluent's outlet treatment, cf. momentum)."""
    ctx, F = st.ctx, st.F
    I, B, c0, c1 = ctx.I, ctx.B, ctx.c0, ctx.c1
    p0 = phi[c0[I]] + np.einsum("ij,ij->i", g[c0[I]], ctx.r0)
    p1 = phi[c1[I]] + np.einsum("ij,ij->i", g[c1[I]], ctx.r1)
    pf = np.where(F[I] >= 0, p0, p1)
    R = np.zeros(ctx.n)
    np.add.at(R, c0[I], F[I] * pf)
    np.add.at(R, c1[I], -F[I] * pf)
    np.add.at(R, c0[B], F[B] * np.where(F[B] >= 0, phi[c0[B]], phib))
    return R


def scalar_diffusion(st, phi, phib, g, gamma, inlet=True, gamma_face="arith", nonorth=True):
    """-sum_f Gamma_f [alpha (phi_N - phi_P) + gbar.(A - alpha d)]; Dirichlet inlets use Gamma_P A/dn.
    Walls, symmetry and outlet carry no diffusive flux (zero-gradient / adiabatic)."""
    ctx = st.ctx
    I, c0, c1 = ctx.I, ctx.c0, ctx.c1
    a, b = gamma[c0[I]], gamma[c1[I]]
    gf = 0.5 * (a + b) if gamma_face == "arith" else 2 * a * b / (a + b)
    fl = st.alpha * (phi[c1[I]] - phi[c0[I]])
    if nonorth:
        fl = fl + np.einsum("ni,ni->n", -ctx.A[I] - st.alpha[:, None] * ctx.d, 0.5 * (g[c0[I]] + g[c1[I]]))
    fl = gf * fl
    R = np.zeros(ctx.n)
    np.add.at(R, c0[I], -fl)
    np.add.at(R, c1[I], fl)
    if inlet:
        s = st.zone(10)
        cc = c0[s]
        np.add.at(R, cc, -gamma[cc] * ctx.area[s] / np.abs(st.dn[s]) * (phib[s] - phi[cc]))
    return R


def outflow_coefficient(st):
    """sum_f max(F_out, 0) per cell (first-order upwind part of a_P)."""
    ctx = st.ctx
    out = np.zeros(ctx.n)
    np.add.at(out, ctx.c0, np.maximum(st.F, 0))
    np.add.at(out, ctx.c1[ctx.I], np.maximum(-st.F[ctx.I], 0))
    return out


def face_sum(st, coef):
    """sum over interior faces of coef_f, accumulated to both cells."""
    ctx = st.ctx
    out = np.zeros(ctx.n)
    np.add.at(out, ctx.c0[ctx.I], coef)
    np.add.at(out, ctx.c1[ctx.I], coef)
    return out


def boundary_cells(st):
    """Bit mask per cell: 1 wall, 2 symmetry, 4 inlet, 8 outlet."""
    tb = np.zeros(st.ctx.n, int)
    for t, bit in ((3, 1), (7, 2), (10, 4), (5, 8)):
        np.bitwise_or.at(tb, st.ctx.c0[st.zone(t)], bit)
    return tb


def pseudo_dt(case_dir, stem):
    """Global pseudo time step, parsed (never evaluated) from the case rpvars."""
    import re

    with h5py.File(os.path.join(case_dir, stem + ".cas.h5"), "r") as f:
        text = f["settings/Rampant Variables"][0].decode(errors="replace")
    return float(re.search(r"\(pseudo-auto-time-step ([-0-9.e+]+)\)", text).group(1))
