"""Sparse cell-gradient, face-reconstruction and limiter operators for the flux/gradient forensics.

A gradient operator is a scipy.sparse matrix G of shape (3n, n + nbf) acting on
phi_all = [phi_cell; phi_boundary_face]; row 3c+i is d(phi)/dx_i in cell c.
Building operators once lets every hypothesis be applied to u, v, w and p alike.
Conventions follow cff.py: boundary faces come first, the right-hand-rule area
A points into c0, and a face's outward area for c0 is -A.
"""
import numpy as np
import scipy.sparse as sp


class Ctx:
    """Mesh connectivity and geometry in the shapes the operators need."""

    def __init__(self, m):
        self.m = m
        self.n = m["n_cells"]
        self.c0, self.c1 = m["c0"], m["c1"]
        inn = self.c1 >= 0
        self.nbf = int(np.nonzero(~inn)[0].max()) + 1
        assert not inn[: self.nbf].any() and inn[self.nbf :].all(), "boundary faces must come first"
        self.I = np.nonzero(inn)[0]
        self.B = np.arange(self.nbf)
        self.N = self.n + self.nbf
        self.cen, self.fc, self.A, self.vol = m["cen"], m["fcent"], m["farea"], m["vol"]
        self.area = np.linalg.norm(self.A, axis=1)
        self.nhat = self.A / self.area[:, None]
        self.zone_type = np.full(self.c0.size, 2)
        for z in m["zones"]:
            self.zone_type[z["lo"] : z["hi"]] = z["type"]
        self.is_wall_face = self.zone_type[: self.nbf] == 3
        # Owner/neighbour vectors to the face centroid, interior faces only.
        self.r0 = self.fc[self.I] - self.cen[self.c0[self.I]]
        self.r1 = self.fc[self.I] - self.cen[self.c1[self.I]]
        self.d = self.cen[self.c1[self.I]] - self.cen[self.c0[self.I]]

    def unit_rows(self, cols):
        return sp.csr_matrix((np.ones(cols.size), (np.arange(cols.size), cols)), shape=(cols.size, self.N))

    def apply(self, G, phi, phib):
        """[n, 3(, m)]: d phi_j / d x_i for vector phi."""
        phi_all = np.concatenate([np.asarray(phi), np.asarray(phib)])
        return (G @ phi_all).reshape((self.n, 3) + phi_all.shape[1:])


def lsq_coef(ctx, cell, dx, w):
    """Per-row weighted least-squares coefficient vectors M_c^-1 w dx."""
    M = np.zeros((ctx.n, 3, 3))
    np.add.at(M, cell, w[:, None, None] * dx[:, :, None] * dx[:, None, :])
    return np.einsum("rij,rj->ri", np.linalg.inv(M)[cell], w[:, None] * dx)


def rows_to_operator(ctx, cell, coef, V):
    """G = sum_rows coef_r (V_r - e_cell)."""
    R = cell.size
    E = sp.csr_matrix((np.ones(R), (np.arange(R), cell)), shape=(R, ctx.N))
    rows = np.concatenate([3 * cell + i for i in range(3)])
    cols = np.concatenate([np.arange(R)] * 3)
    Cm = sp.csr_matrix((coef.T.ravel(), (rows, cols)), shape=(3 * ctx.n, R))
    return (Cm @ (V - E)).tocsr()


def lsq_rows(ctx, wall_rows=True, normal_types=(), omit_types=()):
    """Face-neighbour LSQ rows (Fluent TG 23.3.3: centroid-to-centroid differences).

    Boundary faces act as data points at the face centroid with the boundary value.
    wall_rows=False drops wall faces from the stencil (tested in §1.4 of the report).
    """
    c0, c1, I, B = ctx.c0, ctx.c1, ctx.I, ctx.B
    keep_b = np.ones(ctx.nbf, bool) if wall_rows else ~ctx.is_wall_face
    keep_b &= ~np.isin(ctx.zone_type[B], omit_types)
    Bk = B[keep_b]
    cell = np.concatenate([c0[I], c1[I], c0[Bk]])
    dx = np.concatenate([ctx.d, -ctx.d, ctx.fc[Bk] - ctx.cen[c0[Bk]]])
    normal = np.isin(ctx.zone_type[Bk], normal_types)
    rb = dx[2 * I.size:]
    rb[normal] = np.einsum("ij,ij->i", rb[normal], ctx.nhat[Bk[normal]])[:, None] * ctx.nhat[Bk[normal]]
    V = sp.vstack([ctx.unit_rows(c1[I]), ctx.unit_rows(c0[I]), ctx.unit_rows(ctx.n + Bk)]).tocsr()
    return cell, dx, V


def lsq(ctx, power=2.0, wall_rows=True, normal_types=(), omit_types=()):
    """LSQ cell gradient with weights 1/|dx|^power (Fluent lscb-wf: 0, 1, 2 = 1/d^0, 1/d, 1/d^2)."""
    cell, dx, V = lsq_rows(ctx, wall_rows, normal_types, omit_types)
    return rows_to_operator(ctx, cell, lsq_coef(ctx, cell, dx, np.linalg.norm(dx, axis=1) ** -float(power)), V)


def green_gauss(ctx, Fv):
    """Green-Gauss gradient (1/V) sum phi_f A_out from a face-value operator Fv (nfaces x N)."""
    c0, c1, I, A, n = ctx.c0, ctx.c1, ctx.I, ctx.A, ctx.n
    rows, cols, vals = [], [], []
    for i in range(3):
        rows += [3 * c0 + i, 3 * c1[I] + i]
        cols += [np.arange(c0.size), I]
        vals += [-A[:, i] / ctx.vol[c0], A[I, i] / ctx.vol[c1[I]]]
    Sm = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(3 * n, c0.size))
    return (Sm @ Fv).tocsr()


def arithmetic_faces(ctx):
    """Interior faces: (phi0 + phi1)/2; boundary faces: their own value."""
    I, B = ctx.I, ctx.B
    rows = np.concatenate([I, I, B])
    cols = np.concatenate([ctx.c0[I], ctx.c1[I], ctx.n + B])
    vals = np.concatenate([np.full(I.size, 0.5), np.full(I.size, 0.5), np.ones(ctx.nbf)])
    return sp.csr_matrix((vals, (rows, cols)), shape=(ctx.c0.size, ctx.N))


def _dot_r(ctx, cells, r, rowidx, nrows, G):
    """Rows rowidx: r . (G phi)_cells."""
    out = None
    for i in range(3):
        R = sp.csr_matrix((r[:, i], (rowidx, 3 * cells + i)), shape=(nrows, 3 * ctx.n))
        out = R if out is None else out + R
    return (out @ G).tocsr()


def second_order_faces(ctx, G):
    """Fluent 'Second Order' face value: mean of the two side reconstructions phi_c + G_c . r_c."""
    I = ctx.I
    nf = ctx.c0.size
    return (arithmetic_faces(ctx) + 0.5 * _dot_r(ctx, ctx.c0[I], ctx.r0, I, nf, G) + 0.5 * _dot_r(ctx, ctx.c1[I], ctx.r1, I, nf, G)).tocsr()


def jump_correction(ctx, G, power=2.0):
    """C(G) = per-cell LSQ fit, rows from centroid to face centroid weighted 1/|r_f|^power, of the
    reconstruction jump J_f = (phi_nb + G_nb r_nb) - (phi_c + G_c r_c); boundary: phi_b - (phi_c + G_c r_b)."""
    c0, c1, I, B, n, nbf = ctx.c0, ctx.c1, ctx.I, ctx.B, ctx.n, ctx.nbf
    R = I.size
    rows = np.arange(R)
    own0 = sp.csr_matrix((np.ones(R), (rows, c0[I])), shape=(R, ctx.N)) + _dot_r(ctx, c0[I], ctx.r0, rows, R, G)
    own1 = sp.csr_matrix((np.ones(R), (rows, c1[I])), shape=(R, ctx.N)) + _dot_r(ctx, c1[I], ctx.r1, rows, R, G)
    J0 = own1 - own0
    rb = ctx.fc[B] - ctx.cen[c0[B]]
    Jb = ctx.unit_rows(n + B) - (sp.csr_matrix((np.ones(nbf), (B, c0[B])), shape=(nbf, ctx.N)) + _dot_r(ctx, c0[B], rb, B, nbf, G))
    cell = np.concatenate([c0[I], c1[I], c0[B]])
    rvec = np.concatenate([ctx.r0, ctx.r1, rb])
    coef = lsq_coef(ctx, cell, rvec, np.linalg.norm(rvec, axis=1) ** -float(power))
    Jm = sp.vstack([J0, -J0, Jb]).tocsr()
    rr = np.concatenate([3 * cell + i for i in range(3)])
    cc = np.concatenate([np.arange(cell.size)] * 3)
    Cm = sp.csr_matrix((coef.T.ravel(), (rr, cc)), shape=(3 * n, cell.size))
    return (Cm @ Jm).tocsr()


def barth_jespersen(ctx, phi, phib, g, coefficient=None):
    """Standard (Barth-Jespersen) limiter, cell-to-face, scalar per component.

    min/max over the cell and its face neighbours (boundary faces contribute their face value).
    g: [n, 3] or [n, 3, m]; returns the limited gradient with the same shape.
    coefficient: Kutta's differentiable limiter instead (fv_gradient.kern.h klimit_face):
    Venkatakrishnan's function with threshold coefficient * (global range of each component).
    """
    phi = np.asarray(phi)
    scalar = phi.ndim == 1
    if scalar:
        phi, phib, g = phi[:, None], np.asarray(phib)[:, None], g[:, :, None]
    c0, c1, I, B, n = ctx.c0, ctx.c1, ctx.I, ctx.B, ctx.n
    allv = np.concatenate([phi, phib])
    cell = np.concatenate([c0[I], c1[I], c0[B]])
    oth = np.concatenate([c1[I], c0[I], n + B])
    faces = np.concatenate([I, I, B])
    lo, hi = phi.copy(), phi.copy()
    np.minimum.at(lo, cell, allv[oth])
    np.maximum.at(hi, cell, allv[oth])
    dd = np.einsum("ri,rim->rm", ctx.fc[faces] - ctx.cen[cell], g[cell])
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(dd > 0, (hi[cell] - phi[cell]) / dd, np.where(dd < 0, (lo[cell] - phi[cell]) / dd, np.inf))
        if coefficient is not None:
            eps2 = (coefficient * (phi.max(0) - phi.min(0))) ** 2
            room = np.where(dd > 0, hi[cell] - phi[cell], lo[cell] - phi[cell])
            a, b = room**2 + eps2, dd * room
            ratio = np.where(dd != 0, (a + 2 * b) / (a + 2 * dd**2 + b), np.inf)
    psi = np.ones_like(phi)
    np.minimum.at(psi, cell, np.minimum(ratio, 1.0))
    out = g * np.clip(psi, 0, 1)[:, None, :]
    return out[:, :, 0] if scalar else out


def hops_from_boundary(ctx, maxh=8):
    """Face-neighbour distance (in cells) from any boundary-adjacent cell."""
    h = np.full(ctx.n, 99)
    h[ctx.c0[ctx.B]] = 0
    a, b = ctx.c0[ctx.I], ctx.c1[ctx.I]
    for k in range(1, maxh + 1):
        prev = h == k - 1
        nb = np.zeros(ctx.n, bool)
        nb[b[prev[a]]] = True
        nb[a[prev[b]]] = True
        h[nb & (h == 99)] = k
    return h


def strain_rate(gradU):
    """gradU[n, i, j] = d u_j / d x_i -> S = sqrt(2 S_ij S_ij)."""
    L = np.transpose(gradU, (0, 2, 1))
    Sij = 0.5 * (L + np.transpose(L, (0, 2, 1)))
    return np.sqrt(2.0 * np.einsum("nij,nij->n", Sij, Sij))
