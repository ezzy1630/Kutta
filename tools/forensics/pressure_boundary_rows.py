"""Pressure LSQ boundary rows: Fluent's saved momentum residual and wall-cell conditioning (validation/results/tet_robustness).

Sections:
  1. Momentum residual of Fluent's saved state (SST cases with SV_MU_T, e.g. V6) under pressure-gradient boundary
     treatments, with momentum's boundary pressure extrapolated as p_c + G.r_b. Skipped for other models.
  2. Per-cell conditioning of the 1/d^2 pressure fit (eigenvalue ratio lambda_min/lambda_max) and the gain of its
     extrapolation to each wall face, sum_j |r_b . M^-1 w_j d_j|, with and without boundary rows; and, per boundary
     type, the gain from main's outlet-only fit that decides Kutta's second-order boundary rows (gain > 1).
Usage: python -I pressure_boundary_rows.py [case_dir] [stem] [cache_dir]   (see cff.py)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import cff
import fvops

WALL, SYMMETRY, INLET = 3, 7, 10


def residuals(case_dir, stem, cache):
    import momentum as mo

    st = mo.State(case_dir, stem, cache)
    ctx = st.ctx
    U, Ub, p, pb = st.U, st.Ub, st.p, st.pb
    tb = mo.boundary_cells(st)
    regions = {"wall cells": (tb & 1) > 0, "sym cells": (tb & 2) > 0, "inlet cells": (tb & 4) > 0, "all": np.ones(ctx.n, bool)}
    gu = ctx.apply(mo.lsq_operator(ctx, 2), U, Ub)
    rest = mo.convection(st, mo.bj_limit(ctx, U, Ub, gu)[0]) + mo.interior_viscous(st, gu, gu) + mo.boundary_viscous(st, gu)
    print("1. sum|R| (x+y+z) of Fluent's saved state; boundary p in momentum = p_c + G.r_b; saved p_b = p_c at walls")
    print(f"  {'pressure LSQ boundary rows':52s}" + "".join(f"{k:>12s}" for k in regions))
    for label, G in (
        ("every boundary face, p_b at the face centroid", fvops.lsq(ctx, 2)),
        ("walls at the normal foot, others at the centroid", fvops.lsq(ctx, 2, normal_types=(WALL,))),
        ("walls, symmetry, inlets at the normal foot", fvops.lsq(ctx, 2, normal_types=(WALL, SYMMETRY, INLET))),
        ("walls and symmetry at the centroid, no inlet rows", fvops.lsq(ctx, 2, omit_types=(INLET,))),
        ("no wall rows", fvops.lsq(ctx, 2, wall_rows=False)),
        ("outlet rows only (Kutta second-order before)", fvops.lsq(ctx, 2, wall_rows=False, omit_types=(SYMMETRY, INLET))),
    ):
        a = np.abs(rest + mo.pressure(st, ctx.apply(G, p, pb), boundary="extrap")[0]).sum(1)
        print(f"  {label:52s}" + "".join(f"{a[m].sum():12.3e}" for m in regions.values()))


def conditioning(m):
    ctx = fvops.Ctx(m)
    n = ctx.n
    wall = np.nonzero(ctx.is_wall_face)[0]
    print("\n2. 1/d^2 pressure fit per cell: eigenvalue ratio, and extrapolation gain to its wall faces")
    for label, wall_rows, omit in (("outlet rows only", False, (SYMMETRY, INLET)), ("every boundary face", True, ())):
        cell, dx, _ = fvops.lsq_rows(ctx, wall_rows=wall_rows, omit_types=omit)
        w = np.linalg.norm(dx, axis=1) ** -2.0
        M = np.zeros((n, 3, 3))
        np.add.at(M, cell, w[:, None, None] * dx[:, :, None] * dx[:, None, :])
        ev = np.linalg.eigvalsh(M)
        ratio = ev[:, 0] / ev[:, 2]
        coef = np.einsum("rij,rj->ri", np.linalg.pinv(M)[cell], w[:, None] * dx)
        order = np.argsort(cell, kind="stable")
        starts = np.searchsorted(cell[order], np.arange(n + 1))
        gain = np.zeros(n)
        for f in wall:
            c = ctx.c0[f]
            rows = order[starts[c] : starts[c + 1]]
            gain[c] = max(gain[c], np.abs(coef[rows] @ (ctx.fc[f] - ctx.cen[c])).sum())
        g = gain[gain > 0]
        print(f"  {label}: ratio min {ratio.min():.2e}, cells < 1e-5 {(ratio < 1e-5).sum()}, < 1e-3 {(ratio < 1e-3).sum()},"
              f" < 1e-2 {(ratio < 1e-2).sum()}; wall-cell gain median {np.median(g):.2f}, p99 {np.percentile(g, 99):.2f},"
              f" max {g.max():.1f}, > 2: {(g > 2).sum()}, > 5: {(g > 5).sum()} of {g.size}")
    # Kutta's second-order rule (boundary_extrapolation_gain): from main's fit (outlet rows only), an unprescribed face
    # holds the cell pressure where the gain of the extrapolation to it exceeds 1.
    cell, dx, _ = fvops.lsq_rows(ctx, wall_rows=False, omit_types=(SYMMETRY, INLET))
    w = np.linalg.norm(dx, axis=1) ** -2.0
    M = np.zeros((n, 3, 3))
    np.add.at(M, cell, w[:, None, None] * dx[:, :, None] * dx[:, None, :])
    coef = np.einsum("rij,rj->ri", np.linalg.pinv(M)[cell], w[:, None] * dx)
    order = np.argsort(cell, kind="stable")
    starts = np.searchsorted(cell[order], np.arange(n + 1))
    for label, kind in (("wall", WALL), ("symmetry", SYMMETRY), ("inlet", INLET)):
        faces = np.nonzero(ctx.zone_type[: ctx.nbf] == kind)[0]
        g = np.array([np.abs(coef[order[starts[ctx.c0[f]] : starts[ctx.c0[f] + 1]]] @ (ctx.fc[f] - ctx.cen[ctx.c0[f]])).sum()
                      for f in faces])
        if g.size:
            print(f"  {label} faces, outlet-only fit: gain median {np.median(g):.2f}, max {g.max():.2f},"
                  f" > 0.75: {(g > .75).sum()}, > 1 (rows kept): {(g > 1).sum()}, > 1.5: {(g > 1.5).sum()} of {g.size}")


def main():
    case_dir, stem, cache = cff.args()
    m = cff.load(case_dir, stem, cache)
    if "SV_MU_T" in m["cells"] and "SV_O" in m["cells"]:
        residuals(case_dir, stem, cache)
    else:
        print("1. skipped: not an SST data file")
    conditioning(m)


if __name__ == "__main__":
    main()
