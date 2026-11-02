"""Independent standard k-epsilon budget at Fluent's saved fields.

Usage: python -I ke_budget.py case_dir stem cache_dir diagnostic_dir
diagnostic_dir is bench_elbow_turbulence's `saved 0` output. Uses the actual
generated model's exported production, sinks and cell diffusivities, Fluent's
saved face flux, first-order upwind convection and compact nonorthogonal
diffusion. The last is an independent transport check, NOT the solver's walk
reconstruction; use its monitor for the actual assembled defect. Normalization
is physical destruction, not Fluent's diagonal normalization. No fit is made.
Needs numpy, scipy and h5py, as do the other forensics scripts.
"""
import os
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import cff
import fvops
import momentum as mo


def coefficients(ctx, cell, dx, weight):
    """Same deficient-direction convention as boundary_least_squares()."""
    normal = np.zeros((ctx.n, 3, 3))
    np.add.at(normal, cell, weight[:, None, None] * dx[:, :, None] * dx[:, None, :])
    values, vectors = np.linalg.eigh(normal)
    inverse = np.einsum("nik,nk,njk->nij", vectors, 1 / np.where(values > 1e-10, values, 1), vectors)
    return np.einsum("rij,rj->ri", inverse[cell], weight[:, None] * dx)


def main():
    case_dir, stem, cache, directory = sys.argv[1:]
    mesh = cff.load(case_dir, stem, cache)
    ctx = fvops.Ctx(mesh)
    cell = lambda name: cff.cell_field(mesh["cells"], name, ctx.n)
    face = lambda name: cff.cell_field(mesh["faces"], name, ctx.nbf)
    k, epsilon, mu, mut = (cell(name) for name in ("SV_K", "SV_D", "SV_MU_LAM", "SV_MU_T"))
    u = np.stack([cell("SV_" + name) for name in "UVW"], 1)
    ub = np.stack([face("SV_" + name) for name in "UVW"], 1)
    state = SimpleNamespace(
        ctx=ctx, F=cff.cell_field(mesh["faces"], "SV_FLUX", ctx.c0.size),
        zone=lambda kind: np.where(ctx.zone_type[:ctx.nbf] == kind)[0],
        alpha=ctx.area[ctx.I] ** 2 / np.einsum("ij,ij->i", ctx.A[ctx.I], -ctx.d),
        dn=np.einsum("ij,ij->i", ctx.fc[:ctx.nbf] - ctx.cen[ctx.c0[:ctx.nbf]], ctx.nhat[:ctx.nbf]),
    )
    model = np.genfromtxt(os.path.join(directory, "cells.csv"), names=True, delimiter=",")
    assert len(model) == ctx.n
    assert np.allclose(model["k"], k, rtol=2e-7, atol=1e-14), "requires saved-state diagnostics"
    assert np.allclose(model["epsilon"], epsilon, rtol=2e-7, atol=1e-14)
    wall = np.zeros(ctx.n, bool)
    wall[ctx.c0[state.zone(3)]] = True
    # The base library's inv() cannot represent rank-deficient boundary fits.
    fvops.lsq_coef = coefficients
    for name, q, qb, normal, omitted, source, sink, gamma in (
        ("k", k, face("SV_K"), (3,), (5, 7), "pk", "dk", "gk"),
        ("epsilon", epsilon, face("SV_D"), (), (3, 5, 7), "ps", "ds", "gs"),
    ):
        gradient = ctx.apply(fvops.lsq(ctx, 2, normal_types=normal, omit_types=omitted), q, qb)
        convection = mo.scalar_convection(state, q, qb, np.zeros((ctx.n, 3)))
        diffusion = mo.scalar_diffusion(state, q, qb, gradient, model[gamma])
        production = model[source] * ctx.vol
        destruction = model[sink] * q * ctx.vol
        residual = convection + diffusion - production + destruction
        for region, selected in (("all", np.ones(ctx.n, bool)), ("wall", wall.copy()), ("interior", ~wall)):
            if name == "epsilon":
                selected &= ~wall  # wall rows are fixed, not transport equations
            if not selected.any():
                continue
            scale = destruction[selected].sum()
            values = [a[selected].sum() / scale for a in (convection, diffusion, production, destruction, residual, np.abs(residual))]
            print(name, region, "C D P Y R absR / Y:", " ".join(f"{value:+.9g}" for value in values))
        if name == "k":
            for power in (0, 1, 2):
                base = fvops.lsq(ctx, power, wall_rows=False)
                g = ctx.apply(base, u, ub)
                correction = ctx.apply(fvops.jump_correction(ctx, base), u, ub)
                for coefficient in (0, 0.5, 1):
                    strain = fvops.strain_rate(g + coefficient * correction)
                    candidate = (mut + mu) * strain ** 2 * ctx.vol
                    candidate[wall] = production[wall]
                    r = convection + diffusion - candidate + destruction
                    scale = destruction[~wall].sum()
                    print(f"strain weight_power={power} jump={coefficient}: interior R/Y={r[~wall].sum()/scale:+.9g} absR/Y={np.abs(r[~wall]).sum()/scale:.9g}")


if __name__ == "__main__":
    main()
