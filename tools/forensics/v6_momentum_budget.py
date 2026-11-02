"""Isolate V6 momentum terms at Fluent's saved fields; never advance a field.

Uses the saved conservative SV_FLUX, SV_MU_T and wall shear unless a named
substitution is requested. The Python assembly approximates Kutta's LSQ/SOU
operator; it does not include the boundary extrapolation-gain guard. Its L1
values must not be labelled as SteadySolver monitor values.

Needs numpy, scipy, h5py. Example:
  python -I tools/forensics/v6_momentum_budget.py CASE_DIR --cache CACHE \
    --kutta-flux FLUX.csv --kutta-sst metal-saved-sst.csv
FLUX.csv comes from bench_v6_operators ... flux=FLUX.csv; the SST table comes
from KUTTA_V6_WALL_DIR. These dumps contain reference data; do not commit them.
"""
import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import cff
import fvops
import momentum as mo


def validate_inputs(st):
    for name in ('U', 'Ub', 'p', 'pb', 'rho', 'mu', 'mut', 'k', 'kb', 'F'):
        if not np.isfinite(getattr(st, name)).all():
            raise ValueError('missing or nonfinite saved momentum input: ' + name)
    # State zero-fills non-wall shear entries, which are legitimately absent.
    # Check the original wall rows so a missing wall value cannot become zero.
    shear = cff.cell_field(st.m['faces'], 'SV_WALL_SHEAR', st.ctx.nbf)
    if not np.isfinite(shear[st.zone(3)]).all():
        raise ValueError('missing or nonfinite saved wall shear')


def flux_by_identity(st, path, hashes):
    """Match input point IDs and both adjacent cells, including duplicate pairs.

    V6's CFF node IDs are dense, starting at 1; both readers retain their rank.
    A sorted node set identifies a face; geometry/nearest-neighbour matching is
    deliberately unnecessary. Reject missing, repeated or inconsistent faces.
    """
    ctx, mesh = st.ctx, st.m
    ptr = np.r_[0, np.cumsum(mesh['nnodes'])]
    identities = {}
    for f in range(ctx.c0.size):
        key = tuple(sorted(mesh['fnodes'][ptr[f]:ptr[f + 1]]))
        if key in identities:
            raise ValueError('source mesh has duplicate face node sets')
        identities[key] = f
    result = np.full(ctx.c0.size, np.nan)
    with open(path, newline='') as stream:
        for label, suffix in (('case', 'cas.h5'), ('data', 'dat.h5')):
            if stream.readline().strip() != f'# {label}_sha256={hashes[suffix]}':
                raise ValueError('Kutta flux dump belongs to different case/data inputs')
        for row in csv.DictReader(stream):
            key = tuple(sorted(int(n) for n in row['nodes'].split()))
            f = identities[key]
            if np.isfinite(result[f]):
                raise ValueError('Kutta dump repeats a face')
            owner, neighbour = int(row['owner']), int(row['neighbour'])
            expected = (int(ctx.c0[f]), int(ctx.c1[f]))
            value = float(row['flux'])
            if not np.isfinite(value):
                raise ValueError('nonfinite Kutta flux')
            if (owner, neighbour) == expected:
                result[f] = value
            elif expected[1] >= 0 and (neighbour, owner) == expected:
                result[f] = -value
            else:
                raise ValueError('Kutta face adjacency differs from Fluent')
    if not np.isfinite(result).all():
        raise ValueError('Kutta dump does not cover every source face')
    return result


def face_moments(st):
    """Triangle-fan first area-vector moment about the stored face centroid."""
    mesh, ctx = st.m, st.ctx
    nn, fn, nodes = mesh['nnodes'], mesh['fnodes'], mesh['nodes']
    ptr = np.r_[0, np.cumsum(nn)]
    face = np.repeat(np.arange(nn.size), nn)
    centre = np.zeros_like(ctx.fc)
    np.add.at(centre, face, nodes[fn])
    centre /= nn[:, None]
    nxt = np.arange(fn.size) + 1
    nxt[ptr[1:] - 1] = ptr[:-1]
    a, b, c = nodes[fn], nodes[fn[nxt]], centre[face]
    area = .5 * np.cross(a - c, b - c)
    centroid = (a + b + c) / 3
    moment = np.zeros((nn.size, 3, 3))
    np.add.at(moment, face, (centroid - ctx.fc[face])[:, :, None] * area[:, None, :])
    return moment, (face, area, a, b, c)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case_dir', type=Path)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--kutta-flux', type=Path)
    parser.add_argument('--kutta-sst', type=Path)
    args = parser.parse_args()
    hashes = {suffix: hashlib.sha256((args.case_dir / ('mixing_elbow.' + suffix)).read_bytes()).hexdigest()
              for suffix in ('cas.h5', 'dat.h5')}
    cache = args.cache / hashes['cas.h5']
    cache.mkdir(parents=True, exist_ok=True)
    st = mo.State(str(args.case_dir), 'mixing_elbow', str(cache))
    validate_inputs(st)
    ctx = st.ctx
    I, B, p, n = ctx.I, ctx.B, ctx.c0, ctx.c1
    gu = ctx.apply(fvops.lsq(ctx, 2, normal_types=(7,), omit_types=(5,)), st.U, st.Ub)
    gp = ctx.apply(fvops.lsq(ctx, 2, omit_types=(3, 7, 10)), st.p, st.pb)
    gl = mo.bj_limit(ctx, st.U, st.Ub, gu)[0]
    convection = mo.convection(st, gl)
    pressure = mo.pressure(st, gp, boundary='extrap')[0]
    def viscous():
        return (mo.interior_viscous(st, gu, gu, mu_mode='dist', grad_mode='dist')
                + mo.boundary_viscous(st, gu, sym='none', outlet='transpose'))
    diffusion = viscous()
    residual = convection + pressure + diffusion
    deep = fvops.hops_from_boundary(ctx) >= 3
    rows = {}
    def record(name, value):
        rows[name] = dict(raw_l1_xyz=np.abs(value).sum(0).tolist(),
                          raw_l1_total=float(np.abs(value).sum()),
                          deep_l1_total=float(np.abs(value[deep]).sum()))
    record('saved_flux_viscosity_shear', residual)

    gk = ctx.apply(fvops.lsq(ctx, 2, omit_types=(3, 5), normal_types=(7,)), st.k, st.kb)
    kf = .5 * (st.k[p[I]] + np.einsum('ij,ij->i', gk[p[I]], ctx.r0)
               + st.k[n[I]] + np.einsum('ij,ij->i', gk[n[I]], ctx.r1))
    term = (2 / 3 * st.rho[p[I]] * kf)[:, None] * ctx.A[I]
    isotropic = mo.scatter(ctx, I, -term, term)
    np.add.at(isotropic, p[B], -(2 / 3 * st.rho[p[B]] * st.kb)[:, None] * ctx.A[B])
    record('isotropic_term', isotropic)
    record('add_isotropic', residual + isotropic)
    record('subtract_isotropic', residual - isotropic)
    for mode in ('arith', 'harm'):
        value = mo.interior_viscous(st, gu, gu, mu_mode=mode, grad_mode='dist')
        record('viscosity_' + mode, residual - mo.interior_viscous(st, gu, gu, mu_mode='dist', grad_mode='dist') + value)

    moment, (face, area, a, b, c) = face_moments(st)
    value = -np.einsum('fi,fij->fj', .5 * (gp[p[I]] + gp[n[I]]), moment[I])
    warp = mo.scatter(ctx, I, value, -value)
    np.add.at(warp, p[B], -np.einsum('fi,fij->fj', gp[p[B]], moment[B]))
    record('warped_pressure_moment_term', warp)
    record('add_warped_pressure_moment', residual + warp)

    # Integrate the product of affine mass-flux and transported-velocity
    # variations over the fan, holding the saved face-integrated flux fixed.
    idx = np.nonzero(face >= ctx.nbf)[0]
    ff = face[idx]
    up = np.where(st.F[ff] >= 0, p[ff], n[ff])
    gup, gf = gl[up], .5 * (gu[p[ff]] + gu[n[ff]])
    du = [np.einsum('fi,fij->fj', v[idx] - ctx.fc[ff], gup) for v in (a, b, c)]
    dv = [np.einsum('fi,fij->fj', v[idx] - ctx.fc[ff], gf) for v in (a, b, c)]
    cov = (sum(du) * np.einsum('fi,fi->f', sum(dv), -area[idx])[:, None]
           + sum(x * np.einsum('fi,fi->f', y, -area[idx])[:, None] for x, y in zip(du, dv))) / 12
    cov *= st.rho[up, None]
    value = np.zeros_like(ctx.A)
    np.add.at(value, ff, cov)
    covariance = mo.scatter(ctx, I, value[I], -value[I])
    record('affine_flux_covariance_term', covariance)
    record('add_affine_flux_covariance', residual + covariance)

    if args.kutta_sst:
        data = np.genfromtxt(args.kutta_sst, delimiter=',', names=True)
        if data.size != ctx.n or not np.isfinite(data['mu_t']).all() or (data['mu_t'] < 0).any():
            raise ValueError('Kutta mu_t table must contain one finite nonnegative value per source cell')
        saved_mut = st.mut
        st.mut = data['mu_t']
        record('kutta_mut_saved_flux_shear', convection + pressure + viscous())
        st.mut = saved_mut
    continuity = {}
    def mass_l1(flux):
        div = np.zeros(ctx.n)
        np.add.at(div, p, flux)
        np.add.at(div, n[I], -flux[I])
        return float(np.abs(div).sum())
    continuity['saved_flux'] = mass_l1(st.F)
    if args.kutta_flux:
        saved_flux = st.F
        st.F = flux_by_identity(st, args.kutta_flux, hashes)
        continuity['kutta_flux'] = mass_l1(st.F)
        record('kutta_flux_saved_viscosity_shear', mo.convection(st, gl) + pressure + diffusion)
        st.F = saved_flux
    print(json.dumps(dict(cells=ctx.n, faces=int(p.size), input_sha256=hashes,
                         fluent_raw_l1_xyz=[st.res[k][0] for k in ('x-velocity', 'y-velocity', 'z-velocity')],
                         isotropic_least_squares_coefficient=float(-np.sum(residual * isotropic) / np.sum(isotropic * isotropic)),
                         continuity_l1_kg_s=continuity, momentum=rows), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
