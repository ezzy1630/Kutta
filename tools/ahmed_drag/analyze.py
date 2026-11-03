#!/usr/bin/env python3
"""Aggregate Ahmed diagnostics without exporting restricted reference curves.

Reads each run's result.json (run.py format), walls.csv and optional cells.csv.
Writes only integrated loads and aggregate errors, never reference samples.
Use surface_size, not cell-count^(−1/3), for these locally refined meshes.
"""
import argparse
import json
import math
import pathlib
import re
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.spatial import cKDTree

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools/fsae_aero'))
import report  # noqa: E402


def richardson(spacing, drag):
    """Three-level fit q(h)=q0+C*h**p; reject nonmonotonic/no-positive-order fits.

    A successful algebraic fit does not establish asymptotic convergence.
    Compare adjacent triplets and iterative uncertainty before using q0.
    """
    h, q = np.asarray(spacing, dtype=float), np.asarray(drag, dtype=float)
    if h.shape != (3,) or q.shape != (3,) or not np.all(np.isfinite([h, q])):
        raise ValueError('three finite spacings and drags required')
    if not np.all(h[:-1] > h[1:]) or h[-1] <= 0:
        raise ValueError('spacings must be positive and strictly decreasing')
    d0, d1 = q[0] - q[1], q[1] - q[2]
    if d0 * d1 <= 0:
        return {'accepted': False, 'reason': 'nonmonotonic or zero increments'}
    # expm1 avoids cancellation as p approaches zero.
    a, b = h[0] / h[1], h[1] / h[2]
    def equation(p):
        return b**p * math.expm1(p * math.log(a)) / math.expm1(p * math.log(b)) - d0 / d1
    try:
        p = brentq(equation, 1e-6, 20)
    except ValueError:
        return {'accepted': False, 'reason': 'no positive order in (0,20)'}
    q0 = q[-1] - d1 / math.expm1(p * math.log(b))
    return {'accepted': True, 'order': p, 'extrapolated_cd': q0,
            'fine_correction': q0 - q[-1], 'asymptotic_verified': False}


def table(path):
    """ERCOFTAC ASCII; coordinate VALUES are mm, including BL files labelled m."""
    rows, names = [], []
    for line in path.read_text(encoding='latin-1').splitlines():
        column = re.match(r'%Column\d+Name\s*=\s*(\S+)', line)
        if column:
            names.append(column[1].lower())
        if line.startswith(('#', '%')) or not line.strip():
            continue
        try:
            rows.append([float(v) for v in line.split()])
        except ValueError as e:
            raise ValueError(f'non-numeric data in {path.name}') from e
    data = np.asarray(rows)
    if data.ndim != 2 or data.shape[1] != len(names):
        raise ValueError(f'column metadata does not match data in {path.name}')
    return data, names


def velocity_errors(cells, slant):
    xyz = np.column_stack([cells[k] for k in ('x', 'y', 'z')])
    tree = cKDTree(xyz)
    out = {}
    # Nearest-cell sampling is explicit, bounded and reproducible. These are
    # sampling diagnostics, not interpolated wall values or acceptance gates.
    for source, prefix in [('experiment', 'exp_lienhart_becker/Ahmed_'),
                           ('CFX-5 SST', 'cfd11_cfx5_sst/CFD11_SST_Ahmed_')]:
        for station in ['BL_x=-163', 'BL_x=-123', 'BL_x=-83', 'BL_x=-43', 'BL_x=-3',
                        'x=80', 'x=200', 'x=500']:
            raw, names = table(report.DATA / 'ercoftac_ahmed_body' / f'{prefix}{slant}_{station}.dat')
            points = raw[:, :3] / 1000
            distance, idx = tree.query(points)
            is_bl = station.startswith('BL')
            keep = distance < (.006 if is_bl else .03)
            cols = {name: names.index(name) for name in ('u', 'w')}
            row = {'points': len(raw), 'matched': int(keep.sum()),
                   'max_distance_m': float(distance[keep].max()) if keep.any() else None}
            for name, col in cols.items():
                diff = cells[name][idx[keep]] - raw[keep, col] / 40
                row[f'{name}_rms_over_U'] = float(np.sqrt(np.mean(diff**2))) if keep.any() else None
                row[f'{name}_bias_over_U'] = float(diff.mean()) if keep.any() else None
            out[f'{source} {station}'] = row
    return out


def pressure_groups(xyz, slant, groups):
    """Split the slant at the same Cartesian distance used for force bands."""
    start = -.222 * math.cos(math.radians(slant))
    leading = xyz[:, 0] < start + .04
    return dict(groups, slant_leading_40mm=groups['slant'] & leading,
                slant_remaining=groups['slant'] & ~leading)


def upstream_turbulence(cells, viscosity):
    """Volume-weighted box ahead of the nose, outside ground/body layers.

    U_inlet = rho = 1 in this case. Tu is the mean of sqrt(2k/3)/U_inlet;
    the outlet-gauge Cp is recorded without fitting a reference-pressure shift.
    """
    keep = ((cells['x'] > -1.5) & (cells['x'] < -1.2) &
            (np.abs(cells['y']) < .15) & (cells['z'] > .35) & (cells['z'] < .55))
    values = {'tu_over_U_inlet': np.sqrt(2 * cells['k'] / 3),
              'k': cells['k'], 'omega': cells['omega'],
              'mu_t_over_mu': cells['mu_t'] / viscosity,
              'u': cells['u'], 'cp_outlet_gauge': 2 * cells['p']}
    return dict(cells=int(keep.sum()), **{
        name: float(np.average(value[keep], weights=cells['volume'][keep])) if keep.any() else None
        for name, value in values.items()})


def summarize(path):
    r = json.loads((path / 'result.json').read_text())
    if r.get('configuration_verified') is False:
        raise ValueError(f'unconfirmed benchmark configuration in {path.name}')
    r['dir'] = path
    res = r['result']
    if res['failure']:
        return {'run': path.name, 'failure': res['failure']}
    w = report.walls(path)
    command = r['command']
    area = float(command[command.index('--area') + 1])
    angle = float(r['mesh']['slant_deg'])
    if angle not in (25, 35) or float(r['point']) != angle:
        raise ValueError('reference comparisons require matching 25 or 35 degree mesh/run metadata')
    slant = int(angle)
    regions = report.regions(w, slant)
    regions['other'] = ~np.logical_or.reduce(list(regions.values()))
    # Fluid-outward normals and stress on the fluid: pressure n A, shear -tau A.
    cp_force = w['cp'] * w['nx'] * w['area'] / area
    cf_force = -w['tau_x'] * w['area'] / (.5 * area)
    components = {name: {'pressure': float(cp_force[m].sum()), 'friction': float(cf_force[m].sum()),
                         'area': float(w['area'][m].sum()),
                         'cp_area_mean': float(np.average(w['cp'][m], weights=w['area'][m])) if m.any() else None}
                  for name, m in regions.items()}
    start = -.222 * math.cos(math.radians(slant))
    bottom = .338 - .222 * math.sin(math.radians(slant))
    bands = {'slant_leading_40mm': regions['slant'] & (w['x'] < start + .04),
             'slant_remaining': regions['slant'] & (w['x'] >= start + .04),
             'base_perimeter_30mm': regions['base'] &
             ((np.abs(w['y']) > .389 / 2 - .03) | (w['z'] < .08) | (w['z'] > bottom - .03))}
    bands['base_centre'] = regions['base'] & ~bands['base_perimeter_30mm']
    edge_loads = {name: float(cp_force[mask].sum()) for name, mask in bands.items()}
    cp = {}
    for name, xyz, reference, kutta, distance, groups in report.cp_compare(w, slant):
        groups = pressure_groups(xyz, slant, groups)
        for region, mask in groups.items():
            keep = mask & (distance < .006)
            cp[f'{name} {region}'] = {'points': int(keep.sum()), 'eligible': int(mask.sum()),
                                     'reference_mean': float(reference[keep].mean()) if keep.any() else None,
                                     'kutta_mean': float(kutta[keep].mean()) if keep.any() else None,
                                     'rms': float(np.sqrt(np.mean((kutta[keep] - reference[keep])**2))) if keep.any() else None}
    h = np.genfromtxt(path / 'history.csv', delimiter=',', names=True)
    h = np.atleast_1d(h)
    tail = h[h['iteration'] > h['iteration'][-1] - 2000]
    row = {'run': path.name, 'source': r['source'], 'mesh': r['mesh'],
           'mesh_sha256': r['mesh_sha256'], 'backend': r['backend'], 'command': command,
           'peak_gb': r['peak_footprint_bytes'] / 1e9, 'iterations': res['iterations'],
           'residuals_met': res['residuals_met'], 'settled': res['settled'],
           'settle_lift': res['settle_lift'], 'settle_drag': res['settle_drag'],
           'max_residual': res['max_residual'], 'cd_final': res['drag'],
           'cd_mean_last2000': float(tail['drag'].mean()),
           'pressure_mean_last2000': float(tail['drag_pressure'].mean()),
           'friction_mean_last2000': float((tail['drag'] - tail['drag_pressure']).mean()),
           'cd_range_last2000': [float(tail['drag'].min()), float(tail['drag'].max())],
           'cd_last1000': float(tail['drag'][tail['iteration'] > h['iteration'][-1] - 1000].mean()),
           'pressure_final': res['drag_pressure'], 'friction_final': res['drag'] - res['drag_pressure'],
           'force_closure_pressure': float(cp_force.sum() - res['drag_pressure']),
           'force_closure_total': float(cp_force.sum() + cf_force.sum() - res['drag']),
           'components_final': components, 'edge_pressure_drag_final': edge_loads, 'cp': cp,
           'yplus_mean': res['yplus_mean'], 'yplus_max': res['yplus_max'],
           'topology': report.centreline_topology(w, slant),
           'reversed_slant_fraction': report.separated_fraction(w, slant)}
    for key in ('binary_sha256', 'pressure_scheme', 'configuration_verified'):
        if key in r:
            row[key] = r[key]
    if (path / 'cells.csv').exists():
        cells = np.genfromtxt(path / 'cells.csv', delimiter=',', names=True)
        viscosity = float(command[command.index('--length') + 1]) / float(command[command.index('--re') + 1])
        row['upstream_turbulence'] = upstream_turbulence(cells, viscosity)
        row['velocity'] = velocity_errors(cells, slant)
        if (path / 'owners.csv').exists():
            owners = np.genfromtxt(path / 'owners.csv', delimiter=',', names=True)
            coordinates = lambda data, names: np.column_stack([data[name] for name in names])
            face_distance, face_idx = cKDTree(coordinates(owners, ('x', 'y', 'z'))).query(
                coordinates(w, ('x', 'y', 'z')))
            owner_distance, cell_idx = cKDTree(coordinates(cells, ('x', 'y', 'z'))).query(
                coordinates(owners, ('owner_x', 'owner_y', 'owner_z'))[face_idx])
            # CSV fields have nine significant digits; this is a coordinate
            # identity check, not an adjustable CFD-accuracy threshold.
            if max(face_distance.max(), owner_distance.max()) > 1e-8:
                raise ValueError(f'wall/cell centroid identity mismatch in {path.name}')
            owner_force = 2 * cells['p'][cell_idx] * w['nx'] * w['area'] / area
            row['owner_pressure_diagnostic'] = {
                'max_coordinate_mismatch_m': float(max(face_distance.max(), owner_distance.max())),
                'components_using_owner_pressure': {name: float(owner_force[mask].sum())
                                                     for name, mask in regions.items()},
                'slant_leading_40mm_using_owner_pressure': float(owner_force[bands['slant_leading_40mm']].sum())}
    return row


def common_pressure_comparisons(rows, paths):
    """Use the same spatially supported taps across runs of each slant angle."""
    angles = {r['mesh']['slant_deg'] for r in rows if 'cp' in r}
    for angle in angles:
        population = [(r, p) for r, p in zip(rows, paths) if 'cp' in r and r['mesh']['slant_deg'] == angle]
        comparisons = [report.cp_compare(report.walls(p), int(angle)) for _, p in population]
        for source in range(2):
            groups = [pressure_groups(c[source][1], angle, c[source][5]) for c in comparisons]
            for region in ('slant', 'base', 'slant_leading_40mm', 'slant_remaining'):
                common = np.logical_and.reduce([g[region] & (c[source][4] < .006)
                                                for g, c in zip(groups, comparisons)])
                for (row, _), comparison in zip(population, comparisons):
                    name, _, reference, kutta, _, _ = comparison[source]
                    row.setdefault('cp_common', {})[f'{name} {region}'] = {
                        'points': int(common.sum()),
                        'reference_mean': float(reference[common].mean()) if common.any() else None,
                        'kutta_mean': float(kutta[common].mean()) if common.any() else None,
                        'rms': float(np.sqrt(np.mean((kutta[common] - reference[common])**2))) if common.any() else None}
                    row['cp_common_runs'] = [r['run'] for r, _ in population]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('runs', nargs='+', type=pathlib.Path)
    ap.add_argument('--out', required=True, type=pathlib.Path)
    args = ap.parse_args()
    rows = [summarize(p) for p in args.runs]
    common_pressure_comparisons(rows, args.runs)
    args.out.write_text(json.dumps(rows, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
