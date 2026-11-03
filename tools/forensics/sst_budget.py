"""Physical SST budgets from bench_sst_budget's native scalar assembly.

All output contains reference-derived values; do not commit it. No source term
is inferred from a zero residual: the reported 'needed production' also includes
any transport error and Fluent's finite convergence. Needs numpy, scipy, h5py.

python -I tools/forensics/sst_budget.py CASE_DIR --cache CACHE --dump-flux F.csv
build/bench_sst_budget CASE DATA OUT metal 0 saved F.csv
python -I tools/forensics/sst_budget.py CASE_DIR --cache CACHE --diagnostics OUT
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

import momentum as mo
import mut_limiter as ml
import v6_momentum_budget as mb


def digest(path):
    with open(path, 'rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def native_cells(directory, hashes, n):
    identity = (directory / 'identity.txt').read_text().splitlines()
    if len(identity) != 7 or identity[:2] != hashes:
        raise ValueError('diagnostics case/data identity mismatch')
    cells = np.genfromtxt(directory / 'cells.csv', delimiter=',', names=True)
    ids = cells['cell']
    if len(cells) != n or not np.array_equal(np.sort(ids), np.arange(n)):
        raise ValueError('diagnostics must cover every source cell exactly once')
    cells = cells[np.argsort(ids)]
    if not all(np.isfinite(cells[name]).all() for name in cells.dtype.names):
        raise ValueError('nonfinite diagnostic field')
    return cells, identity


def budget(cells, st, variant):
    """R=C+D-P+Y-X; physical Y excludes Newton's compensating source."""
    rho, vol = st.rho, st.ctx.vol
    f1, k, w = cells['f1'], cells['k'], cells['omega']
    beta = .075 * f1 + .0828 * (1 - f1)
    cross = (1 - f1) * cells['cdkw']
    cross = np.clip(cross, -1e6 * rho * beta * w**2, 1e6 * rho * beta * w**2)
    if variant == 'no-cross':
        cross = np.zeros_like(cross)
    yw = rho * beta * w**2
    # ps=Pw+max(X,0)+Yw, ds*w=2Yw+max(-X,0).
    pw = cells['ps'] - np.maximum(cross, 0) - yw
    result = {}
    for equation, convection, diffusion, production, destruction, x in (
        ('k', -cells['ck'], -cells['dk_transport'], cells['pk'] * vol,
         cells['dk'] * k * vol, np.zeros_like(k)),
        ('omega', -cells['cw'], -cells['dw'], pw * vol, yw * vol, cross * vol),
    ):
        residual = convection + diffusion - production + destruction - x
        actual = -cells['rk' if equation == 'k' else 'rw']
        regions = {'all': np.ones_like(k, bool), 'wall': cells['fixed'] > 0,
                   'interior': cells['fixed'] == 0}
        for region, selected in regions.items():
            if equation == 'omega':
                selected &= cells['fixed'] == 0
            if not selected.any():
                continue
            scale = destruction[selected].sum()
            values = {name: float(value[selected].sum() / scale) for name, value in (
                ('C/Y', convection), ('D/Y', diffusion), ('P/Y', production),
                ('Y/Y', destruction), ('X/Y', x), ('R/Y', residual),
                ('absR/Y', np.abs(residual)), ('native_absR/Y', np.abs(actual)),
                ('split_roundoff/Y', np.abs(residual - actual)))}
            values['Y'] = float(scale)
            result[equation + '/' + region] = values
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case_dir', type=Path)
    parser.add_argument('--stem', default='mixing_elbow')
    parser.add_argument('--cache', required=True)
    parser.add_argument('--dump-flux', type=Path)
    parser.add_argument('--diagnostics', type=Path, action='append', default=[])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    st = mo.State(str(args.case_dir), args.stem, args.cache)
    mb.validate_inputs(st)
    for name in ('w', 'wb', 'y'):
        if not np.isfinite(getattr(st, name)).all():
            raise ValueError('missing or nonfinite saved SST input: ' + name)
    hashes = [digest(args.case_dir / (args.stem + suffix)) for suffix in ('.cas.h5', '.dat.h5')]
    if args.dump_flux:
        ptr = np.r_[0, np.cumsum(st.m['nnodes'])]
        with args.dump_flux.open('w', newline='') as stream:
            stream.write(f'# case_sha256={hashes[0]}\n# data_sha256={hashes[1]}\n')
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(('owner', 'neighbour', 'flux', 'nodes'))
            for face, flux in enumerate(st.F):
                writer.writerow((st.ctx.c0[face], st.ctx.c1[face], format(flux, '.17g'),
                                 ' '.join(map(str, st.m['fnodes'][ptr[face]:ptr[face+1]]))))
    output = {'case_sha256': hashes[0], 'data_sha256': hashes[1], 'diagnostics': {}}
    layer = ml.wall_layers(st.m, 8)
    limited = st.mut * st.w / (st.rho * st.k) < 1 - 1e-9
    for directory in args.diagnostics:
        cells, identity = native_cells(directory, hashes, st.ctx.n)
        entry = {'identity': identity, 'budget': budget(cells, st, identity[5])}
        if identity[6] == '0':
            if not np.allclose(cells['k'], st.k, rtol=2e-7, atol=1e-14) or not np.allclose(cells['omega'], st.w, rtol=2e-7, atol=1e-14):
                raise ValueError('zero-iteration budget differs from saved fields')
            cap_ratio = cells['mu_t'] * cells['strain']**2 / (10 * st.rho * .09 * st.k * st.w)
            entry['production_cap'] = {'max_ratio': float(cap_ratio.max()), 'active': int((cap_ratio > 1).sum())}
            sf = st.rho * .31 * st.k / (st.mut * cells['f2'])
            entry['strain_from_saved_mut'] = {}
            for label, rate in (('strain', cells['strain']), ('vorticity', cells['vorticity'])):
                entry['strain_from_saved_mut'][label] = {
                    'limited_L2': float(np.linalg.norm((rate-sf)[limited]) / np.linalg.norm(sf[limited])),
                    'limited_median_ratio': float(np.median((rate/sf)[limited]))}
            entry['layers'] = {}
            # Raw native residuals preserve the actual projected diffusion.
            for l in (1, 2, 3, 4, 5, 99):
                selected = layer == l
                if selected.any():
                    entry['layers'][str(l)] = {
                        'cells': int(selected.sum()),
                        'k_absR': float(np.abs(cells['rk'][selected]).sum()),
                        'omega_absR': float(np.abs(cells['rw'][selected]).sum())}
        saved = {'k': st.k, 'omega': st.w, 'mu_t': st.mut}
        entry['fields'] = {name: float(np.linalg.norm(cells[name]-reference)/np.linalg.norm(reference))
                           for name, reference in saved.items()}
        output['diagnostics'][str(directory)] = entry
    rendered = json.dumps(output, indent=2, allow_nan=False)
    if args.output:
        args.output.write_text(rendered + '\n')
    print(rendered)


if __name__ == '__main__':
    main()
