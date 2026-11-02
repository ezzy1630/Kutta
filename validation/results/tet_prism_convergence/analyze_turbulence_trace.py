"""Localize the rejected exact-symmetry Coupled control's scalar residuals.

The 250-iteration sampling cannot resolve a temporal mode. These are assembled
steady-equation residuals, not errors against Fluent or physical-time data.
"""
import argparse
import json
import pathlib

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trace', type=pathlib.Path)
    parser.add_argument('output', type=pathlib.Path)
    args = parser.parse_args()
    meta = np.genfromtxt(args.trace / 'cells.csv', delimiter=',', names=True)
    n = len(meta)
    far = ((meta['x'] > 1) | (meta['x'] < -1.5) |
           (np.abs(meta['y']) > .5) | (meta['z'] > .6))
    regions = {'farfield': far, 'body': ~far & (meta['x'] < 0),
               'wake': ~far & (meta['x'] >= 0), 'tets': meta['faces'] == 4,
               'prisms': meta['faces'] == 5, 'poor': meta['poor'] > 0}
    result = {'control': 'exact symmetry constraint, Coupled, Metal',
              'not_retained_option': True, 'sampling_iterations': 250,
              'physical_frequency_unresolved': True, 'equations': {}}
    for equation in ['k', 'omega']:
        records = []
        for iteration in range(1000, 5001, 250):
            values = np.fromfile(args.trace / f'{iteration}_{equation}.bin', dtype=np.float32).reshape(8, n)
            absolute = values[0].astype(np.float64)
            total = absolute.sum()
            monitor = np.fromfile(args.trace / f'{iteration}_monitor.bin', dtype=np.float32).reshape(24, n)
            state = monitor[4 if equation == 'k' else 5]
            order = np.argsort(absolute)[-8:][::-1]
            records.append({'iteration': iteration, 'raw_residual': float(total),
                            'region_fractions': {name: float(absolute[mask].sum() / total)
                                                 for name, mask in regions.items()},
                            'top_cells': [{**{name: float(meta[name][c]) for name in meta.dtype.names},
                                           'absolute_residual': float(absolute[c]),
                                           'signed_algebraic_residual': float(values[1, c]),
                                           'limiter': float(values[2, c]),
                                           'state': float(state[c]),
                                           'source_per_volume': float(values[6, c]),
                                           'sink_coefficient': float(values[7, c])}
                                          for c in order]})
        result['equations'][equation] = records
    args.output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
