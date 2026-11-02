"""Analyze the read-only native-cell traces described in evidence/README.md.

The lowest nonzero FFT bin spans the whole window and may be drift. It is
not evidence of a resolved physical period. Full field dumps stay ignored.
"""
import json
import pathlib
import sys

import numpy as np

FIELDS = ('u v w p k omega limiter_u limiter_v limiter_w '
          'gradpx gradpy gradpz sourcex sourcey sourcez continuity '
          'diagonal correction pressurex pressurey pressurez').split()
SELECTED = 'u v w p limiter_u limiter_v limiter_w gradpz continuity pressurez'.split()


def analyze(root):
    meta = np.genfromtxt(root / 'cells.csv', delimiter=',', names=True)
    count = len(meta)
    files = [root / f'{iteration}_step.bin' for iteration in range(1000, 1100)]
    indices = [FIELDS.index(name) for name in SELECTED]
    # Keep only the ten fields needed for spatial/spectral statistics.
    data = np.stack([np.fromfile(path, dtype=np.float32).reshape(len(FIELDS), count)[indices]
                     for path in files])
    mean, std, delta = data.mean(axis=0), data.std(axis=0), np.diff(data, axis=0)
    continuity = SELECTED.index('continuity')
    limiter = slice(4, 7)
    regions = {
        'all': np.ones(count, dtype=bool),
        'farfield': ((meta['x'] > 1) | (meta['x'] < -1.5) |
                     (np.abs(meta['y']) > .5) | (meta['z'] > .6)),
        'wake': ((meta['x'] >= 0) & (meta['x'] <= 1) &
                 (np.abs(meta['y']) < .5) & (meta['z'] < .6)),
        'body': ((meta['x'] >= -1.5) & (meta['x'] < 0) &
                 (np.abs(meta['y']) < .5) & (meta['z'] < .6)),
        'tets': meta['faces'] == 4,
        'prisms': meta['faces'] == 5,
        'poor': meta['poor'] > 0,
        'boundary': meta['boundary'] > 0,
    }
    result = {'window_iterations': len(files), 'lowest_frequency_bin_may_be_drift': True,
              'n_cells': count, 'iterations': list(range(1000, 1100)), 'regions': {}}
    for name, mask in regions.items():
        result['regions'][name] = {
            'n': int(mask.sum()),
            'mean_raw_continuity': float(mean[continuity, mask].sum()),
            'continuity_fraction': float(mean[continuity, mask].sum() / mean[continuity].sum()),
            'limiter_flips_gt_005_mean': float(
                (np.max(np.abs(delta[:, limiter, mask]), axis=1) > .05).sum(axis=1).mean()),
            'velocity_delta_l1': float(np.linalg.norm(delta[:, :3, mask], axis=1).sum(axis=1).mean()),
            'pressure_std_sum': float(std[3, mask].sum()),
        }
    for name in ['continuity', 'p', 'u', 'gradpz', 'pressurez', 'limiter_u']:
        index = SELECTED.index(name)
        rank = mean[index] if name == 'continuity' else std[index]
        order = np.argsort(rank)[-12:][::-1]
        result[name + '_top'] = [
            {**{key: float(meta[key][cell]) for key in meta.dtype.names},
             'value': float(rank[cell]), 'mean': float(mean[index, cell]),
             'std': float(std[index, cell]), 'first': float(data[0, index, cell]),
             'last': float(data[-1, index, cell])} for cell in order]
    result['series'] = {
        'iteration': list(range(1000, 1100)),
        'continuity': data[:, continuity, :].sum(axis=1).tolist(),
        'limiter_flips_gt_005': (np.max(np.abs(delta[:, limiter, :]), axis=1) > .05).sum(axis=1).tolist(),
    }
    # Sum cell powers so opposite-phase cells cannot cancel each other.
    for region, mask in regions.items():
        if not mask.any():
            continue
        for index in [0, 3, 4, continuity]:
            series = data[:, index, mask] - mean[index, mask]
            power = (np.abs(np.fft.rfft(series, axis=0)) ** 2).sum(axis=1)
            power[0] = 0
            peak = int(power.argmax())
            result['regions'][region][SELECTED[index] + '_dominant_period'] = len(files) / peak if peak else None
    (root / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    print(json.dumps(analyze(pathlib.Path(sys.argv[1]))['regions'], indent=2))
