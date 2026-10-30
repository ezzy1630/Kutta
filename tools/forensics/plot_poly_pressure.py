#!/usr/bin/env python3
"""Plot measured pressure and exact-state face-operator refinement (matplotlib)."""
import argparse
import json
import math
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def rows(path, prefix):
    return [dict(item.split('=', 1) for item in line.split()[1:] if '=' in item)
            for line in path.read_text().splitlines() if line.startswith(prefix)]


def pressure(root, after):
    prefix = 'after_' if after else ''
    out = [float(next(r for r in rows(root / f'{prefix}refine{n}.log' if after else root / f'refine_baseline{n}.log', 'POLY_SOLVE')
                      if r['backend'] == 'cpu')['pressure_L2']) for n in (8, 16)]
    fine = next(r for r in rows(root / ('after_v72.log' if after else 'baseline_v72.log'), 'FSAE V7.2')
                if r['backend'] == 'cpu' and r['poly'] == '1')
    return out + [float(fine['pressure_relative_L2'])]


def operator(records, name, location, poly='1'):
    return next(r for r in records if r['backend'] == 'cpu' and r['poly'] == poly
                and r['wall_data'] == 'analytic' and r['gradient'] == 'lsq'
                and r['operator'] == name and r['location'] == location)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    root = parser.parse_args().directory
    before = json.loads((root / 'operators-baseline.json').read_text())['records']
    after = json.loads((root / 'operators-after.json').read_text())['records']
    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.9), constrained_layout=True)
    red, blue = '#a83b32', '#196581'
    for data, label, color in ((pressure(root, False), 'Before', red), (pressure(root, True), 'After', blue)):
        axes[0].loglog([8, 16, 32], [100*x for x in data], 'o-', color=color, label=label)
    axes[0].axhline(1, color='#555', linestyle='--', label='1% gate (32 rings)')
    axes[0].set(title='Poly Taylor–Couette pressure', ylabel='Relative L2 (%)')
    for records, label, color in ((before, 'Before', red), (after, 'After', blue)):
        r = operator(records, 'viscous_normal_derivative', 'radial_hanging')
        axes[1].loglog(r['rings'], r['rms'], 'o-', color=color, label=f"{label}: p = {r['orders'][-1]:.2f}")
        r = operator(records, 'continuity_divergence', 'interior')
        axes[2].loglog(r['rings'], r['rms'], 'o-', color=color, label=f"{label}: p = {r['orders'][-1]:.2f}")
    wall = operator(after, 'viscous_normal_derivative', 'inner_wall', '0')
    axes[1].loglog(wall['rings'], wall['rms'], ':', color='#777', label='Hex inner wall: p = 0.99')
    axes[1].set(title='Coarse/fine radial viscous flux', ylabel='Area-weighted RMS derivative error')
    axes[2].set(title='Bulk continuity at exact state', ylabel='Volume-weighted RMS divergence')
    for axis in axes:
        axis.set(xlabel='Radial rings', xticks=[8, 16, 32], xticklabels=[8, 16, 32], xlim=(7, 36))
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.legend(fontsize=8)
    fig.savefig(root / 'diagnostics.png', dpi=180)
    fig.savefig(root / 'diagnostics.svg')
    # A readable single-panel copy for the inline report; keep the full study
    # as a separate scientific figure above.
    plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12})
    single, axis = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
    for values, label, color in ((pressure(root, False), 'Before', red), (pressure(root, True), 'After', blue)):
        order = math.log(values[0]/values[2])/math.log(4)
        axis.loglog([8, 16, 32], [100*x for x in values], 'o-', color=color,
                    label=f'{label}: observed order {order:.2f}')
        axis.annotate(f'{100*values[-1]:.3f}%', (32, 100*values[-1]), xytext=(-8, 9),
                      textcoords='offset points', ha='right', color=color)
    axis.axhline(1, color='#555', linestyle='--', label='1% gate at 6,144 cells')
    axis.set(title='Polyhedral Taylor–Couette pressure', xlabel='Radial rings',
             ylabel='Relative L2 (%)', xticks=[8, 16, 32], xticklabels=[8, 16, 32], xlim=(7, 36))
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.legend(fontsize=9)
    single.savefig(root / 'pressure-refinement.png', dpi=180)
    single.savefig(root / 'pressure-refinement.svg')


if __name__ == '__main__':
    main()
