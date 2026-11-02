#!/usr/bin/env python3
"""Plot aggregate Ahmed results; no restricted curves are read or exported."""
import argparse
import json
import pathlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('summary', type=pathlib.Path)
    ap.add_argument('--out', required=True, type=pathlib.Path)
    a = ap.parse_args()
    rows = json.loads(a.summary.read_text())
    family = sorted([r for r in rows if r['run'].startswith('compact_L') and 'diff' not in r['run']],
                    key=lambda r: r['mesh']['surface_size'], reverse=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), constrained_layout=True)
    ax = axes[0]
    ax.axhspan(.26 * .92, .26 * 1.08, color='#319476', alpha=.13, label='unchanged 0.26 ±8% gate')
    ax.axhline(.26, color='#319476', linestyle='--', linewidth=1)
    ax.plot([r['mesh']['surface_size'] * 1000 for r in family],
            [r['cd_final'] if r['settled'] else r['cd_mean_last2000'] for r in family],
            'o-', color='#286abb', label='compact full-model family')
    for r in rows:
        if r in family:
            continue
        ax.plot(r['mesh']['surface_size'] * 1000,
                r['cd_final'] if r['settled'] else r['cd_mean_last2000'],
                's' if 'reference' in r['run'] else 'D', label=r['run'].replace('_', ' '))
    ax.set(xlabel='specified surface spacing, mm', ylabel='drag coefficient Cd',
           title='Fresh runs, at most 5000 iterations')
    ax.invert_xaxis()
    ax.legend(fontsize=7.5)
    ax = axes[1]
    components = [('front', '#8b70b1'), ('slant', '#286abb'), ('base', '#e88732')]
    for i, r in enumerate(family):
        positive = negative = 0
        for name, color in components:
            value = r['components_final'][name]['pressure']
            bottom = positive if value >= 0 else negative
            ax.bar(i, value, bottom=bottom, color=color, label=f'{name} pressure' if i == 0 else None)
            if value >= 0:
                positive += value
            else:
                negative += value
        ax.bar(i, r['friction_final'], bottom=positive, color='#319476',
               label='all friction' if i == 0 else None)
    ax.axhline(0, color='#777777', linewidth=.6)
    ax.set(xticks=range(len(family)), xticklabels=[f"L{r['mesh']['level']:g}" for r in family],
           ylabel='drag coefficient contribution', title='Final-iterate components')
    ax.legend(fontsize=7.5)
    for ax in axes:
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(a.out, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
