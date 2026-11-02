#!/usr/bin/env python3
"""Plot paired solve speedups from summary.json (requires matplotlib 3.11.2)."""
import json
import pathlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent


def plot(rows, filename, title):
    rows = sorted(rows, key=lambda r: (
        ['wing-L2', 'ahmed-L1', 'ahmed-L3'].index(r['case']), r['scheme'] == 'coupled'))
    labels = [f"{r['case'].replace('wing', 'Wing').replace('ahmed', 'Ahmed')} / {r['scheme'].upper() if r['scheme'] == 'simplec' else 'Coupled'}" for r in rows]
    medians = [r['solve_speedup']['median'] for r in rows]
    low = [m-r['solve_speedup']['min'] for m, r in zip(medians, rows)]
    high = [r['solve_speedup']['max']-m for m, r in zip(medians, rows)]
    fig, ax = plt.subplots(figsize=(8.2, 1.5 + .45 * len(rows)), layout='constrained')
    ax.axvline(1, color='#777777', linewidth=1, linestyle='--')
    ax.errorbar(medians, range(len(rows)), xerr=[low, high], fmt='o', capsize=4,
                color='#2469aa', markersize=6)
    ax.set_yticks(range(len(rows)), labels)
    ax.invert_yaxis()
    ax.set_xlabel('Baseline / candidate solve time (higher is faster)')
    ax.set_title(title + '\nMedian and complete paired range')
    ax.grid(axis='x', alpha=.18)
    ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(ROOT / f'{filename}.png', dpi=180)
    svg = ROOT / f'{filename}.svg'
    fig.savefig(svg, metadata={'Date': None})
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
    plt.close(fig)


def main():
    data = json.loads((ROOT / 'summary.json').read_text())
    plot(data['fixed_pairs'], 'speedup',
         'Metal: 200 iterations, four counterbalanced pairs')
    converged = [r for r in data['convergence_pairs']
                 if r['timing_scope'] == 'converged_answer' and r['repeated_pairs']]
    if converged:
        plot(converged, 'converged-speedup',
             'Metal: unchanged convergence criteria, four counterbalanced pairs')


if __name__ == '__main__':
    main()
