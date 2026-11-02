#!/usr/bin/env python3
"""Plot the completed robustness study. Uses benchmark histories, not cell fields.

Usage: python3 plot.py RUN_DIRECTORY OUTPUT_PNG
Requires matplotlib. Native and compact-hex drag statistics are kept separate.
"""
import csv
import pathlib
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from analyze import summarize


def history(root, name):
    with (root / name / 'history.csv').open() as stream:
        return [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(stream)]


def main():
    root, output = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.4), constrained_layout=True)
    curves = [
        ('baseline_cpu_exact', 'CPU before', '#9d3030', '--'),
        ('baseline_metal', 'Metal before', '#d36b45', '--'),
        ('candidate_cpu_simplec', 'CPU after', '#285ba6', '-'),
        ('candidate_metal_simplec', 'Metal after', '#19856a', '-'),
    ]
    for name, label, color, style in curves:
        rows = history(root, name)
        x = [row['iteration'] for row in rows]
        axes[0, 0].semilogy(x, [row['max_velocity'] for row in rows],
                           style, color=color, label=label, linewidth=1.5)
        axes[0, 1].semilogy(x, [max(row[key] for key in
                           ['u', 'v', 'w', 'continuity', 'k', 'omega']) for row in rows],
                           style, color=color, label=label, linewidth=1.5)
    axes[0, 0].set(title='SIMPLEC maximum speed', xlabel='Iteration', ylabel='Speed / U')
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 1].axhline(1e-3, color='#555555', linestyle=':', label='Residual criterion')
    axes[0, 1].set(title='Bounded runs still miss convergence', xlabel='Iteration',
                   ylabel='Maximum scaled residual')
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[1, 0].plot([.327584, .563140, .919620], [.36207, .34275, .33032],
                    'o-', color='#555555')
    axes[1, 0].set(title='Retained compact hex family', xlabel='Cells, millions',
                   ylabel='Last-2000 mean Cd')
    for x, y in zip([.327584, .563140, .919620], [.36207, .34275, .33032]):
        axes[1, 0].annotate(f'{y:.3f}', (x, y), xytext=(0, 7),
                           textcoords='offset points', ha='center', fontsize=9)
    axes[1, 0].set_ylim(.315, .38)
    axes[1, 0].text(.5, .02, 'Unconverged statistics; compact exterior grading',
                    transform=axes[1, 0].transAxes, ha='center', fontsize=8)
    names = ['candidate_cpu_simplec', 'candidate_metal_simplec',
             'candidate_cpu_coupled', 'candidate_metal_coupled']
    labels = ['CPU\nSIMPLEC', 'Metal\nSIMPLEC', 'CPU\nCoupled', 'Metal\nCoupled']
    for i, name in enumerate(names):
        result = summarize(root, name)
        mean = result['cd_mean_last2000']
        lo, hi = result['cd_range_last2000']
        axes[1, 1].errorbar(i, mean, yerr=[[mean - lo], [hi - mean]],
                           fmt='o', color='#285ba6' if i % 2 == 0 else '#19856a',
                           capsize=5)
        axes[1, 1].annotate(f'{mean:.3f}', (i, hi), xytext=(0, 7),
                           textcoords='offset points', ha='center', fontsize=9)
    axes[1, 1].set(title='Native mesh, 130,331 cells, unconverged', ylabel='Cd',
                   xticks=range(4), xticklabels=labels)
    axes[1, 1].margins(x=.2, y=.3)
    axes[1, 1].text(.5, .02, 'Bars show iterative min/max, not uncertainty bounds',
                    transform=axes[1, 1].transAxes, ha='center', fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=.18)
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Native prism pressure correction: blow-up prevention, convergence still open',
                 fontsize=13)
    fig.savefig(output, dpi=160)


if __name__ == '__main__':
    main()
