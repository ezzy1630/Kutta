"""Render the measured iteration histories and the localized period-two mode."""
import argparse
import csv
import pathlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    return {name: np.array([float(row[name]) for row in rows]) for name in rows[0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', type=pathlib.Path)
    parser.add_argument('--candidate', default='final_native_metal_simplec')
    args = parser.parse_args()
    destination = pathlib.Path(__file__).resolve().parent / 'evidence'
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'figure.dpi': 140})
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True, constrained_layout=True)
    for name, label, color in [('base_metal_simplec', 'Unchanged 22e65c7', '#59636d'),
                               (args.candidate, 'Opt-in normal-foot boundary fit', '#126b96')]:
        data = read_csv(args.runs / name / 'history.csv')
        axes[0].semilogy(data['iteration'], data['continuity_raw'], label=label, color=color)
        axes[1].semilogy(data['iteration'], data['k'], color=color)
        axes[2].plot(data['iteration'], data['drag'], color=color)
    axes[0].set_ylabel('Raw continuity')
    axes[0].legend(loc='upper right')
    axes[1].set_ylabel('Scaled k residual')
    axes[1].axhline(1e-3, color='#a63333', linestyle='--', linewidth=1)
    axes[2].set_ylabel('Drag coefficient')
    axes[2].set_xlabel('Steady iteration')
    axes[0].set_title('Native Ahmed tet/prism mesh: neither run meets the full convergence gate', loc='left')
    for axis in axes:
        axis.grid(alpha=.18)
    fig.savefig(destination / 'native_histories.png')
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, constrained_layout=True)
    controls = [('final_native_metal_simplec', 'Standard velocity and scalars', '#59636d'),
                ('diff_velocity_normal_foot_metal_simplec', 'Smooth velocity only', '#126b96'),
                ('diff_turbulence_normal_foot_metal_simplec', 'Smooth velocity and scalars', '#ad6328')]
    for name, label, color in controls:
        data = read_csv(args.runs / name / 'history.csv')
        maximum = np.stack([data[key] for key in ['u', 'v', 'w', 'continuity', 'k', 'omega']]).max(axis=0)
        settling = np.full(len(data['iteration']), np.nan)
        for index, iteration in enumerate(data['iteration']):
            if iteration - data['iteration'][0] < 1000:
                continue
            window = (data['iteration'] <= iteration) & (iteration - data['iteration'] <= 1000)
            settling[index] = max(float(np.max(np.abs(data[key][window] - data[key][index])) /
                                        max(abs(data[key][index]), 1e-12)) for key in ['lift', 'drag'])
        axes[0].semilogy(data['iteration'], maximum, color=color, label=label)
        axes[1].semilogy(data['iteration'], settling, color=color)
    axes[0].set_title('Native mesh limiter controls: no run passes both gates', loc='left')
    axes[0].set_ylabel('Maximum scaled residual')
    axes[0].axhline(1e-3, color='#a63333', linestyle='--', linewidth=1)
    axes[0].legend(loc='upper right')
    axes[1].set_ylabel('Maximum relative force change\nover 1,000 iterations')
    axes[1].axhline(5e-4, color='#a63333', linestyle='--', linewidth=1)
    axes[1].set_xlabel('Steady iteration; all controls use the normal-foot pressure option')
    for axis in axes:
        axis.grid(alpha=.18)
    fig.savefig(destination / 'limiter_controls.png')
    plt.close(fig)

    with (destination / 'cpu_cycle_cells.csv').open() as stream:
        rows = [row for row in csv.DictReader(stream) if row['original'] == '51376']
    data = {name: np.array([float(row[name]) for row in rows]) for name in rows[0]}
    selected = data['iteration'] <= 1020
    fig, axes = plt.subplots(3, 1, figsize=(9, 7), sharex=True, constrained_layout=True)
    for axis, name, label in zip(axes, ['w', 'p', 'pressurez'],
                                 ['Vertical velocity', 'Pressure', 'z pressure source']):
        axis.plot(data['iteration'][selected], data[name][selected], 'o-', color='#126b96', markersize=3)
        axis.set_ylabel(label)
        axis.grid(alpha=.18)
    axes[0].set_title('Unchanged CPU solve, input cell 51376 in the far field', loc='left')
    axes[2].set_xlabel('Steady iteration; source is from the preceding assembly')
    axes[2].set_xticks(range(1000, 1021, 2))
    fig.savefig(destination / 'cpu_cell_cycle.png')
    plt.close(fig)


if __name__ == '__main__':
    main()
