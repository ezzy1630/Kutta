#!/usr/bin/env python3
"""Recreate the mesh-control figure from archived aggregate solver results."""
import json, pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

here = pathlib.Path(__file__).resolve().parent
base = here/'evidence/2026-11-03'
rows = {r['run']: r for r in json.loads((base/'summary.json').read_text())}
bands = {r['run']: r for r in json.loads((base/'cut_band.json').read_text())}
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
fig = plt.figure(figsize=(13, 9.2), layout='constrained')
grid = fig.add_gridspec(2, 2, width_ratios=[.95, 1.5])
ax = fig.add_subplot(grid[:, 0])
names = [['baseline_L1', 'volume_L2_body_L1'], ['body_L2_volume_L1', 'baseline_L2']]
values = np.array([[rows[n]['cd_final'] if rows[n]['settled'] and rows[n]['residuals_met']
                    else rows[n]['cd_mean_last2000'] for n in line] for line in names])
ax.imshow(values, cmap='Blues', vmin=.33, vmax=.375)
for i, line in enumerate(names):
    for j, name in enumerate(line):
        r = rows[name]
        converged = r['settled'] and r['residuals_met']
        label = f"Cd {values[i,j]:.5f}\n{'converged' if converged else 'unconverged'}\n{r['mesh']['cells']:,} cells"
        ax.text(j, i, label, ha='center', va='center', color='white' if values[i,j]>.357 else '#162838')
ax.set_xticks([0, 1], ['39.6 mm', '28 mm'])
ax.set_yticks([0, 1], ['39.6 mm', '28 mm'])
ax.set_xlabel('Prescribed exterior mesh spacing*')
ax.set_ylabel('Prescribed whole-body surface spacing')
ax.set_title('Body surface / exterior controls\nSame first height and wall-layer growth', pad=18)
ax.grid(False)

cp_ax = fig.add_subplot(grid[0, 1])
cp_names = [('baseline_L1', 'L1 baseline'), ('edge10_L1', 'Edge surface 10 mm'),
            ('edge5_L1', 'Edge surface 5 mm'), ('slant10_L1', 'Whole slant 10 mm'),
            ('shear15_L1_pressure05', 'Rear volume 15 mm, pressure 0.5'),
            ('wake20_L1_pressure05', 'Wake 20 mm, pressure 0.5'),
            ('edge10_shear10_L1_pressure05', 'Edge 10 + rear 10 mm, pressure 0.5'),
            ('slant15_shear10_L1_pressure05', 'Slant 15 + rear 10 mm, pressure 0.5'),
            ('baseline_L2_urfu03_urfp05', 'L2 baseline, momentum/pressure 0.3/0.5'),
            ('final_refined_L2_urfu03_urfp05', 'Finer combined, momentum/pressure 0.3/0.5')]
cp_names = [(n, label) for n, label in cp_names if n in rows and 'cp_common' in rows[n]]
errors = [rows[n]['cp_common']['experiment slant_leading_40mm']['rms'] for n, _ in cp_names]
points = rows['baseline_L1']['cp_common']['experiment slant_leading_40mm']['points']
cp_ax.barh(np.arange(len(cp_names)), errors, color='#476f8c')
def gate_mark(name):
    return '' if rows[name]['settled'] and rows[name]['residuals_met'] else ' †'

cp_ax.set_yticks(np.arange(len(cp_names)), [label + gate_mark(n) for n, label in cp_names])
cp_ax.invert_yaxis()
cp_ax.set_xlabel('RMS Cp error')
cp_ax.set_title(f'Leading 40 mm in x: same {points} experimental taps')
cp_ax.set_xlim(left=0)
cp_ax.grid(axis='x', alpha=.2)

load_ax = fig.add_subplot(grid[1, 1])
load_names = [('pressure_urf05_L1', 'Baseline, pressure 0.5'),
              ('edge10_L1', 'Edge surface 10 mm'),
              ('shear15_L1_pressure05', 'Rear 15 mm, pressure 0.5'),
              ('wake20_L1_pressure05', 'Wake 20 mm, pressure 0.5'),
              ('shear15_wake25_L1_pressure05', 'Rear 15 + wake 25 mm, pressure 0.5'),
              ('slant15_shear10_L1_pressure05', 'Slant 15 + rear 10 mm, pressure 0.5'),
              ('final_refined_L2_urfu03_urfp05', 'Finer combined, momentum/pressure 0.3/0.5')]
load_names = [(n, label) for n, label in load_names if n in rows and n in bands and 'cp_common' in rows[n]]
lead = [bands[n]['leading_pressure'] for n, _ in load_names]
rest = [bands[n]['remaining_pressure'] for n, _ in load_names]
y = np.arange(len(load_names))
load_ax.barh(y, lead, color='#cb8544', label='First 40 mm')
load_ax.barh(y, rest, left=lead, color='#476f8c', label='Remaining slant')
def reported_cd(r):
    return r['cd_final'] if r['settled'] and r['residuals_met'] else r['cd_mean_last2000']
load_ax.set_yticks(y, [f'{label}; body Cd {reported_cd(rows[n]):.5f}{gate_mark(n)}' for n, label in load_names])
load_ax.invert_yaxis()
load_ax.set_xlabel('Final-iterate slant pressure Cd')
load_ax.set_title('Equal geometric bands, clipped face polygons', pad=34)
load_ax.grid(axis='x', alpha=.2)
load_ax.legend(frameon=False, loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2)
fig.suptitle('Ahmed 35°: mesh sensitivity; the 0.26 ± 8% gate remains open', fontsize=14)
fig.supxlabel('Cd: final when converged, last-2000 iteration mean otherwise. Cp and bands: final iterate.\n'
              '*Exterior level includes ground sizing. Bands start 40 mm in x (48.8 mm along the slant).\n'
              'Cp retains the outlet pressure gauge; no reference offset is fitted.\n'
              '† Both numerical gates are not met. Selected finite snapshots are shown; the report includes failures.\n'
              'These controls do not establish a mesh-independent limit.', fontsize=9)
out = here/'cause_2026-11-03.png'
fig.savefig(out, dpi=160)
print(out)
