#!/usr/bin/env python3
"""Plot accepted compact-grid controls from the archived own-run aggregates."""
import json
import pathlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = pathlib.Path(__file__).resolve().parent
rows = {r['run']: r for r in json.loads((here/'evidence/2026-11-03/summary.json').read_text())}
selected = [('pressure_urf05_L1', 'Baseline'), ('wake25_L1_pressure05', 'Near wake 25 mm'),
            ('wake20_L1_pressure05', 'Near wake 20 mm'), ('shear15_L1_pressure05', 'Broad rear 15 mm'),
            ('shear10_L1_pressure05', 'Broad rear 10 mm'),
            ('shear15_wake25_L1_pressure05', 'Rear 15 + wake 25 mm')]
accepted = [(rows[n], label) for n, label in selected if n in rows and
            rows[n].get('settled') and rows[n].get('residuals_met')]
if not accepted:
    raise ValueError('no accepted runs in the archived selection')
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
fig, ax = plt.subplots(figsize=(9.2, 2.2 + .55 * len(accepted)), layout='constrained')
ax.axvspan(.26 * .92, .26 * 1.08, color='#dcebd9', label='Experimental gate: 0.26 ± 8%')
ax.axvline(.26, color='#608056', ls='--', lw=1)
for i, (r, _) in enumerate(accepted):
    ax.plot(r['cd_final'], i, 'o', color='#356a91', ms=8)
    ax.annotate(f"{r['cd_final']:.5f}", (r['cd_final'], i), xytext=(9, 0),
                textcoords='offset points', va='center', fontsize=11)
ax.set_yticks(range(len(accepted)),
              [f"{label}\n{r['mesh']['cells']:,} cells; {r['iterations']:,} iterations"
               for r, label in accepted])
ax.invert_yaxis()
ax.set_ylim(len(accepted) - .5, -.8)
ax.set_xlim(.235, .38)
ax.set_xticks([.24, .26, .28, .30, .32, .34, .36, .38])
ax.set_xlabel('Final Cd after both existing gates pass')
ax.grid(axis='x', alpha=.2)
ax.legend(frameon=False, loc='lower left')
ax.set_title('Ahmed 35°: accepted compact-grid controls\nSame body/ground faces; pressure relaxation 0.5', fontsize=13)
fig.supxlabel('Selected runs only; the report includes every failed and unconverged control.\n'
              'These results do not establish a mesh-independent SST limit.', fontsize=10)
out = here/'accepted_2026-11-03.png'
fig.savefig(out, dpi=160)
print(out)
