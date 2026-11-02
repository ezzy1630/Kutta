"""Summarize completed native benchmark records without treating caps as convergence."""
import csv
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
rows = []
for path in sorted(root.glob('*.log')):
    records = [line.removeprefix('RESULT ') for line in path.read_text().splitlines()
               if line.startswith('RESULT ')]
    if not records:
        continue
    if len(records) != 1:
        raise ValueError(f'{path.name}: expected one benchmark record')
    record = json.loads(records[0])
    r = record['result']
    history = root / path.stem / 'history.csv'
    samples = list(csv.DictReader(history.open())) if history.exists() else []
    tail = [s for s in samples if int(s['iteration']) > r['iterations'] - 2000]
    complete = not r['failure'] and not r['budget']
    passed = complete and r['residuals_met'] and r['settled']
    footprint = next((int(line.split()[0]) for line in path.read_text().splitlines()
                      if 'peak memory footprint' in line), None)
    row = {'run': path.stem, 'backend': record['backend'], 'iterations': r['iterations'],
           'failure': r['failure'] or None, 'converged': passed,
           'residuals_met': r['residuals_met'], 'settled': r['settled'],
           'maximum_scaled_residual': r['max_residual'] if samples else None,
           'relative_lift_change': r['settle_lift'] if samples else None, 'relative_drag_change': r['settle_drag'] if samples else None,
           'peak_footprint_bytes': footprint, 'raw_continuity': None,
           'final_drag': r['drag'] if complete else None,
           'mean_drag_last2000': None, 'range_drag_last2000': None,
           'tail_window_requested_iterations': 2000,
           'tail_window_available_iterations': min(2000, r['iterations']),
           'tail_sample_count': len(tail)}
    if samples:
        row['raw_continuity'] = float(samples[-1]['continuity_raw'])
        row['last_residuals'] = {k: float(samples[-1][k]) for k in ['u', 'v', 'w', 'continuity', 'k', 'omega']}
    if complete and tail:
        drag = [float(s['drag']) for s in tail]
        row['mean_drag_last2000'] = sum(drag) / len(drag)
        row['range_drag_last2000'] = [min(drag), max(drag)]
    rows.append(row)
print(json.dumps(rows, indent=2))
