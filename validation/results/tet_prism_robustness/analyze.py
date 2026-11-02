#!/usr/bin/env python3
"""Summarize retained foreground solve logs without accepting unfinished runs.

Usage: python3 validation/results/tet_prism_robustness/analyze.py RUN_DIRECTORY
Outputs JSON to stdout. No meshes, reference curves or cell fields are read.
"""
import csv
import json
import pathlib
import re
import statistics
import sys


def summarize(root, name):
    log = root / (name + '.log')
    if not log.exists():
        return {'run': name, 'available': False}
    match = re.search(r'^RESULT (.*)$', log.read_text(), re.M)
    if not match:
        return {'run': name, 'available': True, 'complete_record': False}
    record = json.loads(match[1])
    result = record['result']
    failed = bool(result['failure'])
    history = root / name / 'history.csv'
    tail = []
    rows = []
    if history.exists():
        with history.open() as stream:
            rows = list(csv.DictReader(stream))
        if rows and not failed:
            last = int(rows[-1]['iteration'])
            tail = [float(row['drag']) for row in rows
                    if int(row['iteration']) > last - 2000]
    footprint = re.search(r'(\d+)\s+peak memory footprint', log.read_text())
    return {
        'run': name,
        'available': True,
        'complete_record': True,
        'backend': record['backend'],
        'source': record['source'],
        'iterations_last_monitor': result['iterations'],
        'failure': result['failure'] or None,
        'converged': not failed and result['residuals_met'] and result['settled'],
        'maximum_residual': result['max_residual'],
        'raw_continuity_last_monitor': float(rows[-1]['continuity_raw']) if rows else None,
        'maximum_velocity_last_monitor': float(rows[-1]['max_velocity']) if rows else None,
        'lift_change': result['settle_lift'],
        'drag_change': result['settle_drag'],
        'cd_final': None if failed else result['drag'],
        'cd_mean_last2000': statistics.mean(tail) if tail else None,
        'cd_range_last2000': [min(tail), max(tail)] if tail else None,
        'peak_footprint_bytes': int(footprint[1]) if footprint else None,
    }


def main():
    root = pathlib.Path(sys.argv[1])
    names = ['baseline_cpu_exact', 'baseline_metal', 'baseline_cpu_equation',
             'baseline_metal_equation', 'incident_metal',
             'candidate_metal_simplec', 'candidate_metal_coupled', 'candidate_cpu_simplec',
             'candidate_cpu_coupled', 'sym_incident_metal', 'gain5_metal',
             'gain5_incident_metal']
    print(json.dumps([summarize(root, name) for name in names], indent=2))


if __name__ == '__main__':
    main()
