"""Run the requested 40-iteration CLI controls serially and compare outputs.

Only the residual timing column and three summary provenance/timing fields
are excluded. No physical field or convergence quantity is normalized.
"""
import argparse
import csv
import hashlib
import json
import os
import pathlib
import re
import subprocess


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def residual_cells(path):
    lines = path.read_text().splitlines()
    comments = [line for line in lines if line.startswith('#')]
    rows = list(csv.reader(line for line in lines if not line.startswith('#')))
    timing = rows[0].index('seconds_per_iteration')
    values = [[value for index, value in enumerate(row) if index != timing] for row in rows]
    canonical = json.dumps([comments, values], separators=(',', ':')).encode()
    return canonical


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=pathlib.Path)
    parser.add_argument('candidate', type=pathlib.Path)
    parser.add_argument('examples', type=pathlib.Path)
    parser.add_argument('output', type=pathlib.Path)
    args = parser.parse_args()
    executables = {'baseline': args.baseline.resolve(), 'candidate': args.candidate.resolve()}
    env = dict(os.environ, MallocSpaceEfficient='1')
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    cases = [('mixing_elbow', 'mixing_elbow/mixing_elbow.cas.h5', 'wfgc'),
             ('elbow_param', 'mixing_elbow/elbow_param.cas.h5', 'compressibility,ke-options')]
    for label, relative, approximations in cases:
        case = args.examples.resolve() / relative
        for backend in ['cpu', 'metal']:
            directories, runs = {}, {}
            for revision, executable in executables.items():
                name = f'{label}_{backend}_{revision}'
                destination = args.output / name
                if destination.exists():
                    raise RuntimeError(f'{name}: output already exists; choose a fresh output root')
                directories[revision] = destination
                log, resources = args.output / f'{name}.log', args.output / f'{name}.resources.log'
                command = ['/usr/bin/time', '-l', str(executable), 'solve', str(case),
                           '--backend', backend, '--iterations', '40', '--approximate', approximations,
                           '--out', str(destination)]
                print(f'START {name}', flush=True)
                with log.open('w') as out, resources.open('w') as err:
                    result = subprocess.run(command, stdout=out, stderr=err, env=env)
                if result.returncode:
                    raise RuntimeError(f'{name}: exit {result.returncode}; inspect its log')
                summary = json.loads((destination / 'summary.json').read_text())
                if summary['iterations'] != 40 or summary['backend'] != backend:
                    raise RuntimeError(f'{name}: runtime did not confirm 40 iterations on {backend}')
                match = re.search(r'(\d+)\s+peak memory footprint', resources.read_text())
                if not match or int(match[1]) > 3_000_000_000:
                    raise RuntimeError(f'{name}: missing or excessive memory readback')
                runs[revision] = {'binary_sha256': digest(executable), 'iterations': summary['iterations'],
                                  'kutta_commit': summary['kutta_commit'],
                                  'peak_footprint_bytes': int(match[1])}
                print(f'END {name}', flush=True)
            baseline, candidate = directories['baseline'], directories['candidate']
            names = sorted(path.name for path in baseline.iterdir() if path.is_file())
            if names != sorted(path.name for path in candidate.iterdir() if path.is_file()):
                raise RuntimeError('Output file populations differ')
            files = []
            for name in names:
                if name in ['summary.json', 'residuals.csv']:
                    continue
                before, after = baseline / name, candidate / name
                files.append({'file': name, 'same_bytes': before.read_bytes() == after.read_bytes(),
                              'baseline_sha256': digest(before), 'candidate_sha256': digest(after),
                              'baseline_bytes': before.stat().st_size, 'candidate_bytes': after.stat().st_size})
            before_res, after_res = residual_cells(baseline / 'residuals.csv'), residual_cells(candidate / 'residuals.csv')
            summaries = {revision: json.loads((directory / 'summary.json').read_text())
                         for revision, directory in directories.items()}
            excluded = ['kutta_commit', 'wall_time_seconds', 'setup_time_seconds']
            omitted = {revision: {key: summary.pop(key) for key in excluded}
                       for revision, summary in summaries.items()}
            record = {'case': relative, 'case_sha256': digest(case), 'approximations': approximations,
                      'backend': backend, 'runs': runs, 'files': files,
                      'residuals_without_timing_same': before_res == after_res,
                      'excluded_residual_column': 'seconds_per_iteration',
                      'residuals_canonical_baseline_sha256': hashlib.sha256(before_res).hexdigest(),
                      'residuals_canonical_candidate_sha256': hashlib.sha256(after_res).hexdigest(),
                      'summary_without_provenance_and_timing_same': summaries['baseline'] == summaries['candidate'],
                      'excluded_summary_values': omitted}
            record['all_physical_values_exact'] = (all(row['same_bytes'] for row in files)
                                                  and before_res == after_res
                                                  and summaries['baseline'] == summaries['candidate'])
            records.append(record)
            (args.output / 'comparison.json').write_text(json.dumps(records, indent=2) + '\n')
    if not all(row['all_physical_values_exact'] for row in records):
        raise SystemExit('Measured differences exist; inspect comparison.json')
    print('COMPLETE: four CPU/Metal case comparisons have exact physical outputs', flush=True)


if __name__ == '__main__':
    main()
