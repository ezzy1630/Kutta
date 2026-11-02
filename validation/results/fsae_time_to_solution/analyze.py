#!/usr/bin/env python3
"""Summarize the retained numeric benchmark receipts; no fixtures are required."""
import collections
import hashlib
import json
import pathlib
import statistics

ROOT = pathlib.Path(__file__).resolve().parent


def rows(name):
    path = ROOT / 'raw' / name
    return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []


def distribution(values):
    return {'n': len(values), 'median': statistics.median(values),
            'min': min(values), 'max': max(values)}


def summarize(name):
    groups = collections.defaultdict(list)
    for row in rows(name):
        groups[(row['case'], row['scheme'], row['label'])].append(row)
    result = []
    for (case, scheme, label), values in sorted(groups.items()):
        result.append({'case': case, 'scheme': scheme, 'label': label,
                       'cells': values[0]['cells'],
                       'seconds_per_iteration': distribution([x['solve_s'] / x['iterations'] for x in values]),
                       'solve_seconds': distribution([x['solve_s'] for x in values]),
                       'process_seconds': distribution([x['process_real_s'] for x in values]),
                       'setup_seconds': distribution([x['setup_s'] for x in values]),
                       'iterations': [x['iterations'] for x in values],
                       'peak_bytes': max(x['peak_bytes'] for x in values),
                       'all_converged': all(x['residuals_met'] and x['settled'] for x in values)})
    return result


def paired(name):
    groups = collections.defaultdict(dict)
    for row in rows(name):
        groups[(row['case'], row['scheme'], row['rep'])][row['label']] = row
    result = collections.defaultdict(list)
    for (case, scheme, rep), pair in sorted(groups.items()):
        if set(pair) != {'baseline', 'candidate'}:
            raise ValueError(f'incomplete pair: {case} {scheme} {rep}')
        before, after = pair['baseline'], pair['candidate']
        prefix = 'measure' if name == 'measurements.jsonl' else 'converge'
        histories = {}
        for label in pair:
            path = ROOT / 'raw' / f'{prefix}-{rep}-{case}-{scheme}-{label}.log'
            lines = [line for line in path.read_text().splitlines() if line.startswith('fsae_aero it=')]
            histories[label] = {'checkpoints': len(lines),
                                'sha256': hashlib.sha256('\n'.join(lines).encode()).hexdigest()}
        if histories['baseline'] != histories['candidate']:
            raise ValueError(f'changed physical history: {case} {scheme} {rep}')
        result[(case, scheme)].append({
            'rep': rep, 'iterations_equal': before['iterations'] == after['iterations'],
            'physical_histories': histories, 'physical_histories_identical': True,
            'solve_speedup': before['solve_s'] / after['solve_s'],
            'process_speedup': before['process_real_s'] / after['process_real_s'],
            'both_converged': all(r['residuals_met'] and r['settled'] for r in pair.values()),
            'final_values_equal': all(before[k] == after[k] for k in
                                      ('lift', 'drag', 'drag_pressure', 'max_residual',
                                       'settle_lift', 'settle_drag', 'yplus_max', 'yplus_mean'))})
    return [{'case': case, 'scheme': scheme,
             'timing_scope': 'converged_answer' if all(x['both_converged'] for x in values) else 'iteration_budget',
             'repeated_pairs': len(values) >= 3,
             'solve_speedup': distribution([x['solve_speedup'] for x in values]),
             'process_speedup': distribution([x['process_speedup'] for x in values]),
             'pairs': values} for (case, scheme), values in sorted(result.items())]


def family(name):
    if name in {'amg_galerkin', 'amg_m_matrix_diagonal', 'krylov_diagonal',
                'coupled_fine_coefficients', 'coupled_fine_gsum', 'coupled_fine_schur',
                'coupled_coarse_gsum', 'coupled_coarse_schur', 'coupled_velocity_block',
                'coupled_rhs', 'coupled_face_coefficient', 'coupled_copy_diagonal',
                'coupled_floor_physical'} or name.startswith('coupled_galerkin'):
        return 'Operator refresh and validation'
    if name == 'sequence':
        return 'Small-operation replay (mixed)'
    if name.startswith(('csr_', 'solver_', 'krylov_', 'reduce')) or name == 'coupled_fine_spmv':
        return 'Krylov/CSR/reductions'
    if name.startswith(('amg_', 'coupled_fine_', 'coupled_coarse_', 'coupled_restrict', 'coupled_prolong')):
        return 'AMG smoothing, transfers and shared residuals'
    if name.startswith('sst_') or name in {'fv_effective_viscosity', 'fv_turbulence_boundary', 'fv_floor'}:
        return 'SST closure'
    if name.startswith('fv_transport') or name in {'fv_projected_transport_faces', 'fv_fixed_rows'}:
        return 'Transport assembly (momentum and turbulence)'
    if any(word in name for word in ('gradient', 'limiter', 'moment_', 'pressure_faces')):
        return 'Gradients, limiter and face reconstruction'
    if name.startswith(('fv_momentum_', 'fv_paired_', 'coupled_paired_')) or name in {
            'fv_backflow_pressure', 'fv_boundary_momentum', 'fv_cell_pressure_boundary',
            'fv_pressure_assemble', 'fv_pressure_boundary', 'fv_rhie_chow',
            'fv_velocity_boundary', 'fv_stress_source', 'fv_transpose_stress'}:
        return 'Momentum, pressure and flux coupling'
    return 'Other FV/runtime'


def profiles():
    result = []
    for stage in ['baseline', 'candidate']:
        for case in ['wing', 'ahmed']:
            for scheme in ['simplec', 'coupled']:
                paths = sorted((ROOT / 'raw').glob(f'profile-quiet-{stage}-{case}-{scheme}-*.log'))
                if not paths:
                    continue
                families, kernels, host = collections.defaultdict(list), collections.defaultdict(list), collections.defaultdict(list)
                for path in paths:
                    groups = collections.defaultdict(float)
                    for line in path.read_text().splitlines():
                        if line.startswith('GPU_PROFILE '):
                            _, kernel, _, seconds = line.split()
                            ms = float(seconds) * 1000 / 50
                            groups[family(kernel)] += ms
                            kernels[kernel].append(ms)
                    for key, ms in groups.items():
                        families[key].append(ms)
                    count_path = path.with_name(path.name.replace('profile-quiet-', 'counts-quiet-'))
                    for line in count_path.read_text().splitlines():
                        if line.startswith('HOST_PROFILE '):
                            for key, value in (x.split('=') for x in line.split()[1:]):
                                host[key].append(float(value) / 50)
                result.append({'stage': stage, 'case': case, 'scheme': scheme,
                               'gpu_ms_per_iteration': {k: distribution(v) for k, v in families.items()},
                               'kernels_ms_per_iteration': {k: distribution(v) for k, v in kernels.items()},
                               'native_host_per_iteration': {k: distribution(v) for k, v in host.items()}})
    return result


if __name__ == '__main__':
    summary = {'fixed_iterations': summarize('measurements.jsonl'),
               'fixed_pairs': paired('measurements.jsonl'),
               'convergence': summarize('convergence.jsonl'),
               'convergence_pairs': paired('convergence.jsonl'), 'profiles': profiles()}
    (ROOT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    for row in summary['fixed_iterations']:
        d = row['seconds_per_iteration']
        print(f"{row['case']} {row['scheme']} {row['label']}: "
              f"{d['median']:.6f} [{d['min']:.6f}, {d['max']:.6f}] s/it, n={d['n']}")
