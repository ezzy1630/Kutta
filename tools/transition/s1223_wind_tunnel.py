#!/usr/bin/env python3
"""Summarize the S1223 investigation without copying restricted curves or meshes.

Reads the fresh paired run.py output below runs/s1223_wind_tunnel.
Reference values are loaded only at the run angles.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / 'runs/s1223_wind_tunnel'
OUT = ROOT / 'validation/results/s1223_wind_tunnel'
POINTS = {2e5: ['0.00', '4.76', '7.93', '10.00', '13.17', '16.00'],
          3e5: ['0.48', '4.00', '8.21', '12.00', '14.74', '16.00']}


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def history(path):
    import gzip
    p = path / 'history.csv'
    if p.exists():
        return np.genfromtxt(p, delimiter=',', names=True)
    with gzip.open(str(p) + '.gz', 'rt') as stream:
        return np.genfromtxt(stream, delimiter=',', names=True)


def values(path):
    r = json.loads((path / 'result.json').read_text())
    h = np.atleast_1d(history(path))
    q = r['result']
    if q['failure']:
        raise ValueError(f'Failed solve has no accuracy coefficient: {path}: {q["failure"]}')
    if not len(h) or h['iteration'][-1] != q['iterations']:
        raise ValueError(f'History and completed receipt disagree: {path}')
    # Preserve the prior report's convention: settled final value, otherwise
    # mean of samples strictly inside the last 2000-iteration window.
    selected = h[h['iteration'] > q['iterations'] - 2000]
    steady = q['residuals_met'] and q['settled'] and not q['budget'] and not q['failure']
    cl, cd = (q['lift'], q['drag']) if steady else (float(selected['lift'].mean()), float(selected['drag'].mean()))
    if not np.isfinite([cl, cd]).all():
        raise ValueError(f'No finite force comparison is available: {path}')
    return r, h, float(cl), float(cd), steady


def freestream(path, re, tu=.001, ratio=1, inlet=-30):
    """Nearest cell samples; reference Tu uses U_infinity=1, not local speed."""
    cells = np.genfromtxt(path / 'cells.csv', delimiter=',', names=True)
    k0 = 1.5 * tu**2
    w0 = k0 * re / ratio
    rows = []
    for x, y in [(-20, 2), (-10, 2), (-2, 2), (-.1, .2), (-.05, .5)]:
        i = np.argmin((cells['x'] - x)**2 + (cells['y'] - y)**2)
        c = cells[i]
        factor = 1 + .0828 * w0 * (c['x'] - inlet)
        expected_k = k0 * factor**(-.09/.0828)
        expected_w = w0 / factor
        rows.append(dict(x=float(c['x']), y=float(c['y']), tu_percent=float(100*np.sqrt(2*c['k']/3)),
                         k_over_analytic=float(c['k']/expected_k), omega_over_analytic=float(c['omega']/expected_w),
                         speed=float(np.hypot(c['u'], c['v']))))
    return rows


def bubbles(path):
    """Reversed shear spans, using the historical report's chord/side convention."""
    w = np.genfromtxt(path / 'walls.csv', delimiter=',', names=True, dtype=None, encoding='ascii')
    w = w[w['patch'] == 'wing']
    p = np.column_stack([w['x'], w['y']])
    te = p[np.argmax(p[:, 0])]
    le = p[np.argmax(np.linalg.norm(p-te, axis=1))]
    chord = te-le
    x = (p-le) @ chord / (chord @ chord)
    upper = -chord[1]*w['nx'] + chord[0]*w['ny'] < 0
    tangent = np.column_stack([-w['ny'], w['nx']])
    tangent *= np.sign(tangent @ chord)[:, None]
    flow = -(w['tau_x']*tangent[:, 0] + w['tau_y']*tangent[:, 1])
    indices = np.flatnonzero(upper & (x > .02))
    indices = indices[np.argsort(x[indices])]
    spans, current = [], []
    for i in indices:
        if flow[i] < 0:
            current.append(i)
        elif current:
            spans.append(current)
            current = []
    if current:
        spans.append(current)
    return [dict(start=float(x[s[0]]), end=float(x[s[-1]]), faces=len(s),
                 max_dx=float(np.diff(x[s]).max()) if len(s)>1 else 0,
                 kind='trailing edge' if x[s[-1]] >= .98 else 'bubble')
            for s in spans if x[s[-1]]-x[s[0]] > .005]


def diagnostics(root):
    rows = []
    for path in sorted(root.glob('*/result.json')):
        path = path.parent
        r, h, cl, cd, steady = values(path)
        q = r['result']
        window = h[h['iteration'] > q['iterations'] - 2000]
        row = dict(run=path.name, lift=cl, drag=cd, settled=steady,
                   failure=q['failure'],
                   iterations=q['iterations'], residual=q['max_residual'],
                   settle_lift=q['settle_lift'], settle_drag=q['settle_drag'],
                   cells=q['cells'], poor_cells=q.get('poor_cells'),
                   peak_footprint_bytes=r.get('peak_footprint_bytes'),
                   drag_range=[float(window['drag'].min()), float(window['drag'].max())],
                   max_velocity_range=[float(window['max_velocity'].min()), float(window['max_velocity'].max())],
                   binary_sha256=r.get('binary_sha256'), mesh_sha256=r.get('mesh_sha256'))
        if (path / 'walls.csv').exists():
            row['reversed_shear'] = bubbles(path)
        if (path / 'cells.csv').exists():
            c = np.genfromtxt(path / 'cells.csv', delimiter=',', names=True)
            speed = np.sqrt(c['u']**2 + c['v']**2 + c['w']**2)
            i = speed.argmax()
            row['fastest_final_cell'] = {k: float(c[i][k]) for k in
                                        ('x', 'y', 'u', 'v', 'p', 'k', 'omega', 'mu_t', 'intermittency')
                                        if k in c.dtype.names}
            row['fastest_final_cell']['speed'] = float(speed[i])
        rows.append(row)
    return rows


def convergence_plots(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    name = 's1223_re2e5_4.76_L2_metal'
    cases = [(root / 'before' / name, 'Original, pressure cutoff .05'),
             (root / 'diagnostics' / (name + '_pressure01'), 'Wall 1x, cutoff .10'),
             (root / 'diagnostics' / (name + '_wall2_pressure01'), 'Wall 2x, cutoff .10'),
             (root / 'diagnostics' / (name + '_wall4_pressure01'), 'Wall 4x, cutoff .10')]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for path, label in cases:
        r, h, _, _, steady = values(path)
        selected = h[h['iteration'] >= 1000]
        for ax, key in zip(axes, ['lift', 'drag']):
            ax.plot(selected['iteration'], selected[key], label=label + ('' if steady else ' *'))
            ax.set_xlabel('Iteration (from 1000)')
            ax.grid(alpha=.2)
    axes[0].set_ylabel('CL')
    axes[1].set_ylabel('CD')
    axes[0].legend(fontsize=7)
    fig.suptitle('Re 200k, 4.76 degrees: local-k transition; * means stopping criteria unmet')
    fig.tight_layout()
    fig.savefig(OUT / 'wall_convergence.png', dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 3.5))
    for suffix, label in [('_wall4', 'Transition, cutoff .05'),
                          ('_wall4_sst', 'Fully turbulent SST, cutoff .05'),
                          ('_wall4_pressure01', 'Transition, cutoff .10')]:
        h = history(root / 'diagnostics' / (name + suffix))
        ax.plot(h['iteration'], h['max_velocity'], label=label)
    ax.set_yscale('log')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Maximum speed / freestream speed')
    ax.set_title('Identical 4x wall mesh: existing pressure fallback control')
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / 'pressure_instability.png', dpi=160)
    plt.close(fig)


def tail_sweep(root):
    """Keep failed tail-mesh runs explicit, without presenting partial means."""
    aero = module('fsae_report', ROOT / 'tools/fsae_aero/report.py')
    rows = []
    lines = ['# Separate eightfold trailing-edge mesh sweep', '',
             'Both models use the molecular omega-production ratio 0.001; gamma also uses harmonic k. Original body spacing, wall-normal resolution and the polygon are retained; outer triangles are regenerated. This is not the source-only before/after comparison.', '',
             '`*` means stopping criteria unmet and a final-2000 mean. A failed solve has no accuracy coefficient; its last finite sample remains in the receipt.', '',
             '| Re | alpha | SST CL | SST CD | gamma CL | gamma CD | gamma dCL | gamma dCD | SST / gamma iterations |',
             '|---|---:|---:|---:|---:|---:|---:|---:|---|']
    for re, points in POINTS.items():
        reference, _ = aero.s1223_reference(re)
        for point in points:
            cr, dr = reference(float(point))
            row = dict(re=re, alpha=float(point), uiuc_cl=cr, uiuc_cd=dr)
            for model, group in [('sst', 'after_sst'), ('gamma', 'after_fluent')]:
                path = root / group / f's1223_re{re/1e5:.0f}e5_{point}_L2_metal'
                r = json.loads((path / 'result.json').read_text())
                q = r['result']
                if q['failure']:
                    cl, cd, steady = None, None, False
                else:
                    _, _, cl, cd, steady = values(path)
                row[model] = dict(cl=cl, cd=cd, settled=steady, failure=q['failure'],
                                  iterations=q['iterations'], residual=q['max_residual'],
                                  settle_lift=q['settle_lift'], settle_drag=q['settle_drag'],
                                  cells=q['cells'], peak_footprint_bytes=r.get('peak_footprint_bytes'),
                                  binary_sha256=r['binary_sha256'], mesh_sha256=r['mesh_sha256'],
                                  dcl=None if cl is None else cl-cr,
                                  dcd_percent=None if cd is None or dr is None else 100*(cd/dr-1))
            if row['sst']['mesh_sha256'] != row['gamma']['mesh_sha256']:
                raise ValueError(f'SST/gamma tail meshes differ at {point}, Re {re}')
            def coefficient(model, key, precision):
                value = row[model][key]
                return 'failed' if value is None else f'{value:.{precision}f}' + ('' if row[model]['settled'] else '*')
            g = row['gamma']
            dcl = 'failed' if g['dcl'] is None else f'{g["dcl"]:+.4f}'
            dcd = 'failed' if g['failure'] else ('n/a' if g['dcd_percent'] is None else f'{g["dcd_percent"]:+.1f}%')
            iterations = lambda m: str(row[m]['iterations']) + (' failed' if row[m]['failure'] else ('' if row[m]['settled'] else '*'))
            lines.append(f'| {re:.0f} | {float(point):g} | {coefficient("sst", "cl", 4)} | {coefficient("sst", "cd", 5)} | '
                         + f'{coefficient("gamma", "cl", 4)} | {coefficient("gamma", "cd", 5)} | {dcl} | {dcd} | '
                         + f'{iterations("sst")} / {iterations("gamma")} |')
            rows.append(row)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'tail_sweep.md').write_text('\n'.join(lines)+'\n')
    (OUT / 'tail_summary.json').write_text(json.dumps(rows, indent=2, allow_nan=False)+'\n')
    print(f'Wrote {len(rows)} paired tail-mesh points')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-root', type=Path, default=RUNS)
    ap.add_argument('--after', default='after_model', help='completed twelve-point transition sweep on the original meshes')
    ap.add_argument('--sst-after', default='after_model_sst', help='completed twelve-point SST control on the original meshes')
    ap.add_argument('--diagnostics-only', action='store_true')
    ap.add_argument('--tail-only', action='store_true', help='summarize the separate paired refined-tail mesh sweep')
    a = ap.parse_args()
    if a.tail_only:
        tail_sweep(a.run_root)
        return
    if a.diagnostics_only:
        OUT.mkdir(parents=True, exist_ok=True)
        rows = diagnostics(a.run_root / 'diagnostics')
        (OUT / 'diagnostics.json').write_text(json.dumps(rows, indent=2) + '\n')
        print(f'Wrote {len(rows)} completed diagnostic runs')
        return
    aero = module('fsae_report', ROOT / 'tools/fsae_aero/report.py')
    lines = ['# S1223 coefficients at every reference angle', '',
             'Lift gate: absolute error <= 0.10 for alpha <= 10 deg. Drag gate: relative error <= 25% wherever UIUC CD exists.', '',
             'After: original meshes, molecular omega-production denominator ratio 0.001; transition also uses harmonic k. Solver defaults are unchanged. The refined-tail sweep is a separate mesh experiment.', '',
             '`*` means forces or residuals did not meet the existing convergence criterion; the value is a last-2000-iteration mean.', '',
             '| Re | alpha | UIUC CL | UIUC CD | SST CL before | SST CD before | SST CL after | SST CD after | gamma CL before | gamma CD before | Fluent-details CL after | Fluent-details CD after | dCL after | dCD after | after gates CL/CD |',
             '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    convergence_lines = ['# Stopping iterations and convergence', '',
                         '`*` means the unchanged residual/force stopping criteria were unmet.', '',
                         '| Re | alpha | SST iterations before / after | gamma iterations before / after | gamma after maximum residual | gamma after CL / CD settling |',
                         '|---|---:|---:|---:|---:|---:|']
    series = []
    for re, points in POINTS.items():
        at, _ = aero.s1223_reference(re)
        for point in points:
            alpha = float(point)
            tag = f's1223_re{re/1e5:.0f}e5_{point}_L2_metal'
            sst = a.run_root / 'before_sst' / tag
            before = a.run_root / 'before' / tag
            for baseline in (sst, before):
                if not (baseline / 'result.json').exists():
                    raise SystemExit(f'Missing required fresh baseline: {baseline}')
            after = a.run_root / a.after / tag
            if not (after / 'result.json').exists():
                raise SystemExit(f'Missing required sweep point: {after}')
            sst_after = a.run_root / a.sst_after / tag
            if not (sst_after / 'result.json').exists():
                raise SystemExit(f'Missing required SST control: {sst_after}')
            sr, _, scl, scd, ss = values(sst)
            ar, _, acl, acd, ast = values(sst_after)
            br, _, bcl, bcd, bs = values(before)
            if not br.get('mesh_sha256') or br.get('mesh_sha256') != sr.get('mesh_sha256'):
                raise SystemExit(f'SST/transition baseline meshes differ at {tag}')
            r, h, cl, cd, settled = values(after)
            if not r.get('mesh_sha256') or r.get('mesh_sha256') != ar.get('mesh_sha256'):
                raise SystemExit(f'SST/transition after meshes differ at {tag}')
            if r['mesh_sha256'] != br['mesh_sha256']:
                raise SystemExit(f'The source-only comparison changed mesh at {tag}')
            cr, dr = at(alpha)
            mark = '' if settled else '*'
            lift_gate = 'trend' if alpha > 10 else ('PASS' if abs(cl-cr) <= .10 else 'FAIL')
            drag_gate = 'n/a' if dr is None else ('PASS' if abs(cd/dr-1) <= .25 else 'FAIL')
            fmt = lambda x, steady, precision: f'{x:.{precision}f}' + ('' if steady else '*')
            lines.append(f'| {re:.0f} | {alpha:g} | {cr:.4f} | ' + ('n/a' if dr is None else f'{dr:.5f}')
                         + f' | {fmt(scl, ss, 4)} | {fmt(scd, ss, 5)} | {fmt(acl, ast, 4)} | {fmt(acd, ast, 5)}'
                         + f' | {fmt(bcl, bs, 4)} | {fmt(bcd, bs, 5)}'
                         + f' | {cl:.4f}{mark} | {cd:.5f}{mark} | {cl-cr:+.4f} | '
                         + ('n/a' if dr is None else f'{100*(cd/dr-1):+.1f}%') + f' | {lift_gate}/{drag_gate} |')
            iters = lambda record, valid: str(record['result']['iterations']) + ('' if valid else '*')
            monitor = lambda v: 'n/a' if v is None else f'{v:.3g}'
            convergence_lines.append(f'| {re:.0f} | {alpha:g} | {iters(sr, ss)} / {iters(ar, ast)} | '
                                     + f'{iters(br, bs)} / {iters(r, settled)} | {monitor(r["result"]["max_residual"])} | '
                                     + f'{monitor(r["result"]["settle_lift"])} / {monitor(r["result"]["settle_drag"])} |')
            convergence_keys = ['iterations', 'max_residual', 'residuals_met', 'settled',
                                'settle_lift', 'settle_drag', 'budget', 'failure']
            convergence = {label: {key: record['result'][key] for key in convergence_keys}
                           for label, record in [('sst_before', sr), ('sst_after', ar),
                                                 ('gamma_before', br), ('gamma_after', r)]}
            series.append(dict(re=re, alpha=alpha, uiuc_cl=cr, uiuc_cd=dr, sst_cl=scl, sst_cd=scd,
                               convergence=convergence, before_cells=br['result']['cells'], after_cells=r['result']['cells'],
                               sst_after_cl=acl, sst_after_cd=acd, sst_after_settled=ast,
                               sst_settled=ss, before_settled=bs,
                               sst_before_source=sr['source'], before_source=br['source'],
                               before_mesh_sha256=br.get('mesh_sha256'),
                               sst_after_source=ar['source'], after_source=r['source'],
                               after_mesh_sha256=r.get('mesh_sha256'),
                               after_binary_sha256=r.get('binary_sha256'),
                               before_cl=bcl, before_cd=bcd, after_cl=cl, after_cd=cd, settled=settled,
                               iterations=r['result']['iterations'], residual=r['result']['max_residual'],
                               settle_lift=r['result']['settle_lift'], settle_drag=r['result']['settle_drag'],
                               last_window_drag_range=[float(h['drag'][h['iteration']>r['result']['iterations']-2000].min()),
                                                       float(h['drag'][h['iteration']>r['result']['iterations']-2000].max())]))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'tables.md').write_text('\n'.join(lines)+'\n')
    (OUT/'convergence.md').write_text('\n'.join(convergence_lines)+'\n')
    (OUT/'summary.json').write_text(json.dumps(series, indent=2)+'\n')
    if series:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex='col')
        for row, re in enumerate(POINTS):
            s = [p for p in series if p['re'] == re]
            for col, key, ylabel in [(0, 'cl', 'CL error'), (1, 'cd', 'CD error (%)')]:
                ax = axes[row, col]
                sub = [p for p in s if p['uiuc_'+key] is not None]
                for model, label in [('sst', 'SST before'), ('sst_after', 'SST after: omega term'),
                                     ('before', 'gamma before: local k'), ('after', 'gamma after: harmonic k + omega term')]:
                    x = [p['alpha'] for p in sub]
                    y = [p[model+'_'+key]-p['uiuc_'+key] if key=='cl' else 100*(p[model+'_'+key]/p['uiuc_'+key]-1) for p in sub]
                    line, = ax.plot(x, y, 'o-', markerfacecolor='white', label=label)
                    state = {'sst': 'sst_settled', 'sst_after': 'sst_after_settled',
                             'before': 'before_settled', 'after': 'settled'}[model]
                    good = [i for i, p in enumerate(sub) if p[state]]
                    ax.scatter([x[i] for i in good], [y[i] for i in good], color=line.get_color(), zorder=3)
                gate = .1 if key=='cl' else 25
                if key=='cl':
                    left, right = ax.get_xlim()
                    ax.fill_between([left, 10], -gate, gate, color='#dcefdc', alpha=.7)
                    ax.set_xlim(left, right)
                else:
                    ax.axhspan(-gate, gate, color='#dcefdc', alpha=.7)
                ax.axhline(0, color='black', linewidth=.6)
                scope = 'CL gate through 10 degrees' if key=='cl' else 'CD gate at all available points'
                ax.set_title(f'Re {re:.0f}: {scope}', fontsize=10)
                ax.set_ylabel(ylabel)
                ax.grid(alpha=.2)
        axes[0, 0].legend(fontsize=9)
        for ax in axes[1]:
            ax.set_xlabel('Angle (degrees)')
        fig.text(.5, .007, 'Filled markers meet stopping rules; hollow markers are unsettled final-2000 means.', ha='center', fontsize=9)
        fig.tight_layout(rect=(0, .03, 1, 1))
        fig.savefig(OUT/'gate_errors.png', dpi=160)
        plt.close(fig)
        convergence_plots(a.run_root)
    print(f'Wrote {len(series)} points to {OUT.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
