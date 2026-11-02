#!/usr/bin/env python3
"""Run the three 40-iteration CLI parity cases on CPU and Metal.

Full fields and converted meshes go to --out (use an ignored runs directory).
Only hashes are printed. Residual normalization removes # lines and the timing
column while retaining every other byte, including whitespace and line endings.
"""
import argparse
import hashlib
import json
import os
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[3]
FILES = ['fluxes.csv', 'forces.csv', 'loads.json', 'result.vtu', 'residuals.csv']


def normalized(path):
    data = path.read_bytes()
    if path.name != 'residuals.csv':
        return data
    lines = data.splitlines(keepends=True)
    header = next(line for line in lines if line.strip() and not line.startswith(b'#'))
    column = header.rstrip(b'\r\n').split(b',').index(b'seconds_per_iteration')
    kept = []
    for line in lines:
        if line.startswith(b'#'):
            continue
        if not line.strip():
            kept.append(line)
            continue
        body = line.rstrip(b'\r\n')
        ending = line[len(body):]
        row = body.split(b',')  # numeric CSV; no quoted comma fields
        kept.append(b','.join(row[:column] + row[column + 1:]) + ending)
    return b''.join(kept)


def execute(command, log):
    result = subprocess.run(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log.write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f'{command[0]} exited {result.returncode}: {result.stdout[-2000:]}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=pathlib.Path, required=True)
    parser.add_argument('--candidate', type=pathlib.Path, required=True)
    parser.add_argument('--out', type=pathlib.Path, required=True)
    args = parser.parse_args()
    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=False)
    examples = pathlib.Path(os.environ.get('KUTTA_FLUENT_EXAMPLES', ROOT / 'validation/fluent_examples')) / 'mixing_elbow'
    setup = out / 'S.json'
    setup.write_text('{"v":1,"overrides":{"methods/pressure-velocity-coupling":"simplec"}}\n')
    converted = out / 'simplec.kutta.json'
    execute([str(args.baseline), 'convert', str(examples / 'mixing_elbow.cas.h5'),
             '--setup', str(setup), '--out', str(converted)], out / 'convert.log')
    cases = {'elbow': (examples / 'mixing_elbow.cas.h5', 'wfgc'),
             'param': (examples / 'elbow_param.cas.h5', 'compressibility,ke-options'),
             'simplec': (converted, 'wfgc')}
    for backend in ['cpu', 'metal']:
        for case, (source, approximate) in cases.items():
            hashes, commits, binaries = {}, {}, {}
            for label, binary in [('baseline', args.baseline), ('candidate', args.candidate)]:
                target = out / f'{case}-{backend}-{label}'
                binaries[label] = hashlib.sha256((ROOT / binary).read_bytes()).hexdigest()
                execute([str(binary), 'solve', str(source), '--approximate', approximate,
                         '--backend', backend, '--iterations', '40', '--out', str(target)],
                        out / f'{case}-{backend}-{label}.log')
                hashes[label] = {name: hashlib.sha256(normalized(target / name)).hexdigest() for name in FILES}
                commits[label] = json.loads((target / 'summary.json').read_text())['kutta_commit']
            same = hashes['baseline'] == hashes['candidate']
            print(json.dumps(dict(case=case, backend=backend, iterations=40, identical=same,
                                  source=commits, binary_sha256=binaries, hashes=hashes)), flush=True)
            if not same:
                raise RuntimeError(f'byte parity failed: {case} {backend}')


if __name__ == '__main__':
    main()
