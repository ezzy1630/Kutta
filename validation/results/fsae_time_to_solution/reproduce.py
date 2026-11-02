#!/usr/bin/env python3
"""Alternate two frozen Metal bench_fsae_aero binaries on the same local meshes.

Run from the checkout root. Meshes and full fields stay out of git. This records
BSD time's process peak footprint as well as the bench's solve-only wall time.
"""
import argparse
import datetime
import hashlib
import json
import os
import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[3]
CASES = {
    'wing-L2': ('zerihan_h0.224_L2.msh', ['--re', '4.62e5', '--downforce', '--moving-ground', '--tu', '0.001', '--ratio', '1']),
    'ahmed-L1': ('ahmed35_L1.msh', ['--re', '2.78e6', '--length', '1.044', '--area', '0.112', '--walls', 'body', '--lift-z', '--tu', '0.01']),
    'ahmed-L3': ('ahmed35_L3.msh', ['--re', '2.78e6', '--length', '1.044', '--area', '0.112', '--walls', 'body', '--lift-z', '--tu', '0.01']),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=pathlib.Path, required=True)
    parser.add_argument('--candidate', type=pathlib.Path, required=True)
    parser.add_argument('--mesh-dir', type=pathlib.Path, default=pathlib.Path('runs/fsae-time'))
    parser.add_argument('--out', type=pathlib.Path, required=True)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--schemes', nargs='+', choices=['simplec', 'coupled'], default=['simplec', 'coupled'])
    parser.add_argument('--iterations', type=int, default=200)
    parser.add_argument('--pairs', type=int, default=4)
    parser.add_argument('--convergence', action='store_true',
                        help='Use the unchanged 6000-iteration limit; do not repeat unconverged pairs')
    args = parser.parse_args()
    if args.iterations <= 0 or args.pairs <= 0:
        parser.error('iterations and pairs must be positive')
    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=False)
    skipped = set()
    iterations = 6000 if args.convergence else args.iterations
    name = 'convergence.jsonl' if args.convergence else 'measurements.jsonl'
    with (out / name).open('w') as records:
        for rep in range(args.pairs):
            order = [('baseline', args.baseline), ('candidate', args.candidate)]
            if rep % 2:
                order.reverse()
            for case in args.cases:
                mesh, settings = CASES[case]
                for scheme in args.schemes:
                    if (case, scheme) in skipped:
                        continue
                    for label, binary in order:
                        tag = f'{rep}-{case}-{scheme}-{label}'
                        cmd = ['/usr/bin/time', '-l', str(binary), '--mesh', str(args.mesh_dir / mesh),
                               *settings, '--iterations', str(iterations), '--check', '50']
                        if scheme == 'coupled':
                            cmd.append('--coupled')
                        started = datetime.datetime.now(datetime.timezone.utc).isoformat()
                        proc = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                              text=True, env=dict(os.environ, MallocSpaceEfficient='1'))
                        (out / f'{tag}.log').write_text(proc.stdout)
                        if proc.returncode:
                            raise RuntimeError(f'{tag} exited {proc.returncode}: {proc.stdout[-2000:]}')
                        result = next(json.loads(line[7:]) for line in proc.stdout.splitlines() if line.startswith('RESULT '))
                        if result['backend'] != 'metal' or result['result']['failure']:
                            raise RuntimeError(f'{tag}: unexpected backend or failed solve: {result}')
                        peak = re.search(r'(\d+)\s+peak memory footprint', proc.stdout)
                        real = re.search(r'^\s*([0-9.]+) real\s+', proc.stdout, re.MULTILINE)
                        row = dict(rep=rep, case=case, scheme=scheme, label=label, command=cmd,
                                   run_started_utc=started, device=result['device'], backend=result['backend'],
                                   source=result['source'], returncode=proc.returncode,
                                   binary_sha256=hashlib.sha256((ROOT / binary).read_bytes()).hexdigest(),
                                   mesh_sha256=hashlib.sha256((ROOT / args.mesh_dir / mesh).read_bytes()).hexdigest(),
                                   process_real_s=float(real[1]) if real else None,
                                   peak_bytes=int(peak[1]) if peak else None, **result['result'])
                        if args.convergence and not (row['residuals_met'] and row['settled']):
                            skipped.add((case, scheme))
                        records.write(json.dumps(row) + '\n')
                        records.flush()
                        print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
