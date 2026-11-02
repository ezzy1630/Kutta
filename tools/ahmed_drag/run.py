#!/usr/bin/env python3
"""One foreground Ahmed run, with mesh hash and a portable reproduction command.

Generate the mesh with tools/fsae_aero/ahmed.py first. Run serially on a shared
machine; only use meshes whose measured/projected footprint fits the budget.
The 5000 iterations, 1000-iteration force window and convergence criteria are
those of the original Ahmed study. There is no wall-clock stopping criterion.
"""
import argparse
import hashlib
import json
import os
import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh', required=True, type=pathlib.Path)
    ap.add_argument('--out', required=True, type=pathlib.Path)
    ap.add_argument('--limiter', choices=['standard', 'differentiable'], default='standard')
    ap.add_argument('--coupled', action='store_true')
    a = ap.parse_args()
    meta = json.loads(pathlib.Path(str(a.mesh) + '.json').read_text())
    a.out.mkdir(parents=True, exist_ok=True)
    command = ['build/bench_fsae_aero', '--mesh', str(a.mesh), '--backend', 'metal',
               '--re', '2.78e6', '--length', '1.044', '--area', '0.112032', '--walls', 'body',
               '--lift-z', '--tu', '0.01', '--window', '1000', '--iterations', '5000', '--cells',
               '--limiter', a.limiter, '--out', str(a.out)]
    if a.coupled:
        command.append('--coupled')
    env = dict(os.environ, MallocSpaceEfficient='1')
    with (a.out / 'solve.log').open('w') as log:
        proc = subprocess.run(['/usr/bin/time', '-l'] + command, cwd=ROOT, env=env,
                              stdout=log, stderr=subprocess.STDOUT)
    log = (a.out / 'solve.log').read_text()
    result = re.search(r'^RESULT (\{.*\})$', log, re.M)
    if not result:
        raise RuntimeError(f'no RESULT record; exit {proc.returncode}; inspect {a.out}/solve.log')
    r = json.loads(result[1])
    footprint = re.search(r'(\d+)\s+peak memory footprint', log)
    if not footprint:
        raise RuntimeError('missing macOS peak memory footprint')
    r.update(command=command, mesh=meta, mesh_sha256=hashlib.sha256(a.mesh.read_bytes()).hexdigest(),
             case='ahmed', point=str(meta['slant_deg']), level=meta['level'],
             backend='metal', limiter=a.limiter, coupled=a.coupled, exit=proc.returncode,
             peak_footprint_bytes=int(footprint[1]))
    (a.out / 'result.json').write_text(json.dumps(r, indent=2) + '\n')
    print(json.dumps(r, indent=2))
    if proc.returncode:
        raise SystemExit(proc.returncode)
    if r['peak_footprint_bytes'] > 3e9:
        raise RuntimeError('run exceeded the 3 GB study budget; do not refine this family further')


if __name__ == '__main__':
    main()
