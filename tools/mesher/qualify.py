#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Run a reproducible mesher-to-Metal qualification on generated public STL.

Input directory comes from generate_geometry.py. Keeps outputs locally; writes
one compact receipt, retaining nonzero exits and unconverged results. It never
accepts a coefficient from a failed or unconverged solve as an accuracy result.
Run with --backend cpu for a portable diagnostic; Metal is the acceptance path.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time


def sha(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()


def run(command, log):
    start=time.monotonic()
    with log.open('w') as f:
        result=subprocess.run([str(a) for a in command],stdout=f,stderr=subprocess.STDOUT)
    return {'exit':result.returncode,'wall_seconds':time.monotonic()-start}


def qualify(args, name):
    out=args.out/name
    out.mkdir(parents=True,exist_ok=True)
    mesh_dir=out/'mesh'
    receipt={'geometry':name,'backend':args.backend,'binary_sha256':sha(args.kutta),
             'stl_sha256':sha(args.geometry/f'{name}.stl')}
    receipt['mesh']=run([args.kutta,'mesh',args.geometry/f'{name}.stl','--setup',args.geometry/f'{name}-coarse.json','--out',mesh_dir,'--case','--json','--force'],out/'mesh.log')
    if receipt['mesh']['exit']:
        return receipt
    report=json.loads((mesh_dir/'mesh-report.json').read_text())
    report.pop('out',None)
    receipt['mesh_report']=report
    receipt['mesh_sha256']=sha(mesh_dir/'mesh.msh')
    # Exercise the separate new-case path, then state the benchmark settings.
    case=mesh_dir/'study.kutta.json'
    speed=40 if name=='ahmed' else 20
    receipt['new_case']=run([args.kutta,'new-case',mesh_dir/'mesh.msh','--speed',speed,'--out',case,'--force','--json'],out/'new-case.log')
    if receipt['new_case']['exit']:
        return receipt
    c=json.loads(case.read_text())
    c['methods']['warped_face_gradient_correction']=False  # all mesh faces planar
    if name=='ahmed':
        c['reference_values'].update(area=.112032,length=1.044)
        c['boundaries']['ground']={'type':'wall'}
    case.write_text(json.dumps(c,indent=2)+'\n')
    receipt['case_settings']=c
    controls={'monitors/check-convergence':True,'monitors/condition':'all'}
    for report in ('drag','downforce'):
        controls.update({f'monitors/report/{report}/active':True,
                         f'monitors/report/{report}/previous-values':50,
                         f'monitors/report/{report}/ignore-before':500,
                         f'monitors/report/{report}/stop-criterion':.0005})
    setup=out/'convergence.json'
    setup.write_text(json.dumps({'v':1,'overrides':controls},indent=2)+'\n')
    receipt['solve']=run([args.kutta,'solve',case,'--backend',args.backend,'--iterations',args.iterations,'--setup',setup,'--out',out/'run','--progress','jsonl'],out/'solve.log')
    events=[]
    text=(out/'solve.log').read_text()
    readings=re.findall(r'kutta-memory .*?current=(\d+) peak=(\d+) heap=\d+ device=(\d+)',text)
    if readings:
        device=max(int(r[2]) for r in readings)
        receipt['memory']={'max_observed_device_bytes':device,
                           'device_bytes_per_cell':device/receipt['mesh_report']['cells'],
                           'max_observed_footprint_bytes':max(int(r[1]) for r in readings),
                           'scope':'observed phase marks, not guaranteed peak; Metal unified memory'}
    for line in text.splitlines():
        if line.startswith('{'):
            events.append(json.loads(line))
    iterations=[e for e in events if e.get('event')=='iteration']
    receipt['last_monitor']=iterations[-1] if iterations else None
    receipt['terminal_event']=next((e for e in reversed(events) if e.get('event') in ('done','error')),None)
    receipt['accepted']=bool(receipt['solve']['exit']==0 and receipt['terminal_event'] and receipt['terminal_event'].get('converged'))
    if (out/'run'/'summary.json').exists():
        summary=json.loads((out/'run'/'summary.json').read_text())
        receipt['summary']={k:v for k,v in summary.items() if k not in ('case_file','mesh_file','data_file')}
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kutta',type=Path,required=True)
    p.add_argument('--geometry',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--kind',choices=['both','ahmed','fsae'],default='both')
    p.add_argument('--backend',choices=['cpu','metal','cuda'],default='metal')
    p.add_argument('--iterations',type=int,default=5000)
    args=p.parse_args()
    args.kutta=args.kutta.resolve()
    if args.iterations<=0:
        p.error('--iterations must be positive')
    results=[qualify(args,n) for n in (('ahmed','fsae') if args.kind=='both' else (args.kind,))]
    (args.out/'receipt.json').write_text(json.dumps({'v':1,'results':results},indent=2)+'\n')
    print(json.dumps({'receipt':str(args.out/'receipt.json'),'accepted':all(r.get('accepted',False) for r in results)}))
    raise SystemExit(0 if all(r.get('accepted',False) for r in results) else 1)


if __name__=='__main__':
    main()
