"""Isolate pivot-scan setup cost after both binaries already have cached face jumps."""
import argparse,datetime,hashlib,importlib.util,json,os,re,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--pairs',type=int,default=3);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
spec=importlib.util.spec_from_file_location('r','validation/results/fsae_time_to_solution/reproduce.py');r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
meshroot=Path(os.environ['KUTTA_RUNS'])/'fsae_aero/mesh'
with (a.out/'receipts.jsonl').open('w') as receipt:
 for rep in range(a.pairs):
  order=[('before',Path('build/speed/candidate/bench_fsae_aero')),('after',Path('build/speed/candidate-final/bench_fsae_aero'))]
  if rep%2:order.reverse()
  for case in ['ahmed-L1','wing-L2']:
   mesh,settings=r.CASES[case]
   for label,binary in order:
    cmd=['/usr/bin/time','-l',str(binary),'--mesh',str(meshroot/mesh),*settings,'--iterations','1','--check','50']
    started=datetime.datetime.now(datetime.timezone.utc).isoformat();log=a.out/f'{rep}-{case}-{label}.log'
    with log.open('w') as stream:proc=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT,env=dict(os.environ,MallocSpaceEfficient='1'))
    text=log.read_text()
    if proc.returncode:raise RuntimeError(text[-2000:])
    result=next(json.loads(l[7:]) for l in text.splitlines() if l.startswith('RESULT '))
    assert not result['result']['failure'] and result['backend']=='metal'
    row=dict(rep=rep,case=case,label=label,started_utc=started,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),**result['result'])
    receipt.write(json.dumps(row)+'\n');receipt.flush();print(case,label,rep,'setup_s',row['setup_s'],flush=True)
