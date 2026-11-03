import argparse,os,subprocess,json,hashlib,datetime,importlib.util
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--pairs',type=int,default=2);p.add_argument('--first',type=int,default=50);p.add_argument('--cases',nargs='+',default=['ahmed-L1','wing-L2']);p.add_argument('--schemes',nargs='+',default=['simplec','coupled']);a=p.parse_args()
a.out.mkdir(parents=True,exist_ok=False)
spec=importlib.util.spec_from_file_location('r','validation/results/fsae_time_to_solution/reproduce.py');r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
meshroot=Path(os.environ['KUTTA_RUNS'])/'fsae_aero/mesh'
with (a.out/'receipts.jsonl').open('w') as receipt:
 for rep in range(a.pairs):
  order=[('baseline',Path('build/speed/profile-baseline')),('candidate',Path('build/speed/profile-candidate'))]
  if rep%2:order.reverse()
  for case in a.cases:
   mesh,settings=r.CASES[case]
   for scheme in a.schemes:
    for label,binary in order:
     tag=f'{rep}-{case}-{scheme}-{label}'
     cmd=['/usr/bin/time','-l',str(binary.resolve()),'--mesh',str(meshroot/mesh),*settings,'--iterations',str(a.first+50),'--check','50']+(['--coupled'] if scheme=='coupled' else [])
     log=a.out/f'{tag}.log';started=datetime.datetime.now(datetime.timezone.utc).isoformat()
     with log.open('w') as stream:
      result=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT,env=dict(os.environ,MallocSpaceEfficient='1',KUTTA_SPEED_PROFILE='1',KUTTA_SPEED_PROFILE_FIRST=str(a.first)))
     text=log.read_text()
     if result.returncode:raise RuntimeError(text[-2000:])
     row=dict(rep=rep,case=case,scheme=scheme,label=label,first=a.first,window=50,start_utc=started,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),result=next(json.loads(l[7:]) for l in text.splitlines() if l.startswith('RESULT ')))
     receipt.write(json.dumps(row)+'\n');receipt.flush();print(tag,'profile done',flush=True)
