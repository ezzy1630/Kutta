"""Resume exact frozen binary pairs, keeping each benchmark's checkpoint log live."""
import argparse,datetime,hashlib,importlib.util,json,os,re,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--pairs',type=int,default=3);a=p.parse_args()
spec=importlib.util.spec_from_file_location('r','validation/results/fsae_time_to_solution/reproduce.py');r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
bins={label:Path('build/speed')/folder/'bench_fsae_aero' for label,folder in [('baseline','baseline'),('candidate','candidate-final')]}
hashes={label:hashlib.sha256(binary.read_bytes()).hexdigest() for label,binary in bins.items()}
meshroot=Path(os.environ['KUTTA_RUNS'])/'fsae_aero/mesh'
a.out.mkdir(parents=True,exist_ok=True);ledger=a.out/'convergence.jsonl'
rows=[json.loads(line) for line in ledger.read_text().splitlines()] if ledger.exists() else []
finished={}
for row in rows:
 assert row['binary_sha256']==hashes[row['label']], 'cannot resume measurements from another binary'
 finished[(row['rep'],row['case'],row['scheme'],row['label'])]=row
skipped=set()
with ledger.open('a') as receipt:
 for rep in range(a.pairs):
  order=['baseline','candidate'] if rep%2==0 else ['candidate','baseline']
  for case in ['wing-L2','ahmed-L1']:
   mesh,settings=r.CASES[case];meshpath=meshroot/mesh;meshhash=hashlib.sha256(meshpath.read_bytes()).hexdigest()
   for scheme in ['simplec','coupled']:
    if (case,scheme) in skipped:continue
    for label in order:
     key=(rep,case,scheme,label)
     if key in finished:
      assert finished[key]['mesh_sha256']==meshhash, 'cannot resume another mesh'
      continue
     binary=bins[label];tag=f'{rep}-{case}-{scheme}-{label}'
     cmd=['/usr/bin/time','-l',str(binary),'--mesh',str(meshpath),*settings,'--iterations','6000','--check','50']+(['--coupled'] if scheme=='coupled' else [])
     log=a.out/f'{tag}.log';started=datetime.datetime.now(datetime.timezone.utc).isoformat()
     print('START',tag,started,flush=True)
     with log.open('w') as stream: proc=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT,env=dict(os.environ,MallocSpaceEfficient='1'))
     text=log.read_text()
     if proc.returncode:raise RuntimeError(f'{tag} exit {proc.returncode}: '+text[-2000:])
     result=next(json.loads(line[7:]) for line in text.splitlines() if line.startswith('RESULT '))
     assert result['backend']=='metal' and not result['result']['failure'], result
     peak=re.search(r'(\d+)\s+peak memory footprint',text);real=re.search(r'^\s*([0-9.]+) real\s+',text,re.MULTILINE)
     row=dict(rep=rep,case=case,scheme=scheme,label=label,command=cmd,run_started_utc=started,device=result['device'],backend=result['backend'],source=result['source'],returncode=proc.returncode,binary_sha256=hashes[label],mesh_sha256=meshhash,process_real_s=float(real[1]) if real else None,peak_bytes=int(peak[1]) if peak else None,**result['result'])
     receipt.write(json.dumps(row)+'\n');receipt.flush();finished[key]=row
     print(json.dumps(row),flush=True)
    pair=[finished[(rep,case,scheme,label)] for label in ['baseline','candidate']]
    if not all(row['residuals_met'] and row['settled'] for row in pair):skipped.add((case,scheme))
