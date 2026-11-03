import hashlib,json,os,subprocess
from pathlib import Path
out=Path('build/speed/overhead-parity');out.mkdir(exist_ok=True)
rows=[]
for backend in ['cpu','metal']:
 hashes={};binaries={}
 for label,folder in [('baseline','baseline'),('candidate','candidate-final')]:
  binary=Path('build/speed')/folder/'bench_overhead';binaries[label]=hashlib.sha256(binary.read_bytes()).hexdigest()
  proc=subprocess.run([str(binary),'parity',backend],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,env=dict(os.environ,MallocSpaceEfficient='1'))
  (out/f'{backend}-{label}.log').write_text(proc.stdout)
  if proc.returncode:raise RuntimeError(proc.stdout[-2000:])
  hashes[label]={}
  for line in proc.stdout.splitlines():
   if line.startswith('PARITY '):
    fields=dict(part.split('=',1) for part in line.split()[1:]);hashes[label][fields['case']]=fields['hash']
  assert len(hashes[label])==4,(backend,label,hashes[label])
 assert hashes['baseline']==hashes['candidate'],(backend,hashes)
 rows.append(dict(backend=backend,identical=True,hashes=hashes,binary_sha256=binaries));print(backend,'4 matching solver hashes',flush=True)
(out/'receipt.json').write_text(json.dumps(rows,indent=2)+'\n')
