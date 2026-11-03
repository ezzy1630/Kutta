import argparse,collections,csv,importlib.util,json,re,statistics
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
spec=importlib.util.spec_from_file_location('old','validation/results/fsae_time_to_solution/analyze.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
def distribution(values):return dict(n=len(values),median=statistics.median(values),min=min(values),max=max(values))
def family(name):
 if name=='runtime_sequence':return 'Small-operation replay (mixed)'
 if name=='fv_face_jump':return 'Gradients, limiter and face reconstruction'
 return old.family(name.removeprefix('runtime_'))
receipts=[json.loads(line) for line in (a.root/'receipts.jsonl').read_text().splitlines()]
groups=collections.defaultdict(list);counters=set()
for receipt in receipts:
 tag=f"{receipt['rep']}-{receipt['case']}-{receipt['scheme']}-{receipt['label']}"
 kernels={};phases=collections.defaultdict(float)
 for line in (a.root/f'{tag}.log').read_text().splitlines():
  if line.startswith('GPU_COUNTER '):counters.add(line[len('GPU_COUNTER '):])
  match=re.match(r'GPU_PROFILE (\S+) calls=(\d+) seconds=([\d.e+-]+) logical_bytes=(\d+)',line)
  if not match:continue
  name=match[1];calls=int(match[2]);seconds=float(match[3]);bytes_=int(match[4]);window=receipt['window']
  ms=1000*seconds/window;phases[family(name)]+=ms
  kernels[name]=dict(calls_per_outer=calls/window,ms_per_outer=ms,argument_bytes_per_outer=bytes_/window,argument_footprint_GB_s=bytes_/seconds/1e9 if bytes_ else None)
 assert kernels,'no profile window captured'
 groups[(receipt['case'],receipt['scheme'],receipt['label'])].append(dict(rep=receipt['rep'],kernels=kernels,phases=dict(phases)))
result=[];csvrows=[]
for (case,scheme,label),samples in sorted(groups.items()):
 names=set().union(*(s['kernels'] for s in samples));phases=set().union(*(s['phases'] for s in samples));summary={}
 for name in sorted(names):
  # Zero is an absent dispatch in a completed window, never missing evidence.
  values=[s['kernels'].get(name,dict(calls_per_outer=0,ms_per_outer=0,argument_bytes_per_outer=0,argument_footprint_GB_s=None)) for s in samples]
  summary[name]={key:distribution([v[key] for v in values]) for key in ['calls_per_outer','ms_per_outer','argument_bytes_per_outer']}
  bw=[v['argument_footprint_GB_s'] for v in values if v['argument_footprint_GB_s'] is not None]
  summary[name]['argument_footprint_GB_s']=distribution(bw) if bw else None
  row=dict(case=case,scheme=scheme,label=label,kernel=name,samples=len(samples))
  for metric in ['calls_per_outer','ms_per_outer']:
   for stat in ['median','min','max']:row[f'{metric}_{stat}']=summary[name][metric][stat]
  for stat in ['median','min','max']:row[f'argument_footprint_GB_s_{stat}']=summary[name]['argument_footprint_GB_s'][stat] if bw else ''
  csvrows.append(row)
 result.append(dict(case=case,scheme=scheme,label=label,samples=len(samples),kernels=summary,phases={name:distribution([s['phases'].get(name,0) for s in samples]) for name in sorted(phases)},total_timestamp_ms_per_outer=distribution([sum(s['phases'].values()) for s in samples])))
a.out.write_text(json.dumps(dict(counter_inventory=sorted(counters),groups=result,bandwidth_note='Whole buffer argument footprint divided by command-buffer GPU timestamp span. Not measured DRAM traffic, achieved DRAM bandwidth, or a traffic bound. Replay traffic is not modeled.'),indent=2)+'\n')
with a.out.with_suffix('.csv').open('w',newline='') as stream:
 writer=csv.DictWriter(stream,fieldnames=list(csvrows[0]));writer.writeheader();writer.writerows(csvrows)
print('Profile groups:',len(result),'counter inventory:',sorted(counters))
