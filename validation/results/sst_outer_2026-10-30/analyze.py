#!/usr/bin/env python3
"""Summarize the preserved FV profiler logs; no solver or acceptance changes."""
import csv
import json
from pathlib import Path
from statistics import median

root=Path(__file__).resolve().parent
records=json.loads((root/'measurements.json').read_text())
summary=[]
for stage in ['baseline','batching','fused','amg2','final','highre-baseline','highre-final']:
    for backend in ['metal','cpu']:
        runs=[r for r in records if r['stage']==stage and r['backend']==backend and r['profile']=='0']
        if not runs:
            continue
        seconds=[float(r['wall_s_per_outer']) for r in runs]
        seconds_median=median(seconds)
        summary.append(dict(stage=stage,backend=backend,samples=len(runs),
            median_ms=1000*seconds_median,min_ms=1000*min(seconds),max_ms=1000*max(seconds),
            outer_per_s=1/seconds_median))

with (root/'throughput.csv').open('w',newline='') as out:
    writer=csv.DictWriter(out,fieldnames=list(summary[0]),lineterminator='\n')
    writer.writeheader();writer.writerows(summary)

phases={}
for row in json.loads((root/'kernels.json').read_text()):
    name=row['kernel']
    group=('SST' if name.startswith('sst_') else 'FV' if name.startswith('fv_') else
        'AMG' if name.startswith('amg_') else 'reductions' if name.startswith('reduce') else
        'runtime' if name.startswith('runtime.') else 'CSR/Krylov/validation')
    key=(row['stage'],row['backend'],group)
    sums=phases.setdefault(key,[0,0])
    sums[0]+=float(row['calls_per_outer'])*float(row['us_per_call'])/1000
    sums[1]+=float(row['logical_bytes_per_cell_per_outer'])
with (root/'kernel_groups.csv').open('w',newline='') as out:
    writer=csv.writer(out,lineterminator="\n");writer.writerow(['stage','backend','kernel_group','profiled_ms','logical_bytes_per_cell'])
    writer.writerows([*key,*values] for key,values in phases.items())

n=983040
final=next(r for r in records if r['stage']=='highre-final' and r['backend']=='metal' and r['profile']=='1')
bytes_outer=n*float(final['all_logical_bytes_per_cell'])
print('stage, backend, samples, median_ms, outer_per_s')
for row in summary:
    print(row['stage'],row['backend'],row['samples'],f"{row['median_ms']:.3f}",f"{row['outer_per_s']:.6f}")
print('Bandwidth-only logical-byte projection; CUDA GPU has not run:')
for card,bandwidth in [('RTX 3080 10GB',760),('RTX 4090',1008)]:
    print(card,f"ideal {1e3*bytes_outer/(bandwidth*1e9):.3f} ms/outer",
        f"{bandwidth*1e9/bytes_outer:.3f} outer/s; 70% bandwidth {0.7*bandwidth*1e9/bytes_outer:.3f} outer/s")
