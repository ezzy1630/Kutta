"""Compare aligned V6 cold-start field exports; no pressure gauge alignment."""
import csv,gzip,json,math,pathlib
root=pathlib.Path(__file__).resolve().parent;result={}
for backend in ['cpu','metal']:
 names=[root/'plateau_hybrid'/f'{backend}-cold-hybrid.csv',root/'plateau_proxy'/f'{backend}-cold-standard-proxy.csv']
 vectors=[]
 for p in names:
  reader=gzip.open(p.with_suffix('.csv.gz'),'rt') if not p.exists() else p.open()
  with reader as f:vectors.append(list(csv.DictReader(f)))
 a,b=vectors;assert len(a)==len(b)==17822
 rows={}
 for name,cols in [('velocity',['u','v','w']),('pressure',['p']),('k',['k']),('omega',['omega'])]:
  diffs=[sum((float(x[c])-float(y[c]))**2 for c in cols) for x,y in zip(a,b)]
  refs=[sum(float(y[c])**2 for c in cols) for y in b]
  rows[name]={'relative_l2':math.sqrt(sum(diffs)/sum(refs)),'max_absolute':math.sqrt(max(diffs)),'relative_max':math.sqrt(max(diffs)/max(refs))}
 dp=[float(x['p'])-float(y['p']) for x,y in zip(a,b)];mean=sum(dp)/len(dp)
 rows['pressure_mean_difference']=mean
 rows['pressure_difference_stddev']=math.sqrt(sum((v-mean)**2 for v in dp)/len(dp))
 result[backend]=rows
(root/'initializer_comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
