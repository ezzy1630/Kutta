"""Plot numeric receipts only; matplotlib 3.11.2, headless, no fixture reads."""
import argparse,json,math,re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--summary',type=Path,required=True);p.add_argument('--logs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
summary=json.loads(a.summary.read_text())
rows=[r for r in summary['groups'] if r['scope']=='converged_answer']
fig,ax=plt.subplots(figsize=(8.5,2.2+.55*len(rows)),layout='constrained')
ax.axvline(1,color='#777',linestyle='--',linewidth=1)
for y,row in enumerate(rows):
 s=row['solve_speedup'];values=[p['solve_speedup'] for p in row['pairs']]
 ax.scatter(values,[y]*len(values),color='#a1a8b1',s=27,zorder=2)
 ax.errorbar(s['median'],y,xerr=[[s['median']-s['min']],[s['max']-s['median']]],fmt='o',color='#2169ad',capsize=5,zorder=3)
ax.set_yticks(range(len(rows)),[f"{r['case']} / {r['scheme'].upper() if r['scheme']=='simplec' else 'Coupled'}\nN={r['baseline']['iterations'][0]}, {len(r['pairs'])} pairs" for r in rows]);ax.invert_yaxis();ax.grid(axis='x',alpha=.18)
ax.set_xlabel('Baseline / candidate solve time; above 1 is faster')
ax.set_title('Time to the unchanged convergence gates\nPaired median and full range on the shared M4 Pro')
ax.spines[['top','right']].set_visible(False)
for suffix in ['png','svg']:fig.savefig(a.out/f'converged-speedup.{suffix}',dpi=180,metadata={'Date':None} if suffix=='svg' else None)
plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(10,6),sharex='col',layout='constrained')
for col,scheme in enumerate(['simplec','coupled']):
 for label,style,color in [('baseline','-','#475569'),('candidate','--','#d05c22')]:
  path=a.logs/f'0-ahmed-L1-{scheme}-{label}.log';values=[]
  for line in path.read_text().splitlines():
   if line.startswith('fsae_aero it='):
    values.append({k:float(v) for k,v in (part.split('=',1) for part in line.split()[1:])})
  for row,key,gate in [(0,'residual',1e-3),(1,'settle',5e-4)]:
   data=[v for v in values if math.isfinite(v[key]) and v[key]>0]
   axes[row,col].plot([v['it'] for v in data],[v[key] for v in data],style,label=label,color=color,linewidth=1.3)
   axes[row,col].set_yscale('log');axes[row,col].grid(alpha=.15)
 axes[0,col].axhline(1e-3,color='#8b2929',linestyle=':',label='unchanged gate')
 axes[1,col].axhline(5e-4,color='#8b2929',linestyle=':')
 axes[0,col].set_title(scheme.upper() if scheme=='simplec' else 'Coupled')
 axes[1,col].set_xlabel('Outer iteration')
axes[0,0].set_ylabel('Maximum scaled residual');axes[1,0].set_ylabel('Maximum force change over 500 iterations')
axes[0,0].legend(frameon=False,fontsize=9)
fig.suptitle('Ahmed L1: residual and force-settling gates are separate')
for suffix in ['png','svg']:fig.savefig(a.out/f'ahmed-gates.{suffix}',dpi=180,metadata={'Date':None} if suffix=='svg' else None)
plt.close(fig)
