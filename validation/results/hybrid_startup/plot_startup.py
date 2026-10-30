"""Render the measured uniform-start pressure RHS trace (matplotlib)."""
import csv,json,pathlib,re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=pathlib.Path(__file__).resolve().parent
def values(name):
 out={}
 for s in (root/name).read_text().splitlines():
  m=re.search(r'STARTUP step=(\d+) kernel=pressure .*rhs_l2=([^\s]+)',s)
  if m: out[int(m[1])]=float(m[2])
 return out
old=values('unbounded_startup_cpu.log');new=values('bounded_startup_cpu.log')
with (root/'startup_rhs.csv').open('w',newline='') as f:
 w=csv.writer(f,lineterminator="\n");w.writerow(['completed_steps','unbounded_rhs_l2','bounded_rhs_l2'])
 for step in sorted(set(old)|set(new)):w.writerow([step,old.get(step,''),new.get(step,'')])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(8.2,4.6),layout='constrained')
ax.semilogy(list(old),list(old.values()),'-o',lw=2,markersize=3,color='#bd3030',label='Original face extrapolation')
ax.semilogy(list(new),list(new.values()),'-o',lw=2,markersize=3,color='#146fa8',label='Bounded face extrapolation')
ax.set(xlabel='Completed outer steps before pressure assembly',ylabel='Pressure-correction RHS L2 norm',xlim=(-.5,23.5),ylim=(1e-5,2e17))
ax.set_title('Uniform bump startup · CPU · 89 × 41 grid',loc='left',weight='bold',pad=12)
ax.grid(axis='y',alpha=.2);ax.legend(loc='upper left',frameon=False)
ax.annotate('Nonfinite solve after step 16',xy=(16,old[16]),xytext=(8,1e12),fontsize=10,color='#a52828',arrowprops={'arrowstyle':'->','color':'#a52828'})
fig.savefig(root/'startup_rhs.png',dpi=180,facecolor='white')
print(json.dumps({'old':{s:old.get(s) for s in [0,1,8,12,16]},'new':{s:new.get(s) for s in [0,1,8,12,16,23]},'matplotlib':matplotlib.__version__},indent=2))
