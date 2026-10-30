#!/usr/bin/env python3
"""Render the measured wall times. Requires matplotlib; see metadata.json."""
import csv
from pathlib import Path
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
rows=list(csv.DictReader((root/'throughput.csv').open()))
metal=[r for r in rows if r['backend']=='metal' and not r['stage'].startswith('highre')]
cpu=[r for r in rows if r['backend']=='cpu' and r['stage'] in ['baseline','final']]
labels={'baseline':'Baseline','batching':'Reuse + adaptive batches','fused':'Reuse fits + kernel fusion',
    'amg2':'SST AMG: 2 + 2 sweeps','final':'Reuse final residual check'}
fig,axes=plt.subplots(1,2,figsize=(10.2,4.1),gridspec_kw={'width_ratios':[1.8,1]})
for ax,data,title,scale,unit in [(axes[0],metal,'Metal',1,'ms / outer'),(axes[1],cpu,'CPU serial oracle',.001,'s / outer')]:
    values=[float(r['median_ms'])*scale for r in data]
    colors=['#0e766c' if r['stage']=='final' else '#8a9cab' for r in data]
    bars=ax.barh(range(len(data)),values,color=colors,height=.6)
    for bar,r,value in zip(bars,data,values):
        ax.text(value+max(values)*.025,bar.get_y()+bar.get_height()/2,
            f"{value:.2f}" if scale!=1 else f"{value:.1f}",va='center',fontsize=9)
    ax.set_yticks(range(len(data)),[labels[r['stage']] if title=='Metal' else r['stage'].capitalize() for r in data],fontsize=9)
    ax.invert_yaxis();ax.set_xlim(0,max(values)*1.18);ax.set_xlabel(unit,fontsize=10)
    ax.set_title(title,loc='left',fontsize=12,fontweight='bold',pad=12)
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y',length=0);ax.grid(axis='x',alpha=.16);ax.set_axisbelow(True)
fig.suptitle('SST tet channel, Re_H = 1,000: incremental measured cost',fontsize=13,x=.02,ha='left',fontweight='bold')
fig.text(.02,.02,'M4 Pro · Release · same transient window (outer 6–8). Metal: median 3 runs except AMG experiment (1).\nCPU: baseline median 2 runs, final 1 run. Setup/warmup excluded. Channel not converged.',fontsize=8,color='#444444')
fig.tight_layout(rect=(0,.13,1,.93),w_pad=2)
fig.savefig(root/'performance.png',dpi=180,facecolor='white')
print('matplotlib',plt.matplotlib.__version__)

# Higher-Re primary benchmark, independently measured after validation finished.
fig,axes=plt.subplots(1,2,figsize=(8,3.6))
for ax,backend,title,scale,unit in [(axes[0],'metal','Metal',1,'ms / outer'),(axes[1],'cpu','CPU serial oracle',.001,'s / outer')]:
    data=[r for r in rows if r['backend']==backend and r['stage'].startswith('highre')]
    values=[float(r['median_ms'])*scale for r in data]
    bars=ax.barh([0,1],values,color=['#8a9cab','#0e766c'],height=.55)
    for bar,value in zip(bars,values):
        ax.text(value+max(values)*.03,bar.get_y()+bar.get_height()/2,
            f"{value:.1f}" if scale==1 else f"{value:.2f}",va='center',fontsize=11)
    ax.set_yticks([0,1],['Baseline','Final'],fontsize=11);ax.invert_yaxis()
    ax.set_xlim(0,max(values)*1.25);ax.set_xlabel(unit,fontsize=11)
    ax.set_title(title,loc='left',fontsize=12,fontweight='bold',pad=10)
    ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0)
    ax.grid(axis='x',alpha=.16);ax.set_axisbelow(True)
fig.suptitle('983,040-cell SST channel · Re_H = 10,000',fontsize=13,x=.03,ha='left',fontweight='bold')
fig.text(.03,.02,'M4 Pro · outer 6–8 · Metal median 3 runs; CPU 1 run.\nSetup and warmup excluded. Transient case; same-channel convergence remains unmet.',fontsize=8,color='#444444')
fig.tight_layout(rect=(0,.15,1,.93),w_pad=2)
fig.savefig(root/'highre-performance.png',dpi=180,facecolor='white')
