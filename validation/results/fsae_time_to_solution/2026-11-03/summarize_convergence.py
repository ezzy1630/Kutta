"""Summarize exact frozen-bin pairs without treating unfinished work as a zero."""
import argparse,collections,hashlib,json,statistics
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--final',action='store_true');a=p.parse_args()
rows=[json.loads(line) for line in (a.root/'convergence.jsonl').read_text().splitlines()]
def dist(v):return {'n':len(v),'median':statistics.median(v),'min':min(v),'max':max(v)}
by_pair=collections.defaultdict(dict);by_group=collections.defaultdict(list)
for row in rows:
 key=(row['case'],row['scheme'],row['rep']);assert row['label'] not in by_pair[key]
 by_pair[key][row['label']]=row;by_group[(row['case'],row['scheme'],row['label'])].append(row)
complete=collections.defaultdict(list);pending=[]
for (case,scheme,rep),pair in sorted(by_pair.items()):
 if set(pair)!={'baseline','candidate'}:pending.append({'case':case,'scheme':scheme,'rep':rep,'present':list(pair)});continue
 b,c=pair['baseline'],pair['candidate'];histories={}
 for label in ['baseline','candidate']:
  path=a.root/f'{rep}-{case}-{scheme}-{label}.log'
  lines=[line for line in path.read_text().splitlines() if line.startswith('fsae_aero it=')]
  assert lines,'empty history'
  histories[label]={'checkpoints':len(lines),'sha256':hashlib.sha256('\n'.join(lines).encode()).hexdigest()}
 keys=['lift','drag','drag_pressure','max_residual','settle_lift','settle_drag','yplus_max','yplus_mean']
 assert b['iterations']==c['iterations'],(case,scheme,'stopping iteration changed')
 assert all(b[key]==c[key] for key in keys),(case,scheme,'physical values changed')
 assert histories['baseline']==histories['candidate'],(case,scheme,'checkpoint history changed')
 both=all(row['residuals_met'] and row['settled'] for row in pair.values())
 if both:
  for row in pair.values():assert row['max_residual']<=1e-3 and row['settle_lift']<=5e-4 and row['settle_drag']<=5e-4
 complete[(case,scheme)].append({'rep':rep,'iterations':b['iterations'],'both_converged':both,'final_values_identical':True,'checkpoint_histories_identical':True,'histories':histories,'solve_speedup':b['solve_s']/c['solve_s'],'process_speedup':b['process_real_s']/c['process_real_s']})
summary=[]
for (case,scheme),pairs in sorted(complete.items()):
 row={'case':case,'scheme':scheme,'scope':'converged_answer' if all(x['both_converged'] for x in pairs) else 'unconverged_budget','pairs':pairs,'solve_speedup':dist([x['solve_speedup'] for x in pairs]),'process_speedup':dist([x['process_speedup'] for x in pairs])}
 for label in ['baseline','candidate']:
  values=by_group[(case,scheme,label)];row[label]={'iterations':sorted(set(x['iterations'] for x in values)),'solve_seconds':dist([x['solve_s'] for x in values]),'process_seconds':dist([x['process_real_s'] for x in values]),'setup_seconds':dist([x['setup_s'] for x in values]),'peak_bytes':max(x['peak_bytes'] for x in values),'final_values':{key:values[0][key] for key in keys}}
 if a.final:assert len(pairs)==(3 if row['scope']=='converged_answer' else 1)
 summary.append(row)
if a.final:assert len(summary)==4 and not pending
result={'groups':summary,'pending_pairs':pending,'completed_runs':len(rows)};a.out.write_text(json.dumps(result,indent=2)+'\n')
for row in summary:print(row['case'],row['scheme'],row['scope'],row['baseline']['iterations'],'pairs',len(row['pairs']),'solve speedup',row['solve_speedup'])
if pending:print('Pending pairs:',pending)
