#!/usr/bin/env python3
"""Preserve failed smoke outputs and resubmit after the portable make preflight."""
import csv,json,shlex
from launch_ibex import B,HERE,R,remote,submit

def main():
 ledger=HERE/'ibex-launch.json';data=json.loads(ledger.read_text());jobs=data['jobs']
 preflight=(HERE/'make-preflight-job.txt').read_text().strip()
 script='''from pathlib import Path
import shutil
r=Path(ROOT)
for arm in ('full','hprc','asian_matched'):
 for name in (('A','native','native_work') if arm=='full' else ('A',)):
  p=r/'arms'/arm/'runs/HG00658'/name
  if not p.exists():continue
  assert not (p/'COMPLETE').exists(), 'Never replace completed output'
  dest=r/'attempts/missing-make'/arm/name
  assert not dest.exists();dest.parent.mkdir(parents=True,exist_ok=True);shutil.move(p,dest)
p=r/'gourraud/t1k4/NA10846'
assert not (p/'COMPLETE').exists()
dest=r/'attempts/t1k-analyzer-segfault/NA10846';assert not dest.exists();dest.parent.mkdir(parents=True,exist_ok=True);shutil.move(p,dest)
'''.replace('ROOT',repr(B+'/hla/dogohla-series20260918'))
 remote('python3 -c '+shlex.quote(script))
 def save():ledger.write_text(json.dumps(data,indent=2)+'\n')
 for arm in ('native','full','hprc','asian_matched'):
  key=arm+'_smoke';jobs[arm+'_missing_make_initial']=jobs[key]
  panel='full' if arm=='native' else arm
  command=['python3',R+'/code/dogohla.py','--root',R+'/arms/'+panel,'--threads','4','--cohort-index','0']
  if arm=='native':command+=['--native-only']
  jobs[key]=submit('dogo-'+arm+'-make-recovery',4,'8G','02:00:00',command,['--partition=debug','--export=ALL,ARM='+panel,'--dependency=afterok:'+preflight]);save()
  dependents=('native','native_long') if arm=='native' else (arm,)
  for dependent in dependents:remote('scontrol update JobId='+jobs[dependent]+' Dependency=afterok:'+jobs[key])
 rows=list(csv.DictReader((HERE/'cohort.gourraud.tsv').open(),delimiter='\t'));i=next(i for i,r in enumerate(rows) if r['donor']=='NA10846')
 jobs['t1k_analyzer_retry']=submit('dogo-t1k-NA10846-retry',4,'8G','00:30:00',['python3',R+'/series/run_t1k4.py',R+'/gourraud',str(i)],['--partition=debug']);save()
 print({k:v for k,v in jobs.items() if k.endswith('_smoke') or k=='t1k_analyzer_retry'})

if __name__=='__main__':main()
