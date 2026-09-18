#!/usr/bin/env python3
"""Read compact intermediate-stage progress and failures, not just final outputs."""
import json,shlex
from launch_ibex import remote

SCRIPT=r'''
from pathlib import Path
import collections,datetime,json,subprocess
r=Path('/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/dogohla-series20260918')
report=dict(observed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),stages={},failures=[])
patterns={'reconstruction':'runs/*/A/manifest.json','native':'runs/*/native/manifest.json','native_long':'runs/*/native-long/manifest.json','phasing':'phase/runs/*/AE-ipd365-pairing/manifest.json','no_graph':'final/runs/*/DogoHLA-no-graph/manifest.json','graph':'final/runs/*/DogoHLA/manifest.json'}
for cohort in ('','gourraud'):
 for arm in ('full','hprc','asian_matched'):
  root=r/cohort/'arms'/arm
  for stage,pattern in patterns.items():
   if stage.startswith('native') and arm!='full':continue
   counts=collections.Counter()
   for p in root.glob(pattern):
    try:
     m=json.loads(p.read_text());counts[m['status']]+=1
     if m['status']=='failed':report['failures'].append(dict(path=str(p),error=m.get('error','UNKNOWN')))
    except (OSError,ValueError):counts['unreadable_manifest']+=1
   if counts:report['stages']['/'.join((cohort or 'matched',arm,stage))]=dict(counts)
q=subprocess.check_output(['squeue','-r','-h','-u','hohndor','-t','R','-o','%C'],text=True)
report['running_allocations']=len(q.splitlines());report['allocated_cpus']=sum(map(int,q.split()))
print(json.dumps(report,indent=2))
'''

if __name__=='__main__':print(remote('python3 -c '+shlex.quote(SCRIPT)))
