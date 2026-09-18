#!/usr/bin/env python3
"""Read compact intermediate-stage progress and failures, not just final outputs."""
import json,shlex
from launch_ibex import remote

SCRIPT=r'''
from pathlib import Path
import collections,datetime,hashlib,json,subprocess
r=Path('/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/dogohla-series20260918')
report=dict(observed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),stages={},failures=[],final_failures=[])
patterns={'reconstruction':'runs/*/A/manifest.json','native':'runs/*/native/manifest.json','native_long':'runs/*/native-long/manifest.json','phasing':'phase/runs/*/AE-ipd365-pairing/manifest.json','no_graph':'final/runs/*/DogoHLA-no-graph/manifest.json','graph':'final/runs/*/DogoHLA/manifest.json'}
for cohort in ('','gourraud'):
 for arm in ('full','hprc','asian_matched'):
  root=r/cohort/'arms'/arm
  for stage,pattern in patterns.items():
   if stage.startswith('native') and arm!='full':continue
   counts=collections.Counter()
   for p in root.glob(pattern):
    try:
     m=json.loads(p.read_text());state=m['status'];marker=p.parent/'TERMINAL_FAILURE.json'
     if state=='failed' and marker.exists():
      audit=json.loads(marker.read_text())
      if audit.get('manifest_sha256')==hashlib.sha256(p.read_bytes()).hexdigest():state='failed_final'
     counts[state]+=1
     if state=='failed':report['failures'].append(dict(path=str(p),error=m.get('error','UNKNOWN')))
     if state=='failed_final':report['final_failures'].append(dict(path=str(p),error=m.get('error','UNKNOWN')))
    except (OSError,ValueError):counts['unreadable_manifest']+=1
   if counts:report['stages']['/'.join((cohort or 'matched',arm,stage))]=dict(counts)
prep=r.parent/'t1k-pangenome'
for stage,pattern in [('reference_graphs','graphs/build-v1/*/*/manifest.json'),
                      ('reference_graphs_repaired','graphs/build-v2/*/*/manifest.json'),
                      ('graph_pilot','graphs/pilot-v2/*/*/manifest.json'),
                      ('personalization_smoke','smoke/personalization-v[23]/manifest.json'),
                      ('reserved_reads','validation-reads/*/manifest.json'),
                      ('linear_hprc','development/linear-v1/hprc/*/manifest.json'),
                      ('linear_hprc_asian','development/linear-v1/hprc_asian/*/manifest.json')]:
 counts=collections.Counter()
 for p in prep.glob(pattern):
  try:
   m=json.loads(p.read_text());state=m['status']
   if state=='failed' and stage=='reference_graphs':
    repaired=prep/'graphs/build-v2'/p.parent.parent.name/p.parent.name/'COMPLETE.json'
    if repaired.exists():
     audit=json.loads(repaired.read_text())
     if audit.get('status')=='complete' and audit.get('previous_manifest_sha256')==hashlib.sha256(p.read_bytes()).hexdigest():state='repaired'
   counts[state]+=1
   if state=='failed':report['failures'].append(dict(path=str(p),error=m.get('error','UNKNOWN')))
  except (OSError,ValueError):counts['unreadable_manifest']+=1
 if counts:report['stages']['t1k_preparation/'+stage]=dict(counts)
q=subprocess.check_output(['squeue','-r','-h','-u','hohndor','-t','R','-o','%C'],text=True)
report['running_allocations']=len(q.splitlines());report['allocated_cpus']=sum(map(int,q.split()))
print(json.dumps(report,indent=2))
'''

if __name__=='__main__':print(remote('python3 -c '+shlex.quote(SCRIPT)))
