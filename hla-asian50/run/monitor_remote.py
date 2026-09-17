#!/usr/bin/env python3
"""Concise live status plus output presence checks for the submitted jobs."""
import csv,collections,datetime,json,re,subprocess
from pathlib import Path
R=Path(__file__).resolve().parent
jobs={'T1K':20632380,'SpecHLA':20632470,'PanGenie':20632571}
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'jobs':{},'outputs':{}}
s=subprocess.check_output(['sacct','-X','-n','-P','-j',','.join(map(str,[*jobs.values(),20632534])),'--format=JobID%40,State%30,ExitCode'],text=True)
for name,job in jobs.items():
 states=collections.Counter();fail=[]
 for line in s.splitlines():
  f=line.split('|');jid,state,code=f[:3]
  if re.fullmatch(str(job)+r'_\d+',jid):
   states[state]+=1
   if state not in ['COMPLETED','RUNNING','PENDING','CONFIGURING','COMPLETING']:fail.append([jid,state,code])
 out['jobs'][name]={'states':dict(states),'failed':fail,'expected':40}
for line in s.splitlines():
 f=line.split('|')
 if f[0]=='20632534':out['panel_index']={'state':f[1],'exit_code':f[2]}
rows=list(csv.DictReader(open(R/'source/donors.tsv'),delimiter='\t'))
for tool in jobs:
 good=[];missing=[]
 for r in rows:
  d=r['donor']
  files=([R/'t1k'/d/(d+'_genotype.tsv')] if tool=='T1K' else [R/'spechla'/d/'hla.result.txt'] if tool=='SpecHLA' else [R/'pangenie'/d/a/'calls_genotyping.vcf' for a in ['full','hprc']])
  (good if all(f.is_file() and f.stat().st_size>0 for f in files) else missing).append(d)
 out['outputs'][tool]={'donors_with_files':len(good),'missing':missing}
out['index_logs_complete']=sum('total wallclock time:' in p.read_text(errors='replace') for p in (R/'panels').glob('fold*/*/index.log'))
out['all_complete']=all(out['jobs'][t]['states'].get('COMPLETED',0)==40 and out['outputs'][t]['donors_with_files']==40 for t in jobs)
print(json.dumps(out,indent=2))
