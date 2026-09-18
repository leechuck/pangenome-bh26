#!/usr/bin/env python3
"""Persist Slurm state every minute, including dependency and failure details."""
import collections,datetime,json,subprocess,time,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
HOST='hohndor@ilogin.ibex.kaust.edu.sa'
B='/ibex/scratch/projects/c2014/rob/dogohla-benchmark'
cycle=0
while True:
 try:
  launch=json.loads((HERE/'ibex-launch.json').read_text())
  jobs={k:v for k,v in launch['jobs'].items() if not k.endswith('_initial')}
  jobs['recruit']=(HERE/'ibex-recruit-job.txt').read_text().strip()
  jobs['assets']=(HERE/'ibex-assets-job.txt').read_text().strip()
  ids=','.join(jobs.values())
  cmd='sacct -X -P -n -j '+ids+' --format=JobID,JobName,State,ExitCode,Elapsed,NodeList; squeue -h -j '+ids+' -o "QUEUE|%i|%T|%R"'
  r=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20',HOST,cmd],capture_output=True,text=True,timeout=55)
  if r.returncode:raise RuntimeError(r.stderr[-2000:])
  rows=[];queue=[]
  for line in r.stdout.splitlines():
   fields=line.split('|')
   if fields[0]=='QUEUE':queue.append(fields[1:]);continue
   if len(fields)>=6:rows.append(dict(zip(['job','name','state','exit','elapsed','nodes'],fields[:6])))
  counts=dict(collections.Counter(x['state'] for x in rows))
  record=dict(observed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),jobs=jobs,states=counts,records=rows,queue=queue)
  tmp=HERE/'ibex-status.tmp';tmp.write_text(json.dumps(record,indent=2)+'\n');tmp.replace(HERE/'ibex-status.json')
  print(record['observed_at'],counts,flush=True)
  if cycle%3==0:
   subprocess.run([sys.executable,str(HERE/'fetch_ibex.py')],check=True,timeout=110)
   subprocess.run([sys.executable,str(HERE/'score_series.py')],check=True,timeout=120)
  if cycle%10==0:
   subprocess.run([sys.executable,str(HERE/'fetch_ibex.py'),'--gourraud'],check=True,timeout=110)
   subprocess.run([sys.executable,str(HERE/'score_gourraud.py')],check=True,timeout=180)
  cycle+=1
  if rows and not queue and all(x['state'] in ('COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL') for x in rows):break
 except Exception as e:
  print(datetime.datetime.now(datetime.timezone.utc).isoformat(),'MONITOR_ERROR',repr(e),flush=True)
 time.sleep(60)
