#!/usr/bin/env python3
"""Persist Slurm state every minute, including dependency and failure details."""
import collections,datetime,json,subprocess,time,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
HOST='hohndor@ilogin.ibex.kaust.edu.sa'
B='/ibex/scratch/projects/c2014/rob/dogohla-benchmark'
T1K=HERE.parent.parent/'t1k-pangenome'

def development_jobs():
 result={}
 def add(prefix,value):
  if isinstance(value,dict):
   for key,item in value.items():add(prefix+'/'+key,item)
  elif isinstance(value,str) and value.isdigit():
   result['t1k/'+prefix]=value
 for name in ('linear-control-launch.json','genome-control-launch.json',
              'development-map-launch.json','development-map-all-loci-launch.json',
              'native-evidence-development-launch.json','development-evidence-launch.json'):
  path=T1K/name
  if not path.exists():continue
  record=json.loads(path.read_text())
  for key in ('jobs','job','smoke','pilot_job','array_job','reference_job'):
   if key in record:add(name+'/'+key,record[key])
 return result

cycle=0
previous_failures=set()
while True:
 try:
  launch=json.loads((HERE/'ibex-launch.json').read_text())
  jobs={k:v for k,v in launch['jobs'].items() if not k.endswith('_initial')}
  jobs['recruit']=(HERE/'ibex-recruit-job.txt').read_text().strip()
  jobs['assets']=(HERE/'ibex-assets-job.txt').read_text().strip()
  jobs.update(development_jobs())
  ids=','.join(sorted(set(jobs.values())))
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
  progress=subprocess.run([sys.executable,str(HERE/'inspect_ibex.py')],capture_output=True,text=True,check=True,timeout=55)
  pipeline=json.loads(progress.stdout)
  tmp=HERE/'pipeline-status.tmp';tmp.write_text(json.dumps(pipeline,indent=2)+'\n');tmp.replace(HERE/'pipeline-status.json')
  failures={json.dumps(f,sort_keys=True) for f in pipeline['failures']}
  for failure in sorted(failures-previous_failures):print('NEW_PIPELINE_FAILURE',failure,flush=True)
  previous_failures=failures
  terminal=bool(rows) and not queue and all(x['state'].split()[0] in ('COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL') for x in rows)
  if cycle%3==0 or terminal:
   subprocess.run([sys.executable,str(HERE/'fetch_ibex.py')],check=True,timeout=110)
   subprocess.run([sys.executable,str(HERE/'score_series.py')],check=True,timeout=120)
  if cycle%3==0 or terminal:
   subprocess.run([sys.executable,str(HERE/'fetch_ibex.py'),'--gourraud'],check=True,timeout=110)
   subprocess.run([sys.executable,str(HERE/'score_gourraud.py')],check=True,timeout=180)
   if any(key.startswith('t1k_preparation/linear_') for key in pipeline['stages']):
    subprocess.run([sys.executable,str(T1K/'fetch_linear_controls.py')],check=True,timeout=180)
  cycle+=1
  if terminal:break
 except Exception as e:
  print(datetime.datetime.now(datetime.timezone.utc).isoformat(),'MONITOR_ERROR',repr(e),flush=True)
 time.sleep(60)
