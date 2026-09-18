#!/usr/bin/env python3
"""Replace only pending tasks with eight-sample allocations; preserve running work."""
import csv,hashlib,io,json,shlex,subprocess,sys
from launch_ibex import B,HERE,HOST,R,remote,submit

def upload(path,data):
    subprocess.run(['ssh','-o','BatchMode=yes',HOST,'python3 -c '+shlex.quote('import sys;from pathlib import Path;p=Path('+repr(path)+');p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(sys.stdin.buffer.read())')],input=data,check=True)

def queued(job):
    out=remote('squeue -r -h -j '+job+' -o "%i|%T"')
    return {int(i.rsplit('_',1)[1]):s for i,s in (line.split('|') for line in out.splitlines())}

def submit_cohort(key,method,root,arm='full'):
    ledger=HERE/'ibex-launch.json';data=json.loads(ledger.read_text());jobs=data['jobs']
    assert key not in jobs,'Already submitted'
    cohort_text=remote('cat '+shlex.quote(root.split('/arms/')[0].replace('/home/leechuck/hla',B+'/hla')+'/cohort.tsv'))+'\n'
    rows=list(csv.DictReader(io.StringIO(cohort_text),delimiter='\t'))
    tasks=[dict(index=i,donor=row['donor']) for i,row in enumerate(rows)]
    plan=dict(root=root,method=method,cohort_sha256=hashlib.sha256(cohort_text.encode()).hexdigest(),chunks=[tasks[i:i+8] for i in range(0,len(tasks),8)],threads_per_sample=4,cpus_per_allocation=32)
    path=HERE/'packed'/f'{key}.json';path.parent.mkdir(exist_ok=True);path.write_text(json.dumps(plan,indent=2)+'\n')
    upload(B+'/hla/dogohla-series20260918/packed/'+path.name,path.read_bytes())
    job=submit('dogo-packed-'+key,32,'64G','02:00:00',['python3',R+'/series/run_packed.py',R+'/packed/'+path.name],['--array=0-'+str(len(plan['chunks'])-1),'--export=ALL,ARM='+arm])
    jobs[key]=job;ledger.write_text(json.dumps(data,indent=2)+'\n');print(key,job,'samples',len(tasks),'allocations',len(plan['chunks']),flush=True)

def pack(key,method,root,arm='full',exclude=(),dependencies=()):
    ledger=HERE/'ibex-launch.json';data=json.loads(ledger.read_text());jobs=data['jobs']
    old=jobs[key];assert key+'_unpacked' not in jobs,'Already repacked'
    indices=[i for i,s in queued(old).items() if s=='PENDING' and i not in exclude]
    if not indices:print('No pending work',key);return
    spec=old+'_['+','.join(map(str,indices))+']'
    remote('scontrol hold '+shlex.quote(spec))
    # Check again after holding: never cancel a task that started in the race.
    states=queued(old);indices=sorted(i for i in indices if states.get(i)=='PENDING')
    cohort_text=remote('cat '+shlex.quote(root.split('/arms/')[0].replace('/home/leechuck/hla',B+'/hla')+'/cohort.tsv'))+'\n'
    rows=list(csv.DictReader(io.StringIO(cohort_text),delimiter='\t'))
    tasks=[dict(index=i,donor=rows[i]['donor']) for i in indices]
    plan=dict(root=root,method=method,cohort_sha256=hashlib.sha256(cohort_text.encode()).hexdigest(),chunks=[tasks[i:i+8] for i in range(0,len(tasks),8)],replaced_job=old,threads_per_sample=4,cpus_per_allocation=32)
    if method=='recruit':plan['recruit_manifest']=R+'/recruit.gourraud.tsv'
    path=HERE/'packed'/f'{key}.json';path.parent.mkdir(exist_ok=True);path.write_text(json.dumps(plan,indent=2)+'\n')
    upload(B+'/hla/dogohla-series20260918/packed/'+path.name,path.read_bytes())
    # Cancel only the held tasks after persisting their exact replacement plan.
    remote('scancel '+shlex.quote(old+'_['+','.join(map(str,indices))+']'))
    extra=['--array=0-'+str(len(plan['chunks'])-1),'--export=ALL,ARM='+arm]
    if dependencies:extra.append('--dependency=afterok:'+':'.join(dependencies))
    job=submit('dogo-packed-'+key,32,'64G','02:00:00',['python3',R+'/series/run_packed.py',R+'/packed/'+path.name],extra)
    jobs[key+'_unpacked']=old;jobs[key]=job
    if exclude:jobs[key+'_smoke']=old+'_'+str(exclude[0])
    data['scope']='Matched 64 and graph/family-excluded Gourraud 946; no deadline; pending tasks packed eight per 32-CPU allocation.'
    ledger.write_text(json.dumps(data,indent=2)+'\n');print(key,job,'samples',len(tasks),'allocations',len(plan['chunks']),flush=True)

if __name__=='__main__':
    key,method,root=sys.argv[1:4]
    pack(key,method,root)
