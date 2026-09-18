#!/usr/bin/env python3
"""Eight independent four-thread samples per 32-CPU Slurm allocation."""
import concurrent.futures,csv,hashlib,json,os,subprocess,sys,time
from pathlib import Path

def command(plan,index):
    base=Path('/home/leechuck/hla/dogohla-series20260918')
    root=Path(plan['root']);method=plan['method']
    if method=='recruit':return ['bash',str(base/'code/recruit.sh')]
    if method=='t1k4':return [sys.executable,str(base/'series/run_t1k4.py'),str(root),str(index)]
    driver='native_long_control.py' if method=='native_long' else 'dogohla.py'
    args=[sys.executable,str(base/'code'/driver),'--root',str(root),'--threads','4','--cohort-index',str(index)]
    if method=='native':args.append('--native-only')
    return args

def main():
    p=Path(sys.argv[1]);plan=json.loads(p.read_text());chunk=int(os.environ['SLURM_ARRAY_TASK_ID'])
    cohort=Path(plan['root'])/'cohort.tsv'
    assert hashlib.sha256(cohort.read_bytes()).hexdigest()==plan['cohort_sha256'],'Cohort changed'
    tasks=plan['chunks'][chunk];workers=min(8,len(tasks))
    assert int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=4*workers,'Insufficient allocation'
    rows=list(csv.DictReader(cohort.open(),delimiter='\t'))
    for task in tasks:assert rows[task['index']]['donor']==task['donor']
    out=p.parent/'logs'/p.stem;out.mkdir(parents=True,exist_ok=True)
    def run(task):
        started=time.time();log=out/(str(task['index'])+'-'+task['donor']+'.log')
        env=dict(os.environ,SLURM_ARRAY_TASK_ID=str(task['index']))
        if plan['method']=='recruit':env['DONORS']=plan['recruit_manifest']
        with log.open('a') as f:result=subprocess.run(command(plan,task['index']),env=env,stdout=f,stderr=subprocess.STDOUT)
        record=dict(**task,returncode=result.returncode,seconds=time.time()-started,log=str(log))
        print(json.dumps(record),flush=True);return record
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:results=list(pool.map(run,tasks))
    (out/(str(chunk)+'.json')).write_text(json.dumps(results,indent=2)+'\n')
    if any(r['returncode'] for r in results):sys.exit(1)

if __name__=='__main__':main()
