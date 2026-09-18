"""Submit remaining frozen stages as queue slots become free; never read calls."""
import fcntl
import json
import subprocess
from launch_development_evidence import HERE,remote,submit


def completed(job,start,end):
    output=remote('sacct -X -n -P -j '+job+' --format=JobID,State,ExitCode')
    observed={}
    for line in output.splitlines():
        identifier,state,code,*_=line.split('|')
        if identifier.startswith(job+'_') and identifier[len(job)+1:].isdigit():
            observed[int(identifier[len(job)+1:])]=(state,code)
    return all(observed.get(i)==('COMPLETED','0:0') for i in range(start,end+1))


def advance():
    ledger=HERE/'validation-graph-launch.json'
    if not ledger.exists():return
    with (HERE/'work/validation-dispatch.lock').open('w') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return
        record=json.loads(ledger.read_text());jobs=record['jobs']['cohort']
        if 'refine' in jobs:return
        # Unknown submission outcomes require reconciliation, never blind retries.
        if record.get('submission_in_progress'):
            raise RuntimeError('Reconcile uncertain submission before retry: '+str(record['submission_in_progress']))
        stages=(('mapping','calibration',1,83,16,1343,4,'12G'),
                ('evidence','mapping',16,1343,16,1343,1,'8G'),
                ('join','evidence',16,1343,2,167,1,'8G'),
                ('refine','join',2,167,2,167,1,'12G'))
        for name,parent,first,last,start,end,cpus,mem in stages:
            if name in jobs:continue
            if parent not in jobs or not completed(jobs[parent],first,last):return
            if name=='join' and not completed(jobs['native'],1,83):return
            code='/home/leechuck/hla/t1k-pangenome/validation-graph-v1-code'
            command=['python3',code+'/run_validation_graph.py','--assets',code,
                     '--freeze-sha256',record['freeze_sha256'],'--index','$SLURM_ARRAY_TASK_ID','--stage',name]
            record['submission_in_progress']=name
            ledger.write_text(json.dumps(record,indent=2)+'\n')
            try:
                job=submit('t1k-val-cohort-'+name,cpus,mem,'02:00:00',command,
                           ['--array='+str(start)+'-'+str(end)+'%300','--partition=batch,debug'])
            except subprocess.CalledProcessError as error:
                if 'QOSMaxSubmitJobPerUserLimit' in (error.stderr or ''):
                    record.pop('submission_in_progress')
                    record['last_submission_delay']='Slurm QOSMaxSubmitJobPerUserLimit; retry after queue capacity frees'
                    ledger.write_text(json.dumps(record,indent=2)+'\n')
                    return
                raise
            jobs[name]=job;record.pop('submission_in_progress')
            record['dispatch_policy']='Remaining stages submitted after verified predecessor completion to respect per-user queued-job cap'
            ledger.write_text(json.dumps(record,indent=2)+'\n')
            print(name,job,flush=True)
            return


if __name__=='__main__':
    (HERE/'work').mkdir(exist_ok=True)
    advance()
