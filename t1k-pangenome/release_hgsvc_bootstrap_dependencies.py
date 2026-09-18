"""Release only projection tasks whose exact bootstrap parent succeeded."""
import json
import shlex
from pathlib import Path
from launch_development_evidence import B,remote

HERE=Path(__file__).resolve().parent
ledger=json.loads((HERE/'validation-next/EXECUTION_LAUNCH.json').read_text())
parent=ledger['jobs']['bootstrap'];child=ledger['jobs']['bootstrap_evidence']
script='''import subprocess,json,pathlib,datetime
parent=PARENT_LITERAL;child=CHILD_LITERAL
output=pathlib.Path(OUTPUT_LITERAL)
if output.exists():raise RuntimeError('Audit exists; do not repeat blindly')
record=dict(parent=parent,child=child,started=str(datetime.datetime.utcnow()),actions=[])
def save():output.write_text(json.dumps(record,indent=2)+'\\n')
save()
raw=subprocess.check_output(['sacct','-X','-n','-P','-j',parent,'--format=JobID,State,ExitCode'],universal_newlines=True)
parents={r.split('|')[0]:r.split('|')[1:3] for r in raw.splitlines()}
queue=subprocess.check_output(['squeue','-r','-j',child,'-h','-o','%i|%T|%R|%E'],universal_newlines=True)
for row in queue.splitlines():
 job,state,reason,dependency=row.split('|')
 if state!='PENDING' or reason!='(Dependency)':continue
 index=job.rsplit('_',1)[1];source=parent+'_'+index
 if not index.isdigit() or parents.get(source)!=['COMPLETED','0:0']:continue
 current=subprocess.check_output(['scontrol','show','job',job,'-o'],universal_newlines=True)
 values=dict(x.split('=',1) for x in current.split() if '=' in x)
 if values.get('JobState')!='PENDING' or values.get('Dependency')!='aftercorr:'+parent+'_*(unfulfilled)':continue
 item=dict(job=job,parent=source,parent_state='COMPLETED',parent_exit='0:0',previous_dependency=values['Dependency'],status='updating')
 record['actions'].append(item);save()
 result=subprocess.run(['scontrol','update','JobId='+job,'Dependency=0'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
 item.update(status='released' if result.returncode==0 else 'update_failed',returncode=result.returncode,stderr=result.stderr);save()
record['finished']=str(datetime.datetime.utcnow());save()
print(json.dumps(record))
'''.replace('PARENT_LITERAL',repr(parent)).replace('CHILD_LITERAL',repr(child)).replace('OUTPUT_LITERAL',repr(B+'/hla/t1k-pangenome/hgsvc-validation-v1/bootstrap-dependency-release.json'))
record=json.loads(remote('python3 -c '+shlex.quote(script)))
(HERE/'validation-next/BOOTSTRAP_DEPENDENCY_RELEASE.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(dict(checked_successful_parents=len(record['actions']),released=sum(r['status']=='released' for r in record['actions']))))
