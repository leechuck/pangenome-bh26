#!/usr/bin/env python3
"""Remove sequence-free IPD tombstones from binning indices before inference starts."""
import concurrent.futures,json,os,subprocess,sys
from pathlib import Path
R=Path('/home/leechuck/hla/dogohla-series20260918');sys.path.insert(0,str(R/'code'))
from build_refs import read_fasta,write_fasta,V2
from experiment import sha,json_write
paths=[R/'env/share/spechla/db/ref'/V2]+[R/'arms'/a/'db/fold0/ref'/V2 for a in ('full','hprc','asian_matched')]
def repair(p):
 records=read_fasta(p);kept=[(n,s) for n,s in records if s];removed=[n for n,s in records if not s]
 assert kept and removed, 'Unexpected repair input; inspect before rerunning'
 backup=p.with_name(p.name+'.with-empty-records');p.rename(backup)
 write_fasta(p,kept)
 assert all(s for _,s in read_fasta(p))
 with p.with_suffix('.repair.log').open('w') as log:
  subprocess.run(['bowtie2-build','--threads','2','-q',str(p),str(p)],stdout=log,stderr=subprocess.STDOUT,check=True)
 return dict(path=str(p),before_sha256=sha(backup),after_sha256=sha(p),removed_empty_names=removed,retained_records=len(kept))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(repair,paths))
prepared=json.loads((R/'PREPARED.json').read_text());prepared['ignored_empty_records']=records[0]['removed_empty_names'];prepared['database_files']['ref/'+V2]=sha(paths[0]);json_write(R/'PREPARED.json',prepared)
for arm,record in zip(('full','hprc','asian_matched'),records[1:]):
 root=R/'arms'/arm
 p=root/'db/build_log.json';data=json.loads(p.read_text())
 for row in data:
  row['binning_db']['imgt_lite']-=len(record['removed_empty_names']);row['binning_db']['total']-=len(record['removed_empty_names'])
  row['binning_db']['ignored_empty_ipd_records']=len(record['removed_empty_names'])
 json_write(p,data)
 phase=root/'phase/PREPARED.json';data=json.loads(phase.read_text());data['exclusion_manifest_sha256']=sha(p);json_write(phase,data)
json_write(R/'EMPTY_RECORD_REPAIR.json',dict(reason='Sequence-free IPD records generate invalid zero-length SAM @SQ headers; no sequence-bearing candidates removed',records=records))
print('EMPTY_RECORD_REPAIR_COMPLETE',flush=True)
