#!/usr/bin/env python3
"""Transfer the compact custom bundle in independent resumable, checked chunks."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from remote import remote
HERE=Path(__file__).resolve().parent
STAGE=HERE/'transfer-cache';STAGE.mkdir(exist_ok=True)
HOST='hohndor@ilogin.ibex.kaust.edu.sa'
DEST='/ibex/scratch/projects/c2014/rob/dogohla-benchmark/custom-parts'
source='/home/leechuck/hla/deploy-ibex-20260918/custom.tar.gz'
expected=remote('sha256sum '+source,capture_output=True,text=True).stdout.split()[0]
bundle=STAGE/'custom.tar.gz'
if not bundle.exists() or hashlib.sha256(bundle.read_bytes()).hexdigest()!=expected:
 with bundle.open('wb') as f:remote('cat '+source,stdout=f)
assert hashlib.sha256(bundle.read_bytes()).hexdigest()==expected
subprocess.run(['ssh','-o','BatchMode=yes',HOST,'mkdir -p '+DEST],check=True)
blob=bundle.read_bytes();size=4*1024*1024
chunks=[(f'part{i//size:04d}',blob[i:i+size]) for i in range(0,len(blob),size)]
def upload(item):
 name,data=item;digest=hashlib.sha256(data).hexdigest();p=DEST+'/'+name
 check=f'test -f {p} && sha256sum {p}'
 r=subprocess.run(['ssh','-o','BatchMode=yes',HOST,check],capture_output=True,text=True)
 if r.returncode==0 and r.stdout.split()[0]==digest:return name,'already verified'
 for attempt in range(3):
  r=subprocess.run(['ssh','-o','BatchMode=yes',HOST,f'cat > {p}.tmp && mv {p}.tmp {p} && sha256sum {p}'],input=data,capture_output=True)
  if r.returncode==0 and r.stdout.decode().split()[0]==digest:return name,'verified'
 raise RuntimeError(name+' failed checksum/transfer')
with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:
 for result in pool.map(upload,chunks):print(*result,flush=True)
manifest={'sha256':expected,'bytes':len(blob),'parts':len(chunks),'destination':DEST}
(HERE/'custom-transfer-complete.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('ALL_CHUNKS_VERIFIED',flush=True)
