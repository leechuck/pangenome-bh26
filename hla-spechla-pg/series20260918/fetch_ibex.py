#!/usr/bin/env python3
"""Incrementally fetch immutable completed outputs, retaining hash validation."""
import io,json,subprocess,tarfile,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from experiment import completed,sha


def verified_cache(out):
    """Only reuse snapshots whose configuration and every output still validate."""
    known={}
    for marker in out.rglob('COMPLETE'):
        p=marker.parent
        if out.name=='snapshot' and p.relative_to(out).parts[0]=='gourraud':continue
        try:
            m=json.loads((p/'manifest.json').read_text())
            if m['status']=='complete' and completed(p,m['configuration']):known[str(p.relative_to(out))]=sha(marker)
        except (OSError,ValueError,KeyError):pass
    return known


SCRIPT=r'''
import csv,hashlib,json,sys,tarfile
from pathlib import Path
root=Path('/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/dogohla-series20260918')
COHORT_ROOT
known=KNOWN_CACHE
donors=list(csv.DictReader((root/'cohort.tsv').open(),delimiter='\t'))
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as tar:
 for name in ('cohort.tsv','PREPARED.json','code/MANIFEST.json'):
  if (root/name).exists():tar.add(root/name,arcname=name)
 for r in donors:
  d=r['donor'];paths=[root/'t1k4'/d,root/'arms/full/runs'/d/'native',root/'arms/full/runs'/d/'native-long']
  paths += [root/'arms'/a/'final/runs'/d/m for a in ('full','hprc','asian_matched') for m in ('DogoHLA','DogoHLA-no-graph')]
  for p in paths:
   marker=p/'COMPLETE'
   # COMPLETE binds the configuration, output hashes and completion time. Reuse
   # the already verified immutable snapshot when this content address agrees.
   if marker.exists() and str(p.relative_to(root)) in known:
    if hashlib.sha256(marker.read_bytes()).hexdigest()==known[str(p.relative_to(root))] and json.loads((p/'manifest.json').read_text())['status']=='complete':continue
   for n in ('manifest.json','COMPLETE','TERMINAL_FAILURE.json'):
    if (p/n).exists():tar.add(p/n,arcname=str((p/n).relative_to(root)))
   if marker.exists():
    m=json.loads(marker.read_text())
    for n in m['outputs']:
     if Path(n).is_absolute() or '..' in Path(n).parts:raise ValueError('Unsafe output path')
     tar.add(p/n,arcname=str((p/n).relative_to(root)))
'''


def main():
    extended='--gourraud' in sys.argv
    out=HERE/'snapshot'/('gourraud' if extended else '');out.mkdir(parents=True,exist_ok=True)
    known=verified_cache(out)
    script=SCRIPT.replace('COHORT_ROOT',"root=root/'gourraud'" if extended else '').replace('KNOWN_CACHE',repr(known))
    r=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','hohndor@ilogin.ibex.kaust.edu.sa','python3 -'],input=script.encode(),capture_output=True,check=True,timeout=90)
    with tarfile.open(fileobj=io.BytesIO(r.stdout),mode='r:gz') as tar:tar.extractall(out,filter='data')
    print('Fetched',len(r.stdout),'bytes; reused',len(known),'verified completed runs')


if __name__=='__main__':main()
