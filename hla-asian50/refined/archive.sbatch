#!/bin/bash
#SBATCH --job-name=asian_refined_archive
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/refined/archive-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run/refined
python3 - <<'PY'
import csv,hashlib,json,tarfile
from pathlib import Path
r=Path('.');ds=list(csv.DictReader(open('source/donors.tsv'),delimiter='\t'));manifest=[]
for d in ds:
 for category in ['pangenie','hla_calls']:
  for arm in ['full','hprc']:
   p=r/category/d['donor']/arm/'calls_genotyping.vcf';h=hashlib.sha256()
   with open(p,'rb') as f:
    for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
   manifest.append(dict(donor=d['donor'],category=category,arm=arm,path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest()))
Path('results/output_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
with tarfile.open('results/hla_native_calls.tar.gz','w:gz') as t:
 for x in manifest:
  if x['category']=='hla_calls':t.add(x['path'])
with tarfile.open('results/panel_metadata.tar.gz','w:gz') as t:
 for pattern in ['panels/*/*/audit.json','panels/*/*/samples.txt','hla_panels/*/*/allele_labels.tsv.gz','hla_panels/*/*/samples.txt','hla_panels/audit.json','hla_panels/reference_checks.json']:
  for p in r.glob(pattern):t.add(p)
Path('results/input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [r/'hla_sources.tar.gz',r/'prepare_panels.py',r/'prepare_hla.py',r/'score_variants.py',r/'score_hla.py',r/'DESIGN.md']},indent=2)+'\n')
print('Archived and hashed',len(manifest),'output VCFs')
PY
