#!/bin/bash
#SBATCH --job-name=asian_named_HLA_archive
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/refined/locityper-archive-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run/refined
python3 - <<'PY'
import csv,hashlib,json,tarfile
from pathlib import Path
r=Path('.');ds=list(csv.DictReader(open('source/donors.tsv'),delimiter='\t'))
paths=[]
for d in ds:
 for arm in ['full','hprc']:
  for gene in ['A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1']:
   p=r/f"locityper/calls/{d['donor']}/{arm}/loci/HLA-{gene}/res.json.gz"
   assert p.is_file(),p
   paths.append(p)
assert len(paths)==640
for pattern in ['locityper/panels/*/*/labels.tsv.gz','locityper/panels/*/*/targets.bed','locityper/panels/*/*/target.log','locityper/panel_audit.json','locityper/preproc-*.log','locityper/genotype-*.log','locityper/preproc/*/*.json*','locityper/GRCh38.jf.COMPLETE']:
 paths.extend(sorted(r.glob(pattern)))
paths.extend(sorted(r.glob('*.py')));paths.extend(sorted(r.glob('*.sbatch')));paths.append(r/'DESIGN.md')
manifest=[]
with tarfile.open('results/locityper_native.tar.gz','w:gz') as t:
 for p in paths:
  t.add(p)
  manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
Path('results/locityper_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Archived',len(paths),'files, including 640 native locus results')
PY
