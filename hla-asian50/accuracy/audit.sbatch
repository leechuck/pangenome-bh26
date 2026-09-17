#!/bin/bash
#SBATCH --job-name=asian_accuracy_audit
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/accuracy-audit-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run
mkdir -p accuracy/source
sha256sum graph/full.top.vcf.gz graph/reference.fa source/donors.tsv source/excluded_paths.tsv accuracy/score_panels.py accuracy/score_linear.py > accuracy/source/evaluation_inputs.sha256
cp /home/leechuck/hla/mhc_all/GRCh38_0.mhc.tsv accuracy/source/reference_extraction.tsv
python3 - <<'PY'
import subprocess,hashlib,json
from pathlib import Path
ref=''.join(Path('graph/reference.fa').read_text().splitlines()[1:]).upper()
linear=subprocess.check_output(['/home/leechuck/hla/mm/envs/asian50-genotyping/bin/samtools','faidx','/home/leechuck/hla/kg/chr6.fa','chr6:28410701-33580488'],text=True)
linear=''.join(linear.splitlines()[1:]).upper()
assert ref==linear
Path('accuracy/source/reference_validation.json').write_text(json.dumps(dict(equal=True,bases=len(ref),sha256=hashlib.sha256(ref.encode()).hexdigest(),chr6_interval_1based='28410701-33580488',offset=28410700),indent=2)+'\n')
PY
