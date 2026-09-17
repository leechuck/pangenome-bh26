#!/bin/bash
#SBATCH --job-name=asian50_fullgraph
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/graph-%j.log
set -euo pipefail
export PATH=/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin:/home/leechuck/hla/mm/envs/hla/bin:$PATH
cd /home/leechuck/hla/codex-asian50/run
mkdir -p graph
vg version > source/vg_version.txt
vg paths -M -x /home/asianhla/data/upload/HLA/mhc_graph/MHC.full.gbz > graph/paths.tsv
vg deconstruct -t 8 -C -p 'GRCh38#0#ctg1' /home/asianhla/data/upload/HLA/mhc_graph/MHC.full.gbz | bgzip -@2 > graph/full.top.vcf.gz
bcftools index -t graph/full.top.vcf.gz
bcftools stats graph/full.top.vcf.gz > graph/full.top.stats.txt
python3 - <<'PY'
from pathlib import Path
s=Path('/home/leechuck/hla/mhc_all/GRCh38_0.mhc.fa').read_text().splitlines()
assert sum(x.startswith('>') for x in s)==1
Path('graph/reference.fa').write_text('>ctg1\n'+'\n'.join(s[1:])+'\n')
PY
samtools faidx graph/reference.fa
bcftools norm -f graph/reference.fa -c e -Ou graph/full.top.vcf.gz > /dev/null
sha256sum graph/reference.fa graph/full.top.vcf.gz > source/fullgraph_input.sha256
echo COMPLETE
