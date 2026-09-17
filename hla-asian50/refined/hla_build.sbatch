#!/bin/bash
#SBATCH --job-name=asian_hla_index
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/refined/hla-build-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run/refined
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
mkdir -p hla_source
tar xzf hla_sources.tar.gz -C hla_source
python3 prepare_hla.py --run .. --source hla_source --out hla_panels
for fold in 0 1 2 3 4; do
 for arm in full hprc; do
  p=hla_panels/fold$fold/$arm
  PanGenie-index -v "$p/bubbles.vcf" -r hla_panels/reference.fa -t 4 -e 20000000 -o "$p/index" > "$p/index.log" 2>&1
 done
done
echo COMPLETE
