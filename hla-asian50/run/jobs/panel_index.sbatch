#!/bin/bash
#SBATCH --job-name=asian50_panel_index
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/panel-index-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
python3 prepare_panels.py
for fold in 0 1 2 3 4; do
 for arm in full hprc; do
  p=panels/fold$fold/$arm
  PanGenie-index -v "$p/bubbles.vcf" -r graph/reference.fa -t 4 -e 20000000 -o "$p/index" > "$p/index.log" 2>&1
 done
done
echo COMPLETE
