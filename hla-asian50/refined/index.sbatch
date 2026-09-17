#!/bin/bash
#SBATCH --job-name=asian_missing_index_parallel
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/refined/index-%A_%a.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run/refined
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
fold=$((SLURM_ARRAY_TASK_ID / 2))
if (( SLURM_ARRAY_TASK_ID % 2 == 0 )); then arm=full; else arm=hprc; fi
p=panels/fold$fold/$arm
PanGenie-index -v "$p/bubbles.vcf" -r graph/reference.fa -t 4 -e 20000000 -o "$p/index_v2"
echo COMPLETE
