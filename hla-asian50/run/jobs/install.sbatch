#!/bin/bash
#SBATCH --job-name=asian50_install
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/install-%A_%a.log
set -euo pipefail
mm=/home/leechuck/hla/bin/micromamba
case "$SLURM_ARRAY_TASK_ID" in
0) envname=asian50-genotyping; packages=(pangenie locityper jellyfish bwa samtools bcftools minimap2) ;;
1) envname=asian50-spechla; packages=(spechla) ;;
esac
"$mm" create -y -r /home/leechuck/hla/mm -p "/home/leechuck/hla/mm/envs/$envname" -c conda-forge -c bioconda "${packages[@]}"
"$mm" list -p "/home/leechuck/hla/mm/envs/$envname" --json > "/home/leechuck/hla/codex-asian50/run/source/$envname.packages.json"
echo COMPLETE
