#!/bin/bash
#SBATCH --job-name=c4_target_env
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=/home/leechuck/hla/codex-targeted/install-%j.log
set -euo pipefail
/home/leechuck/hla/bin/micromamba create -y -r /home/leechuck/hla/mm -p /home/leechuck/hla/mm/envs/c4-targeted -c conda-forge -c bioconda r-base=4.3 r-data.table r-stringr r-plotly r-argparser bowtie2 samtools
echo COMPLETE
