#!/bin/bash
#SBATCH --job-name=asian50_spechla
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/spechla-%A_%a.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run
export CONDA_PREFIX=/home/leechuck/hla/mm/envs/asian50-spechla
export PATH=$CONDA_PREFIX/bin:$PATH
donor=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" source/donors.tsv | cut -f1)
test -s "reads/$donor/r1.fq.gz"
test -e "reads/$donor/COMPLETE"
mkdir -p spechla
spechla -n "$donor" -1 "$PWD/reads/$donor/r1.fq.gz" -2 "$PWD/reads/$donor/r2.fq.gz" -o "$PWD/spechla" -j 4 -u 0 -p Unknown
find "spechla/$donor" -maxdepth 1 -name '*result*' -size +0c > "spechla/$donor/results.list"
test -s "spechla/$donor/results.list"
echo COMPLETE "$donor"
