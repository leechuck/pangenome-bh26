#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --job-name=s4-run
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114/logs/%x_%A_%a.out
#SBATCH --error=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114/logs/%x_%A_%a.err
# usage: sbatch --job-name=s4-xxx --array=0-N%12 s4.sbatch <script.py> <tasks.txt> [extra args]
#        (tasks.txt: one argument line per task)
set -euo pipefail
N=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114
export TMPDIR=$N/tmp PYTHONNOUSERSITE=1
umask 0002
module load singularity/3.9.7
export SINGULARITYENV_TMPDIR=$N/tmp
PY=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin/python
script=$1; tasks=$2; shift 2
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$tasks")
echo "host $(hostname) task $SLURM_ARRAY_TASK_ID: $line $*"
exec $PY $N/src/$script $line "$@"
