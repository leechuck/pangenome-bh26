#!/bin/bash
# usage: sbatch --array=0-N%12 --job-name=s5-xxx gp.sbatch SCRIPT TASKFILE
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/graph-paths-20261114/logs/%x_%A_%a.out
set -euo pipefail
S=/encrypted/e3008/hohndor-dogohla-20261020/graph-paths-20261114
export TMPDIR=$S/tmp PYTHONNOUSERSITE=1 SINGULARITYENV_TMPDIR=$S/tmp
umask 0002
module load singularity/3.9.7
PY=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin/python
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$2")
echo "host $(hostname) task $SLURM_ARRAY_TASK_ID: $line"
exec $PY $S/src/$1 $line
