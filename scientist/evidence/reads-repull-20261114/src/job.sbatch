#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --job-name=s2-job
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%A_%a.out
#SBATCH --error=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%A_%a.err
# usage: sbatch --job-name=s2-xxx --array=0-N%12 job.sbatch <script.py> <tasks.txt>  (tasks: args per line)
set -euo pipefail
ST=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114
export TMPDIR=$ST/tmp PYTHONNOUSERSITE=1
umask 0002
module load singularity/3.9.7
export SINGULARITYENV_TMPDIR=$ST/tmp
PY=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin/python
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$2")
echo "host $(hostname) task $SLURM_ARRAY_TASK_ID: $line"
exec $PY $ST/src/$1 $line
