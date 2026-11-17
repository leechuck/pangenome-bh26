#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --job-name=s1-rescore
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/rescore-partial-20261114/logs/%x_%A_%a.out
#SBATCH --error=/encrypted/e3008/hohndor-dogohla-20261020/rescore-partial-20261114/logs/%x_%A_%a.err
# usage: sbatch [--job-name=s1-xxx] run.sbatch <script.py> [args...]
set -euo pipefail
S1=/encrypted/e3008/hohndor-dogohla-20261020/rescore-partial-20261114
export TMPDIR=$S1/tmp PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
umask 0002
PY=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin/python
script=$1; shift
echo "host $(hostname) $script $*"
exec $PY $S1/src/$script "$@"
