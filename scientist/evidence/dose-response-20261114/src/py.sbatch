#!/bin/bash
# generic python job: sbatch ... py.sbatch <script.py> [args]
set -euo pipefail
S=/encrypted/e3008/hohndor-dogohla-20261020/dose-response-20261114
export TMPDIR=$S/tmp PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
umask 0002
module load singularity/3.9.7
export SINGULARITYENV_TMPDIR=$S/tmp
PY=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin/python
echo "host $(hostname) job $SLURM_JOB_ID args $*"
script=$1; shift
exec $PY $S/src/$script "$@"
