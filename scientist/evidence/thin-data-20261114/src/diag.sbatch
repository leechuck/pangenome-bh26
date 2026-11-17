#!/bin/bash
set -euo pipefail
TS=/encrypted/e3008/hohndor-dogohla-20261020/thin-data-20261114
export TMPDIR=$TS/tmp PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
umask 0002
module load singularity/3.9.7
export SINGULARITYENV_TMPDIR=$TS/tmp
PY=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin/python
exec $PY $TS/src/"$@"
