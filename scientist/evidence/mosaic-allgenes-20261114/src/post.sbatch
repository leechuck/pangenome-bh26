#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --job-name=s4-post
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114/logs/%x_%j.out
#SBATCH --error=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114/logs/%x_%j.err
# usage: sbatch --job-name=s4-xxx post.sbatch "<script.py args>" ["<script.py args>" ...]   (run in order)
set -euo pipefail
N=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114
export TMPDIR=$N/tmp PYTHONNOUSERSITE=1
umask 0002
PY=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin/python
for c in "$@"; do echo "== $c"; $PY $N/src/$c; done
