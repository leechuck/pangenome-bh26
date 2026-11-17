#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --job-name=s2b-post
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%j.out
#SBATCH --error=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%j.err
# usage: sbatch [--job-name=s2b-xxx] post3.sbatch "<script.py args>" ["<script.py args>" ...]   (run in order)
set -euo pipefail
ST=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114
export TMPDIR=$ST/tmp PYTHONNOUSERSITE=1
umask 0002
PY=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin/python
for c in "$@"; do echo "== $c"; $PY $ST/src/$c; done
