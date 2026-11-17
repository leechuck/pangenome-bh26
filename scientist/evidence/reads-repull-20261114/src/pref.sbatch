#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --job-name=s2b-pref
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%A_%a.out
#SBATCH --error=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%A_%a.err
# usage: sbatch --array=0-1 pref.sbatch donors.txt
set -euo pipefail
ST=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114
Z=/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hprc-zenodo-20261004
export TMPDIR=$ST/tmp PYTHONNOUSERSITE=1
umask 0002
PY=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin/python
donor=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$1")
W=$ST/work/pref/$donor; mkdir -p $W; rm -f $W/COMPLETE.json $W/FAILURE.json
trap 'echo "{\"status\":\"failed\",\"donor\":\"'$donor'\"}" > $W/FAILURE.json' ERR
[ -s $ST/work/pref/bait_both.txt ] || { echo "missing patterns"; exit 1; }
for m in 1 2; do
  ( set +o pipefail; zcat $Z/baselines/prepared/wgs.$donor/r$m.fq.gz | paste - - - - | cut -f1,2 | LC_ALL=C grep -F -f $ST/work/pref/bait_both.txt > $W/hits_r$m.tsv; st=("${PIPESTATUS[@]}"); echo "r$m pipestatus ${st[*]}"; [ "${st[0]}" -eq 0 ] && [ "${st[3]}" -le 1 ] ) &
done
wait %1; wait %2
$PY $ST/src/pref_check.py evaluate $donor
rm -f $W/hits_r1.tsv $W/hits_r2.tsv
echo '{"status":"complete"}' > $W/COMPLETE.json
