#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --job-name=s2-scan
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=10:00:00
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%A_%a.out
#SBATCH --error=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%A_%a.err
# usage: sbatch --array=0-N scan.sbatch donors.txt     (one donor per line)
set -euo pipefail
ST=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114
B=/ibex/scratch/projects/c2014/rob/dogohla-benchmark
export TMPDIR=$ST/tmp PYTHONNOUSERSITE=1
umask 0002
module load singularity/3.9.7
PY=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin/python
donor=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$1")
W=$ST/work/${SCAN_SUB:-scan}/$donor
mkdir -p $W
rm -f $W/COMPLETE.json $W/FAILURE.json
trap 'echo "{\"status\":\"failed\",\"donor\":\"'$donor'\",\"job\":\"'$SLURM_JOB_ID'\"}" > $W/FAILURE.json' ERR
echo "host $(hostname) donor $donor"
[ -s $W/patterns.raw.txt ] || $PY $ST/src/scan.py kmers $donor
SRC=$B/hprc-zenodo-20261004/reads/$donor/$donor.final.cram
REF=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
case $donor in NA18944|NA18945|NA18970)
  SRC=$B/hprc-zenodo-20261004/japanese-read-fingerprints-extension-v1/$donor/outputs/alignment.cram;;
esac
IDX=$B/hprc-zenodo-20261004/wgs-read-fingerprints-v1/outputs/$donor/input.cram.crai
[ -s $IDX ] || IDX=$(ls $B/hprc-zenodo-20261004/MHC-recruitment-*/$donor/input.cram.crai $B/hprc-zenodo-20261004/MHC-recruitment-control-$donor-*/input.cram.crai 2>/dev/null | head -1)
case $donor in NA18944|NA18945|NA18970) IDX=$SRC.crai;; esac
test -s $SRC
test -s $IDX
SAM="singularity exec --cleanenv --no-home --bind $B/hla:/home/leechuck/hla,$B/hprc-zenodo-20261004:$B/hprc-zenodo-20261004 $B/containers/ubuntu24.sif /home/leechuck/hla/mm/envs/asian50-spechla/bin/samtools"
# repeat filter (truth-free): drop k-mers seen in reads aligned to chr1:1-60,000,000 (far from the MHC)
if [ ! -s $W/patterns.raw.txt ]; then
  cp $W/patterns.txt $W/patterns.raw.txt
  set +o pipefail
  $SAM view -@ 6 -T $REF -X $SRC $IDX chr1:1-60000000 | cut -f10 | LC_ALL=C grep -o -F -f $W/patterns.raw.txt | sort -u > $W/repeat_kmers.txt
  st=("${PIPESTATUS[@]}")
  set -o pipefail
  echo "filter pipestatus ${st[*]}"
  [ "${st[0]}" -eq 0 ]
  $PY $ST/src/scan.py filter $donor
fi
if [ ! -s $W/hits.sam.done ]; then
  set +o pipefail
  $SAM view -@ 6 -T $REF $SRC \
    | LC_ALL=C grep -F -f $W/patterns.txt > $W/hits.sam
  st=("${PIPESTATUS[@]}")
  set -o pipefail
  echo "pipestatus ${st[*]}"
  [ "${st[0]}" -eq 0 ] && [ "${st[1]}" -le 1 ]
  echo "$SRC $(stat -c %s $SRC)" > $W/hits.sam.done
fi
$PY $ST/src/scan.py analyze $donor
echo "{\"status\":\"complete\",\"donor\":\"$donor\",\"job\":\"$SLURM_JOB_ID\",\"source\":\"$SRC\"}" > $W/COMPLETE.json
