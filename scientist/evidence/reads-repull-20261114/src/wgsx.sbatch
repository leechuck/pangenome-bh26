#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --job-name=s2b-wgsx
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%A_%a.out
#SBATCH --error=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%A_%a.err
# Whole-genome T1K candidate extraction for donors without a baseline full-WGS T1K run:
# all primary pairs of the CRAM (samtools collate + fastq -F 0x900, as baselines/prepared/wgs.<donor>), then
# T1K fastq-extractor -s 0.97 against genomic.fa (as run-t1k --preset hla-wgs). Temporaries in this stage only.
# usage: sbatch --array=0-N%12 wgsx.sbatch donors.txt [OUTSUB]     (OUTSUB default wgsx)
set -euo pipefail
ST=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114
B=/ibex/scratch/projects/c2014/rob/dogohla-benchmark
Z=$B/hprc-zenodo-20261004
export TMPDIR=$ST/tmp PYTHONNOUSERSITE=1
umask 0002
module load singularity/3.9.7
donor=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$1")
SUB=${2:-wgsx}
O=$ST/work/$SUB/$donor; T=$ST/tmp/wgsx_$SUB/$donor
mkdir -p $O; rm -f $O/COMPLETE.json $O/FAILURE.json
trap 'echo "{\"status\":\"failed\",\"donor\":\"'$donor'\",\"job\":\"'$SLURM_JOB_ID'\"}" > $O/FAILURE.json' ERR
rm -rf $T; mkdir -p $T
SRC=$Z/reads/$donor/$donor.final.cram
REC=$Z/reads/$donor/$donor.final.cram.receipt.json
SIF="singularity exec --cleanenv --no-home --bind $B/hla:/home/leechuck/hla,$Z:/benchmark,$ST:$ST $B/containers/ubuntu24.sif"
SAMT=/home/leechuck/hla/mm/envs/asian50-spechla/bin/samtools
REF=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
want=$(grep -o '"sha256": "[0-9a-f]*"' $REC | tail -1 | grep -o '[0-9a-f]\{64\}')
got=$(sha256sum $SRC | cut -d' ' -f1)
[ "$want" = "$got" ] || { echo "CRAM sha mismatch"; false; }
CS=/benchmark/reads/$donor/$donor.final.cram
$SIF $SAMT collate -@ 3 -u -O --reference $REF $CS $T/collate-spill | \
  $SIF $SAMT fastq -@ 3 -F 0x900 -n -1 $T/r1.fq.gz -2 $T/r2.fq.gz -0 $T/other.fq.gz -s $T/singletons.fq.gz - 2> $O/fastq.log
cat $O/fastq.log
$SIF /benchmark/baselines/assets/t1k/fastq-extractor -t 8 -f /benchmark/baselines/assets/references/genomic.fa \
  -o $O/calls_candidate -s 0.97 -1 $T/r1.fq.gz -2 $T/r2.fq.gz > $O/extract.log 2>&1
n1=$(wc -l < $O/calls_candidate_1.fq); n2=$(wc -l < $O/calls_candidate_2.fq)
[ "$n1" = "$n2" ] && [ $((n1 % 4)) -eq 0 ]
rm -rf $T
echo "{\"status\":\"complete\",\"donor\":\"$donor\",\"job\":\"$SLURM_JOB_ID\",\"cram_sha256\":\"$got\",\"candidate_pairs\":$((n1/4))}" > $O/COMPLETE.json
