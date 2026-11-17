#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --job-name=s2-hdr
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%j.out
#SBATCH --error=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%j.err
# CRAM header @PG lines and idxstats for one donor (no reads printed)
set -euo pipefail
ST=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114
B=/ibex/scratch/projects/c2014/rob/dogohla-benchmark
module load singularity/3.9.7
donor=$1
SRC=$B/hprc-zenodo-20261004/reads/$donor/$donor.final.cram
case $donor in NA18944|NA18945|NA18970) SRC=$B/hprc-zenodo-20261004/japanese-read-fingerprints-extension-v1/$donor/outputs/alignment.cram;; esac
IDX=$B/hprc-zenodo-20261004/wgs-read-fingerprints-v1/outputs/$donor/input.cram.crai
[ -s $IDX ] || IDX=$SRC.crai
SAM="singularity exec --cleanenv --no-home --bind $B/hla:/home/leechuck/hla,$B/hprc-zenodo-20261004:$B/hprc-zenodo-20261004 $B/containers/ubuntu24.sif /home/leechuck/hla/mm/envs/asian50-spechla/bin/samtools"
mkdir -p $ST/work/hdr
$SAM view -H $SRC | grep -v '^@SQ' > $ST/work/hdr/$donor.header_noSQ.txt
$SAM view -H $SRC | grep -c '^@SQ' > $ST/work/hdr/$donor.nSQ.txt
$SAM idxstats -X $SRC $IDX > $ST/work/hdr/$donor.idxstats.txt
