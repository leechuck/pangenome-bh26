#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --job-name=s2-bait
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%j.out
#SBATCH --error=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/logs/%x_%j.err
set -euo pipefail
ST=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114
B=/ibex/scratch/projects/c2014/rob/dogohla-benchmark
export TMPDIR=$ST/tmp PYTHONNOUSERSITE=1
umask 0002
module load singularity/3.9.7
PY=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin/python
W=$ST/work/bait
$PY $ST/src/bait.py prep
JF="singularity exec --cleanenv --no-home --bind $B/hla:/home/leechuck/hla,$ST:$ST $B/containers/ubuntu24.sif /home/leechuck/hla/mm/envs/asian50-spechla/bin/jellyfish"
n=$(grep -c '^>' $W/ipd_kmers.fa)
$JF count -m 31 -C -t 8 -s $((n*2)) --if $W/ipd_kmers.fa -o $W/found.jf $W/ref_outside_mhc.fa
$JF dump -c $W/found.jf > $W/found_outside_mhc.txt
$PY $ST/src/bait.py final
rm -f $W/ref_outside_mhc.fa $W/found.jf
echo done
