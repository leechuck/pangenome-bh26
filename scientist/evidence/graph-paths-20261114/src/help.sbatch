#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --job-name=s5-help
#SBATCH --partition=batch
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/encrypted/e3008/hohndor-dogohla-20261020/graph-paths-20261114/logs/%x_%j.out
set -uo pipefail
S=/encrypted/e3008/hohndor-dogohla-20261020/graph-paths-20261114
export TMPDIR=$S/tmp
module load singularity/3.9.7
B=/ibex/scratch/projects/c2014/rob/dogohla-benchmark
VG=/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin/vg
R="singularity exec --cleanenv --no-home --bind $B/hla:/home/leechuck/hla $B/containers/ubuntu24.sif"
$R $VG version
for c in haplotypes gbwt paths index giraffe pack call minimizer convert stats; do echo "=== $c"; $R $VG $c --help 2>&1; done
echo "=== kmc"; $R /home/leechuck/hla/t1k-pangenome/tools/kmc-3.2.4/kmc 2>&1 | head -40
