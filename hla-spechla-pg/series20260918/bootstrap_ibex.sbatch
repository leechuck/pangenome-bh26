#!/bin/bash
#SBATCH --job-name=dogohla-runtime
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --chdir=/ibex/scratch/projects/c2014/rob/dogohla-benchmark
#SBATCH --output=/ibex/scratch/projects/c2014/rob/dogohla-benchmark/runtime-%j.log
set -euo pipefail
module load singularity/3.9.7
mkdir -p containers cache tmp
export SINGULARITY_CACHEDIR="$PWD/cache"
export SINGULARITY_TMPDIR="$PWD/tmp"
if [[ ! -s containers/ubuntu24.sif ]]; then
    singularity pull containers/ubuntu24.sif docker://ubuntu:24.04
fi
singularity exec containers/ubuntu24.sif cat /etc/os-release
sha256sum containers/ubuntu24.sif > containers/SHA256SUMS
