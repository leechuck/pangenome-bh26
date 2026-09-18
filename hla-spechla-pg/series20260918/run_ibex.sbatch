#!/bin/bash
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
#SBATCH --chdir=/ibex/scratch/projects/c2014/rob/dogohla-benchmark
set -euo pipefail
module load singularity/3.9.7
base=/home/leechuck/hla/dogohla-series20260918
workdir="$base/workdirs/${SLURM_JOB_ID:-manual}"
mkdir -p "$PWD/hla/dogohla-series20260918/workdirs/${SLURM_JOB_ID:-manual}"
export SINGULARITYENV_DOGOHLA_PANEL="$base/panels/${ARM:-full}"
export SINGULARITYENV_PATH=/home/leechuck/hla/mm/envs/asian50-spechla/bin:/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_SLURM_ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
export SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-4}"
export SINGULARITYENV_DONORS="${DONORS:-}"
export SINGULARITYENV_LC_ALL=C.UTF-8
export SINGULARITYENV_LANG=C.UTF-8
export SINGULARITYENV_BASH_ENV=/dev/null
export SINGULARITYENV_PYTHONNOUSERSITE=1
exec singularity exec --cleanenv --no-home --pwd "$workdir" --bind "$PWD/hla:/home/leechuck/hla" --bind /etc/pki/tls/certs/ca-bundle.crt:/etc/ssl/certs/ca-certificates.crt containers/ubuntu24.sif "$@"
