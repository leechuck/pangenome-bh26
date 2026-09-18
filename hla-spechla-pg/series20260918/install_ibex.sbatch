#!/bin/bash
#SBATCH --job-name=dogohla-install
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --array=0-2
#SBATCH --chdir=/ibex/scratch/projects/c2014/rob/dogohla-benchmark
#SBATCH --output=/ibex/scratch/projects/c2014/rob/dogohla-benchmark/install-%A_%a.log
set -euo pipefail
module load singularity/3.9.7
names=(asian50-spechla typing pggb053)
name=${names[$SLURM_ARRAY_TASK_ID]}
mkdir -p "hla/bootstrap/$name" "hla/mm/envs"
curl --fail --location --retry 3 https://micro.mamba.pm/api/micromamba/linux-64/latest -o "hla/bootstrap/$name/micromamba.tar.bz2"
tar xjf "hla/bootstrap/$name/micromamba.tar.bz2" -C "hla/bootstrap/$name" bin/micromamba
singularity exec --bind "$PWD/hla:/home/leechuck/hla" --bind /etc/pki/tls/certs/ca-bundle.crt:/etc/ssl/certs/ca-certificates.crt containers/ubuntu24.sif "/home/leechuck/hla/bootstrap/$name/bin/micromamba" create -y \
 --root-prefix "/home/leechuck/hla/bootstrap/$name/cache" \
 --prefix "/home/leechuck/hla/mm/envs/$name" \
 --file "/home/leechuck/hla/bootstrap/explicit/$name.txt"
touch "hla/bootstrap/$name/COMPLETE"
