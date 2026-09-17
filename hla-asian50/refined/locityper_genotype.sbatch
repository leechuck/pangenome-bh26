#!/bin/bash
#SBATCH --job-name=asian_locityper_HLA
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/refined/locityper-%A_%a.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run/refined
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
read -r donor stratum population family paired truth fold url selhash < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" source/donors.tsv)
mkdir -p locityper/preproc "locityper/calls/$donor"
locityper preproc -i "reads/$donor/r1.fq.gz" "reads/$donor/r2.fq.gz" -r graph/reference.fa -j locityper/GRCh38.jf -b ctg1:100000-3100000 -o "locityper/preproc/$donor" -@ 2 --seed 20260917 > "locityper/preproc-$donor.log" 2>&1
for arm in full hprc; do
 locityper genotype -i "reads/$donor/r1.fq.gz" "reads/$donor/r2.fq.gz" -p "locityper/preproc/$donor" -d "locityper/panels/fold$fold/$arm/db" -o "locityper/calls/$donor/$arm" -@ 2 --seed 20260917 > "locityper/genotype-$donor-$arm.log" 2>&1
done
echo COMPLETE "$donor"
