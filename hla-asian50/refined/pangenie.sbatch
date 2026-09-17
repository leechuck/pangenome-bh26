#!/bin/bash
#SBATCH --job-name=asian_missing_genotype
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/refined/pangenie-%A_%a.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run/refined
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
read -r donor stratum population family paired truth fold url selhash < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" source/donors.tsv)
test -e "reads/$donor/COMPLETE"
for arm in full hprc; do
 mkdir -p "pangenie/$donor/$arm"
 PanGenie -f "panels/fold$fold/$arm/index_v2" -i <(gzip -dc "reads/$donor/r1.fq.gz" "reads/$donor/r2.fq.gz") -s "$donor" -j 2 -t 1 -e 20000000 -o "pangenie/$donor/$arm/calls"
 test -s "pangenie/$donor/$arm/calls_genotyping.vcf"
done
echo COMPLETE "$donor"
