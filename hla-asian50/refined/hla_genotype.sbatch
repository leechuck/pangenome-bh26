#!/bin/bash
#SBATCH --job-name=asian_whole_hla
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/refined/hla-%A_%a.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run/refined
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
read -r donor stratum population family paired truth fold url selhash < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" source/donors.tsv)
for arm in full hprc; do
 mkdir -p "hla_calls/$donor/$arm"
 PanGenie -f "hla_panels/fold$fold/$arm/index" -i <(gzip -dc "reads/$donor/r1.fq.gz" "reads/$donor/r2.fq.gz") -s "$donor" -j 2 -t 1 -e 20000000 -o "hla_calls/$donor/$arm/calls"
 test -s "hla_calls/$donor/$arm/calls_genotyping.vcf"
done
echo COMPLETE "$donor"
