#!/bin/bash
#SBATCH --job-name=hlatyper_db
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/leechuck/hla/hla-typer/logs/db-%x-%j.log
# Build one Locityper database for panel arm $ARM; augment loci listed in $LOCI (default: all).
set -euo pipefail
cd /home/leechuck/hla/hla-typer
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
: "${ARM:?set ARM}"
db=db/$ARM; mkdir -p db
if [ ! -e "$db/TARGET_COMPLETE" ]; then
 locityper target -d "$db" -r ref/reference.fa -j ref/GRCh38.jf -L "panels/$ARM/targets.bed" -@ 16 --force
 touch "$db/TARGET_COMPLETE"
fi
for locus in ${LOCI:-HLA-A HLA-C HLA-B HLA-DRB1 HLA-DQ HLA-DP HLA-DRB345}; do
 [ -e "$db/AUGMENT_$locus" ] && continue
 start=$(date +%s)
 locityper augment -d "$db" --subset-loci "$locus" -@ 16
 echo -e "$locus\t$(( $(date +%s) - start ))" >> "$db/augment_seconds.tsv"
 touch "$db/AUGMENT_$locus"
done
echo COMPLETE "$ARM"
