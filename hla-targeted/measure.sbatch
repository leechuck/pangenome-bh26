#!/bin/bash
#SBATCH --job-name=c4_targeted_markers
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/home/leechuck/hla/codex-targeted/measure-%A_%a.log
set -euo pipefail
export PATH=/home/leechuck/hla/mm/envs/hla/bin:$PATH
cd /home/leechuck/hla/codex-targeted
cohort=${COHORT:-development}
donor=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${DONOR_LIST:-all_development_donors.txt}")
cram=/home/asianhla/data/upload/1000G_MHC/cram/${donor}.mhc.cram
if [ "$cohort" = validation ]; then cram=fresh_reads/${donor}.mhc.cram; fi
mkdir -p "counts/$cohort"
ref=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools collate -@1 -u -O --reference "$ref" "$cram" "counts/$cohort/${donor}.tmp" \
 | samtools fastq -n -F 0xF04 - \
 | ./count_groups marker_groups.tsv "counts/$cohort/${donor}.bin" "counts/$cohort/${donor}.groups.tsv"
samtools view -T "$ref" -u -F 0xF04 "$cram" chr6:31960000-32110000 \
 | samtools depth -a -q20 -Q0 -b rccx_depth_targets.bed - \
 | gzip > "counts/$cohort/${donor}.depth.tsv.gz"
echo COMPLETE "$donor"
