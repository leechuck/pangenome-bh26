#!/bin/bash
#SBATCH --job-name=c4_fresh_reads
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/home/leechuck/hla/codex-targeted/recruit-%A_%a.log
set -euo pipefail
export PATH=/home/leechuck/hla/mm/envs/hla/bin:$PATH
cd /home/leechuck/hla/codex-targeted
mkdir -p fresh_reads
read -r donor family cohort url selection < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" validation_donors.tsv)
test -n "$donor"
ref=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools view -@1 -T "$ref" -M -L mhc_regions.bed -C -o "fresh_reads/${donor}.partial.cram" "$url"
samtools quickcheck "fresh_reads/${donor}.partial.cram"
mv "fresh_reads/${donor}.partial.cram" "fresh_reads/${donor}.mhc.cram"
samtools index "fresh_reads/${donor}.mhc.cram"
echo COMPLETE "$donor"
