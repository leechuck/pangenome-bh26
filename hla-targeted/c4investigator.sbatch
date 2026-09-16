#!/bin/bash
#SBATCH --job-name=c4_dedicated_baseline
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=/home/leechuck/hla/codex-targeted/c4investigator-%A_%a.log
set -euo pipefail
export PATH=/home/leechuck/hla/mm/envs/c4-targeted/bin:$PATH
work=/home/leechuck/hla/codex-targeted
donor=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$work/${DONOR_LIST:-development_donors.txt}")
test -n "$donor"
ref=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
cram=/home/asianhla/data/upload/1000G_MHC/cram/${donor}.mhc.cram
if [ "${COHORT:-development}" = validation ]; then cram=$work/fresh_reads/${donor}.mhc.cram; fi
out=$work/baseline/${COHORT:-development}/$donor
mkdir -p "$out/fastq" "$out/results"
samtools collate -@1 -u -O --reference "$ref" "$cram" "$out/collate" \
 | samtools fastq -F 0xF04 -n -1 "$out/fastq/${donor}_C4_1.fq" -2 "$out/fastq/${donor}_C4_2.fq" -0 /dev/null -s /dev/null -
cd "$work/C4Investigator"
Rscript C4Investigator_run.R --fqDirectory "$out/fastq" --resultsDirectory "$out/results" --fastqPattern fq --threads 2
test -s "$out/results/C4Investigator_c4_summary.csv"
echo COMPLETE "$donor"
