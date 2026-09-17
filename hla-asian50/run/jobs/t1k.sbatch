#!/bin/bash
#SBATCH --job-name=asian50_t1k
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/t1k-%A_%a.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run
export PATH=/home/leechuck/hla/mm/envs/typing/bin:$PATH
read -r donor stratum population family paired truth fold url selhash < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" source/donors.tsv)
test -n "$donor"
mkdir -p "reads/$donor" "t1k/$donor"
ref=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
cram=/home/asianhla/data/upload/1000G_MHC/cram/$donor.mhc.cram
if [ ! -s "$cram" ]; then cram=/home/leechuck/hla/codex-targeted/fresh_reads/$donor.mhc.cram; fi
if [ ! -s "$cram" ]; then
 cram="reads/$donor/mhc.cram"
 if [ ! -s "$cram" ]; then
  samtools view -@2 -T "$ref" -M -L source/mhc_regions.bed -C -o "reads/$donor/partial.cram" "$url"
  samtools quickcheck "reads/$donor/partial.cram"
  mv "reads/$donor/partial.cram" "$cram"
  samtools index "$cram"
 fi
fi
printf '%s\n' "$cram" > "reads/$donor/input_path.txt"
samtools quickcheck "$cram"
if [ ! -e "reads/$donor/COMPLETE" ]; then
 samtools collate -@2 -u -O --reference "$ref" "$cram" "reads/$donor/collate" |
 samtools fastq -@2 -F 0x900 -n -1 "reads/$donor/r1.fq.gz" -2 "reads/$donor/r2.fq.gz" -0 /dev/null -s "reads/$donor/singletons.fq.gz" -
 touch "reads/$donor/COMPLETE"
fi
run-t1k -1 "reads/$donor/r1.fq.gz" -2 "reads/$donor/r2.fq.gz" --preset hla-wgs -f source/t1k_dna.fa -t 4 -o "t1k/$donor/$donor"
test -s "t1k/$donor/${donor}_genotype.tsv"
echo COMPLETE "$donor"
