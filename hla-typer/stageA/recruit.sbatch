#!/bin/bash
#SBATCH --job-name=hlatyper_recruit
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/home/leechuck/hla/hla-typer/logs/recruit-%A_%a.log
# Recruit HLA reads from a full 1000G CRAM: chr6 MHC (+/-100 kb, identical to graph ctg1),
# all chr6 alts, all 525 HLA-* contigs, plus unplaced unmapped pairs. Pairs are kept when both
# mates fall in this set. Mate fetching (samtools view -P) was dropped: on NA18620 only 1.5% of
# region reads have a mate elsewhere, scattered genome-wide, and fetching them read >5.8 GB of the
# remote CRAM (>30 min per donor).
# DONORS (tsv: donor, cram path or URL) selects the cohort; array index = data line number - 1.
set -euo pipefail
cd /home/leechuck/hla/hla-typer
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
ref=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
read -r donor cram < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "${DONORS:-source/donors_dev.tsv}")
test -n "$donor"
out=reads/$donor; mkdir -p "$out"
if [ -e "$out/COMPLETE" ]; then echo ALREADY "$donor"; exit 0; fi
bed=reads/regions.bed
if [ ! -s "$bed" ]; then
 { printf 'chr6\t28410700\t33580488\n'; awk -v OFS='\t' '$1~/^chr6_.*_alt$/ || $1~/^HLA-/{print $1,0,$2}' "$ref.fai"; } > "$bed.tmp"
 mv "$bed.tmp" "$bed"
fi
samtools view -@4 -T "$ref" -M -L "$bed" -C -o "$out/regions.cram" "$cram"
samtools view -@4 -T "$ref" -C -o "$out/unmapped.cram" "$cram" '*'
samtools cat -o "$out/recruited.cram" "$out/regions.cram" "$out/unmapped.cram"
rm "$out/regions.cram" "$out/unmapped.cram"
samtools collate -@4 -u -O --reference "$ref" "$out/recruited.cram" "$out/collate" |
 samtools fastq -@4 -F 0x900 -n -1 "$out/r1.fq.gz" -2 "$out/r2.fq.gz" -0 /dev/null -s "$out/singletons.fq.gz" -
samtools flagstat -@4 "$out/recruited.cram" > "$out/flagstat.txt"
samtools view -@4 -F 0x900 "$out/recruited.cram" | awk -F'\t' '{c=($3=="*")?"unmapped":($3~/^HLA-/)?"HLA_contig":($3~/_alt$/)?"chr6_alt":($3=="chr6")?"chr6":"other_mate";n[c]++} END{for(k in n)print k"\t"n[k]}' | sort > "$out/contig_classes.tsv"
touch "$out/COMPLETE"
echo COMPLETE "$donor"
