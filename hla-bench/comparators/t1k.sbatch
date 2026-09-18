#!/bin/bash
#SBATCH --job-name=bench_t1k
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/home/leechuck/hla/hla-typer/logs/t1k-%A_%a.log
# T1K 1.0.6, WGS preset, database built from IPD-IMGT/HLA 3.65.0 (t1kdb/SHA256SUMS), on the shared recruited reads.
set -euo pipefail
cd /home/leechuck/hla/hla-typer
export PATH=/home/leechuck/hla/mm/envs/typing/bin:$PATH
read -r donor cram < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "${DONORS:-source/donors_dev.tsv}")
test -e "reads/$donor/COMPLETE"
out=comparators/t1k/$donor; mkdir -p "$out"
[ -e "$out/COMPLETE" ] && exit 0
run-t1k -1 "reads/$donor/r1.fq.gz" -2 "reads/$donor/r2.fq.gz" --preset hla-wgs -f t1kdb/hla365_dna_seq.fa -t 4 -o "$out/$donor"
test -s "$out/${donor}_genotype.tsv"
touch "$out/COMPLETE"
