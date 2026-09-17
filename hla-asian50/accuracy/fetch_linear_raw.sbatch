#!/bin/bash
#SBATCH --job-name=asian_linear
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/linear-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
mkdir -p accuracy/source
cut -f1 source/donors.tsv | tail -n +2 > accuracy/source/samples.txt
url=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.recalibrated_variants.vcf.gz
printf '%s\n' "$url" > accuracy/source/linear_raw_url.txt
bcftools view -r chr6:28510121-33480577 -S accuracy/source/samples.txt -Oz -o accuracy/source/linear_raw.vcf.gz "$url"
bcftools index -t accuracy/source/linear_raw.vcf.gz
sha256sum accuracy/source/linear_raw.vcf.gz > accuracy/source/linear_raw.sha256
bcftools --version > accuracy/source/bcftools_version.txt
