#!/bin/bash
#SBATCH --job-name=dogohla-assets
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --chdir=/ibex/scratch/projects/c2014/rob/dogohla-benchmark
#SBATCH --output=/ibex/scratch/projects/c2014/rob/dogohla-benchmark/assets-%j.log
set -euo pipefail
mkdir -p hla/cactus hla/ref1kg/local
curl --fail --location --retry 3 https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v3.3.0/cactus-bin-v3.3.0.tar.gz -o hla/cactus/release.tar.gz
tar xzf hla/cactus/release.tar.gz -C hla/cactus
base=https://s3.amazonaws.com/1000genomes/technical/reference/GRCh38_reference_genome
ref=GRCh38_full_analysis_set_plus_decoy_hla.fa
for file in "$ref" "$ref.fai"; do
 curl --fail --location --retry 3 "$base/$file" -o "hla/ref1kg/local/$file.tmp"
 mv "hla/ref1kg/local/$file.tmp" "hla/ref1kg/local/$file"
done
sha256sum hla/cactus/release.tar.gz hla/ref1kg/local/* > public-assets.sha256
touch PUBLIC_ASSETS_COMPLETE
