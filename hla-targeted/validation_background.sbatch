#!/bin/bash
#SBATCH --job-name=c4_validation_specificity
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=/home/leechuck/hla/codex-targeted/validation-background-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-targeted
mkdir -p validation_background
while read -r entry; do
 /home/leechuck/hla/codex-structural/count_markers fasta markers.txt "/home/leechuck/hla/mhc_all/${entry}.mhc.fa" "validation_background/${entry}.bin"
done < validation_entries.txt
echo COMPLETE
