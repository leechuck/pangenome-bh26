#!/bin/bash
#SBATCH --job-name=c4_target_background
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --output=/home/leechuck/hla/codex-targeted/background-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-targeted
g++ -O3 -std=c++17 count_groups.cpp -o count_groups
mkdir -p background
while read -r entry; do
 /home/leechuck/hla/codex-structural/count_markers fasta markers.txt "/home/leechuck/hla/mhc_all/${entry}.mhc.fa" "background/${entry}.bin"
done < training_entries.txt
echo COMPLETE
