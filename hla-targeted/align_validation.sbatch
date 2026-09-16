#!/bin/bash
#SBATCH --job-name=c4_diagnostic_targets
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/home/leechuck/hla/codex-targeted/targets-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-targeted
/home/leechuck/hla/mm/envs/hla/bin/minimap2 -t2 -x asm5 -c --cs=long --secondary=no C4_reference.fa validation_C4_genes.fa > validation_C4_genes.paf
echo COMPLETE
