#!/bin/bash
#SBATCH --job-name=hlatyper_loci
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/leechuck/hla/hla-typer/logs/loci-%j.log
set -euo pipefail
cd /home/leechuck/hla/hla-typer
export PATH=/home/leechuck/hla/mm/envs/hla/bin:$PATH
python3 loci/extract_loci.py
echo COMPLETE
