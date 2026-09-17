#!/bin/bash
#SBATCH --job-name=asian_named_HLA_score
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/refined/named-HLA-score-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run/refined
python3 score_locityper.py --run . --out results
