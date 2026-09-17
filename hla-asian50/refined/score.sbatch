#!/bin/bash
#SBATCH --job-name=asian_refined_score
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/refined/score-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run/refined
python3 score_variants.py --run . --original .. --out results
