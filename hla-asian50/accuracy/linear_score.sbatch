#!/bin/bash
#SBATCH --job-name=asian_linear_score
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/linear-score-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run
python3 accuracy/score_linear.py --run . --out accuracy/results
