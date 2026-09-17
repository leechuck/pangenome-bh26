#!/bin/bash
#SBATCH --job-name=asian50_verify
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/verify-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run
python3 verify_outputs.py
