#!/bin/bash
#SBATCH --job-name=hlatyper_scan
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/home/leechuck/hla/hla-typer/logs/scan-%j.log
set -euo pipefail
cd /home/leechuck/hla/hla-typer
export PATH=/home/leechuck/hla/mm/envs/hla/bin:$PATH
python3 loci/scan_graph.py --out results/graph_scan
echo COMPLETE
