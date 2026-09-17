#!/bin/bash
#SBATCH --job-name=asian_locityper_setup
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/home/leechuck/hla/codex-asian50/run/refined/locityper-setup-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-asian50/run/refined
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
mkdir -p locityper
if ! test -e locityper/GRCh38.jf.COMPLETE; then
primary_chroms=()
for chr in {1..22} X Y M; do primary_chroms+=("chr$chr"); done
samtools faidx /home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa "${primary_chroms[@]}" | jellyfish count --canonical --lower-count 2 --out-counter-len 2 --mer-len 25 --threads 8 --size 3G --output locityper/GRCh38.jf /dev/fd/0
touch locityper/GRCh38.jf.COMPLETE
fi
python3 prepare_locityper.py
for fold in 0 1 2 3 4; do
 for arm in full hprc; do
  p=locityper/panels/fold$fold/$arm
  locityper target -d "$p/db" -r graph/reference.fa -j locityper/GRCh38.jf -L "$p/targets.bed" -@ 8 --force > "$p/target.log" 2>&1
 done
done
echo COMPLETE
