#!/bin/bash
#SBATCH --job-name=hlatyper_c1
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/home/leechuck/hla/hla-typer/logs/c1-%A_%a.log
# Stage C1 for one donor: candidate haplotypes (Stage B top options x IPD CDS grafts [+ T1K alleles]),
# per-donor Locityper database, exon-weighted genotyping with the donor's Stage B preprocessing.
# VARIANT=t1k (ensemble) or ipd (no T1K candidates); ARM = Stage B panel arm.
set -euo pipefail
cd /home/leechuck/hla/hla-typer
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
py=/home/leechuck/hla/mm/envs/hlatyper/bin/python3
read -r donor cram < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "${DONORS:-source/donors_dev.tsv}")
arm=${ARM:-full}; variant=${VARIANT:-t1k}; partial=${PARTIAL:-0}
if [ "$partial" = 1 ]; then variant="$variant-p"; fi
if [ "${PRIORS:-0}" = 1 ]; then variant="$variant-pr"; fi
out=c1/$arm-$variant/$donor
if [ -e "$out/COMPLETE" ]; then echo ALREADY "$donor"; exit 0; fi
test -e "calls/$arm/$donor/COMPLETE"
extra=()
if [ "$partial" = 1 ]; then extra+=(--partial); fi
case "$variant" in
 t1k*) test -e "comparators/t1k/$donor/COMPLETE"; extra+=(--t1k "comparators/t1k/$donor/${donor}_genotype.tsv") ;;
 *)    extra+=(--no-t1k) ;;
esac
# PRIORS=1 adds panel-frequency genotype priors (variant suffix -pr); REUSE=<c1 dir> reuses that run's
# candidates and database instead of rebuilding them (same candidate set, different genotyping options).
mkdir -p "$out"
if [ -n "${REUSE:-}" ] && [ -e "$REUSE/$donor/COMPLETE" ]; then
 ln -sfn "$PWD/$REUSE/$donor/candidates" "$out/candidates"; ln -sfn "$PWD/$REUSE/$donor/db" "$out/db"
else
 (cd stageC1 && $py build_candidates.py --donor "$donor" --arm "$arm" --out "../$out/candidates" "${extra[@]}")
 locityper target -d "$out/db" -r ref/reference.fa -j ref/GRCh38.jf -L "$out/candidates/targets.bed" -@ 8 --force
 locityper augment -d "$out/db" -@ 8
fi
popt=()
if [ "${PRIORS:-0}" = 1 ]; then
 (cd stageC1 && $py make_priors.py --donor "$donor" --candidates "../$out/candidates" --db "../$out/db" --out "../$out/priors.tsv")
 popt=(--priors "$out/priors.tsv")
fi
locityper genotype -i "reads/$donor/r1.fq.gz" "reads/$donor/r2.fq.gz" -p "preproc/$donor" -d "$out/db" -o "$out/calls" \
  --reg-weights "$out/candidates/reg_weights.tsv" "${popt[@]}" -@ 8 --seed 20260917
touch "$out/COMPLETE"
echo COMPLETE "$donor"
