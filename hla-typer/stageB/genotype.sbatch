#!/bin/bash
#SBATCH --job-name=hlatyper_gt
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/home/leechuck/hla/hla-typer/logs/gt-%A_%a.log
# Stage B for one donor: whole-genome background calibration (chr17 slice of the donor's own CRAM),
# then leave-one-donor-out genotyping against every database in $ARMS. TAG suffixes the output
# directory (e.g. TAG=-w with EXTRA="--reg-weights weights/reg_weights.tsv").
set -euo pipefail
cd /home/leechuck/hla/hla-typer
export PATH=/home/leechuck/hla/mm/envs/asian50-genotyping/bin:$PATH
grch38=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
read -r donor cram < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "${DONORS:-source/donors_dev.tsv}")
test -e "reads/$donor/COMPLETE"
bg=background/$donor; mkdir -p "$bg" preproc calls
if [ ! -e "$bg/COMPLETE" ]; then
 samtools view -@8 -T "$grch38" -C -o "$bg/chr17.cram" "$cram" chr17:70000000-78000000
 samtools index "$bg/chr17.cram"
 touch "$bg/COMPLETE"
fi
if [ ! -e "preproc/$donor/COMPLETE" ]; then
 locityper preproc -a "$bg/chr17.cram" -r "$grch38" -j ref/GRCh38.jf -o "preproc/$donor" -@ 8 --seed 20260917
 touch "preproc/$donor/COMPLETE"
fi
# every panel sample name that belongs to this donor (includes JaSaPaGe re-assemblies)
leave=$(awk -F'\t' -v d="$donor" 'NR>1 && $3==d {split($2,a,"."); print a[1]}' source/haplotype_donors.tsv | sort -u | tr '\n' ' ')
# WEIGHTS=1: exon/intron/intergenic region weights. After --leave-out Locityper resolves the left-out
# names and their identical substitutes to the same haplotype, so the donor's own rows must be removed
# from the BED files or the coverage check fails ("not fully covered").
wopt=()
if [ "${WEIGHTS:-0}" = 1 ]; then
 wd=weights/perdonor/$donor; mkdir -p "$wd"; : > "$wd/reg_weights.tsv"
 pat=$(printf '%s\n' $leave | sed 's/[.]/\\./g' | paste -sd'|')
 while read -r locus bed; do
  awk -v re="^(${pat})[.][0-9]+$" '$1 !~ re' "$bed" > "$wd/$(basename "$bed")"
  echo -e "$locus\t$PWD/$wd/$(basename "$bed")" >> "$wd/reg_weights.tsv"
 done < weights/reg_weights.tsv
 wopt=(--reg-weights "$wd/reg_weights.tsv")
fi
for arm in ${ARMS:-full}; do
 out=calls/$arm${TAG:-}/$donor; mkdir -p "calls/$arm${TAG:-}"
 [ -e "$out/COMPLETE" ] && continue
 start=$(date +%s)
 locityper genotype -i "reads/$donor/r1.fq.gz" "reads/$donor/r2.fq.gz" -p "preproc/$donor" -d "db/$arm" -o "$out" \
   ${leave:+--leave-out $leave} "${wopt[@]}" -@ 8 --seed 20260917 ${EXTRA:-}
 echo -e "$donor\t$arm\t$(( $(date +%s) - start ))\t$leave" > "$out/run.tsv"
 touch "$out/COMPLETE"
done
echo COMPLETE "$donor"
