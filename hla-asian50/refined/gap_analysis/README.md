# Why the graph panel does not yet beat SpecHLA (2026-09-17)

Post-hoc diagnosis of the frozen comparison in [../REPORT.md](../REPORT.md). No predictions
were rerun. Inputs: `../results/named_HLA_scores.tsv`, `../results/locityper_native.tar.gz`,
and from DDBJ `run/refined/hla_source/` and `run/refined/locityper/panels/*/full/labels.tsv.gz`.
`diag.py` and `oracle.py` expect those copied next to them as `hla_source/`, `panels/`, `loc/`
and need `edlib`.

## 1. The benchmark itself is compromised for 21 of 40 donors

Read recruitment used three CRAM sources (`run/jobs/t1k.sbatch`). The 19 donors found in
`/home/asianhla/data/upload/1000G_MHC/cram` include reads placed on the 525 `HLA-*` contigs.
The 21 freshly downloaded donors used `run/source/mhc_regions.bed` (chr6 MHC plus chr6 alts
only), so reads that bwa placed on `HLA-*` contigs are missing. Locityper read counts at
HLA-A/B/C/DRB1 are about 35% lower in these donors; DPA1/DPB1 are unaffected.

Assembly endpoint, exact two-field genotype accuracy:

| Method | 19 intact donors (151) | 21 depleted donors (168) |
|---|---:|---:|
| Locityper full | 137 (90.7%) | 143 (85.1%) |
| Locityper HPRC-only | 132 (87.4%) | 141 (83.9%) |
| SpecHLA | 141 (93.4%) | 95 (56.5%) |
| T1K | 147 (97.4%) | 107 (63.7%) |

The reported assembly-endpoint lead over SpecHLA and T1K is an artefact of this defect.
All eight experimental donors are intact, so the experimental endpoint is unaffected.
On intact donors the order is T1K > SpecHLA > full panel on both endpoints.

## 2. The experimental endpoint cannot show superiority

39 genotypes. Two DRB1 truths are historical artefacts (NA18608 `14:01;14:44`,
NA18952 `14:01`; assemblies and every method say 14:05/14:07 and 14:54). The ceiling is
37/39 and SpecHLA has 36. One genotype separates SpecHLA from a perfect method.

## 3. Locityper errors (full panel): the panel is mostly not the problem

`locityper_error_diagnosis.txt` lists, for every wrong genotype, the edit distance from each
held-out truth haplotype to the nearest panel haplotypes and to the picked haplotypes.

- In most wrong picks a panel haplotype at distance 0 to 5 with the correct label exists,
  while the pick is 100 to 1,000 edits away, with low Locityper quality (3 to 20).
- One unlabelled DRB1 haplotype (`HIFI032302D#1`, protein-exact DRB1*08:03) is picked in six
  donors at distance 600 to 1,100 from truth and causes most DRB1 no-calls. Untested
  hypothesis: it absorbs paralog reads (DRB3 and DRB4 are absent from the GRCh38 primary
  assembly, so their k-mers are not masked as repeated).
- A*02:06 replaces 02:01/02:03/02:07/02:11 in six donors, all of them read-depleted.
- Setup weaknesses visible in the logs: haplotypes shorter than 10 kb with non-homologous
  ends (`target.log`: 140/158 HLA-A and 411/414 DRB1 haplotypes do not match the reference
  boundary; lengths range 6.7 to 8.1 kb at HLA-B, 9.6 to 13.8 kb at DPA1); alignment recovery
  disabled (no pairwise haplotype alignments); error and depth calibration on the MHC itself
  (mismatch rate 0.29%, depth 23.6 for 30x samples).
- 9 no-calls come from the exact-CDS-or-nothing label rule.

## 4. A closed panel has a ceiling near SpecHLA

`oracle_ceiling.txt`: if a perfect genotyper always chose the panel haplotype nearest to each
held-out haplotype, the two-field genotype accuracy would be 299/319 (93.7%); Locityper reaches
280/319. Only 223/319 genotypes have both haplotypes exactly in the panel; 291/319 have both
within 5 edits. Alleles absent from the panel (B*15:27, B*15:07, A*02:09, A*74:02, DRB1*14:07,
DRB1*15:06, DQA1*05:09, DQB1*02:07, DPB1*16:01/414:01) are often one CDS substitution from
the nearest panel haplotype. SpecHLA and T1K see all IPD-IMGT/HLA alleles.

## 5. What is required

1. Repair the benchmark: recruit chr6 MHC, chr6 alts, all `HLA-*` contigs, unmapped reads and
   mates for all 40 donors; rerun all four methods; same IPD release for all; score
   DRB1*14:01/14:54 era labels at the resolution of the historical assay. Enlarge the endpoint
   (leave-one-donor-out over all panel donors with 1000G reads, plus the Gourraud 2014 and
   Lai 2024 labels); 39 genotypes cannot separate 92% from 95%.
2. Repair the locus genotyping: cut loci from the graph between conserved anchor nodes with
   at least 10 kb homologous flanks; one joint DRB target (DRB1 with DRB3/4/5 and pseudogenes)
   or paralog decoys; pairwise alignments so alignment recovery works; background calibration
   on a non-polymorphic region from the whole-genome CRAM.
3. Open the panel: after haplotype-pair selection, realign reads to the two selected
   haplotypes, call residual variants per haplotype, apply them to the CDS and look the CDS up
   in IPD-IMGT/HLA. Alternatively add IPD genomic alleles as extra panel haplotypes.
4. Replace exact-CDS-or-no-call labels with nearest-IPD labels (exact CDS, else ARS exons
   / G group).
