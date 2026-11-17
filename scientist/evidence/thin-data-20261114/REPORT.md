# Step 6: thin-data test (stage thin-data-20261114), 2026-11-14

Question: when reads are few, do a matched-population reference panel (arm B) and matched-population allele
frequencies (rule P_matched) help more than at full depth? All rules were fixed in advance; nothing was tuned.

## Answer in short
No. Neither the panel nor the population prior gives a clear gain when reads are thin.
- Panel (B minus A, P0f rule, two-field): full depth real -0.2 pp, 1/16 real -1.2 pp (95% -2.1 to -0.2); Middle East simulated: full -0.1, 1/16 -0.5 (-2.0 to 0.8). Four-field is lower with the panel at every depth (real -1.8 pp at full, -0.8 at 1/16; Middle East -1.6 and -1.1). Interaction (is B minus A larger at 1/16 than at full), two-field: real -0.9 pp (-1.9 to 0.0), Middle East -0.4 (-2.0 to 1.0). No evidence that it grows at low depth.
- The panel does recover sequences absent from IPD, but less as depth falls: real 22.3% of absent copies at full depth, 5.7% at 1/16; Middle East 11.7% to 3.7%.
- Population prior (P_matched minus P0f): zero at full depth and 1/4 depth for two-field; at 1/16 depth real +0.3 pp (-0.1 to 0.7), Middle East +1.2 pp (0.3 to 2.1). Four-field at 1/16: real +1.2 (0.2 to 2.2), Middle East +0.9 (-0.7 to 2.2). P_matched minus P_other (does matching matter): almost all intervals include 0 (only borderline exception seen: arm B, real, four-field, 1/16: +1.0, 0.1 to 2.0). So the effect is small, and mostly not specific to the matching population.
- What does collapse at low depth is the share of ambiguous calls (T1K reports several alleles): 0.8% of slots at full depth, 51% at 1/16 (real). Counting an ambiguous slot as correct only if all listed alleles agree (the frozen scorer) gives 51% two-field at 1/16; picking any single listed allele gives 75%. That choice (+24 pp) is far larger than any panel or prior effect.

## What was run
Scripts: `src/` (identical to IBEX `$E/thin-data-20261114/src`). Jobs (all `--account=pi-hohndor`, prefix `s6-`):
tabulation 55582731; pilot array 55583252 (4 donors); scorer check 55583268, 55583740; nested check 55583996;
main array 55584512 (143 donors, %12); scoring 55592263; summaries 55592264, then 55592412 (rerun after fixing one metric, see below).
- Donors: 147 = 86 real (D-real 20 + T-real 66; EAS 50, SAS 36) + 61 simulated Arab/Saudi (T-sim; Middle East, ME). Real and simulated never pooled. D-real and T-real are pooled as "real". Simulated copies of HPRC donors were not used.
- Reads: the full recruited pair files of the arm A run (r1/r2.fq.gz, same files in all arms; checked). 482k to 1.12M pairs per donor.
- Subsampling: pair kept if sha256(seed 20261114 : read name without /1,/2) mod 16 < 8, 4, 2, 1. Pilot check on 4 donors: mates match, files equal the rule exactly, levels nested (16 inside 8 inside 4 inside 2), fractions 0.4994-0.5003, 0.2494-0.2504, 0.1247-0.1252, 0.0623-0.0625 (`results/pilot_nested_check.json`). Subsampled FASTQ of the 143 main donors were deleted after typing (sha256 and pair counts kept in `work/reads/*/*/SUBSAMPLE.json`); the 4 pilot donors keep theirs.
- T1K: exactly `run_donor.t1k` (hla-wgs, alleleDigitUnits 4, delimiter :), 4 threads, container ubuntu24.sif. Arms: A = IPD genomic.fa; B = existing `genomic_core_hprc_asian` reference of the donor's fold; N = `genomic_core_nonasian` reference of the fold, at 1/4 and 1/16 (exists for all 147 donors); M (added, see deviations) = `genomic_core_arab_saudi`, ME donors only. Full depth = existing runs (table sha verified against receipts).
- Counts: 1714 thin T1K runs complete, 0 failed (86 x 10 + 61 x 14); 147 of 147 donor tasks COMPLETE.json, 0 FAILURE.json. No table missing at scoring (`results/rows_meta.json`: missing_tables = []). Task time: 30.0 core-task hours in total (mean 735 s per donor task; the 4 pilot tasks 1105-2625 s because they also reran full depth). The pilot projected about 4 h wall for the array at %12; it took about 3.3 h.

## Sanity checks
1. Reproduction of existing full-depth tables with my wrapper (pilot: 2 real, 2 simulated donors; arms A, B, N, M, 14 reruns, 4 threads, including simulated donors whose originals used 2 or 8 threads): 14 of 14 genotype tables byte-identical. T1K was deterministic here; not tested beyond these 4 donors.
2. Scorer: my slot handling with rule P0s reproduces the frozen scorer (`score.score`) on all 4702 full-depth copy rows of arms A and B: 0 differences in 2-, 3-, 4-field correctness (`results/scorer_reproduction.json`).
3. Pair counts / nested files: see above.

## Deviations and choices (please read)
- P0. The brief says P0 = first listed = the frozen scorer's rule. The frozen scorer does not do that: a field level counts as correct only if ALL listed alleles agree at that level (`part0.slot_value`). I therefore report both: P0f = first listed allele (primary comparison for the prior, as the brief states), and P0s = frozen unanimity rule. At full depth they differ by under 1 pp.
- Accuracy definitions are the frozen ones: four-field over IPD-present copies with a single four-field name; three-field and two-field "decoded" (a panel allele counts through its CDS/protein-identical IPD name), over copies with an IPD allele of identical CDS (3-field) or protein (2-field). Slots are assigned to copies by the frozen best-permutation rule (uses truth), separately for each rule.
- Prior choice is among ALL listed alleles, panel (PG) alleles included. Frequencies: truth `genomic_labels` of the group's donors at 4, then 3, then 2 fields; ties and all-zero go to first listed. For a panel allele the 4-field count is the count of copies with the identical sequence. Own family component excluded (FOLDS groups; the 9 Saudi donors form one block, so Saudi donors use Arab donors only). Dev donors' truth counts are included in the tables.
- Arm M added (not in the brief): for Middle East donors arm B (HPRC Asian panel) is not the matched population, so I also ran the existing Arab/Saudi leave-family-out panel at all depths. For ME donors, "B" below is therefore the HPRC-Asian panel as specified.
- Ambiguity share, no-call share are measured on the slot assigned to each truth copy (about 2 per gene and donor); "mean list length" is the mean number of listed alleles among ambiguous slots (can be long because T1K lists all tied alleles).
- "Fully correct sequence" = every allele of the (chosen) call has a sequence byte-identical to the truth core sequence. "Recovered exactly" = same, restricted to IPD-absent copies, for the chosen call (arm A is 0 by construction).
- Depth estimate (truth-free): median over the 8 genes of (sum of T1K abundance of the called slots)/2, and median T1K candidate and assigned fragments. T1K abundance units are not verified to be read depth; use the ratios between levels, not absolute values (`results/depth_summary.tsv`).
- Intervals: cells use a family-clustered bootstrap (4000 draws, seed 20261020, numpy, same procedure as the frozen scorer but not the same random stream); contrasts use the frozen `uncertainty()` function itself (paired by copy, family-clustered bootstrap and sign-flip p). Interaction = per-copy difference at 1/16 minus per-copy difference at full, same function.
- A first summary run had a bug in the "recovered absent" metric (counted all copies in the numerator); fixed and rerun (job 55592412); all tables here come from the fixed run.
- Arm N and M tables at 1/2 and 1/8 are not available by design (N only at 1/4 and 1/16).

## Estimated depth (arm A, median over donors)
| read type, group | level | median pairs | T1K assigned fragments | T1K abundance depth |
|---|---|---|---|---|
| real EAS | full / 1/2 / 1/4 / 1/8 / 1/16 | 675,736 / 337,758 / 168,865 / 84,303 / 42,214 | 23,408 / 11,691 / 5,853 / 2,856 / 1,421 | 54.1 / 26.8 / 12.4 / 6.0 / 3.2 |
| real SAS | same | 679,636 / 339,272 / 169,492 / 84,813 / 42,365 | 24,114 / 12,097 / 5,948 / 2,898 / 1,441 | 57.0 / 26.4 / 12.7 / 5.8 / 2.9 |
| sim ME | same | 554,612 / 277,505 / 138,348 / 69,217 / 34,553 | 18,556 / 9,241 / 4,565 / 2,201 / 1,085 | 34.1 / 15.6 / 7.5 / 3.3 / 2.1 |

## Accuracy by depth (percent, 95% family-clustered interval)
Rules: P0f first listed; PM matched population; PO other populations; P0s frozen unanimity. EAS and SAS separately, and all arms/metrics, are in `results/figure_table.tsv`; per-gene two-field tables in `results/pergene_acc2.tsv`.

real reads (REAL), two-field accuracy %, 95% interval

| arm, rule | full | 1/2 | 1/4 | 1/8 | 1/16 |
|---|---|---|---|---|---|
| A, P0f | 99.0 (98.5-99.4) | 98.5 (97.7-99.1) | 95.8 (94.8-96.8) | 83.6 (81.7-85.5) | 74.9 (72.7-77.0) |
| A, PM | 99.0 (98.5-99.4) | 98.5 (97.7-99.1) | 95.8 (94.8-96.8) | 83.6 (81.7-85.5) | 75.2 (73.0-77.3) |
| A, PO | 99.0 (98.5-99.4) | 98.5 (97.7-99.1) | 95.8 (94.8-96.8) | 83.5 (81.6-85.4) | 75.0 (72.8-77.1) |
| A, P0s | 99.0 (98.5-99.4) | 98.3 (97.6-99.0) | 93.5 (92.1-94.7) | 72.0 (69.7-74.5) | 51.1 (48.3-54.1) |
| B, P0f | 98.8 (98.2-99.3) | 98.5 (97.9-99.1) | 95.7 (94.7-96.7) | 83.6 (81.8-85.6) | 73.7 (71.5-76.0) |
| B, PM | 98.8 (98.2-99.3) | 98.5 (97.9-99.1) | 95.7 (94.7-96.7) | 83.6 (81.7-85.5) | 74.2 (71.9-76.4) |
| B, PO | 98.8 (98.2-99.3) | 98.5 (97.9-99.1) | 95.7 (94.7-96.7) | 83.5 (81.6-85.5) | 74.0 (71.7-76.2) |
| B, P0s | 98.8 (98.2-99.3) | 98.4 (97.7-99.1) | 93.5 (92.2-94.8) | 72.2 (69.9-74.7) | 49.7 (46.7-52.9) |

real reads (REAL), four-field accuracy %, 95% interval

| arm, rule | full | 1/2 | 1/4 | 1/8 | 1/16 |
|---|---|---|---|---|---|
| A, P0f | 94.2 (92.9-95.4) | 92.4 (90.7-93.9) | 85.7 (83.6-87.8) | 69.5 (66.6-72.2) | 55.0 (51.8-58.1) |
| A, PM | 94.2 (92.9-95.4) | 92.4 (90.7-93.9) | 85.9 (83.8-88.0) | 69.9 (67.1-72.6) | 56.2 (53.1-59.2) |
| A, PO | 94.2 (92.9-95.4) | 92.4 (90.7-93.9) | 85.8 (83.7-87.9) | 69.7 (66.9-72.4) | 55.5 (52.4-58.5) |
| A, P0s | 93.6 (92.2-95.0) | 90.5 (88.5-92.5) | 76.9 (74.3-79.5) | 44.2 (40.9-47.6) | 19.7 (16.9-22.5) |
| B, P0f | 92.4 (91.0-93.6) | 91.1 (89.3-92.8) | 84.3 (82.0-86.5) | 68.8 (66.0-71.6) | 54.2 (51.1-57.3) |
| B, PM | 92.4 (91.0-93.6) | 91.0 (89.2-92.7) | 84.5 (82.2-86.7) | 69.2 (66.4-71.9) | 55.5 (52.6-58.5) |
| B, PO | 92.4 (91.0-93.6) | 91.1 (89.3-92.8) | 84.4 (82.2-86.6) | 68.9 (66.1-71.7) | 54.5 (51.4-57.7) |
| B, P0s | 91.9 (90.5-93.3) | 89.3 (87.2-91.4) | 75.7 (73.1-78.4) | 44.1 (40.9-47.5) | 19.1 (16.3-21.9) |

sim reads (ME), two-field accuracy %, 95% interval

| arm, rule | full | 1/2 | 1/4 | 1/8 | 1/16 |
|---|---|---|---|---|---|
| A, P0f | 97.1 (95.7-98.2) | 96.5 (94.9-97.8) | 85.5 (83.1-87.8) | 67.3 (64.3-70.0) | 60.3 (57.5-63.6) |
| A, PM | 97.1 (95.7-98.2) | 96.5 (94.9-97.8) | 85.5 (83.1-87.8) | 67.4 (64.5-70.0) | 61.5 (58.7-64.8) |
| A, PO | 97.1 (95.7-98.2) | 96.5 (94.9-97.8) | 85.5 (83.1-87.8) | 67.6 (64.6-70.2) | 61.0 (58.2-64.4) |
| A, P0s | 97.1 (95.7-98.2) | 96.2 (94.5-97.5) | 81.4 (78.5-83.9) | 52.7 (49.3-55.5) | 37.0 (34.4-39.5) |
| B, P0f | 97.0 (95.4-98.2) | 96.3 (94.8-97.4) | 85.9 (83.4-88.4) | 68.9 (65.9-71.5) | 59.8 (57.1-62.7) |
| B, PM | 97.0 (95.4-98.2) | 96.3 (94.8-97.4) | 85.9 (83.4-88.4) | 69.2 (66.3-71.7) | 60.9 (58.0-64.0) |
| B, PO | 97.0 (95.4-98.2) | 96.3 (94.8-97.4) | 85.8 (83.1-88.4) | 69.3 (66.3-71.9) | 60.7 (57.8-63.9) |
| B, P0s | 97.0 (95.4-98.2) | 96.0 (94.5-97.2) | 81.5 (78.4-84.2) | 54.1 (50.8-56.7) | 37.3 (34.2-40.0) |

sim reads (ME), four-field accuracy %, 95% interval

| arm, rule | full | 1/2 | 1/4 | 1/8 | 1/16 |
|---|---|---|---|---|---|
| A, P0f | 88.2 (86.0-90.0) | 86.4 (83.8-88.6) | 68.2 (65.1-71.4) | 44.4 (41.4-47.3) | 39.3 (35.8-43.4) |
| A, PM | 88.2 (86.0-90.0) | 86.4 (83.8-88.6) | 68.2 (65.3-71.2) | 45.3 (42.2-48.1) | 40.2 (37.0-44.0) |
| A, PO | 88.2 (86.0-90.0) | 86.4 (83.8-88.6) | 68.2 (65.2-71.3) | 45.5 (42.5-48.3) | 40.3 (37.0-44.3) |
| A, P0s | 87.9 (85.6-89.8) | 84.6 (81.6-87.2) | 57.1 (53.8-60.3) | 26.2 (23.5-29.0) | 13.2 (10.9-15.9) |
| B, P0f | 86.6 (84.0-88.6) | 84.3 (81.6-86.5) | 66.2 (63.0-69.4) | 44.7 (41.7-47.6) | 38.2 (34.9-42.1) |
| B, PM | 86.4 (83.9-88.5) | 84.1 (81.4-86.3) | 65.8 (62.8-68.8) | 45.3 (42.2-48.0) | 38.8 (35.7-42.5) |
| B, PO | 86.6 (84.0-88.6) | 84.3 (81.6-86.5) | 65.9 (62.9-69.0) | 45.6 (42.6-48.5) | 39.6 (36.3-43.4) |
| B, P0s | 86.3 (83.8-88.4) | 82.7 (79.9-85.1) | 55.3 (52.0-58.6) | 26.0 (23.3-28.7) | 12.7 (10.5-15.3) |

Recovery of IPD-absent copies, exact sequence (rule P0f), %:

| read type, arm | full | 1/2 | 1/4 | 1/8 | 1/16 |
|---|---|---|---|---|---|
| real, B | 22.3 (15.8-28.9) | 19.4 (13.4-25.8) | 16.6 (10.7-22.8) | 12.0 (7.0-17.5) | 5.7 (2.4-9.6) |
| real, N | 2.3 (0.0-5.3) | NA | 2.3 (0.0-5.3) | NA | 0.0 (0.0-0.0) |
| sim, B | 11.7 (7.2-15.8) | 11.7 (7.2-15.8) | 10.5 (5.9-14.8) | 8.0 (4.2-12.2) | 3.7 (1.4-6.2) |
| sim, M | 8.6 (5.1-12.4) | 8.0 (4.5-11.7) | 8.6 (5.0-12.5) | 5.6 (2.1-10.1) | 1.9 (0.0-4.0) |
| sim, N | 4.3 (1.4-7.8) | NA | 2.5 (0.6-4.7) | NA | 0.6 (0.0-1.6) |

No-call share and ambiguity share (A, P0f), %:

| read type, metric | full | 1/2 | 1/4 | 1/8 | 1/16 |
|---|---|---|---|---|---|
| real, A nocall | 0.0 (0.0-0.0) | 0.0 (0.0-0.0) | 0.2 (0.0-0.5) | 3.4 (2.3-4.6) | 7.3 (6.0-8.8) |
| real, A ambig | 0.8 (0.2-1.5) | 2.4 (1.4-3.5) | 10.5 (8.7-12.3) | 31.2 (28.5-34.0) | 50.8 (47.5-54.1) |
| real, A listlen | 12.0 (9.3-16.4) | 11.6 (9.7-14.1) | 9.1 (8.2-10.1) | 6.1 (5.6-6.7) | 8.7 (7.8-9.7) |
| real, A seqok_all | 78.1 (76.5-79.6) | 76.5 (74.9-78.2) | 71.1 (69.0-73.2) | 57.4 (54.9-59.9) | 45.6 (42.8-48.3) |
| sim, A nocall | 0.2 (0.0-0.5) | 0.2 (0.0-0.6) | 2.2 (1.3-3.2) | 13.0 (11.2-14.9) | 16.1 (13.7-18.8) |
| sim, A ambig | 1.1 (0.5-1.9) | 3.4 (2.3-4.6) | 13.6 (11.8-15.8) | 30.2 (26.8-33.5) | 41.6 (38.2-45.5) |
| sim, A listlen | 27.5 (13.1-58.1) | 17.5 (13.0-24.5) | 8.7 (7.8-9.6) | 6.7 (6.2-7.4) | 16.7 (13.2-21.3) |
| sim, A seqok_all | 64.5 (61.7-66.8) | 62.6 (59.4-65.1) | 49.7 (47.0-52.4) | 32.8 (30.3-35.1) | 28.8 (26.4-31.9) |

## Contrasts (percentage points, paired by copy, 95% interval; full list in `results/contrasts.tsv`)
Two-field accuracy, real (EAS+SAS) / Middle East simulated:
| contrast | full | 1/4 | 1/16 |
|---|---|---|---|
| B - A (P0f) real | -0.2 (-0.6, 0.1) | -0.1 (-0.7, 0.4) | -1.2 (-2.1, -0.2) |
| B - A (P0f) ME | -0.1 (-0.9, 0.5) | +0.4 (-0.4, 1.4) | -0.5 (-2.0, 0.8) |
| M - A (P0f) ME | 0.0 (-0.6, 0.5) | +0.6 (-0.2, 1.6) | -0.5 (-1.6, 0.5) |
| A: PM - P0f, real | 0.0 | 0.0 | +0.3 (-0.1, 0.7) |
| A: PM - P0f, ME | 0.0 | 0.0 | +1.2 (0.3, 2.1) |
| B: PM - P0f, real | 0.0 | 0.0 | +0.4 (0.0, 0.9) |
| B: PM - P0f, ME | 0.0 | 0.0 | +1.1 (0.1, 2.1) |
| A: PM - PO, real | 0.0 | 0.0 | +0.2 (-0.1, 0.6) |
| A: PM - PO, ME | 0.0 | 0.0 | +0.5 (-0.2, 1.2) |
Four-field accuracy: B - A real -1.8 (-2.6, -1.1) full, -1.4 (-2.3, -0.7) at 1/4, -0.8 (-1.8, 0.1) at 1/16; ME -1.6 / -2.0 / -1.1. A: PM - P0f real +1.2 (0.2, 2.2) at 1/16, ME +0.9 (-0.7, 2.2). PM - PO real +0.8 (-0.2, 1.8) at 1/16.
Interaction (1/16 minus full), two-field: B - A real -0.9 (-1.9, 0.0), ME -0.4 (-2.0, 1.0); PM - P0f (arm A) real +0.3 (-0.1, 0.7), ME +1.2 (0.3, 2.1); arm B real +0.4 (0.0, 0.9), ME +1.1 (0.1, 2.1). Four-field: B - A real +1.0 (-0.1, 2.1), ME +0.5 (-1.1, 1.9); PM - P0f arm B real +1.3 (0.3, 2.4), ME +0.7 (-0.8, 2.0).
Exact recovery of IPD-absent copies, B - A at 1/16 minus at full: real -16.6 pp (-22.2, -11.5), ME -8.0 (-12.0, -4.1) (B recovers fewer absent copies when reads are thin); panel N recovers almost none (real 2.3% to 0%; ME 4.3% to 0.6%).
Many contrasts were examined (about 2,500 rows); none has a multiplicity correction, so isolated intervals at the edge of 0 should not be over-read. The ME prior effect at 1/16 (+1.2 pp, 61 donors) is the only one that holds in both arms; its non-matching-population control (PM - PO, +0.5, -0.2 to 1.2) is not clearly different from 0.

## Files
Local (`scientist/evidence/thin-data-20261114/`): `REPORT.md`, `src/`, `results/rows.tsv` (141,752 row-level lines: donor, group, read type, depth, arm, rule, gene, haplotype, slot, listed-allele count, chosen allele(s), correct at 2/3/4 fields, no-call, seq-correct; names and flags only), `figure_table.tsv` (read type x group x depth x arm x rule x metric with interval), `contrasts.tsv`, `pergene_acc2.tsv`, `depth_summary.tsv`, `depth_rows.tsv`, `donor_table.json`, `tables_headline.md`, `tables_accuracy.md`, check files. IBEX: `/encrypted/e3008/hohndor-dogohla-20261020/thin-data-20261114/` (work/ has all T1K outputs, 428 MB).

## Not verified / open
- Determinism of T1K only checked on 4 donors (14 runs). Thin-run tables for other donors were not rerun.
- T1K abundance as depth unit is an assumption; no independent coverage measurement (e.g. read mapping) was made.
- The ambiguity-list order (what "first listed" means) is whatever T1K writes; P0f depends on it.
- Per-gene tables exist but were not examined for gene-specific stories (e.g. DPA1, DRB1).
- ME reads are simulated; any real-read statement does not carry over.
