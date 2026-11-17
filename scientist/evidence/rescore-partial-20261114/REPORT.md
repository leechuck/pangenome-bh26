# Step 1: partial-credit re-scoring (stage rescore-partial-20261114)

Scoring only. No typing was run. All computation ran on IBEX through sbatch (account pi-hohndor, jobs `s1-*`); the stage is `$E/rescore-partial-20261114`, mirrored here. Existing stages were only read.
Full tables (every arm, subset, gene, contrast): `results/TABLES_confirm.md`, `TABLES_typer_T-sim.md`, `TABLES_typer_T-real.md`, `TABLES_drb1_{T-sim,T-real,D-sim,D-real}.md`. Row-level files (names, flags, numbers only): `results/rows_{confirm,typer_*,drb1_*}.tsv`. Code: `src/` (identical to IBEX).

## What was run
Final jobs (all COMPLETED, 0 failed; these are single jobs, not arrays, so there are no COMPLETE/FAILURE.json files; success = exit 0 and row counts below):
confirmatory rows 55584369; mosaic-typer rows 55584275 (T-sim), 55584276 (T-real); DRB1 rows 55584277 (T-sim), 55584278 (T-real), 55584279 (D-sim), 55584280 (D-real); tables 55584466; spot check 55584333. Earlier runs of the same code (55582747, 55582776-81, 55583789-94, 55583988, 55584274) were superseded after fixing the pairing tie-break and clearing gene-body values for the ipd arm (see "Changes during the work").
Row counts: confirmatory 2104 = 526 true copies x 4 arms; typer T-sim 4875 = 975 x 5 arms, T-real 5280 = 1056 x 5; DRB1 T-sim 488 = 122 x 4 (D4 of APR051 missing), T-real 528, D-sim 80, D-real 160.

## Metric as implemented (choices where the brief was silent)
- Unit: one true gene copy. Confirmatory: the 262 eligible donor-genes (2 copies each) = 526 copies in set A, 254 in set P (16 donors); the single truth-ineligible donor-gene is excluded, as in the frozen scorer. Typer: all 8 genes, 2 copies per donor-gene (T-sim 975, T-real 1056; one donor-gene has one copy in T-sim). DRB1: DRB1 only.
- Gene body = first coding-exon start to last coding-exon end (T1K header coordinates, **UTRs excluded**), as in the DRB1 stage's "gene body". Coding = concatenated coding exons. Distance = edlib NW edit distance. Non-coding = max(0, gene body - coding); the clip was never needed (0 negatives in every table).
- Exon coordinates: IPD alleles and cohort (PG) alleles: header coordinates. Novel/patched (MOS_, PG, patched) alleles: header coordinates of the `additions.fa` / reference the stage wrote. Mosaic-typer and DRB1 truth: `core_exons` of SOURCE_RECORDS. Confirmatory truth: Immuannot CDS + stop_codon features of that copy mapped into the core (TRUTH.json has no exons). Check: for all 416 truth copies whose core equals an IPD genomic allele the Immuannot-derived exons equal the IPD header exons exactly; the other 29 present copies (core differs from the allele ends) cannot be checked this way.
- Slot choice: calls are decoded with each frozen scorer's own code (T1K genotype table, copies==1 duplicates the slot, quality 0 = no call; typer/DRB1: acceptance rule and DRB1 overrides applied exactly as the frozen scorers; a slot holding named and novel alleles counts as the named one). **Tied alleles in one slot: the slot is represented by its worst allele** (largest gene-body distance), so a copy is "correct" only if all tied alleles are right. This is the strict rule of the frozen four-field scorers; the frozen novel-recovery code is optimistic ("any tied"). The best-allele variant is in every TABLES file ("Sensitivity to the tie rule"); it changes the gene-body-exact share by at most 1.9 pp (confirmatory nonasian arm, set A) and by 0-1.1 pp elsewhere.
- Pairing: the two slots are paired with the two true copies to minimise the total gene-body distance (coding distance for exon-only calls); ties are broken by more gene-body-identical copies, then smaller coding distance.
- Missing call (no allele in the slot): not correct in every share (denominator = all copies), excluded from distance means. Missing: confirmatory 0 in all arms; typer T-real 0; typer T-sim 2 copies in every arm (1 IPD-absent); DRB1 T-sim D4 2 copies (donor APR051, pass-2 run absent), all other DRB1 sets 0.
- ipd arm (confirmatory): `ipd.original.fa` mixes exon-only and truncated-intron sequences (for example a DPB1 allele with 2.3 kb against a 10 kb true gene body), so gene-body values are **not available** and cleared in the row file; only coding distance is reported.
- Contrasts: pairs by copy, family-clustered bootstrap (4000 draws, seed 20261020) and family sign-flip p, copied from the frozen `uncertainty()` (confirmatory families from PLAN.json/relatedness rule as in `confirm_score.py`, typer/DRB1 families from FOLDS). Share contrasts in percentage points over all copies of the subset (missing = not correct). Mean contrasts use copies called in both arms; negative = first arm closer to truth. No multiplicity adjustment.

## Validation (distance == 0 against the frozen reports)
Reproduced exactly:
- Confirmatory IPD-absent copies recovered (whole-sequence identical): set P genomic 0/35, genomic_core_hprc_asian **14/35**, nonasian 2/35, ipd 0/35; set A 0/81, **21/81**, 2/81, 0/81.
- Confirmatory four-field named (donor-genes with both copies gene-body identical vs frozen correct): set A genomic 185/197 (frozen 185), asian 184 (frozen 182), nonasian 179 (frozen 173); set P 93/93/90 (frozen 93/93/88). Frozen right and mine wrong: 0 donor-genes. All excess (2 asian, 6 nonasian in A; 2 nonasian in P) are gene C donor-genes where T1K lists two tied alleles with identical gene body; the frozen scorer calls these "not called". Whole-sequence identity is lower (163/156/152) because many truth cores are shorter than the IPD allele at the ends; that is why the comparison uses the gene body.
- Mosaic typer, IPD-absent recovery (whole-sequence identical in the final call): T-sim A 0, B 19, C 2, D 20, E 104 of 162; T-real A 0, B 28, C 2, D 29, E 85 of 134; equal to the frozen `recovered_final` in every arm (one exception below).
- DRB1: T-sim exact(final) A 0, B 16, C4 10, D4 21 (frozen 22); T-real A 0, B 14, C4 2, D4 15; D-sim C4 6; D-real B 4, C4 2, D4 6.
Mismatches left (cannot be removed without changing the brief's rules):
1. Four-field named on IPD-present copies (typer/DRB1): mine (gene body identical) is higher than the frozen name count by 11-26 copies per arm in the typer sets (out of 804 and 914): tied alleles with identical gene body (frozen: not correct) and different 4-field names with identical gene body (a 4th-field difference outside the coding-exon span, i.e. in the UTRs). Reverse direction (frozen right, mine wrong): 0 in T-sim typer, 1 per arm in T-real typer, 0 in DRB1 (the frozen scorer matches a name; I find a 1-2 base gene-body difference between the called IPD allele and the truth). Not further traced.
2. DRB1 gene-body "exact" (frozen counts a candidate in ANY slot): T-sim B frozen 17 vs mine 16, D4 26 vs 24; T-real B 15 vs 16 (mine higher, whole-core copies whose gene body equals an IPD allele), D4 17 vs 17. The T-sim B/D4 difference is copy APR013 copy 2 (slot with four tied PG alleles, one right: frozen credits "any", mine "all") and APR038 copy 2 in D4 (a tied PG|MOS slot that my pairing assigned to the other copy). 
3. Typer T-real arm E: 1 copy (HG02155 DPA1) frozen recovered_final but mine not exact; pairing difference between the frozen score-based and my distance-based assignment, not traced further.
4. IPD-absent copies whose gene body equals an IPD allele: A has 15/81 (confirmatory A), 34/162 (T-sim), 23/134 (T-real), 2/85 (DRB1 T-real) gene-body-exact IPD-absent copies without any novel call. These are copies whose core (allele ends/UTR) is absent from IPD but whose coding-exon span is not. This is expected under the new metric and is not a mismatch; it also means "IPD-absent" is a statement about the whole core.
- Spot check (`results/spot_check.txt`, lengths and distances only, recomputed from the sequences with an edlib alignment path): 5 copies, all 5 file distances equal the recomputed ones (e.g. HG03995 DRB1 nonasian: true gene body 10648, call 10650, distance 2 = 2; NA21133 DPB1 genomic: 10197 vs 10197, distance 0). The spot-check ipd copy (HG03788 DRB1, call body 2806 vs true 13009) shows why ipd gene-body values were removed (coding distance 3 = 3).

## Main results
Shares: denominator = all true copies of the subset. gb = gene body, cod = coding. "wrong bases" = mean (median) over called copies.

### Confirmatory test, 33 donors, real reads (set A: 526 copies; IPD-absent 81)
| arm | gb exact, all | cod exact, all | gb <=1 / <=5 / <=20, all | mean (median) wrong bases gb / cod / non-coding, all |
|---|---|---|---|---|
| ipd | NA | 465 (88.4%) | NA | NA / 4.1 (0) / NA |
| genomic | 443 (84.2%) | 517 (98.3%) | 89.7 / 96.4 / 99.6 % | 0.9 (0) / 0.0 (0) / 0.8 (0) |
| genomic_core_hprc_asian | 457 (86.9%) | 519 (98.7%) | 92.0 / 97.0 / 99.6 % | 0.8 (0) / 0.0 (0) / 0.8 (0) |
| genomic_core_nonasian | 431 (81.9%) | 513 (97.5%) | 89.4 / 96.4 / 99.6 % | 0.9 (0) / 0.0 (0) / 0.8 (0) |

IPD-absent copies (81) / IPD-present copies (445):
| arm | absent: gb exact | absent: cod exact | absent: gb <=1 / <=5 / <=20 | absent: wrong bases gb / cod | present: gb exact | present: wrong bases gb / cod |
|---|---|---|---|---|---|---|
| ipd | NA | 69 (85.2%) | NA | NA / 2.7 | NA | NA / 4.4 |
| genomic | 15 (18.5%) | 76 (93.8%) | 42.0 / 76.5 / 97.5 % | 5.2 (2) / 0.1 | 428 (96.2%) | 0.1 / 0.0 |
| genomic_core_hprc_asian | 33 (40.7%) | 79 (97.5%) | 58.0 / 80.2 / 97.5 % | 4.7 (1) / 0.0 | 424 (95.3%) | 0.1 / 0.0 |
| genomic_core_nonasian | 15 (18.5%) | 74 (91.4%) | 39.5 / 76.5 / 97.5 % | 5.0 (2) / 0.1 | 416 (93.5%) | 0.1 / 0.0 |

Gene with most error: DRB1 (all arms: 14-28 of 64 copies gene-body exact, mean 4.1-4.4 wrong bases); DQA1 next (42-52 of 66). Per-gene table in TABLES_confirm.md.

Contrasts, set A (set P in TABLES_confirm.md; same direction, wider intervals):
| contrast | subset | share (pp) | 95% CI | mean wrong bases difference | 95% CI |
|---|---|---|---|---|---|
| genomic - ipd (coding) | all | +9.89 | [+6.56, +13.31] | -4.11 | [-5.88, -2.74] |
| genomic - ipd (coding) | IPD-absent | +8.64 | [+1.35, +15.38] | -2.64 | [-4.81, -0.93] |
| asian - genomic (gb) | all | +2.66 | [+0.76, +5.32] | -0.06 | [-0.21, +0.08] |
| asian - genomic (gb) | IPD-absent | +22.22 | [+12.09, +34.29] | -0.49 | [-1.30, +0.48] |
| asian - genomic (gb) | IPD-present | -0.90 | [-1.64, 0.00] | +0.02 | [0.00, +0.04] |
| asian - nonasian (gb) | all | +4.94 | [+3.14, +7.14] | -0.05 | [-0.20, +0.10] |
| asian - nonasian (gb) | IPD-absent | +22.22 | [+12.09, +34.29] | -0.30 | [-1.15, +0.70] |
Set P: asian - genomic (gb, all) +4.33 pp [+1.57, +7.87]; IPD-absent +31.43 pp [+13.33, +48.65]; IPD-present +0.00. asian - nonasian: +6.69 pp [+3.57, +10.24]. genomic - ipd coding: +10.24 pp [+5.12, +15.60].
Reading: the cohort panel makes IPD-absent copies fully right more often, but the typical error of the other copies is unchanged: mean wrong bases per IPD-absent copy fall only 5.2 -> 4.7 and the interval includes 0. Coding errors are rare (<2 % of copies) in all genomic arms. The ipd arm is wrong in the coding sequence for 11.6 % of copies (mean 4.1 bases).

### Mosaic-typer test (all 8 genes; gb exact, all copies | IPD-absent: gb exact, mean wrong bases gb)
| arm | T-sim all (975; 2 missing) | T-sim IPD-absent (162) | T-real all (1056) | T-real IPD-absent (134) |
|---|---|---|---|---|
| A | 757 (77.6%) | 34 (21.0%), 30.7 | 902 (85.4%) | 23 (17.2%), 4.6 |
| B | 784 (80.4%) | 62 (38.3%), 30.2 | 923 (87.4%) | 46 (34.3%), 3.5 |
| C | 762 (78.2%) | 40 (24.7%), 30.6 | 905 (85.7%) | 26 (19.4%), 4.5 |
| D | 790 (81.0%) | 67 (41.4%), 30.2 | 925 (87.6%) | 48 (35.8%), 3.5 |
| E | 849 (87.1%) | 129 (79.6%), 0.8 | 968 (91.7%) | 89 (66.4%), 2.7 |
Coding exact (all): T-sim 94.8 / 94.7 / 94.8 / 94.8 / 96.4 %; T-real 98.3 / 98.4 / 98.5 / 98.5 / 99.0 % (A..E). T-sim IPD-absent copies carry large errors (mean 30 bases, a long tail from simulated Arab alleles; 10 copies of 162 are off by more than 20); T-real errors are small (97.8 % within 20 bases).
Contrast B - A: T-sim all +2.77 pp [+1.49, +3.90] (gb exact), mean wrong bases -0.06 [-0.13, +0.03]; IPD-absent +17.28 pp [+11.46, +22.93], mean -0.42 [-0.74, -0.01]. T-real all +1.99 pp [+1.04, +3.03], mean -0.15 [-0.21, -0.09]; IPD-absent +17.16 pp [+10.07, +24.62], mean -1.07 [-1.52, -0.63]. IPD-present copies: -0.12 pp [-0.87, +0.63] (T-sim), -0.22 pp [-0.54, 0.00] (T-real): no cost in gene-body terms. D - B and E - A in the TABLES files (E - A IPD-absent: +58.64 pp T-sim, +49.25 pp T-real).

### DRB1 patch builder (DRB1 only; gb exact | IPD-absent gb exact | mean wrong bases gb: all / IPD-absent)
| arm | T-sim all (122) | T-sim absent (76) | T-sim mean all / absent | T-real all (132) | T-real absent (85) | T-real mean all / absent |
|---|---|---|---|---|---|---|
| A | 35 (28.7%) | 0 | 31.1 / 49.4 | 33 (25.0%) | 2 (2.4%) | 4.6 / 6.5 |
| B | 51 (41.8%) | 16 (21.1%) | 30.6 / 48.6 | 46 (34.8%) | 16 (18.8%) | 3.6 / 5.0 |
| C4 | 54 (44.3%) | 14 (18.4%) | 30.3 / 48.2 | 44 (33.3%) | 6 (7.1%) | 3.5 / 5.0 |
| D4 | 63 (51.6%) | 24 (31.6%) | 30.3 / 48.3 (2 missing) | 56 (42.4%) | 17 (20.0%) | 2.9 / 4.1 |
Development sets (labelled development, not test): D-sim (40 copies): A 10 -> C4 27 gene-body exact; D-real (40): A 8, B 12, C4 15, D4 19.
Contrasts, gene-body exact share (pp) and mean wrong bases, all copies: 
| contrast | T-sim | T-real |
|---|---|---|
| C4 - A | +15.57 [+9.33, +23.15]; -0.85 [-1.13, -0.54] | +8.33 [+3.79, +12.88]; -1.07 [-1.55, -0.62] |
| D4 - B | +9.84 [+4.17, +16.36]; -0.70 [-1.08, -0.39] | +7.58 [+3.79, +12.12]; -0.72 [-1.14, -0.36] |
IPD-absent only: C4 - A +18.42 pp [+10.26, +28.38] (T-sim), +4.71 pp [+1.12, +9.52] (T-real, sign-flip p 0.125); D4 - B +10.53 pp [+3.08, +18.31] (T-sim), +1.18 pp [0.00, +3.75] (T-real).
Coding is almost unchanged (+1 to +3 copies); the gain is in the non-coding part.

Novel (patched) final calls vs the arm they replace (gene-body distance to truth; change = new minus old, negative = closer):
| set | arm pair | copy type | novel final calls | closer | equal | farther | not exact: closer / equal / farther |
|---|---|---|---|---|---|---|---|
| T-sim | C4 vs A | IPD-absent | 28 | 28 | 0 | 0 | 14 / 0 / 0 |
| T-sim | C4 vs A | IPD-present | 1 | 1 | 0 | 0 | 1 / 0 / 0 |
| T-sim | D4 vs B | IPD-absent | 45 | 17 | 27 | 1 | 8 / 12 / 1 |
| T-sim | D4 vs B | IPD-present | 4 | 1 | 3 | 0 | 1 / 2 / 0 |
| T-real | C4 vs A | IPD-absent | 36 | 35 | 0 | 1 | 31 / 0 / 1 |
| T-real | D4 vs B | IPD-absent | 53 | 24 | 28 | 1 | 23 / 14 / 1 |
| T-real | D4 vs B | IPD-present | 2 | 0 | 2 | 0 | 0 / 1 / 0 |
| D-sim | C4 vs A | IPD-absent | 17 | 17 | 0 | 0 | 3 / 0 / 0 |
| D-real | C4 vs A | IPD-absent | 6 | 6 | 0 | 0 | 3 / 0 / 0 |
| D-real | D4 vs B | IPD-absent | 12 | 6 | 6 | 0 | 2 / 4 / 0 |
(Last column: the subset where the new call is not exactly right, distance > 0.) Size of the change, new minus old, gene-body wrong bases: C4 vs A, T-sim IPD-absent median -2 (range -11 to -1); T-real median -2 (range -14 to +12; the one worse copy is +12). D4 vs B: median 0 (equal when the D4 call is the cohort allele arm B already had), range -14 to +12. Full distribution bins are in TABLES_drb1_*.md.
Summary for novel calls that are not exactly right (test sets): against A (C4) 46 of 47 are closer (T-sim 15 of 15 incl. 1 IPD-present copy, T-real 31 of 32); against B (D4) 32 are closer, 29 equal and 2 farther, of 63 (T-sim 9 / 14 / 1; T-real 23 / 15 / 1).

## Answers to the three questions
(a) Wrong bases per copy: for real reads almost all genomic-arm copies are right in the coding sequence (97.5-99 %), and the mean error is below 1 base per copy (confirmatory) with the whole error in non-coding positions. The mean is carried by IPD-absent copies (about 5 bases per copy; 80 % within 5 bases for the best arm). The simulated Arab set is worse for IPD-absent copies (mean about 30, in a few large errors).
(b) With all copies in the denominator, exact gene-body shares are 78-92 % over all copies (confirmatory 82-87 %, typer T-sim 78-87 %, T-real 85-92 %), so the old exactness numbers on IPD-present copies were not hiding larger errors; the IPD-absent copies are where exactness (17-41 % for arms A-D) and distances differ.
(c) Novel calls are closer to the truth than the standard call in almost every case, also when they are not exactly right (see table above).

## Changes during the work (all recorded)
- First run used whole-sequence exactness as a tie-break; a copy pair where two assignments had the same total distance (DRB1 APR049) was assigned differently from the frozen scorer. Tie-break changed to gene-body-identical first; everything rerun.
- ipd-arm gene-body values were first stored; cleared after the spot check showed truncated introns.
- Gene-body exactness in confirmatory validation first compared whole sequences; switched to gene body because many truth cores are shorter than the IPD allele at the ends (documented above).

## Not verified / open
- The unexplained single-copy mismatches listed in "Mismatches left" (items 1 and 3).
- The ipd arm: whether its truncated sequences could be completed from `genomic.fa` for a gene-body comparison was not done (not requested).
- The Immuannot exon mapping was checked only on copies identical to an IPD allele (416 of 445 present copies).
- D and E arms of the typer were scored for T-sim and T-real only; D-sim and D-real typer sets were not run (not requested).
- Contrast p-values are exploratory and unadjusted (about 100 contrasts).
- Simulated and real read sets are never pooled; development (D-) sets are shown only in TABLES_drb1_D-*.md and the one line above.
