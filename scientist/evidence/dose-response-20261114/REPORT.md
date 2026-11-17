# Step 3: dose-response of a same-population haplotype panel (stage `dose-response-20261114`)

Date 2026-11-14/15. IBEX stage `/encrypted/e3008/hohndor-dogohla-20261020/dose-response-20261114`; local mirror here (no sequences).

## Short answer

Typing tool and settings fixed (T1K, exactly as arms A/B); only the reference content changed (IPD genomic alleles plus a panel of
cohort gene sequences drawn from whole family components, never the test donor's own). Test donors: 147 (ME 61 simulated reads, EAS 50 and SAS 36 real reads).
One random draw per (donor, cell). IPD-absent copies: ME 162, EAS 99, SAS 76. IPD only recovers 0 of them in every group.

Byte-exact recovery of IPD-absent copies by panel size (matched population; family-clustered 95% interval):

| n requested | ME (sim) | EAS (real) | SAS (real) |
|---|---|---|---|
| 0 (IPD only, existing arm A) | 0/162 = 0.0% | 0/99 = 0.0% | 0/76 = 0.0% |
| 5 | 1.9% [0.0, 5.1] | 0.0% [0.0, 0.0] | 14.5% [6.0, 24.0] |
| 10 | 3.1% [0.5, 6.7] | 2.0% [0.0, 5.0] | 11.8% [3.9, 21.2] |
| 20 | 4.3% [1.3, 8.2] | 5.1% [1.0, 10.9] | 28.9% [16.9, 40.7] |
| 40 | 9.3% [5.1, 14.0] | 7.1% [2.1, 13.3] | 32.9% [20.6, 45.2] |
| ceiling: copies whose exact sequence is in the n=40 panel | 23.5% | 21.2% | 46.0% |
| conversion at n=40 (recovered / supplied) | 15/38 = 39% | 7/21 = 33% | 25/35 = 71% |
| existing arm B (whole fold training panel; context only) | 12.3% | 11.1% | 38.2% |

Same n, matched vs other-population panel (byte-exact recovery; for n=40 the paired difference in percentage points [95%]):

| test group | matched n=40 | from ME | from EAS | from SAS | from non-Asian HPRC |
|---|---|---|---|---|---|
| ME (sim) | 9.3% | - | 1.9% (diff +7.4 [3.6, 11.7]) | 4.9% (+4.3 [-0.7, 10.4]) | 3.7% (+5.6 [0.6, 11.2]) |
| EAS (real) | 7.1% | 0.0% (+7.1 [2.1, 13.3]) | - | 3.0% (+4.0 [-1.9, 10.7]) | 0.0% (+7.1 [2.1, 13.3]) |
| SAS (real) | 32.9% | 14.5% (+18.4 [5.4, 31.6]) | 10.5% (+22.4 [9.5, 35.1]) | - | 2.6% (+30.3 [17.6, 43.7]) |

(n=40 source SAS has only 36 donors, so SAS-sourced cross panels have 36; SAS matched n=40 has 35. See deviations.)

Reading, in plain terms:
- A panel from the right population recovers unseen alleles; IPD alone recovers none. The curve rises with n in all three groups, steeply in SAS (14.5% to 32.9%), more slowly in ME and EAS (reaching 9.3% and 7.1% at n=40). The rise from n=5 to n=40 is significant in ME (+7.4 pp [3.6, 11.5]) and SAS (+18.4 [7.9, 29.9]) and EAS (+7.1 [2.1, 13.3]); single steps (5 to 10, 10 to 20) are within noise and in SAS n=10 is below n=5 (different random draws; the ceiling also dips, 21.1% to 19.7%).
- Matched is at or above a same-size panel from another population in nearly every comparison by exact sequence (the one exception: EAS donors, n=10, SAS panel 3 vs matched 2 copies, not significant), strongly for SAS (all three comparisons significant) and for ME vs EAS, EAS vs ME / non-Asian. Matched vs a neighbouring population is not significant for EAS vs SAS and ME vs SAS (exact sequence). The non-Asian HPRC panels supply the least.
- At gene-body level (ignoring allele ends) a SAS panel does as well as or better than the matched ME panel for ME donors (16.7% vs 9.9% at n=40; difference -6.8 pp [-14.9, 0.0]); so the ME advantage at exact-sequence level is partly about allele ends and about which alleles are supplied. Gene-body tables T2, T3b.
- Conversion (recovered/supplied) is 33% to 71% at n=40; supply is the main limit, then T1K's choice (T1K misses supplied copies, notably ME DQA1 and DPA1 in simulated reads: ME DQA1 supplied 19 recovered 2; DPA1 supplied 4 recovered 0).
- Safety: IPD-present four-field named accuracy is unchanged within about 0.4 pp for all cells except ME matched n=40: -1.0 pp [-1.76, -0.14] (2 gains, 10 losses of 804). Whole-set "all copies fully correct" rises for SAS (82.3% to 86.4% at n=40) and EAS (81.1% to 81.8%), is flat for ME (73.2% to 73.7%).
- Real reads only: split by read set (T-real vs D-real) is in T11; the direction is the same in both.

## What was run

| step | job | outcome |
|---|---|---|
| inventory, draws, truth-only supply (steps 1, 2) | 55582812 (first), 55583020 (superseded), 55583109 (final) | COMPLETED; see "Deviations" for the component fix |
| pilot: 6 donors (2 per group chosen by sha256 of donor id, ME: APR015, KSA007; EAS: NA18971, HG00544; SAS: HG03831, NA21144), 12 runs each = 72 T1K runs | 55583251 (array 0-5%6; a first launch 55583088 was cancelled by me, see deviations) | 72/72 COMPLETE, 0 FAILURE |
| pilot sanity and timing | 55587138 | see below |
| pilot scoring (code check) | 55587139 | 1152 rows |
| full run: 147 donors, 10 cells each = 1470 T1K runs (60 of them already done in the pilot and reused as is) | 55587218 (array 0-146%12) | 147/147 tasks COMPLETED; 1482 COMPLETE.json (72 pilot + 1410 new), 0 FAILURE.json |
| scoring | 55601319 | 28,212 rows = 147 donors x 16 copies x 12 cell labels (n0, 10 panel cells, existing arm B) |

Counts recomputed from `results/dose_rows_full.tsv`: 28,212 rows; incomplete (missing/failed) donor-cells: 0 (`dose_avail_full.json` is empty).

## Sanity checks (step 3)

(i) n=0: rerun of T1K with the IPD reference on the same read files reproduces the existing arm A `calls_genotype.tsv` byte for byte in 6/6 pilot donors.
(ii) The existing arm B reference (`genomic_core_hprc_asian`, per fold, sha verified) rerun on the same reads reproduces the existing arm B table byte for byte in 6/6 pilot donors.
T1K is deterministic in these 12 reruns (0 genes with different calls). `results/pilot_sanity_pilot.tsv`.
Scoring self-check: n0 and existing arm B rows from my scorer reproduce YIELD.md exactly: arm A 0/162, 0/134 (T-real), 0/41 (D-real); arm B 20/162 (T-sim), 29/134, 11/41 (real).

## Timing and projection (step 4)

Pilot, 72 runs: mean 431 s per T1K run on 4 CPUs (range about 370 to 533 s; panel size has no visible effect: n=5 432 s, n=40 431 s). Projected full run: 1470 runs x 431 s = 176 task-hours of 4 CPUs (705 CPU-hours), far below the 3,000 task-hour limit, so the n=5 cell was kept. Actual: the 147-task array consumed 103.7 task-hours of 4 CPUs (1410 runs; sacct elapsed sum), plus 8.6 task-hours for the pilot. Wall time about 15 hours at %12 with other users' jobs in the queue.

## Design as implemented (and choices where the brief was silent)

- Donors and groups: 147 donors of FOLDS.json (ME 61: Arab 52 + Saudi 9; EAS 50; SAS 36). Reads: exactly the read files of the donor's existing arm A run (`-1/-2` of its COMPLETE.json; arm B uses the same files, asserted). ME = T-sim simulated reads; EAS/SAS = real reads (T-real 66 donors, D-real 20 donors). Simulated HPRC copies not used. Simulated and real never pooled.
- Components: connected components of FOLDS `groups` (union-find, via the frozen `common.exclusion_map()`); donors in no multi-member group are their own component. The test donor's component is never drawn. A component contributes only its members that belong to the source group.
- Draw: seed = sha256("donor|cell|20261114"); components shuffled, taken whole until the drawn donor count >= n (n_drawn recorded in `cells.tsv`, can exceed n); if the source runs out, all of it is used.
- Cells per donor: matched n=5,10,20,40; each other group n=10,40; non-Asian n=10,40 = 10 panel cells. n=0 reuses the existing arm A tables (verified by sha).
- Panel content: all truth `core_sequence` records of the drawn donors (8 genes, both copies), de-duplicated per gene, minus sequences identical to any IPD genomic sequence, minus sequences with N or no exon annotation. Reference = `genomic.fa` bytes unchanged + additions, alias `HLA-<gene>*PG<sha256[:12]>`, header via `common.t1k_header` (same logic as `t1k_reference.header`), as `t1k_reference.build` does for the existing arms. T1K command, container, preset and options are the ones of `run_donor.t1k` (copied into `src/s3_lib.py`), `-t 4`.
- Non-Asian source: rebuilt per donor from the `nonasian` panel sidecars of all five folds of `nonasian-panel-control-20261018` (the panels behind `genomic_core_nonasian`; core = sequence[2000:-2000], exons shifted by 2000, as `build_core.py`): 140 donors, all with 16 copies, disjoint from the 147 test donors. Each non-Asian donor is its own component (relatedness among them is unknown; not screened). The confirmatory stage's `inputs/panel_donors.txt` was listed only and not used (no per-donor sequences there; nothing under confirmatory-test was read for truth or reads).
- Scoring: frozen scorer rules imported from `mosaic-typer-20261113/src/score.py` (slot assignment by best permutation, `exact` = byte equality of the final call, `recovered` = candidate with that sequence in any called slot, `named4`, bootstrap/sign-flip function `uncertainty`). The primary tables use `recovered` (any called slot) as in YIELD Y2; the final-call version is T2b and `*_final` metrics in the TSV. Gene-body exactness = sha of seq[first exon start : last exon end] of the truth copy equals that of a novel candidate (any slot; `_final` variants use the assigned slot's final call). "Supplied by panel" = the truth copy's sequence (or gene body) is among the panel additions of that (donor, cell) for that gene. Conversion = recovered / supplied. "Fully correct": IPD-absent copy -> exact in the final call; IPD-present copy with one four-field truth label -> named 4-field; other copies (ambiguous truth label) are excluded from this share only. Clusters = FOLDS components, 4000 bootstrap resamples, seed 20261020; contrasts are paired by copy.
- Existing arm B (`B_existing`) is scored beside the cells for context only; it uses the full fold panel and no per-cell supply is computed for it.

## Deviations, problems, what was not verified

1. Component bug, caught before any result was used. My first component map let one FOLDS group overwrite another for donors listed in two overlapping groups (APR-F/APR-M/APR-S). I cancelled my own pilot array 55583088 (about 1 minute in), deleted its work folders, replaced the map by the frozen union-find exclusion map, reran the inventory (55583109) and relaunched the pilot (55583251). Nothing from the cancelled run was kept.
2. SAS source has 36 donors: all SAS-sourced n=40 panels have 36 donors (35 for SAS test donors, whose own component is excluded). 147 cells are short of n requested (all SAS n=40 cells); recorded per cell in `cells.tsv`.
3. n_drawn exceeds n by up to a few donors because whole components are drawn.
4. One random draw per cell, so panel-to-panel variability is not separated from the donor sampling; the bootstrap resamples families, not draws. Non-monotonic steps (SAS n=5 above n=10) are consistent with that.
5. Many comparisons were examined; sign-flip p-values are not corrected and not pre-registered. Treat only the large contrasts (SAS matched vs all; n=40 vs n=5) as firm.
6. T1K candidate-read and aligned-read files were deleted after every run to save space (tables, VCF, logs and `additions.fa` kept). Runs cannot be re-scored from reads, only from tables.
7. Not verified: that my panel builder yields a reference byte-identical to what `t1k_reference.build` + `build_core.py` would produce for the same donor set (it follows the same rules; the existing arm B reference itself was rerun and reproduced, but I did not rebuild it with my code and diff). T1K determinism was tested on 6 donors only (12 reruns).
8. ME results use simulated reads only; the read-type control showed simulated reads are not a safe stand-in for real reads (DPA1 especially), so ME curves must not be read as real-read curves.
9. D-real donors are the development donors of the mosaic stage; no tuning happened here (T1K fixed), and results for T-real and D-real are shown separately in T11.

## Output files (local mirror `scientist/evidence/dose-response-20261114/results/`, same on IBEX)

- `dose_rows_full.tsv`: row-level (donor, group, read type, read set, cell, source, n requested, n drawn, gene, copy, IPD-absent flag, supplied by panel (byte, body), call, exact (final), recovered_any, body_exact, body_any, named4, no_call). No sequences.
- `figure_ready_dose_response.tsv` (= `dose_summary_full.tsv`): group x cell x metric, columns k, n, pct, lo, hi.
- `dose_contrasts_full.tsv`, `dose_by_gene_full.tsv`, `truth_ceiling_summary.tsv`, `supply_rows.tsv`, `cells.tsv` (draws, n_drawn), `donors.tsv`, `pilot_sanity_pilot.tsv`, `pilot_timing_pilot.tsv`, `inventory_info.json`.
- Code: `src/` (identical to IBEX): `s3_lib.py`, `s3_inventory.py`, `s3_worker.py`, `s3_score.py`, `s3_sanity.py`, `s3_manifest.py`, `s3_report.py`, `*.sbatch`.

## Result tables (generated by src/s3_report.py from results/*.tsv)

Intervals are family-clustered 95% bootstrap intervals (4000 resamples, FOLDS components as clusters, seed 20261020, as the frozen scorer). Each shared denominator is the same copies in every cell of a group.

#### T1. IPD-absent copies recovered, byte-exact (sequence appears in any called slot, the frozen `recovered` definition)

| cell | ME (Arab+Saudi), simulated reads | East Asian, real reads | South Asian, real reads |
|---|---|---|---|
| n0 | 0/162 = 0.0% [0.0, 0.0] | 0/99 = 0.0% [0.0, 0.0] | 0/76 = 0.0% [0.0, 0.0] |
| match_n05 | 3/162 = 1.9% [0.0, 5.1] | 0/99 = 0.0% [0.0, 0.0] | 11/76 = 14.5% [6.0, 24.0] |
| match_n10 | 5/162 = 3.1% [0.5, 6.7] | 2/99 = 2.0% [0.0, 5.0] | 9/76 = 11.8% [3.9, 21.2] |
| match_n20 | 7/162 = 4.3% [1.3, 8.2] | 5/99 = 5.0% [1.0, 10.9] | 22/76 = 28.9% [16.9, 40.7] |
| match_n40 | 15/162 = 9.3% [5.1, 14.0] | 7/99 = 7.1% [2.1, 13.3] | 25/76 = 32.9% [20.6, 45.2] |
| B_existing | 20/162 = 12.3% [7.8, 16.8] | 11/99 = 11.1% [4.7, 18.4] | 29/76 = 38.2% [28.4, 48.6] |

#### T2. IPD-absent copies recovered, gene-body exact (first exon start to last exon end)

| cell | ME (Arab+Saudi), simulated reads | East Asian, real reads | South Asian, real reads |
|---|---|---|---|
| n0 | 0/162 = 0.0% [0.0, 0.0] | 0/99 = 0.0% [0.0, 0.0] | 0/76 = 0.0% [0.0, 0.0] |
| match_n05 | 3/162 = 1.9% [0.0, 5.1] | 0/99 = 0.0% [0.0, 0.0] | 11/76 = 14.5% [6.0, 24.0] |
| match_n10 | 5/162 = 3.1% [0.5, 6.7] | 2/99 = 2.0% [0.0, 5.0] | 9/76 = 11.8% [3.9, 21.2] |
| match_n20 | 8/162 = 4.9% [1.4, 8.7] | 5/99 = 5.0% [1.0, 10.9] | 22/76 = 28.9% [16.9, 40.7] |
| match_n40 | 16/162 = 9.9% [5.3, 14.6] | 7/99 = 7.1% [2.1, 13.3] | 25/76 = 32.9% [20.6, 45.2] |
| B_existing | 39/162 = 24.1% [17.8, 30.9] | 13/99 = 13.1% [6.5, 20.6] | 29/76 = 38.2% [28.4, 48.6] |

#### T2b. IPD-absent copies recovered, byte-exact AND it is the slot's final call

| cell | ME (Arab+Saudi), simulated reads | East Asian, real reads | South Asian, real reads |
|---|---|---|---|
| n0 | 0/162 = 0.0% [0.0, 0.0] | 0/99 = 0.0% [0.0, 0.0] | 0/76 = 0.0% [0.0, 0.0] |
| match_n05 | 2/162 = 1.2% [0.0, 3.3] | 0/99 = 0.0% [0.0, 0.0] | 11/76 = 14.5% [6.0, 24.0] |
| match_n10 | 4/162 = 2.5% [0.5, 5.1] | 2/99 = 2.0% [0.0, 5.0] | 9/76 = 11.8% [3.9, 21.2] |
| match_n20 | 6/162 = 3.7% [1.3, 6.8] | 5/99 = 5.0% [1.0, 10.9] | 22/76 = 28.9% [16.9, 40.7] |
| match_n40 | 13/162 = 8.0% [4.5, 11.7] | 7/99 = 7.1% [2.1, 13.3] | 25/76 = 32.9% [20.6, 45.2] |
| B_existing | 19/162 = 11.7% [7.4, 15.7] | 11/99 = 11.1% [4.7, 18.4] | 28/76 = 36.8% [27.0, 47.5] |

#### T3. Matched vs cross-population panel at the same n, byte-exact recovery of IPD-absent copies

| group | n | matched | from ME | from EAS | from SAS | from non-Asian HPRC |
|---|---|---|---|---|---|---|
| ME | 10 | 5/162 = 3.1% [0.5, 6.7] | - | 1/162 = 0.6% [0.0, 1.9] | 2/162 = 1.2% [0.0, 3.3] | 1/162 = 0.6% [0.0, 2.1] |
| ME | 40 | 15/162 = 9.3% [5.1, 14.0] | - | 3/162 = 1.9% [0.0, 5.0] | 8/162 = 4.9% [2.0, 7.8] | 6/162 = 3.7% [1.3, 6.3] |
| EAS | 10 | 2/99 = 2.0% [0.0, 5.0] | 0/99 = 0.0% [0.0, 0.0] | - | 3/99 = 3.0% [0.0, 6.7] | 0/99 = 0.0% [0.0, 0.0] |
| EAS | 40 | 7/99 = 7.1% [2.1, 13.3] | 0/99 = 0.0% [0.0, 0.0] | - | 3/99 = 3.0% [0.0, 6.7] | 0/99 = 0.0% [0.0, 0.0] |
| SAS | 10 | 9/76 = 11.8% [3.9, 21.2] | 4/76 = 5.3% [0.0, 11.2] | 0/76 = 0.0% [0.0, 0.0] | - | 0/76 = 0.0% [0.0, 0.0] |
| SAS | 40 | 25/76 = 32.9% [20.6, 45.2] | 11/76 = 14.5% [7.5, 21.1] | 8/76 = 10.5% [3.8, 18.7] | - | 2/76 = 2.6% [0.0, 6.4] |

#### T3b. Same, gene-body exact

| group | n | matched | from ME | from EAS | from SAS | from non-Asian HPRC |
|---|---|---|---|---|---|---|
| ME | 10 | 5/162 = 3.1% [0.5, 6.7] | - | 7/162 = 4.3% [1.3, 8.2] | 15/162 = 9.3% [5.1, 14.6] | 1/162 = 0.6% [0.0, 2.1] |
| ME | 40 | 16/162 = 9.9% [5.3, 14.6] | - | 18/162 = 11.1% [6.3, 17.2] | 27/162 = 16.7% [10.9, 23.2] | 11/162 = 6.8% [3.3, 11.1] |
| EAS | 10 | 2/99 = 2.0% [0.0, 5.0] | 0/99 = 0.0% [0.0, 0.0] | - | 3/99 = 3.0% [0.0, 6.7] | 0/99 = 0.0% [0.0, 0.0] |
| EAS | 40 | 7/99 = 7.1% [2.1, 13.3] | 2/99 = 2.0% [0.0, 5.0] | - | 3/99 = 3.0% [0.0, 6.7] | 0/99 = 0.0% [0.0, 0.0] |
| SAS | 10 | 9/76 = 11.8% [3.9, 21.2] | 4/76 = 5.3% [0.0, 11.2] | 0/76 = 0.0% [0.0, 0.0] | - | 0/76 = 0.0% [0.0, 0.0] |
| SAS | 40 | 25/76 = 32.9% [20.6, 45.2] | 11/76 = 14.5% [7.5, 21.1] | 8/76 = 10.5% [3.8, 18.7] | - | 2/76 = 2.6% [0.0, 6.4] |

#### T4. Ceiling: IPD-absent copies whose exact sequence is in the drawn panel (truth only, no typing)

| cell | ME (Arab+Saudi), simulated reads | East Asian, real reads | South Asian, real reads |
|---|---|---|---|
| n0 | 0/162 = 0.0% [0.0, 0.0] | 0/99 = 0.0% [0.0, 0.0] | 0/76 = 0.0% [0.0, 0.0] |
| match_n05 | 17/162 = 10.5% [5.7, 16.6] | 10/99 = 10.1% [4.8, 16.7] | 16/76 = 21.1% [10.6, 31.7] |
| match_n10 | 21/162 = 13.0% [7.9, 19.3] | 13/99 = 13.1% [7.3, 20.0] | 15/76 = 19.7% [10.3, 30.6] |
| match_n20 | 24/162 = 14.8% [9.6, 21.3] | 18/99 = 18.2% [10.8, 26.8] | 32/76 = 42.1% [30.3, 52.9] |
| match_n40 | 38/162 = 23.5% [17.3, 30.9] | 21/99 = 21.2% [13.1, 31.2] | 35/76 = 46.0% [33.8, 57.5] |
| B_existing | n/a | n/a | n/a |

#### T4b. Ceiling, gene-body

| cell | ME (Arab+Saudi), simulated reads | East Asian, real reads | South Asian, real reads |
|---|---|---|---|
| n0 | 0/162 = 0.0% [0.0, 0.0] | 0/99 = 0.0% [0.0, 0.0] | 0/76 = 0.0% [0.0, 0.0] |
| match_n05 | 18/162 = 11.1% [6.2, 17.3] | 10/99 = 10.1% [4.8, 16.7] | 17/76 = 22.4% [11.0, 34.2] |
| match_n10 | 22/162 = 13.6% [8.5, 20.1] | 14/99 = 14.1% [8.2, 21.0] | 17/76 = 22.4% [12.2, 33.8] |
| match_n20 | 29/162 = 17.9% [12.6, 24.1] | 21/99 = 21.2% [14.0, 29.6] | 35/76 = 46.0% [33.8, 57.3] |
| match_n40 | 43/162 = 26.5% [20.3, 33.9] | 25/99 = 25.2% [16.8, 34.8] | 39/76 = 51.3% [38.5, 62.9] |
| B_existing | n/a | n/a | n/a |

#### T5. Conversion = recovered / supplied by the panel (byte-exact; denominator = supplied copies)

| cell | ME (Arab+Saudi), simulated reads | East Asian, real reads | South Asian, real reads |
|---|---|---|---|
| n0 | n/a | n/a | n/a |
| match_n05 | 3/17 = 17.6% [0.0, 42.1] | 0/10 = 0.0% [0.0, 0.0] | 11/16 = 68.8% [41.2, 93.3] |
| match_n10 | 5/21 = 23.8% [4.8, 45.5] | 2/13 = 15.4% [0.0, 38.5] | 9/15 = 60.0% [30.8, 86.7] |
| match_n20 | 7/24 = 29.2% [9.5, 48.3] | 5/18 = 27.8% [5.9, 52.4] | 22/32 = 68.8% [48.3, 85.3] |
| match_n40 | 15/38 = 39.5% [23.1, 55.6] | 7/21 = 33.3% [14.3, 52.4] | 25/35 = 71.4% [53.1, 86.1] |
| B_existing | n/a | n/a | n/a |

#### T5b. Conversion, gene-body

| cell | ME (Arab+Saudi), simulated reads | East Asian, real reads | South Asian, real reads |
|---|---|---|---|
| n0 | n/a | n/a | n/a |
| match_n05 | 3/18 = 16.7% [0.0, 40.0] | 0/10 = 0.0% [0.0, 0.0] | 11/17 = 64.7% [37.5, 92.3] |
| match_n10 | 5/22 = 22.7% [4.5, 44.0] | 2/14 = 14.3% [0.0, 35.7] | 9/17 = 52.9% [25.0, 81.2] |
| match_n20 | 8/29 = 27.6% [8.0, 46.4] | 5/21 = 23.8% [5.0, 45.5] | 22/35 = 62.9% [42.4, 81.2] |
| match_n40 | 16/43 = 37.2% [21.6, 53.1] | 7/25 = 28.0% [9.1, 44.4] | 25/39 = 64.1% [46.0, 81.8] |
| B_existing | n/a | n/a | n/a |

#### T6. Four-field named correct among IPD-present copies (safety)

| cell | ME (Arab+Saudi), simulated reads | East Asian, real reads | South Asian, real reads |
|---|---|---|---|
| n0 | 707/804 = 87.9% [85.6, 89.9] | 643/694 = 92.7% [90.7, 94.5] | 471/496 = 95.0% [92.8, 96.8] |
| match_n05 | 705/804 = 87.7% [85.4, 89.6] | 643/694 = 92.7% [90.7, 94.5] | 472/496 = 95.2% [93.2, 97.0] |
| match_n10 | 706/804 = 87.8% [85.7, 89.7] | 642/694 = 92.5% [90.5, 94.3] | 472/496 = 95.2% [93.2, 97.0] |
| match_n20 | 705/804 = 87.7% [85.3, 89.6] | 643/694 = 92.7% [90.7, 94.5] | 470/496 = 94.8% [92.8, 96.5] |
| match_n40 | 699/804 = 86.9% [84.8, 88.9] | 642/694 = 92.5% [90.6, 94.3] | 469/496 = 94.6% [92.5, 96.4] |
| B_existing | 694/804 = 86.3% [83.8, 88.4] | 631/694 = 90.9% [89.0, 92.7] | 463/496 = 93.3% [91.0, 95.3] |

#### T7. Share of ALL scored copies fully correct (IPD-absent: exact in final call; IPD-present with one four-field truth label: named 4-field)

| cell | ME (Arab+Saudi), simulated reads | East Asian, real reads | South Asian, real reads |
|---|---|---|---|
| n0 | 707/966 = 73.2% [70.6, 75.4] | 643/793 = 81.1% [79.0, 83.1] | 471/572 = 82.3% [79.4, 85.1] |
| match_n05 | 707/966 = 73.2% [70.6, 75.4] | 643/793 = 81.1% [79.0, 83.1] | 483/572 = 84.4% [81.5, 87.1] |
| match_n10 | 710/966 = 73.5% [71.1, 75.6] | 644/793 = 81.2% [79.2, 83.1] | 481/572 = 84.1% [81.0, 86.9] |
| match_n20 | 711/966 = 73.6% [71.0, 75.8] | 648/793 = 81.7% [79.7, 83.7] | 492/572 = 86.0% [83.2, 88.6] |
| match_n40 | 712/966 = 73.7% [71.1, 75.8] | 649/793 = 81.8% [79.8, 83.8] | 494/572 = 86.4% [83.6, 88.9] |
| B_existing | 713/966 = 73.8% [70.9, 76.2] | 642/793 = 81.0% [78.7, 83.1] | 491/572 = 85.8% [83.2, 88.4] |

Cross-population cells (same table, other sources) are in `results/dose_summary_full.tsv` / `figure_ready_dose_response.tsv`.

### T8. Paired contrasts (same copies; percentage points of IPD-absent copies; sign-flip p from the frozen scorer)

| group | kind | A vs B | metric | gains/losses | diff pp [95%] | p |
|---|---|---|---|---|---|---|
| EAS | matched vs IPD only | match_n05 vs n0 | absent_exact_any | 0/0 | +0.00 [+0.00, +0.00] | 1.0 |
| EAS | matched vs IPD only | match_n10 vs n0 | absent_exact_any | 2/0 | +2.02 [+0.00, +4.95] | 0.5 |
| EAS | matched vs IPD only | match_n20 vs n0 | absent_exact_any | 5/0 | +5.05 [+0.96, +10.89] | 0.125 |
| EAS | matched vs IPD only | match_n40 vs n0 | absent_exact_any | 7/0 | +7.07 [+2.06, +13.33] | 0.0312 |
| EAS | trend within matched | match_n10 vs match_n05 | absent_exact_any | 2/0 | +2.02 [+0.00, +4.95] | 0.5 |
| EAS | trend within matched | match_n10 vs match_n05 | absent_body_any | 2/0 | +2.02 [+0.00, +4.95] | 0.5 |
| EAS | trend within matched | match_n20 vs match_n10 | absent_exact_any | 4/1 | +3.03 [-1.82, +8.91] | 0.5 |
| EAS | trend within matched | match_n20 vs match_n10 | absent_body_any | 4/1 | +3.03 [-1.82, +8.91] | 0.5 |
| EAS | trend within matched | match_n40 vs match_n20 | absent_exact_any | 3/1 | +2.02 [-1.89, +6.02] | 0.625 |
| EAS | trend within matched | match_n40 vs match_n20 | absent_body_any | 3/1 | +2.02 [-1.89, +6.02] | 0.625 |
| EAS | trend within matched | match_n40 vs match_n05 | absent_exact_any | 7/0 | +7.07 [+2.06, +13.33] | 0.0312 |
| EAS | trend within matched | match_n40 vs match_n05 | absent_body_any | 7/0 | +7.07 [+2.06, +13.33] | 0.0312 |
| EAS | trend within matched | match_n40 vs match_n10 | absent_exact_any | 5/0 | +5.05 [+0.95, +10.87] | 0.125 |
| EAS | trend within matched | match_n40 vs match_n10 | absent_body_any | 5/0 | +5.05 [+0.95, +10.87] | 0.125 |
| EAS | matched vs cross-population, same n | match_n10 vs x_ME_n10 | absent_exact_any | 2/0 | +2.02 [+0.00, +4.95] | 0.5 |
| EAS | matched vs cross-population, same n | match_n10 vs x_ME_n10 | absent_body_any | 2/0 | +2.02 [+0.00, +4.95] | 0.5 |
| EAS | matched vs cross-population, same n | match_n40 vs x_ME_n40 | absent_exact_any | 7/0 | +7.07 [+2.06, +13.33] | 0.0312 |
| EAS | matched vs cross-population, same n | match_n40 vs x_ME_n40 | absent_body_any | 7/2 | +5.05 [-1.00, +12.00] | 0.2344 |
| EAS | matched vs cross-population, same n | match_n10 vs x_NONASIAN_n10 | absent_exact_any | 2/0 | +2.02 [+0.00, +4.95] | 0.5 |
| EAS | matched vs cross-population, same n | match_n10 vs x_NONASIAN_n10 | absent_body_any | 2/0 | +2.02 [+0.00, +4.95] | 0.5 |
| EAS | matched vs cross-population, same n | match_n40 vs x_NONASIAN_n40 | absent_exact_any | 7/0 | +7.07 [+2.06, +13.33] | 0.0312 |
| EAS | matched vs cross-population, same n | match_n40 vs x_NONASIAN_n40 | absent_body_any | 7/0 | +7.07 [+2.06, +13.33] | 0.0312 |
| EAS | matched vs cross-population, same n | match_n10 vs x_SAS_n10 | absent_exact_any | 2/3 | -1.01 [-4.44, +2.13] | 1.0 |
| EAS | matched vs cross-population, same n | match_n10 vs x_SAS_n10 | absent_body_any | 2/3 | -1.01 [-4.44, +2.13] | 1.0 |
| EAS | matched vs cross-population, same n | match_n40 vs x_SAS_n40 | absent_exact_any | 7/3 | +4.04 [-1.90, +10.68] | 0.3594 |
| EAS | matched vs cross-population, same n | match_n40 vs x_SAS_n40 | absent_body_any | 7/3 | +4.04 [-1.90, +10.68] | 0.3594 |
| ME | matched vs IPD only | match_n05 vs n0 | absent_exact_any | 3/0 | +1.85 [+0.00, +5.07] | 0.5 |
| ME | matched vs IPD only | match_n10 vs n0 | absent_exact_any | 5/0 | +3.09 [+0.53, +6.74] | 0.125 |
| ME | matched vs IPD only | match_n20 vs n0 | absent_exact_any | 7/0 | +4.32 [+1.33, +8.23] | 0.0312 |
| ME | matched vs IPD only | match_n40 vs n0 | absent_exact_any | 15/0 | +9.26 [+5.07, +13.97] | 0.0005 |
| ME | trend within matched | match_n10 vs match_n05 | absent_exact_any | 2/0 | +1.23 [+0.00, +3.14] | 0.5 |
| ME | trend within matched | match_n10 vs match_n05 | absent_body_any | 2/0 | +1.23 [+0.00, +3.14] | 0.5 |
| ME | trend within matched | match_n20 vs match_n10 | absent_exact_any | 4/2 | +1.23 [-1.80, +4.14] | 0.6875 |
| ME | trend within matched | match_n20 vs match_n10 | absent_body_any | 5/2 | +1.85 [-1.49, +5.10] | 0.5312 |
| ME | trend within matched | match_n40 vs match_n20 | absent_exact_any | 8/0 | +4.94 [+1.90, +8.59] | 0.0156 |
| ME | trend within matched | match_n40 vs match_n20 | absent_body_any | 8/0 | +4.94 [+1.90, +8.59] | 0.0156 |
| ME | trend within matched | match_n40 vs match_n05 | absent_exact_any | 12/0 | +7.41 [+3.61, +11.49] | 0.002 |
| ME | trend within matched | match_n40 vs match_n05 | absent_body_any | 13/0 | +8.02 [+3.77, +12.14] | 0.002 |
| ME | trend within matched | match_n40 vs match_n10 | absent_exact_any | 11/1 | +6.17 [+2.21, +10.39] | 0.0176 |
| ME | trend within matched | match_n40 vs match_n10 | absent_body_any | 12/1 | +6.79 [+2.34, +11.11] | 0.0176 |
| ME | matched vs cross-population, same n | match_n10 vs x_EAS_n10 | absent_exact_any | 5/1 | +2.47 [-0.59, +6.21] | 0.3125 |
| ME | matched vs cross-population, same n | match_n10 vs x_EAS_n10 | absent_body_any | 3/5 | -1.23 [-4.83, +2.07] | 0.7266 |
| ME | matched vs cross-population, same n | match_n40 vs x_EAS_n40 | absent_exact_any | 12/0 | +7.41 [+3.55, +11.70] | 0.002 |
| ME | matched vs cross-population, same n | match_n40 vs x_EAS_n40 | absent_body_any | 11/13 | -1.23 [-8.28, +4.66] | 0.8441 |
| ME | matched vs cross-population, same n | match_n10 vs x_NONASIAN_n10 | absent_exact_any | 5/1 | +2.47 [-0.62, +6.33] | 0.3125 |
| ME | matched vs cross-population, same n | match_n10 vs x_NONASIAN_n10 | absent_body_any | 5/1 | +2.47 [-0.62, +6.33] | 0.3125 |
| ME | matched vs cross-population, same n | match_n40 vs x_NONASIAN_n40 | absent_exact_any | 15/6 | +5.56 [+0.60, +11.19] | 0.078 |
| ME | matched vs cross-population, same n | match_n40 vs x_NONASIAN_n40 | absent_body_any | 16/11 | +3.09 [-3.85, +9.43] | 0.4821 |
| ME | matched vs cross-population, same n | match_n10 vs x_SAS_n10 | absent_exact_any | 5/2 | +1.85 [-1.44, +5.81] | 0.5312 |
| ME | matched vs cross-population, same n | match_n10 vs x_SAS_n10 | absent_body_any | 5/15 | -6.17 [-12.32, -0.60] | 0.0633 |
| ME | matched vs cross-population, same n | match_n40 vs x_SAS_n40 | absent_exact_any | 14/7 | +4.32 [-0.70, +10.37] | 0.1923 |
| ME | matched vs cross-population, same n | match_n40 vs x_SAS_n40 | absent_body_any | 13/24 | -6.79 [-14.86, +0.00] | 0.091 |
| SAS | matched vs IPD only | match_n05 vs n0 | absent_exact_any | 11/0 | +14.47 [+5.95, +24.00] | 0.0078 |
| SAS | matched vs IPD only | match_n10 vs n0 | absent_exact_any | 9/0 | +11.84 [+3.90, +21.21] | 0.0156 |
| SAS | matched vs IPD only | match_n20 vs n0 | absent_exact_any | 22/0 | +28.95 [+16.92, +40.74] | 0.0001 |
| SAS | matched vs IPD only | match_n40 vs n0 | absent_exact_any | 25/0 | +32.89 [+20.59, +45.21] | 0.0 |
| SAS | trend within matched | match_n10 vs match_n05 | absent_exact_any | 4/6 | -2.63 [-12.50, +7.69] | 0.8125 |
| SAS | trend within matched | match_n10 vs match_n05 | absent_body_any | 4/6 | -2.63 [-12.50, +7.69] | 0.8125 |
| SAS | trend within matched | match_n20 vs match_n10 | absent_exact_any | 15/2 | +17.11 [+4.05, +30.23] | 0.043 |
| SAS | trend within matched | match_n20 vs match_n10 | absent_body_any | 15/2 | +17.11 [+4.05, +30.23] | 0.043 |
| SAS | trend within matched | match_n40 vs match_n20 | absent_exact_any | 3/0 | +3.95 [+0.00, +10.67] | 0.5 |
| SAS | trend within matched | match_n40 vs match_n20 | absent_body_any | 3/0 | +3.95 [+0.00, +10.67] | 0.5 |
| SAS | trend within matched | match_n40 vs match_n05 | absent_exact_any | 14/0 | +18.42 [+7.94, +29.87] | 0.0039 |
| SAS | trend within matched | match_n40 vs match_n05 | absent_body_any | 14/0 | +18.42 [+7.94, +29.87] | 0.0039 |
| SAS | trend within matched | match_n40 vs match_n10 | absent_exact_any | 16/0 | +21.05 [+9.52, +32.89] | 0.0039 |
| SAS | trend within matched | match_n40 vs match_n10 | absent_body_any | 16/0 | +21.05 [+9.52, +32.89] | 0.0039 |
| SAS | matched vs cross-population, same n | match_n10 vs x_EAS_n10 | absent_exact_any | 9/0 | +11.84 [+3.90, +21.21] | 0.0156 |
| SAS | matched vs cross-population, same n | match_n10 vs x_EAS_n10 | absent_body_any | 9/0 | +11.84 [+3.90, +21.21] | 0.0156 |
| SAS | matched vs cross-population, same n | match_n40 vs x_EAS_n40 | absent_exact_any | 18/1 | +22.37 [+9.52, +35.14] | 0.0059 |
| SAS | matched vs cross-population, same n | match_n40 vs x_EAS_n40 | absent_body_any | 18/1 | +22.37 [+9.52, +35.14] | 0.0059 |
| SAS | matched vs cross-population, same n | match_n10 vs x_ME_n10 | absent_exact_any | 8/3 | +6.58 [-3.49, +17.81] | 0.3438 |
| SAS | matched vs cross-population, same n | match_n10 vs x_ME_n10 | absent_body_any | 8/3 | +6.58 [-3.49, +17.81] | 0.3438 |
| SAS | matched vs cross-population, same n | match_n40 vs x_ME_n40 | absent_exact_any | 20/6 | +18.42 [+5.41, +31.58] | 0.0179 |
| SAS | matched vs cross-population, same n | match_n40 vs x_ME_n40 | absent_body_any | 20/6 | +18.42 [+5.41, +31.58] | 0.0179 |
| SAS | matched vs cross-population, same n | match_n10 vs x_NONASIAN_n10 | absent_exact_any | 9/0 | +11.84 [+3.90, +21.21] | 0.0156 |
| SAS | matched vs cross-population, same n | match_n10 vs x_NONASIAN_n10 | absent_body_any | 9/0 | +11.84 [+3.90, +21.21] | 0.0156 |
| SAS | matched vs cross-population, same n | match_n40 vs x_NONASIAN_n40 | absent_exact_any | 24/1 | +30.26 [+17.57, +43.66] | 0.0001 |
| SAS | matched vs cross-population, same n | match_n40 vs x_NONASIAN_n40 | absent_body_any | 24/1 | +30.26 [+17.57, +43.66] | 0.0001 |

### T9. Safety contrasts: IPD-present four-field named correct, matched cell minus IPD only (pp)

| group | cell | gains/losses | diff pp [95%] |
|---|---|---|---|
| EAS | match_n05 | 0/0 | +0.00 [+0.00, +0.00] |
| EAS | match_n10 | 0/1 | -0.14 [-0.44, +0.00] |
| EAS | match_n20 | 0/0 | +0.00 [+0.00, +0.00] |
| EAS | match_n40 | 0/1 | -0.14 [-0.44, +0.00] |
| ME | match_n05 | 2/4 | -0.25 [-0.90, +0.32] |
| ME | match_n10 | 3/4 | -0.12 [-0.62, +0.42] |
| ME | match_n20 | 4/6 | -0.25 [-1.15, +0.60] |
| ME | match_n40 | 2/10 | -1.00 [-1.76, -0.14] |
| SAS | match_n05 | 2/1 | +0.20 [+0.00, +0.61] |
| SAS | match_n10 | 2/1 | +0.20 [-0.41, +0.82] |
| SAS | match_n20 | 3/4 | -0.20 [-1.03, +0.61] |
| SAS | match_n40 | 3/5 | -0.40 [-1.22, +0.40] |

### T10. By gene: IPD-absent copies recovered (byte-exact), matched n=40 | supplied by that panel | IPD only (0 by construction), per group

| group | A | B | C | DPA1 | DPB1 | DQA1 | DQB1 | DRB1 |
|---|---|---|---|---|---|---|---|---|
| ME | 0 rec, 0 sup, of 2 | 0 rec, 0 sup, of 12 | 0 rec, 1 sup, of 6 | 0 rec, 4 sup, of 16 | 0 rec, 0 sup, of 15 | 2 rec, 19 sup, of 32 | 0 rec, 0 sup, of 3 | 13 rec, 14 sup, of 76 |
| EAS | 0 rec, 0 sup, of 2 | 0 rec, 0 sup, of 4 | 0 rec, 0 sup, of 2 | 0 rec, 0 sup, of 6 | 1 rec, 1 sup, of 10 | 0 rec, 0 sup, of 3 | 0 rec, 0 sup, of 1 | 6 rec, 20 sup, of 71 |
| SAS | - | 0 rec, 0 sup, of 4 | 0 rec, 0 sup, of 1 | 2 rec, 3 sup, of 9 | 0 rec, 0 sup, of 5 | 14 rec, 14 sup, of 18 | - | 9 rec, 18 sup, of 39 |

### T11. Real reads split by read set (matched cells, byte-exact recovery of IPD-absent copies)

| group | read set | donors | n0 | n=5 | n=10 | n=20 | n=40 |
|---|---|---|---|---|---|---|---|
| EAS | T-real | 33 | 0/76 | 0/76 | 2/76 | 5/76 | 5/76 |
| EAS | D-real | 13 | 0/23 | 0/23 | 0/23 | 0/23 | 2/23 |
| SAS | T-real | 27 | 0/58 | 7/58 | 6/58 | 15/58 | 16/58 |
| SAS | D-real | 6 | 0/18 | 4/18 | 3/18 | 7/18 | 9/18 |
