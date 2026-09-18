# Frozen validation: 84 previously unexamined donors

The frozen Asian-graph candidate improves four-field genotype accuracy over original T1K, but the improvement is explained by the genomic-IPD reference. This validation does **not** establish an incremental Asian-graph benefit, and the complete predeclared acceptance gate did **not** pass.

All 84 donors and all four methods completed successfully (336 results). No failed donor was removed. Inference, endpoint code and the cohort were frozen before outcome inspection. The graph exclusion audit removed 400 haplotypes covering reserved donor/family groups and found zero overlap among 3,515 retained source occurrences.

| Method | Two-field correct / eligible | Four-field correct / eligible |
|---|---:|---:|
| Original T1K | 646/670 (96.42%) | 199/387 (51.42%) |
| T1K with full genomic IPD | 645/670 (96.27%) | 335/387 (86.56%) |
| Genomic IPD + HPRC graph | 645/670 (96.27%) | 336/387 (86.82%) |
| Genomic IPD + HPRC+Asian graph | 645/670 (96.27%) | 333/387 (86.05%) |

The two-field totals are effectively unchanged: 646 versus 645 correct. This is not evidence of a meaningful decline. However, there were 20 gains and 21 losses across individual genotypes, rather than only one changed genotype. The paired 95% interval is −2.09 to +1.79 percentage points, so the predeclared one-percentage-point non-inferiority requirement is not established.

The primary four-field difference is **+34.63 percentage points** (95% family-bootstrap interval **+29.07 to +40.31**), with 152 gains and 18 losses. This establishes higher four-field accuracy on the eligible sequence-derived truth subset, not superiority at every resolution or a graph-specific gain.

Compared with genomic IPD alone, the Asian-graph candidate has one gain and three losses: **−0.52 percentage points**, 95% interval **−1.57 to +0.50**. Compared with the HPRC-only graph, it has one gain and four losses: **−0.78 percentage points**, interval **−2.02 to +0.26**. Neither contrast establishes a benefit or a statistically resolved deterioration.

## Scope and limitations

- Truth comes from the existing sequence catalogue. Four-field eligibility requires genomic-resolution labels for both haplotypes; 387 of 672 donor–gene pairs qualify. The remaining 285 are reported as ineligible, not as correct or incorrect. Two-field eligibility is 670/672.
- This is independent evaluation of a candidate frozen before this cohort was inspected. The 84 donors are now exposed; tuning on them and scoring them again would be post hoc, not a new independent validation.
- The cohort comprises 37 AFR, 16 AMR, 13 EAS, 2 EUR and 16 SAS donors. Subpopulation estimates, especially EUR, have limited precision. The cohort is not a matched Asian/non-Asian population comparison.
- Confidence intervals use 10,000 family-cluster bootstrap replicates, seed 20260918. Original T1K is the prespecified primary comparator; genomic-IPD and graph-panel contrasts are exploratory.
- The literal acceptance checks are retained: four-field lower bound positive = true; two-field point not decreased = false; two-field interval excludes a loss greater than one point = false. `protocol_pass` remains false. No threshold was changed after seeing these results.

## Graph changes requiring diagnosis

| Donor | Population | Gene | Versus genomic IPD |
|---|---|---|---|
| HG02080 | EAS | C | loss |
| HG02257 | AFR | A | loss |
| HG02258 | AFR | DQA1 | rescue |
| HG03704 | SAS | DQA1 | loss |

These four cases are diagnostic observations, not a rule for selecting favourable outcomes. Next work should address the fragment likelihood and coarse-genotype preservation using development data, then use a newly frozen independent evaluation. Increasing the score threshold on these observed losses and reusing this cohort would not establish improvement.

## Reproducibility

`../UNBLINDED.json` records first evaluation and frozen-input hashes. `../execution-handoff.json` records all 336 verified runs and output hashes. `report.json`, `summary.tsv` and `gene_scores.tsv` preserve the frozen endpoint output; `graph-change-diagnostics.json` records the four accuracy-changing graph decisions. The complete small-output snapshot remains at `t1k-pangenome/work/validation-snapshot/`, and the source runs remain on IBEX.

```bash
python3 t1k-pangenome/score_reserved_snapshot.py \
  --snapshot t1k-pangenome/work/validation-snapshot \
  --output /path/to/new-output-directory
```
