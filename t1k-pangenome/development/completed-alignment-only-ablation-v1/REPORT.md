# Alignment-only ablation on exposed donors

All 184 runs completed: eight development donors and 84 previously evaluated
donors, each with both graph panels. This is an exploratory, post-unblinding
ablation. The original frozen validation and its failed combined acceptance gate
remain unchanged.

| Cohort | Original T1K | Genomic IPD | Alignment-only HPRC | Alignment-only HPRC+Asian |
|---|---:|---:|---:|---:|
| Development, four-field | 21/47 | 41/47 | 41/47 | 42/47 |
| Exposed validation, four-field | 199/387 | 335/387 | 335/387 | 335/387 |

The Asian candidate improves from the frozen 333/387 to 335/387 by preventing
the HG02080 HLA-C and HG03704 DQA1 losses diagnosed as length-bias reversals.
The HG02257 HLA-A loss and HG02258 DQA1 rescue remain, balancing to zero net
gain over genomic IPD. On development data, removing length-derived support
also loses one of the frozen candidate's two genotype rescues: 43/47 becomes
42/47. The HPRC-only frozen 336/387 becomes 335/387.

Two-field calls remain identical to genomic IPD by construction and verified
output comparison: 62/63 on development and 645/670 on exposed validation,
versus original T1K's 63/63 and 646/670. This ablation does not address differences
introduced by changing the original T1K reference.

The comparison verifies the normalization mechanism, not an independent graph
benefit. The remaining HLA-A loss requires read/variant-level investigation;
simply tuning the score threshold against this exposed cohort would not provide
new validation. A future candidate must preserve the intended coarse-genotype
contract and be evaluated on newly reserved donors.

`gene_scores.tsv` contains every donor/gene/method row. `provenance.json` records
all 184 completion-manifest hashes and the plan/scorer hashes. The scorer checks
reference-native provenance, output hashes, model hashes, complete donor counts,
and unchanged two-field calls. The small-output snapshot is retained under
`t1k-pangenome/work/alignment-only-ablation-v1/`; source outputs remain on IBEX.
