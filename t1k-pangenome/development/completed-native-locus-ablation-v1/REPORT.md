# Native-locus guarded graph ablation

All 184 runs passed: eight development donors and 84 already-exposed validation
donors, each with HPRC-only and HPRC+Asian graph panels. This is a post-unblinding
exploratory result, not a second independent validation.

| Cohort | Original T1K, four-field | Genomic IPD | Guarded HPRC | Guarded HPRC+Asian |
|---|---:|---:|---:|---:|
| Development | 21/47 | 41/47 | 41/47 | 41/47 |
| Exposed validation | 199/387 | 335/387 | 336/387 | 339/387 |

The Asian graph adds four correct genotypes relative to genomic IPD on the
exposed cohort, with no new eligible genotype errors. All four are DQA1:
HG02258, HG02451, HG03579 and NA19043. The HLA-A error in HG02257 is removed.
The native-locus guard also removes the remaining development rescue seen in
the alignment-only ablation (42/47 becomes 41/47). Both effects are retained in
the report; the model is not selected solely on the favourable cohort.

Two-field outputs remain identical to genomic IPD: 62/63 on development and
645/670 on exposed validation. They do not yet preserve original T1K's calls,
which score 63/63 and 646/670 respectively. No original acceptance gate is
claimed to pass for this post-hoc candidate.

The mechanism is documented in
`../../validation/postvalidation-read-audit-v1/REPORT.md`: the erroneous HLA-A
support came from fragments with positive native assignments exclusively to
HLA-E, a locus missing from the eight-locus graph comparison. The candidate
rejects graph-locus evidence contradicted by exclusively other native loci.
It retains mixed native assignments and graph-only fragments; it is not a
calibrated joint likelihood model.

The scorer verified complete counts, output and code hashes, baseline linkage,
and unchanged two-field calls for all 184 outputs. Per-gene rows, summary counts,
and all completion-manifest hashes are included alongside this report. The
small-output snapshot remains in `work/native-locus-ablation-v1`; source runs
remain on IBEX. Five focused tests cover length normalization, locus ambiguity,
exclusive paralog evidence, graph-only evidence, and mixed native assignments.

A metadata-only HGSVC audit found 28 candidate donors outside indexed prior
cohorts/graph donors and their pedigree families, with indexed short reads.
These are **not yet a frozen validation cohort**. Additional exposure/alias,
assembly-quality and read-availability checks are required. Candidate selection
and endpoint freezing must precede any evaluation on them.
