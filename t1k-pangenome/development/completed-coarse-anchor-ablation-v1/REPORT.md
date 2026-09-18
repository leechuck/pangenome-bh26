# Coarse-call preserving composition on exposed donors

This optional candidate keeps original T1K two-field predictions exactly. It
accepts a resolved four-field pair from genomic IPD or a graph refinement only
when that pair projects to the same resolved original two-field pair. Otherwise
it retains the complete original call. The decision rule does not read truth.

| Cohort | Original T1K, 4 fields | Anchored genomic IPD | Anchored HPRC | Anchored HPRC+Asian |
|---|---:|---:|---:|---:|
| Development | 21/47 | 41/47 | 41/47 | 41/47 |
| Exposed validation | 199/387 | 335/387 | 336/387 | 339/387 |

All anchored methods retain original T1K's two-field performance: 63/63 on
development and 646/670 on exposed validation. The scorer additionally verifies
per-gene prediction identity, not merely equal totals. Four-field genotype totals
are unchanged from the unanchored candidates, but individual calls and allele
matches differ; the composition is not being described as identical inference.

Original ambiguity and missing-call handling are preserved. Three focused tests
cover allele-slot reordering, conflicting or ambiguous coarse calls, and absent
candidate genes. This is post-unblinding development, not independent proof.
The original frozen validation and its failed combined acceptance gate remain
unchanged. All upstream variants remain available for comparison.

`provenance.json` binds original T1K tables, genomic and graph completion markers,
and composition/scorer code. `decisions.json` records every fallback or accepted
refinement. The prospective 28-donor HGSVC reservation is separate, and its
execution and evaluation freeze remain outstanding.
