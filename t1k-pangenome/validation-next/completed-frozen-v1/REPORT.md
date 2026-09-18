> **Superseded evaluation — two-field scoring error.** Retained for audit only.
> Use [the corrected report](../completed-scoring-correction-v2/REPORT.md) and
> [the correction record](../SCORING_CORRECTION.json). Original numeric outputs
> have not been modified.

# Independent HGSVC28 validation

All 28 reserved donors; four methods; eight loci. Four-field truth requires exact whole-genomic matches. Two-field truth is also genomic-derived. Original two-field prediction identity was verified for each anchored method.

| Method | Fields | Correct | Eligible | Failed donors |
|---|---:|---:|---:|---:|
| T1K | 2 | 0 | 144 | 0 |
| T1K | 4 | 71 | 138 | 0 |
| ipd_genome | 2 | 0 | 144 | 0 |
| ipd_genome | 4 | 120 | 138 | 0 |
| graph_hprc | 2 | 0 | 144 | 0 |
| graph_hprc | 4 | 120 | 138 | 0 |
| graph_hprc_asian | 2 | 0 | 144 | 0 |
| graph_hprc_asian | 4 | 120 | 138 | 0 |

Primary combined protocol pass: **True**.

The original-T1K contrast measures the complete pipeline change. Incremental genomic-IPD and HPRC contrasts are exploratory; reference-only gains do not prove an Asian graph advantage. Per-ancestry intervals, failure accounting and paired contrasts are in report.json. The earlier 84-donor failed combined gate remains a separate result.
