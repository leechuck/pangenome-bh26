# Completed reference controls (development only)

All six methods completed all 64 donors, with matched reads and IPD 3.65.
These references use native T1K alignment; they are not graph inference results.

| Reference | Two-field | Four-field |
|---|---:|---:|
| Original T1K | 495/509 | 163/323 |
| Original + HPRC contexts | 389/509 | 165/323 |
| Original + HPRC + Asian contexts | 465/509 | 232/323 |
| Full genomic IPD | 494/509 | 285/323 |
| Full genomic IPD + HPRC contexts | 431/509 | 203/323 |
| Full genomic IPD + HPRC + Asian contexts | 470/509 | 259/323 |

Full genomic IPD gives the strongest four-field result. Its two-field point
accuracy is effectively unchanged: a net one-call difference (13 gains, 14 losses).
Paired family bootstrap, 10,000 replicates, gives a two-field difference of
-0.20 percentage points (95% interval -2.17 to +1.76) and a four-field difference
of +37.77 points (+32.38 to +43.26). This does not prove equivalence or independent
improvement. The interval source hash binds the complete six-method score table
in ../genomic-ipd-paired-uncertainty.json.

Adding observed sequences as separate native T1K candidates reduces accuracy
relative to full genomic IPD alone. Asian contexts improve over HPRC contexts,
but that comparison does not establish incremental benefit over genomic IPD.
The graph inference experiments remain separate. The 84 reserved validation
donors remain unscored. This archive supersedes v1's four-method snapshot.
