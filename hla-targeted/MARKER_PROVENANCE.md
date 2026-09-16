# Marker vocabulary and family exclusion

The 299 new targeted probe groups were derived solely from the 543 retained training haplotypes. All 24 validation families are absent from those haplotypes and candidate reference paths.

The assay also reuses the older 11,917-key measurement vocabulary. That vocabulary was a global union of canonical 31-mers passing the fixed sequence-only rule `splitmix64(key) & 63 == 0`; the vocabulary itself was not constructed from a family-excluded panel. This distinction matters for provenance.

The caller determines feature eligibility, normalization anchors and conserved-C4 markers from the retained training rows only. The post-evaluation audit verifies that all **1,085 inherited features actually used** pass that fixed hash rule and occur in at least three retained training families. All inherited training-profile columns match the original sequence profiles exactly. Recomputing the hash sketch from the retained training sequences would therefore include every inherited feature used by this caller; vocabulary entries present only in held-out paths do not enter scoring or normalization.

Thus the family-exclusion claim applies to newly discovered targeted probes, training statistics and candidate paths. It should not be read as a claim that the inherited global vocabulary or original graph construction never contained validation samples. The original graph includes their assemblies, as disclosed throughout this experiment.

Evidence: `results/marker_provenance_audit.json`, `audit_marker_provenance.py`, the frozen `infer.py`, and the earlier `hla-structural/prepare_markers.py`. This audit changes no predictions or thresholds.
