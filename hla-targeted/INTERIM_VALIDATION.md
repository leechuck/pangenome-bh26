# Historical interim validation result

The comparator has since completed. See [the final report](REPORT.md) and [calibration supplement](COMPARATOR_CALIBRATION.md). The text below preserves the earlier interim state.

The targeted assay and depth comparisons are complete for all 24 held-out donors. The dedicated C4Investigator runs are still pending; placeholder rows in intermediate metric files must not be interpreted as its performance. The final report will include the completed native and separately calibrated total-C4 comparator.

- Total C4 dosage: 24/24 correct, tied with calibrated ordinary reference depth.
- A/B marginal dosage: 22/24 correct; one incorrect call and one abstention.
- Long/short marginal dosage: 22/24 correct; two abstentions.
- Full RCCX signature pairs: 21/24 correct, versus 22/24 for the old sketch.
- Module-start context ablation: 21/24, unchanged.
- The post-freeze assembly audit covers 91 C4 genes; all diagnostic motifs and long/short labels agree with the inherited annotations.

There is no held-out accuracy gain from the targeted constraints. The development gain did not transfer. All donors retain multiple dosage-compatible structure pairs; see `STRUCTURAL_ERRORS.md` for concrete phase counterexamples. One error is an absent reference structure, another is an incorrect A/B constraint, and another is a ranking error despite correct marginal dosage.

This is new-donor read validation within the original assembly panel, with all selected families excluded from marker discovery and candidate references. It is not independent-cohort or clinical validation. The copy-number distribution is one two-copy, six three-copy, fourteen four-copy and three five-copy donors; no six-copy case was selected. Validation populations are AFR=9, AMR=7, EAS=4 and SAS=4, with no EUR donors.

The model remains frozen; no changes were made to improve these observed outcomes.
