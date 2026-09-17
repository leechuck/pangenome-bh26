# HLA pangenome analysis results

- [HLA recovery and candidate discovery](hla-analysis/FINAL_REPORT.md)
- [Structural typing benchmark](hla-structural/REPORT.md)
- [Structural typing limitations and next experiments](hla-structural/FAILURE_ANALYSIS.md)
- [Targeted C4/RCCX validation](hla-targeted/REPORT.md)
- [Dedicated C4 comparator calibration](hla-targeted/COMPARATOR_CALIBRATION.md)
- [Structural phase counterexamples](hla-targeted/STRUCTURAL_ERRORS.md)
- [Targeted assay reproduction](hla-targeted/README.md)
- [Mixed-Asian pilot execution](hla-asian50/run/README.md) (40 1000 Genomes donors submitted; Arabian validation deferred)
- [Revised mixed-Asian benchmark design](hla-asian50/MIXED_ASIAN_DESIGN.md) (stratified analysis; East/South Asian current run; Arabian extension deferred)
- [Asian 50-donor feasibility experiment and literature assessment](hla-asian50/README.md) (accuracy benchmark not yet run)
- [Extended HLA panel](hla-structural/results/extended_HLA_panel.tsv)
- [HLA reproduction instructions](hla-analysis/REPRODUCIBILITY.md)
- [Structural reproduction instructions](hla-structural/README.md)

The repository includes code, reports, result tables, literature, provenance and
compact structural measurement archives. Large reference inputs, reconstructed
FASTA sequences, graph traversal dumps, matrices and compiled/cache files remain
in the working dataset rather than Git. Their identities and hashes are preserved
in the analysis manifests; source locations and cluster job scripts are documented
in the reproduction instructions. A Git checkout alone is not a complete input
dataset. Restore the recorded inputs before running the full workflows.

For the structural workflow, extract `source/background_counts.tar.gz` and
`source/measurements.tar.gz` within `hla-structural/source/`. Restore or regenerate
`locus_sequences.fa.gz` using the documented DDBJ extraction job, then run
`prepare_markers.py` to rebuild the sequence-profile matrix. Do not interpret the
presence of a manifest entry as evidence that its corresponding large file is
included in Git. Historical manifests describe the complete working dataset.
