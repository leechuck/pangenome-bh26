# Targeted C4/RCCX experiments — in progress

This directory is a development checkpoint, not a completed benchmark. Final validation results are not yet available.

We are testing C4 A/B diagnostic sites, long/short HERV insertion junctions, and module-start context probes, alongside the published C4Investigator caller. Module-start context probes are not established structural-variant breakpoints. Dosage evidence must be distinguished from chromosome phase and module order.

The previous 106 donors are development data only. A fixed set of 24 previously untested donors was selected before examining their read-based typing results; all their recorded families are excluded from marker discovery and reference paths. These donors are drawn from the existing assembly panel, not an independent external cohort. Public MHC-recruited WGS reads have been retrieved for all 24. Selection and source URLs are recorded in `source/validation_selection.json` and `source/validation_donors.tsv`.

Development references contain 543 haplotypes and 1,090 annotated C4 genes. Diagnostic-site auditing found three noncanonical A-like motifs; these are candidates for further investigation, not validated novel alleles. All 299 targeted probes passed a full-training-MHC specificity audit, and the fragment counter passed synthetic contract tests. Development read measurements and C4Investigator runs are underway. Dosage calibration, model freezing and fresh-donor evaluation remain to be completed. See `PROTOCOL.md` for endpoints and limitations.

## Reproduction

The scripts use the local cached inputs documented in `../HLA_RESULTS.md` and DDBJ paths in the Slurm scripts. Run `select_validation.py`, `prepare_targets.py`, `align_targets.sbatch`, `build_markers.py`, `prepare_profiles.py`, and `compile_background.sbatch`; validate marker specificity before running `measure.sbatch`. Install and run the dedicated comparator using `install.sbatch` and `c4investigator.sbatch`. Freeze the final model before inspecting fresh-donor results.

The third-party C4Investigator checkout is excluded from Git. Clone https://github.com/Hollenbach-lab/C4Investigator into `source/C4Investigator` and use the commit in `source/tool_version.json`. Large reconstructed FASTA, PAF and NumPy files, binaries and read caches are excluded and can be regenerated. The archived preprint is `literature/C4Investigator_preprint.xml`; the journal article is DOI 10.1111/tan.15273.

Canonical repository: https://github.com/leechuck/pangenome-bh26
