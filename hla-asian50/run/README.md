# Mixed Asian pilot execution

**Completed:** all 40 T1K, 40 SpecHLA and 40 paired-arm PanGenie jobs finished successfully; output verification passed. See the [completion record and saved call outputs](COMPLETION.md). Accuracy evaluation remains outstanding.

Execution authorised by Robert: “Run this now.” Initial jobs submitted 17 September 2026 on DDBJ, account `asianhla-group`, partition `asianhla-c32`. This supersedes earlier decisions not to launch an accuracy benchmark. The submitted calling run is complete; this does not mean the accuracy benchmark is complete.

## Current cohort

Forty donors are frozen in `source/donors.tsv`: 20 East Asian and 20 South Asian, with five eight-donor folds (four per stratum). Eight East Asian donors have experimental classical HLA labels in the current resource. Robert subsequently limited the current validation to 1000 Genomes and deferred Saudi/Arabian samples. The current cohort is therefore these 40 donors; no additional Arabian slots are required for this run. Keep East and South Asian results separate.

Selection balances populations before using experimental-label availability and a fixed hash; no prediction outcomes are used. Both haplotypes and known family/alias paths are excluded from each corresponding inference panel. Earlier analyses have examined some donors, so this is not a fresh blind cohort. `cohort_freeze.json` records selection time and file hashes. The original East-Asian-only feasibility inputs in the parent directory remain historical.

DDBJ's APR README states that raw reads were not downloaded; the inspected `data/fastq` directory is empty. The public APR ENA project PRJNA1108179 has 18 Illumina runs, all with Hi-C run titles. Thirteen have misleading `WGS` library-strategy metadata. `APR_ena_runs.tsv`, `APR_illumina_audit.tsv` and the individual run XML records preserve that audit. We do not substitute Hi-C for ordinary short-read WGS. This source audit is historical. Locating or downloading Saudi/Arabian reads is deferred and does not block the current benchmark.

## Submitted stages

1. Install PanGenie/Locityper and SpecHLA into dedicated environments. Native resolved package inventories are retained on DDBJ. Locityper is installed but not yet run.
2. Extract paired reads from existing MHC CRAMs, fetching the same regions from public CRAMs when absent. Run T1K 1.0.6 with `--preset hla-wgs` and a copied, hashed genomic reference database. Retain singleton reads separately; both initial typers receive the same paired reads.
3. Run SpecHLA full-length mode (`-u 0`), default `Unknown` population prior uniformly across strata, on the same pairs. These are the two required specialist HLA comparators. HLA*LA remains an optional third comparator, not a completed/submitted run.
4. Deconstruct top-level bubbles from **MHC.full.gbz**, not the pre-existing clipped VCF, and verify reference allele identity. The first attempt failed because vg 1.63.1 could not load GBZ v2. It was replaced with the graph-building Cactus bundle's vg 1.76.1 development build, without modifying the graph.
5. Build five family-excluded panels for each of full and HPRC-only arms. Keep canonical donors with two input haplotypes. Filter sites with missing/unphased training calls, non-ACGT alleles or overlapping reference spans; never change missing calls to reference. Remove alleles seen only in the held-out fold. Save allele remapping, retained sample lists and all filter counts. These conservative filters may remove important SVs and must be included in coverage accounting.
6. Index those panels and run PanGenie 4.2.1 for both arms on each donor, dependent on successful graph preparation and read jobs.

The global graph was originally constructed with validation assemblies. These are family-excluded inference panels on a **transductively constructed graph**, not strict training-only graph reconstruction. Full versus HPRC-only uses the same global graph topology; it tests available panel haplotypes under that representation, not independently rebuilt graph topologies.

## Limits and work still required

This initial run uses reference-recruited MHC reads, excluding unmapped/off-region fragments absent from the existing CRAMs. It cannot measure whole-WGS recruitment rescue. Ordinary whole-WGS recruitment, independently audited SNV/SV truth and callable intervals, appropriate linear variant-calling/SV-genotyping baselines, translation of graph-inferred sequences to classical HLA labels, scoring and discordance review remain required before the proposed benchmark is complete. T1K/SpecHLA native databases differ; their version/provenance must accompany comparisons.

Top-level graph bubbles are not automatically individual SNVs or SVs. Length-change site counts are preliminary bookkeeping, not a complete SV classification. Do not report graph-bubble agreement as independent SNV/SV accuracy or claim improved SV detection from these jobs alone. Test-only alleles filtered from inference remain part of the truth/coverage denominator when scoring, rather than being silently removed.

## Operation

- Remote work: `/home/leechuck/hla/codex-asian50/run`.
- Transport: `python3 remote.py 'COMMAND'`, using the established DDBJ SSH aliases and keys without inspecting credentials.
- Job IDs and initial states: `source/jobs.json`. Query current Slurm state before interpreting a stored snapshot.
- `select_cohort.py` generates the initial cohort. Do not rerun it to change a frozen cohort after reading outcomes; make a new version if expansion requires changes.
- `test_prepare_panels.py` verifies family-path exclusion, removal of held-out-only alleles and handling of missing training calls using a synthetic VCF.
- Raw reads, graphs, native indices and intermediate outputs remain on DDBJ; only scripts, manifests and audit metadata are tracked here.
