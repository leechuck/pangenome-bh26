# Asian HLA/MHC benchmark: literature and 50-donor feasibility experiment

**Execution update:** Robert explicitly authorised execution. [The initial mixed-cohort run](run/README.md) has started on DDBJ (20 East Asian + 20 South Asian; Arabian validation deferred by Robert). The feasibility status and decision below are historical.

**Scope update:** Robert has requested separate East Asian, South Asian and West Asian/Arabian analyses. See [the revised mixed-Asian design](MIXED_ASIAN_DESIGN.md). The East-Asian-only selection and results below are historical; the current run is now limited to 1000 Genomes; West Asian/Arabian validation is deferred.

Status (17 September 2026): literature review, deterministic cohort selection, remote input audit and held-out structural-representation experiment complete. **The 50-donor short-read typing / SNV / SV accuracy benchmark has not run.** There are no new read-derived accuracy results in this directory.

The general proposed approach already exists. The useful remaining hypothesis is whether **this Asian-enriched panel** improves inference beyond an HPRC-only panel and existing specialist tools. Robert's execution condition was “If not done and can improve, do it”; the literature does not satisfy the “not done” condition. We therefore completed a bounded feasibility experiment and specified the comparison below, without launching a redundant large benchmark or claiming novelty. The absence of a coarse-structure coverage gain does not rule out a sequence-level gain.

## Completed small experiment

We selected 50 East Asian 1000 Genomes donors with two eligible RCCX haplotype annotations, from 55 eligible donors. East Asian scope was the default proposed to Robert; South Asian donors are not included. Selection uses metadata and experimental-label availability, never read predictions or their correctness. It retains all 23 donors with experimental HLA labels, then fills by population balance and a fixed hash. The resulting populations are CDX 5, CHB 5, CHS 10, JPT 19, KHV 11. Five folds contain ten distinct families each. All duplicate graph paths for an excluded donor/family are listed for exclusion.

The cohort deliberately prioritizes paired assembly availability for SV evaluation. It is not population-representative; JPT enrichment partly reflects experimental-label availability. An alternative classical-HLA-only cohort could have 50 experimental labels, but would have less paired assembly coverage.

For each test fold we removed all fold families from the existing eligible structural catalogue and checked whether each test haplotype's annotated signature exists in the remaining panel:

| Catalogue-level endpoint | Full panel | HPRC-only |
|---|---:|---:|
| DRB signature represented | 100/100 | 100/100 |
| RCCX signature represented | 95/100 | 95/100 |

These are **representation ceilings for exact signature matching**, not genotyping accuracy. They are derived from existing graph-path/assembly annotations, which have their own QC limitations. RCCX signatures include provisional CYP21/TNX labels. DRB signatures are coarse gene-content/order classes. Neither endpoint captures all nucleotide variation, breakpoints or regulatory differences. A method that constructs previously unobserved haplotypes is not bounded by this catalogue ceiling. “HPRC-only” includes Asian HPRC donors: it tests the added cohorts, not Asian ancestry versus no Asian ancestry.

The immediate inference is narrow: the added non-HPRC assemblies contribute no additional represented *signature classes for these test haplotypes*. Any advantage must come from finer sequence, better allele evidence, frequency estimates, or performance on another cohort. It cannot be assumed from panel size.

## Truth and input audit

- 50/50 donors have both MHC assembly FASTAs on DDBJ and eligible paired RCCX annotations.
- Only 23/50 have experimental classical HLA labels in the current Gourraud 2014 resource: at most 115 donor–gene comparisons for A, B, C, DRB1 and DQB1. Retain missing/ambiguous labels at scoring; do not silently turn ambiguity into a single allele. These data do not establish whole-gene SNV or SV truth.
- 33/50 donors were previously evaluated in our structural work. This is a new prespecified comparison, **not a fresh blind validation cohort**. Report the 17 previously unused donors separately; avoid tuning on either set's new predictions.
- On a001, 30/50 full WGS CRAM paths are readable; 33/50 have an existing MHC-recruited CRAM (including our prior cache). Public CRAM URLs were recovered for all 50; this records locations, not successful new downloads.
- 30/50 have legacy T1K results. Those used a different preset and an incompletely frozen database; do not silently reuse them as the final controlled comparator.
- Existing MHC CRAMs omit unmapped and some off-region reads. Equal restricted reads permit a conditional comparison but do not prove reduced genome-wide recruitment bias. Whole-WGS recruitment is required for that claim.
- `MHC.raw.vcf.gz` was generated from the **clipped** graph (the header records `0.clip.raw.vcf.gz`); its chromosome is local `ctg1`, not genomic chr6. Do not feed it straight into an unqualified full-MHC benchmark. Full graph paths, reference sequence and coordinate conversion need auditing.

## Prespecified experiment to implement

Primary question: does the full Asian-enriched panel improve held-out sequence/SV genotype accuracy relative to an HPRC-only panel, with the same inference method, read input, coverage, QC and reference coordinates?

1. Use the fixed five family-separated folds. Exclude both haplotypes, donor aliases and known relatives before constructing each inference panel. Retain test-only alleles in truth but not training; score them as out-of-panel, not silently remove them. Derive frequencies, marker selection and confidence thresholds from training/development only. Any use of a globally constructed graph must be explicitly labelled transductive; strict inductive claims require a training-only reconstruction.
2. Use a validated PanGenie-compatible phased bubble panel from the full graph for SNV/indel/SV re-genotyping. PanGenie is an existing method, not our invention. Audit non-overlap, sequence-resolved alleles, missing genotypes, copy-number representation, reference coordinates and haplotype reconstruction before inference. For repeated RCCX structures that cannot be represented faithfully in that panel, evaluate complete locus haplotypes with Locityper or COSIGT as a separate endpoint. Do not force duplicated C4 copies into ordinary diploid per-copy SNP calls.
3. Run **T1K, SpecHLA and HLA*LA** as the three classical HLA comparators. All are established WGS-capable methods; this is not an assertion that one is universally best. Use T1K's documented `hla-wgs` preset and SpecHLA's whole-gene WGS mode. Record software commits, native database releases and any unavoidable database differences. HLA*LA is itself graph-based; beating it would not establish “graph versus linear reference.”
4. Add appropriate non-HLA baselines: single-reference small-variant calling (a frozen DeepVariant WGS pipeline), SV calling (a frozen supported SV pipeline, selected before predictions), and known-site SV genotyping if the claim concerns genotyping rather than discovery. Reuse C4Investigator and calibrated depth as C4 dosage controls. A classical HLA typer's failure to output an SV is **not an SV false negative**.
5. Use identical whole-WGS-derived inputs for primary comparisons. Include alternate/decoy and unmapped reads when recruiting. Record reads/fragments retained, depth and runtime. A restricted-MHC fast pilot may be reported separately, with its recruitment limitation.
6. Freeze assembly-derived truth and callable masks independently of predictions, verify important discordances against original long-read evidence, and harmonize variants by sequence rather than graph node IDs. Graph-deconstructed assembly genotypes alone are a secondary, representation-dependent comparison. Original assembly availability does not mean a validated SNV/SV truth VCF already exists.
7. Keep three separate scorecards: experimental classical two-field HLA calls (23 labelled donors); sequence-resolved SNV/indel/SV calls within audited callable intervals (eligible subset of 50, denominators explicit); DRB/C4 gene content, dosage, order and phased locus-pair accuracy. Report noncoding and coding small variants separately, SVs by class/size, seen versus absent alleles, unresolved calls and method failures.
8. Pair comparisons by donor. Report genotype precision/recall/F1 and non-reference genotype concordance rather than raw reference-dominated accuracy. Use donor-level confidence intervals/bootstrap; variants and loci within a donor are not independent replicates. For HLA show both allele-level and complete diploid-genotype accuracy. Count abstentions in the all-attempted denominator while reporting call coverage separately. With 50 donors this is a pilot, not a conclusive rare-variant or disease study.

The informative ablations are full panel versus HPRC-only, and (if resources permit) graph inference versus the same haplotypes represented as a sequence panel. This distinguishes reference diversity from graph topology or storage format. Fifty donors cannot establish disease associations or clinical utility.

## Completion statement for a subsequent accuracy benchmark

“On the frozen 50-donor East Asian cohort, complete and publish family-excluded full-panel versus HPRC-only inference, classical HLA comparison with T1K/SpecHLA/HLA*LA, and SNV/SV comparison with appropriate single-reference baselines; use explicit independent experimental and assembly-derived truth tiers, quantify uncertainty and abstentions, audit discordances, and report whether an improvement exists without requiring a positive result.”

This statement is a proposed acceptance criterion, not an active automatically created goal.

## Reproduction and files

Run `python3 hla-asian50/select_panel.py` from this repository to regenerate the cohort, fold exclusions, experimental labels and structural-representation results. Run `python3 hla-asian50/check_setup.py` to validate them. The remote path audit is `preflight_remote.py`; it only checks file availability and must run on a001 with the same relative `source/donors.tsv` path.

- [Literature assessment](literature/README.md), full-text XML/TXT articles and download hashes.
- [Selected donors](source/donors.tsv), [fold exclusions](source/excluded_paths.tsv), [summary](source/cohort_summary.json).
- [Representation results](source/structural_representation_oracle.tsv).
- [Experimental truth](source/experimental_HLA_truth.tsv), [read locations](source/public_reads.tsv), [DDBJ availability](source/remote_availability.tsv).

No 50-donor prediction jobs were submitted by this feasibility study. All existing completed HLA/C4 studies remain unchanged.
