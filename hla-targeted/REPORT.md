# Targeted RCCX/C4 typing: held-out read experiment

<!-- completed-controls:start -->
**Main finding:** targeted diagnostic constraints did not improve held-out structural-pair accuracy: 21/24 versus 22/24 with the old sketch. Total C4 dosage was 24/24, tied with calibrated ordinary reference depth. Module-start context probes added no accuracy in this set. These are new-donor results within the original assembly panel, not external-cohort or clinical validation.

The completed native C4Investigator total-copy result is 14/24; its separately declared pilot-calibrated total-copy result is **24/24**. Native component comparisons can inherit total-depth bias and are not evidence of general superiority of the targeted method. See [the calibration supplement](COMPARATOR_CALIBRATION.md) and [the configuration check](COMPARATOR_METHOD_CHECK.md).

The three structural errors illustrate incorrect A/B dosage, unresolved linkage of A/B with long/short forms, and an absent reference structure. Every donor retains multiple compatible structures. See [the actual phase counterexamples](STRUCTURAL_ERRORS.md). All 91 validation genes contain their full diagnostic probe sets, and the whole-MHC control finds no extra diagnostic-probe matches outside those genes in any of the 48 validation haplotypes. This excludes those specific explanations for the A/B error, not all sampling or mapping effects.
<!-- completed-controls:end -->

This experiment evaluates 24 newly examined donors after using all 106 previously examined donors for development. All 24 selected families were excluded from new targeted-probe discovery, training statistics and candidate reference paths. The inherited global marker vocabulary is distinguished from its training-only eligibility filters in [the marker provenance audit](MARKER_PROVENANCE.md). These samples remain part of the original assembly graph: this is new-donor read validation within the panel, not an independent external-cohort validation.

## Marginal dosage

| Method | Feature | Correct / attempted | Calls made |
|---|---|---:|---:|
| Targeted fragment probes | total | 24/24 | 24 |
| Targeted fragment probes | A | 22/24 | 23 |
| Targeted fragment probes | B | 22/24 | 23 |
| Targeted fragment probes | OTHER | 23/24 | 23 |
| Targeted fragment probes | L | 22/24 | 22 |
| Targeted fragment probes | S | 22/24 | 22 |
| C4Investigator | total | 14/24 | 24 |
| C4Investigator, pilot-calibrated total | total | 24/24 | 24 |
| C4Investigator | A | 17/24 | 24 |
| C4Investigator | B | 21/24 | 24 |
| C4Investigator | L | 18/24 | 24 |
| C4Investigator | S | 20/24 | 24 |
| Reference depth / WHR1 | total | 24/24 | 24 |
| Reference depth / TNXB | total | 24/24 | 24 |

Abstentions count as failures in all-attempted accuracy. A/B and long/short are separate marginal dosages; they do not directly determine AL/AS/BL/BS combinations or physical phase. OTHER denotes a noncanonical diagnostic motif, not a confirmed functional allele class.

## Structural-pair imputation

| Method | Exact signature pair / attempted | Truth retained in compatible set | Median compatible pairs |
|---|---:|---:|---:|
| Old sketch + total CN | 22/24 | 23/24 | 213 |
| Targeted + marginal dosage | 21/24 | 22/24 | 17 |
| Targeted without module contexts | 21/24 | 22/24 | 17 |

The validation set contains nine AFR, seven AMR, four EAS and four SAS donors. There are no EUR donors in this selected set, limiting direct extrapolation to predominantly European UK Biobank samples.
Compared with the old sketch, 23 donors are unchanged, 0 improve and 1 worsen. The three errors comprise one incorrect marginal dosage, one wrong ranking despite a compatible truth, and one absent reference structure. Full donor-level dispositions are in `results/validation_case_analysis.tsv`.

The targeted model changes exact signature-pair accuracy from 22/24 to 21/24. Omitting module-start context probes gives 21/24. These are small-sample exploratory comparisons; ranked imputation does not establish physical phase or a graph-specific advantage.

All compatible reference structure pairs are retained in the prediction table. A unique numerical optimum is labelled imputation, not direct molecular phase. The sets use heuristic dosage calls and are not calibrated confidence sets; dosage errors and absent reference structures can exclude the truth. Full signatures include provisional CYP21/TNX prototype annotations, not established functional alleles.

## Assay and development

The reference contains 543 training haplotypes and 1,090 annotated C4 genes. All 299 targeted probes were confined to RCCX in the full training MHC sequences. Whole-genome specificity was not established. Diagnostic groups are counted once per fragment, with canonical exact 31-mers and Q20 filtering. A/B relative efficiency was 1; long/short relative efficiency was 0.9423. These were fitted using development donors only. Total dosage retains the earlier pilot-derived correction of 1.15263.

The assay was frozen at 2026-09-16T04:22:36.798348+00:00; the file hashes are in `source/validation_freeze.json`. The 106 development outcomes are exploratory resubstitution results, not new validation. Development C4Investigator covers a prespecified subset of 18; validation C4Investigator covers all 24 donors. Its upstream version and command were fixed independently of the targeted model.

## Assembly labels and diagnostic audit

The held-out assembly audit assessed 91 C4 genes. 0 differed from the inherited annotation or had incomplete diagnostic alignment. Original annotation comparisons are retained in the main accuracy tables. The full audit is in `results/validation_C4_diagnostic_audit.tsv`.

Three training genes had a noncanonical A-like five-site motif (GGTGC), rather than the canonical A motif (CGTGC). Their observations are candidates for further sequence investigation, not validated novel HLA types or functional C4 alleles. NA19700 has 16 deduplicated fragments supporting the noncanonical motif in the development data; the other two candidates were not read-tested here. See [the candidate record](NONCANONICAL_CANDIDATES.md). Read support does not establish novelty, function or physical phase.

## Limits and biological interpretation

The read source is public WGS recruited to MHC intervals through existing linear-reference mappings. This tests typing conditional on that recruitment, not whole-WGS recruitment sensitivity. The assembly annotations are comparator labels rather than independent clinical truth. Recorded family exclusion does not rule out all cryptic relatedness.

The long/short HERV feature and A/B diagnostic sites are too far apart for ordinary short fragments to provide their physical linkage. Module-start probes interrogate sequence contexts around WHR1-family annotations; they are not proven SV breakpoints. New structures outside the reference cannot be safely rejected as unknown by the current ranking procedure.

These features can extend classical HLA typing with C4 dosage, HERV state and an explicitly uncertain RCCX structural interpretation. The earlier study identified RCCX variation among haplotypes sharing complete classical two-field HLA labels. Disease association, expression consequences, and added prediction beyond classical HLA require later validation; they are not established here. UK Biobank WGS/WES and phenotype access is available but remains outside this hackathon analysis.

## Prior work

C4Investigator already performs comprehensive short-read C4 dosage and sequence analysis; this experiment is not the first C4 structural typing method. Its published limitation on phasing A/B with long/short is directly relevant here. See [Marin et al., HLA, DOI 10.1111/tan.15273](https://doi.org/10.1111/tan.15273), the archived preprint in `literature/`, and [the upstream software](https://github.com/Hollenbach-lab/C4Investigator). Broader graph and RCCX prior art is archived in the preceding HLA studies.

## Reproduction and deliverables

Selection, probe construction, counting, calibration, inference, evaluation, post-freeze diagnostic auditing and completion checks are implemented in the scripts in this directory. `PROTOCOL.md` records the split and endpoint definitions. `source/job_manifest.json` records cluster jobs. Raw sequence inputs and reads remain outside Git. Compact derived profiles and read measurements are bundled in `source/inference_inputs.tar.gz`, with a verified byte-identical replay recorded in `results/replay_check.json`. The separate held-out MHC background archive supports the specificity control; full sequence-level checks require the original sequence inputs. Paths and regeneration procedures are in `README.md` and the Slurm scripts.

- `results/validation_dosage_predictions.tsv`: estimates, calls and abstentions.
- `results/validation_path_predictions.tsv`: ranked pairs and complete dosage-compatible sets.
- `results/validation_metrics.tsv` and `validation_errors.tsv`: all-attempted metrics and disagreements.
- `results/validation_depth_baseline.tsv`: independently calibrated depth comparator.
- `results/validation_C4_diagnostic_audit.tsv`: post-freeze assembly audit.
- `results/completion_checks.json`: provenance and denominator checks.
- `results/validation_summary.png` and `.pdf`: comparison figure.

Repository: https://github.com/leechuck/pangenome-bh26
