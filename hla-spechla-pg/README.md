# DōgoHLA — population-enriched HLA reconstruction

**DōgoHLA** (`dogohla`, version **0.1.0**) combines SpecHLA's read-based phasing
with the AsianHLA pangenome. Named after [Dōgo Onsen in Matsuyama](https://dogo.jp/en/),
where this BioHackathon work came together.

DōgoHLA is an experimental **extension of SpecHLA**, not an independently
validated replacement. Please cite Wang et al. (2023),
[SpecHLA](https://doi.org/10.1016/j.crmeth.2023.100589), and acknowledge the
underlying pangenome projects. The original `hla-spechla-pg` directory and old
experiment names remain to preserve reproducibility.

## What is frozen in 0.1.0

- AsianHLA-panel-augmented read collection with donor/family exclusions.
- SpecHLA alignment and small-variant calling.
- Corrected phase-block scoring and complementary-haplotype pairing; IPD 3.65
  phase references. Adding panel sequences at this stage is not enabled because
  it did not improve the pilot.
- Path-preserving DRB1 graph alignment and represented-allele genotyping.
- Guarded reconstruction of confident homozygous noncoding long indels, keeping
  read-phased coding differences. Low-confidence, heterozygous and ambiguous
  boundary calls are withheld. Thresholds are frozen in `structural_overlay.py`.
- Native allele naming with the fixed DRB1 query-length restriction removed.

Eight-donor development results: **36 → 57 exact whole genes out of 128**,
**61 → 62 correct two-field genotypes out of 63**, and **22.8% fewer total
sequence edits**, compared with the saved native SpecHLA baseline. These are
selected development results. The earlier 43/128 reference was an already
panel-augmented control, not vanilla SpecHLA.

[Results and limitations](RESULTS-2026-09-18.md) ·
[Implementation and reproducibility](IMPLEMENTATION.md) ·
[Workflow review](review-2026-09-18/REVIEW.md) ·
[Earlier experimental workflow](LEGACY-WORKFLOW.md)

## Validation

The next run freezes the method before evaluating the other 32 1000 Genomes
samples (16 EAS, 16 SAS), with family-excluded references and a native SpecHLA
comparator. All donors and outcomes, including failures, must be reported.
This is a prospective extension of the development cohort, not the untouched
228/946 locked validation. See `validation/20260918/PROTOCOL.md` for the
prespecified design, deadlines and scoring rules.

## Tests

```sh
python3 -m unittest discover -s hla-spechla-pg -p 'test_*.py' -v
```

The implementation uses SpecHLA 1.0.12, vg 1.76.1, IPD-IMGT/HLA 3.65 genomic
phase references and the pinned AsianHLA panel. Alignment/phasing runs belong
on Slurm; the supplied deployment and job scripts target DDBJ. Bulk read and
sequence data are kept outside Git. This release is for research evaluation;
heterozygous structural phasing and generalisation remain open work.
