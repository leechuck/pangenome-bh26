# Independent HGSVC28 follow-up

The original 84-donor validation remains a failed combined protocol gate. Model
changes on those exposed donors are development results. This follow-up retains
all 28 donors in RESERVATION.json; no selection by genotype, prediction, eligibility
or accuracy is permitted. It is ancestry-imbalanced, not an ancestry-matched test.

## Selected candidate and controls

The selected candidate is genomic-IPD T1K followed by unsampled HPRC+Asian graph
refinement, using alignment-only fragment emissions and the native-locus veto,
then the original-T1K coarse anchor. Parameters and fallback rules are those in
`graph_pair_refinement_native_locus.py`, `run_graph_pair_native_locus.py` and
`anchor_t1k_coarse.py`, as hashed in PRIMARY_CANDIDATE.json. No outcome-dependent
tuning is permitted on these 28 donors.

Compare original T1K, coarse-anchored genomic-IPD T1K, coarse-anchored HPRC
refinement, and coarse-anchored HPRC+Asian refinement, with identical reads and
IPD 3.65. Coarse anchoring may add a resolved four-field call only when its
projection agrees with the original two-field call. Report fallback counts and
verify per-gene two-field prediction identity, including ambiguity and no-calls.
An execution failure is not silently converted into a successful fallback.

## Fixed truth derivation

Use HGSVC3's published sense-oriented HLA allele archive. Ignore its annotated
allele labels. Parse donor, haplotype 1/2 and gene identifiers only; require one
record per donor/haplotype/gene. Remove exactly 1,000 bases from each end, as
specified in the archive README. Compare the entire remaining gene sequence to
exact genomic sequences in the pinned IPD 3.65 FASTAs. No local alignment,
substring matching, boundary adjustment, CDS substitution or name padding.

Keep all exact genomic name alternatives. Empty cores, ambiguous/non-ACGT cores,
missing or duplicate haplotype-gene records and unmatched cores are unresolved.
A diploid genotype is eligible only if both haplotypes pass. Four-field eligibility
additionally requires every exact name to have four numeric fields. Two-field
truth is the two-field projection of these exact genomic names (including exact
names shorter than four fields when they have at least two). This two-field
eligibility differs from the original cohort's CDS-derived truth and must be
reported explicitly. No invented CDS annotation is used.

Use all eight loci: A, B, C, DPA1, DPB1, DQA1, DQB1 and DRB1. Preserve all donors
and report eligibility reasons for every haplotype. Novel/unmatched truth is not
credited as correctly typed. Every method has the same eligible denominator.
Numeric endpoints omit expression suffixes, matching the earlier evaluator.
Ambiguous predictions are correct only if all alternatives are compatible with
truth under a valid unphased diploid assignment. Missing/failed eligible calls
score zero. Published assembly-derived truth is not an independent clinical
four-field typing assay; boundary/assembly errors remain a limitation.

## Endpoints and interpretation

Primary: candidate minus original T1K exact diploid four-field accuracy, paired
family-cluster bootstrap with 10,000 replicates, seed 20260918. Require a positive
lower 95% bound; also require nondecreasing two-field point accuracy and a paired
lower 95% bound greater than minus one percentage point. Retain all failures.
Report the full table, per-ancestry results, discordant genotypes and intervals.

Also report candidate minus anchored genomic-IPD and minus anchored HPRC under
the same scoring rules. These incremental contrasts remain exploratory and
cannot be replaced by the larger original-T1K contrast when describing graph
benefit. The reference-only gain is not evidence that Asian graph paths help.
Do not pool this cohort with exposed donors as an independent validation result.
Do not claim improvement if the gate fails; any subsequent tuning makes this
cohort development data and requires another independent test.

## Release conditions

Truth-matching code has synthetic tests, but must not process held-out sequences
until the full inference/evaluation execution freeze is recorded. Freeze all
code, assets, read hashes and metadata, and verify reference graph manifests
against the audited source panels before typing. Evaluation requires the complete
terminal donor/method grid and run provenance, not just successful cases.
PRIMARY_CANDIDATE.json selects the model and truth rules; it does not certify
that the execution freeze or validation has completed.
