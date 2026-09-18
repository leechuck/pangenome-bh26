# T1K plus an Asian-enriched personalized graph

Status: preparation. No improvement established. Existing benchmark runs remain
frozen and finish before new improvement experiments are launched.

## Scope and controls

Complete the matched-64 and Gourraud-946 benchmark and its existing comparators,
including final failure accounting, per-ancestry output, and panel ablations.
HLA-HD is postponed indefinitely at Robert's request on 18 September 2026.
It is outside the current benchmark completion requirements and has no results.

The new inference pipeline retains T1K's joint allele competition. It uses read
k-mers to select local paths from HPRC plus Asian haplotypes, aligns reads to the
resulting reduced graph, and refines allele pairs with genomic-context evidence.
All-IPD alternatives remain available; graph absence is not allele absence.
Both mates and competing loci/pseudogenes must be considered. Uncertain four-field
calls remain ambiguous rather than becoming confident population imputations.

Compare these stages on identical read recruitment and IPD 3.65:

1. Frozen T1K 1.0.6, explicit four-field output (current baseline).
2. Graph-derived sequence references with T1K's existing alignment/inference.
   Include a full-genomic-IPD control: substitute available, valid genomic
   entries while retaining every original candidate and all partial-reference
   fallbacks. This tests whether any gain comes from restoring intronic sequence
   that T1K's default DNA database truncates, independent of Asian graph content.
3. T1K candidates plus graph-supported candidates, personalized graph alignment
   and joint candidate-pair refinement.
4. The same refinement with HPRC alone versus HPRC plus Asian haplotypes.
   Compare selection enabled/disabled to isolate personalization.

Do not replace HPRC with a size-matched Asian subset in the primary comparison.
Retain that older experiment as a secondary panel-composition control. Adding
sequence information and changing alignment are separate effects.

## Development and locked evaluation

Already examined matched-64 and original development samples are development
data. Gourraud results have also been observed and cannot become an untouched
confirmatory test merely by splitting the file now. Do not select samples on
where T1K is wrong.

reserve_validation.py reserves every remaining assembly-cohort family with no
indexed previous exposure, consuming only donor IDs and metadata. The entire
reserved cohort and its known relatives are excluded from every inference panel,
including alternate assembly aliases. Development donors/families are excluded
too. A reservation is immutable. Verify read availability before freezing the
final analysis set; any exclusion must be outcome-independent and documented.
The reserved cohort is ancestry-imbalanced; report the strata and do not call it
an ancestry-matched experiment.

Freeze code, references, parameters, candidate-selection rule, failure handling,
and one primary candidate before scoring reserved outcomes. Reference preparation
may use training-path sequence labels only. Inference must not read held-out
assemblies, catalogue rows or truth calls. The evaluator alone receives them.

## Evidence required to claim improvement

Primary endpoint: exact diploid four-numeric-field genotype accuracy among loci
whose two truth haplotypes have exact genomic labels. Use the existing ambiguity-
aware scoring rule, identical eligibility for all methods, no padding of short
names, and zero credit for failed/no-call eligible genotypes. Novel/unlabelled
truth is reported separately, never silently counted as successfully typed.

Require a positive lower bound of a prespecified 95% paired family-cluster
bootstrap interval (10,000 replicates, seed 20260918) for the selected candidate
minus frozen T1K on the complete reserved cohort. Also require that the two-field
point accuracy does not decrease, and its paired 95% interval excludes a loss
larger than one percentage point. This supports a four-field improvement claim,
not an unqualified improvement at every resolution or ancestry.

To attribute improvement to the Asian graph, report the additive HPRC-plus-Asian
versus HPRC-only contrast separately under the same inference procedure. Do not
attribute a change caused only by additional IPD sequence to graph alignment.
Secondary method/ancestry comparisons are exploratory unless separately frozen
with multiplicity handling before evaluation.

If the primary test fails, report that result. Further development must use a
new independent validation set or explicitly account for repeated testing; never
relabel the failed validation cohort as untouched. Improvement is an experimental
outcome, not a guaranteed property of the implementation.

## Technical checks before graph refinement

Confirm that the existing graph satisfies the chain and path requirements of
vg haplotype sampling. The personalized-pangenome paper assumes adequate coverage
(at least about 20x) and graph-unique informative k-mers. Measure effective HLA
coverage and specificity against paralog/decoy sequences; do not assume that
nominal WGS coverage or locus-local uniqueness suffices. Count from recruited
reads before hard per-gene binning. Preserve multiple candidate paths when k-mer
evidence is inconclusive, and report fallback frequency.

T1K's default DNA references truncate long introns and may fill unknown introns
with modal sequence. Keep observed assembly sequence, uncertain sequence and
imputed sequence distinct. Deduplicate identical paths so representation counts
do not masquerade as independent read evidence or measured allele frequencies.

## Implementation checkpoint, 18 September

The experimental joint model consumes paired-fragment graph placements and fits
diploid pairs under a shared mixture over competing genes. Duplicate placements
do not add evidence; tied candidate pairs remain explicit. Alignment scores are
converted to relative emissions, not calibrated probabilities. Coordinate ascent
can find a local optimum; its score gap is not a confidence measure.

The synthetic training-path integration recovered A*02:01:01:01 plus
A*02:06:01:01 from 2,004 pairs (IBEX job 52044705; joint-smoke-result.json).
This single-locus, error-free test verifies integration only. All-IPD fallback,
real library calibration, cross-locus/decoy evidence and independent validation
are still required before this becomes the proposed refinement method.

run_linear_control.py prepares the separate native-T1K development experiment
with added observed sequences. It retains raw genotype tables, decodes internal
context identifiers to explicit two-/four-field alternatives, hashes inputs and
records failures. These controls use native T1K alignment, not graph alignment.

The matched-64 controls are queued behind the remaining legacy benchmark jobs
(linear-control-launch.json). Each panel has a one-donor execution pilot followed
by 63 samples conditional on pilot success. The scorer reproduces all 1,024
baseline donor/locus/resolution rows, including eligibility, called genotypes,
correct genotypes and allele matches. Baseline totals remain 495/509 at two fields
and 163/323 at four fields. All predicted alternatives must agree with truth;
unknown alternatives cannot be silently removed. These samples remain development
data, irrespective of the resulting improvements or regressions.

The isolated native-evidence exporter modifies only T1K 1.0.6's two assignment
printing sites. It appends the existing weight, qual and adjustWeight fields
before read assignments are coalesced. Those values encode T1K's own heuristics;
they are not posterior probabilities and must not be multiplied by graph scores
as though the same reads were independent observations. Native filtering still
applies, including the maximum assignments per fragment. No baseline executable
or queued linear-control executable is replaced. The build compares original
and instrumented genotype tables and the first four assignment columns on the
synthetic library. Real-data preservation and graph integration remain required.

The synthetic export preservation check passed (evidence-build-result.json):
883,433 assignment rows preserved, with identical genotype tables. The parser
retains 708 assigned fragments from the 2,004-pair training library; the graph
adapter supports all 2,004. This is a difference in available reference-context
evidence on synthetic data, not a demonstrated typing improvement.

kmer_coverage.py prepares canonical 29-mer markers conserved in at least 90% of
the distinct additive-panel sequences, excludes markers repeated within any
path, and removes cross-gene matches in all supplied panel and IPD sequences.
Coverage uses the median across all markers, including zero-count markers. The
initial development gate requires 100 markers, median k-mer depth at least 20,
and at most 25% zero-count markers. These thresholds are uncalibrated development
parameters. Low support must trigger unsampled-graph/native fallback rather than
force a path selection. Marker specificity has not yet been screened against
the entire genome. No validation outcome is used to choose markers or thresholds.
The first synthetic coverage check measured HLA-A median k-mer depth 69 from
1,321 markers (nominal expectation approximately 65); all seven absent loci had
zero depth and failed the personalization gate. A development-donor coverage
pilot is queued separately from genotype inference.

The HG00658 development coverage pilot completed with usable marker depths of
27–49 across all eight loci (1,477,842 unbinned reads). map_personalized.py checks
that the coverage report hashes match the exact input reads, uses that measured
depth for path sampling, and retains all alignment records including unmapped
reads. Selection-disabled or low-coverage operation uses the unsampled graph.
It builds a full distance index for mapping; the distanceless index used during
haplotype preprocessing is insufficient for giraffe. Both sampled and unsampled
synthetic mapping retained and mapped all 4,008 reads. HG00658 HLA-A mapping pilots
for both panels and both selection settings are queued after the remaining
legacy benchmark. These outputs are evidence inputs, not complete HLA calls.

join_evidence.py verifies native/graph input read hashes and the mapping-to-path-
support-to-fragment evidence chain. A disk-backed join retains one row per
fragment across native assignments and multiple graph loci. Mate suffix handling
matches T1K's reader; collisions within an evidence source fail explicitly.
The synthetic integration retains exactly 2,004 fragments: 708 with both sources
and 1,296 with graph-only support. Scores are not combined as independent reads.
The first real native-evidence pilot reuses frozen T1K candidate FASTQs and must
reproduce its genotype table byte-for-byte before its evidence is accepted.

The full-genomic-IPD control reference contains the same 29,429 candidate IDs as
the frozen database: 26,745 genomic substitutions and 2,684 retained fallback
entries. The input audit found 27 inherited baseline exon-coordinate exceptions;
they are recorded rather than silently corrected. All genomic replacements pass
coordinate validation. Its matched-64 pilot and array are queued behind the
legacy benchmark. Marker specificity is also being rechecked against the fuller
IPD sequences, which expose intronic paralog matches absent from truncated entries.

The full-genomic-IPD marker screen completed (jobs 52045734/52045735;
coverage-genome-development-result.json). DRB1 retained only 49 specific markers,
down from 189 with the truncated reference, and therefore fails the prespecified
100-marker development gate despite a median depth of 52. It must use unsampled
graph evidence; high apparent coverage alone does not justify personalization.
The other seven loci retain usable marker sets, with median depths of 27–36.
These are HLA/IPD-relative specificity results, not genome-wide uniqueness tests.

All eight loci are queued for HG00658 under both panels and both selection
settings using this fuller screen (development-map-all-loci-launch.json).
Each array depends on its original HLA-A mapping pilot and the remaining legacy
benchmark. Code is copied into a separately hashed directory. HLA-A is repeated
because its marker set and measured coverage changed; original pilot outputs
remain available. DRB1 in the selection-enabled arm must record its fallback.
These mappings still provide evidence only; joint inference and calibration are
not completed by this experiment.

All 84 reserved donors now have verified read preparation records, archived in
validation/READ_PREPARATION.json (65,467,081 paired fragments). This archives the
preparation manifests and read hashes; it does not genotype or unblind them.

## References

- Song et al. (2023), https://doi.org/10.1101/gr.277585.122
- Siren et al. (2024), https://doi.org/10.1038/s41592-024-02407-2
- Source/paper review: ../research/t1k-review/REVIEW.md
