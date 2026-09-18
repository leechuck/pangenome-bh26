# T1K plus an Asian-enriched personalized graph

Status: development experiments running. Independent improvement is not yet
established. The frozen matched-64 and Gourraud-946 benchmarks have finished,
including reproducible-failure accounting.

## Scope and controls

Complete the matched-64 and Gourraud-946 benchmark and its existing comparators,
including final failure accounting, per-ancestry output, and panel ablations.
HLA-HD is postponed indefinitely at Robert's request on 18 September 2026.
It is outside the current benchmark completion requirements and has no results.

The completed Gourraud benchmark scored full/additive DogoHLA at 3255/4723,
HPRC-only DogoHLA at 3246/4723, and frozen T1K at 4437/4723 eligible two-field
genotypes. Full minus HPRC was +0.19 percentage points (paired family-bootstrap
95% interval −0.32 to +0.71), so this benchmark does not establish an additive
graph advantage. Full minus T1K was −25.03 percentage points (−26.11 to −23.98).
The final full-panel NA20790 phasing failure reproduced with identical log hashes
and unchanged configuration; its eligible genotypes remain in the denominator.
See ../hla-spechla-pg/series20260918/analysis/gourraud/ for complete per-ancestry
tables, intervals and provenance. This experimental truth supports two fields;
the separate matched-64 results also evaluate exact genomic four-field truth.

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
   Also add each observed panel to this full-genomic-IPD control. Comparing long
   observed contexts only against truncated IPD entries confounds panel content
   with the sequence available to competing alleles. These additional controls
   retain every full-genomic/fallback IPD entry and use unchanged native T1K;
   reference/pilot/array jobs are in genomic-context-control-launch.json.
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

The scorer archives IBEX's execution-cohort.tsv separately from the original
selection table: only CRAM locations differ (public URLs replace DDBJ-local
paths). It verifies identical donor order and every non-location metadata field,
requires each control to hash to that exact execution table, and still checks
paired-read hashes against the frozen baseline. A changed donor, family, fold or
selection record is rejected. This resolves the execution-table byte-hash mismatch
without accepting a change of samples or input reads.

The HG00658 execution pilot completed for both observed-sequence controls
(development-linear-pilot-result.json). Frozen T1K scored 8/8 two-field and 2/5
four-field eligible genotypes; HPRC additions scored 7/8 and 3/5; HPRC-plus-Asian
additions scored 8/8 and 4/5. Full-genomic-IPD substitution scored 8/8 and 5/5.
Thus the stronger genomic control explains more of this donor's four-field gain;
the Asian graph must also be assessed against that control. This is one already-observed development donor and
uses native T1K alignment. It cannot establish independent improvement or an
effect of graph alignment. The remaining 63 development donors are still required
to interpret these preliminary gains.

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

The real HLA-A pilots completed in 68–79 seconds each and retained all 1,477,842
input reads (development-map-pilot-result.json). HPRC mapped 12,478 reads without
selection and 8,734 with selection; HPRC-plus-Asian mapped 12,950 and 8,639,
respectively. Mapping counts do not establish specificity or genotype accuracy.
The selection effect must be assessed using concordant fragment support and
development genotypes, while keeping all-IPD alternatives available.

join_evidence.py verifies native/graph input read hashes and the mapping-to-path-
support-to-fragment evidence chain. A disk-backed join retains one row per
fragment across native assignments and multiple graph loci. Mate suffix handling
matches T1K's reader; collisions within an evidence source fail explicitly.
The synthetic integration retains exactly 2,004 fragments: 708 with both sources
and 1,296 with graph-only support. Scores are not combined as independent reads.
The first real native-evidence pilot reuses frozen T1K candidate FASTQs and must
reproduce its genotype table byte-for-byte before its evidence is accepted.
It passed for HG00658 in 156.6 seconds: native-evidence-development-result.json
records identical genotypes and hashes for the instrumented assignment export.

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

The development evidence jobs project the 32 locus/panel/selection mappings onto
observed paths, retain concordant mate pairs and join each eight-locus set with
native T1K fragment assignments. They depend on mapping success and a synthetic
byte-parity check of disk-backed fragment grouping. Grouping uses SQLite so that
arbitrary mate order does not require the entire sample in memory; duplicate and
ambiguous placement semantics are unchanged. Large evidence files are hashed in
chunks. Code snapshots and dependencies are recorded in
development-evidence-launch.json. Native weights and graph scores remain separate
in these joined outputs: calibration and genotype inference are still required.
The disk-grouping parity check passed on IBEX (52045886): the output SHA-256
matches the original 2,004-fragment file exactly (fragment-disk-parity-result.json).

All 84 reserved donors now have verified read preparation records, archived in
validation/READ_PREPARATION.json (65,467,081 paired fragments). This archives the
preparation manifests and read hashes; it does not genotype or unblind them.

## Completed development controls and integration repairs

Both observed-context controls completed all 64 development donors. Original
T1K scored 495/509 at two fields and 163/323 at four fields; HPRC additions scored
389/509 and 165/323; HPRC-plus-Asian additions scored 465/509 and 232/323. The
four-field gain from Asian context is accompanied by a two-field regression and
does not satisfy the success criterion. These are exploratory development results,
not held-out validation. Full-genomic-IPD and genomic-IPD-plus-context controls
are separate experiments needed to distinguish reference completeness from panel
content. They keep the same reads, candidates, caller and scoring rules.

Full-genomic-IPD substitution subsequently completed all 64 donors at 494/509
two-field and 285/323 four-field genotypes, versus original T1K's 495/509 and
163/323. This is the strongest reference control. The net one-call two-field difference
is effectively unchanged in practical terms, not evidence of meaningful deterioration.
It does not meet the literal no-decrease point rule; that rule has not been
changed retrospectively. Paired development uncertainty is reported separately. Its gains cannot be attributed to
the Asian graph. Completed donor/gene results and ancestry summaries are archived
in development/completed-reference-controls-v1/. Reserved outcomes remain unscored.

Three initial fragment projections failed because GBZ import splits long source
GFA nodes: mapping coordinates therefore cannot always index the original GFA.
The repair exports P-line paths directly from the full mapping GBZ using
`vg convert -f -W --vg-algorithm`, verifies every complete path sequence against
the observed reference and normalizes only the previously documented path-name
aliases. Revised outputs use fragments-v2/join-v2; earlier attempts are retained.
The GBZ-repair pilot gates all 32 revised projections. Unit tests cover split-node
coordinates, reverse path sequence and rejection of changed sequences.
The repaired HPRC DRB1 pilot completed with all 63 full path sequences verified,
738,921 fragments retained and 5,396 concordant supported fragments
(projection-repair-pilot-result.json). This validates projection, not genotyping.

Mapping logs separately show failed insert-distribution learning at some loci
when supplied unbinned reads, followed by single-end mapping. The coordinate fix
does not solve this calibration issue. It remains a requirement before the final
graph inference method is frozen. No validation truth has been used.

The six completed reference controls are archived in
`development/completed-reference-controls-v2/`. Full genomic IPD remains strongest;
naively adding observed contexts gives 431/509 and 203/323 (HPRC), or 470/509
and 259/323 (HPRC plus Asian). These findings motivate separate graph evidence
rather than treating each observed context as an independent allele candidate.

Development library calibration now uses one unique fragment once across loci,
requires high alignment identity and mapping quality, and excludes competing-locus
and insert-span ambiguities. It uses no genotype truth. HG00658 yielded 10,036
retained fragments across eight loci, mean 452.45 bp and standard deviation
105.60 bp. Mapping-v3 supplies these same parameters to both graph panels and
selection conditions; read hashes must match calibration. Earlier mapping outputs
remain preserved. This controls failed per-locus insert learning; it is not a
claim of improved genotype accuracy.

Calibrated alignments feed the same sequence-verified projection and paired-fragment
pipeline under `fragments-v3` and `join-v3`. Each completed condition can proceed
independently; any unfinished mapping condition gates its extraction through a
live Slurm dependency. The initial launch waited for complete manifests and found
two Asian DRB1 mappings still running; neither was failed or restarted.

Native assignment export now also accepts completed linear-reference controls.
The genomic-IPD development export verifies reference, reads, completion records
and output hashes, then requires byte-identical genotypes to that control. It
retains native alternatives for subsequent refinement; exported weights and graph
alignment scores are not treated as independent read observations. The original
native-evidence export remains preserved for the original-baseline comparison.

### Graph-assisted recruitment development experiment

Before replacing T1K's inference machinery, test whether graph-supported paired
reads recover useful evidence missed by its candidate-read extraction. Retain the
union of all native genomic-IPD candidate pairs and every pair with a concordant
placement on at least one graph path. Duplicate paths or loci never duplicate a
pair. T1K then aligns and jointly genotypes those selected reads against the same
full genomic IPD reference, with unchanged inference parameters. Graph alignment
scores are not multiplied into T1K weights or interpreted as allele probabilities.
This tests graph-assisted recruitment, not the unfinished joint path-score model
or discovery of alleles absent from IPD.

The HG00658 pilot includes native-only read selection (mandatory byte-identical
genotype parity), all initially recruited reads without graph filtering, HPRC
sampled/unsampled, and HPRC-plus-Asian sampled/unsampled. All graph conditions retain
the native candidate reads; no graph absence can remove an IPD alternative. The
all-read arm distinguishes a graph-specific selection effect from simply bypassing
native candidate extraction. Input, reference, executable and graph provenance are
checked; raw/native paired sequences and qualities must agree. The read-selection
unit tests reject foreign fragments, changed sequences, duplicates and broken
pairs. No held-out outcomes are consumed. Pilot outputs are not confirmatory.

The full-genomic-IPD native assignment export completed with byte-identical
baseline genotypes; see native-genome-evidence-result.json. It remains available
for the separate graph/native evidence refinement work.

The native-only recruitment pilot passed byte-identical genotype parity. It
retained all 35,192 native candidate pairs from 738,921 input pairs. Graph unions
added 17,120 (HPRC sampled), 17,917 (HPRC unsampled), 13,069 (additive sampled),
and 14,259 (additive unsampled) pairs. These counts measure candidate-read recovery,
not genotype accuracy. All four calibrated graph evidence joins are complete;
`calibrated-joined-evidence-result.json` records their provenance.

The all-read control is expanded to all 64 fixed development donors, with the
existing HG00658 pilot retained as index zero. It uses the same code as the pilot,
with a pinned cohort hash, and the native parity job gates submission execution.
Scoring was checked against every one of the 2,048 baseline and genomic-IPD
archived donor/gene/resolution rows. Pending/failed runs retain eligibility and
are marked explicitly; partial cohort totals must not be interpreted as accuracy
contrasts. This control does not consume the reserved validation cohort.

The first all-read cohort array used an incorrect `series/cohort.tsv` path and
all 63 tasks exited before read processing. The corrected retry uses `cohort.tsv`
at the run root; its exact SHA was verified before submission. The initial job
and repair are retained in all-read-control-launch.json. The pilot was unaffected.

### Balanced graph recruitment expansion

All four graph recruitment conditions completed the HG00658 pilot with 8/8
eligible two-field genotypes and 5/5 exact four-field genotypes, matching the
full-genomic-IPD control. Thus the pilot establishes execution and preservation,
not incremental accuracy benefit from graph recruitment or Asian paths.

The next development panel contains the first two entries per ancestry stratum
in the existing execution-cohort order (EAS, SAS, EUR, AFR), yielding eight
independent known families. The existing HG00658 pilot is retained; seven new
donors are processed in graph-development-batch.json. Selection reads metadata
only and does not depend on observed errors. Per donor, the pipeline measures
coverage, exports native genomic-IPD evidence, bootstraps unsampled additive
mapping for insert calibration, maps four graph conditions with the same measured
insert distribution, projects paired evidence, joins competing-locus evidence,
and runs the graph-assisted native T1K call. All native candidate pairs remain.
Dependency arrays prevent downstream work from consuming incomplete upstream
outputs. Reserved validation donors remain excluded and unscored.

### Reserved baseline execution without outcome inspection

The original T1K reference, preset and explicit four-field recipe are frozen in
validation-baseline-plan.json. Baseline typing for all 84 reserved donors can
proceed while development finishes: output genotype files remain remote and are
not fetched or scored before the primary candidate is frozen. The runner consumes
only cohort metadata and verified read-preparation hashes, never truth or assembly
catalogue rows. It rehashes reads under Slurm, checks the frozen wrapper/reference,
and verifies the exact reservation and preparation manifests. An execution pilot
gates the remaining 83 jobs. Monitoring observes job state only. This advances
baseline computation without changing the locked evaluation protocol.

The all-read HG00658 pilot also completed at 8/8 and 5/5, matching every other
genomic-IPD condition. Its complete eight-method snapshot is archived in
`development/completed-graph-recruitment-pilot-v2/`. There is no incremental
recruitment benefit demonstrated on this one donor.

### Endpoint implementation

`paired_evaluation.py` requires the complete fixed donor × eight-locus × two-
resolution table for both candidate and baseline. It rejects pending results,
duplicate/missing rows, altered donor metadata, unequal eligibility and credit
assigned to failed runs. Eligible failed calls stay in the denominator. It reports
paired gains/losses and ancestry strata and resamples complete known families.
Families are sorted by identifier before the seeded draws, making results invariant
to score-table order. It uses 10,000 replicates and seed 20260918. This stable
ordering explains small Monte Carlo differences from the earlier development
interval, which used insertion order; the underlying point estimates are unchanged.

The module reports the existing protocol gates literally, without interpreting a
one-call point difference as meaningful deterioration. It cannot by itself certify
independence: candidate freeze, untouched outcomes, provenance, exclusions and
final failure audits remain separate requirements. Its test on the completed
reference controls is explicitly development-only and does not unblind validation.

### Scheduling correction

The development graph stages initially inherited the pilot's debug partition.
A live scheduler audit found debug restricted to three nodes / 204 CPUs, while
batch exposed 367 nodes / 32,000 CPUs, with 12,191 idle at observation time. Pending
stages and pending bootstrap tasks were therefore moved to batch without restarting
running tasks or changing inputs/inference. Future cohort launchers default to
batch. `graph-batch-partition-update.json` records scheduler update results.

The reserved baseline execution pilot completed successfully in about 158 seconds;
only its completion metadata was fetched, not genotype contents. Its execution
manifest is validation-baseline-pilot-execution.json. No reserved outcome was scored.

### Completed all-read control and residual-error direction

All 64 all-read controls completed at 493/509 two-field and 285/323 four-field
correct genotypes, versus 494/509 and 285/323 for genomic IPD with native candidate
extraction. The four-field contrast comprises one gain and one loss, with paired
95% interval approximately −0.92 to +0.92 percentage points. This does not support
bypassing read extraction as the source of further four-field gains. Complete
results are in development/completed-all-read-control-v1/.

An audit of all 38 residual genomic-IPD four-field errors found 19 at DQA1;
33 already had correct two-field calls. Both truth haplotypes have exact genomic
labels in the additive training graph for 32 errors, versus 18 in HPRC alone.
This is label availability, not evidence that short reads distinguish the paths.
The audit reads exposed development truth only and does not filter the independent
validation cohort or redefine its endpoint. See development/residual-error-audit/.

`graph_pair_refinement.py` implements an exploratory refinement core, not yet a
validated caller. It retains the native two-field pair, assigns each qualifying
fragment to at most one locus, takes maximum support across redundant contexts,
and considers diploid allele pairs. Native alleles absent from the graph retain
IPD fallback; ties or an unlabelled winning path also retain the native call.
Refinement requires at least 20 fragments with candidate-dependent emissions and score gaps of at
least 10 against both the native and next-best pair, using noise 0.01 and score
temperature 10. These are development thresholds, not calibrated probabilities.
Unlabelled and incompatible context paths must remain nuisance alternatives in
the real-data adapter, rather than being discarded before normalization. The core
passes synthetic tests for false-heterozygosity correction, fallback, ambiguity,
coarse-label preservation and competing loci. Real-data integration and evaluation
remain necessary before any candidate freeze or improvement claim.

### Real-data graph pair adapter

The real-data adapter now consumes all-locus joined paired evidence, native
full-genomic-IPD calls, verified observed reference labels and the same measured
insert calibration used for mapping. It checks every provenance chain, exact path
sequence/label hashes and graph-panel identity. Every path remains represented:
unknown or incompatible labels become nuisance alternatives. Only genomic labels
can support four fields; CDS identity cannot create a four-field call. Duplicate
contexts contribute through maximum support, not summed path counts. Fragment IDs
must be sorted and unique, and native two-field calls are checked unchanged.

Four HG00658 development refinement jobs are submitted, with the additive
unsampled pilot gating the other conditions. Calls and per-gene decisions are
scored separately from native T1K tables; no artificial T1K abundance or quality
columns are generated. The development scorer checks unchanged native coarse calls
and requires every changed gene to have a recorded refinement decision. Independent
validation remains unscored. Real-data performance is pending; passing adapter
unit tests is not an improvement claim.

The real-data pilot completed all four conditions without changing any native
call (8/8 two-field, 5/5 four-field). Complete decisions and scores are archived in
`development/completed-graph-pair-pilot-v1/`. The same code and parameters are
submitted for the fixed seven additional development donors, after their evidence
joins. This remains characterization; no primary candidate has been frozen.
The current informative-fragment counter measures variation across all candidate
columns, not specifically between the native and proposed pair. It is not a count
of independently diagnostic variant observations or calibrated confidence.

`audit_graph_refinement.py` revalidates the development outputs and joins each
four-field decision to its genomic-IPD baseline outcome. It reports rescues,
losses, changed calls and reasons for remaining errors for each graph condition,
with completed-donor counts for partial batches. The monitor refreshes this audit
after fetching refinement calls. It neither reads reserved outcomes nor selects
parameters from them.

The residual-error audit additionally checks eligibility for the current native
fallback rule using the same genomic alias resolver as inference. Of 38 errors,
12 have an unresolved native four-field pair. Another 13 contain a native allele
absent from the additive graph (15 for HPRC alone). Thus only 13 additive-panel
errors, or 11 HPRC-panel errors, can currently enter pair refinement. Eligibility
does not establish that either truth pair is represented or distinguishable.
This identifies a coverage limitation: improving the remaining cases would
require comparable evidence for missing IPD candidates or explicit ambiguity
handling, rather than simply removing the fallback guard.

### Scheduling follow-up (18 September)

The eight-donor pair-refinement batch subsequently completed all 32 graph
conditions. HPRC sampled and unsampled retain 41/47 correct four-field calls;
both additive conditions have 40/47, versus genomic IPD 41/47 and original T1K
21/47. Each additive condition introduces one error, with no rescues: HG02132
DPB1 under sampling and HG01530 A without sampling. These small development
differences provide no evidence of incremental graph benefit. Complete results
and decisions are in `development/completed-graph-pair-batch-v1/`.

`diagnose_refinement_changes.py` reproduces the recorded changed-pair score gains
and partitions contributions by gene-span overlap and native-candidate support.
It is a read-only development diagnostic, not a revised caller. The competing
paths at each changed allele have equal lengths, so a simple path-length
difference is not sufficient to explain these errors. Fragment-level causes
remain under investigation; reserved validation outcomes remain unexamined.

After 138 of 224 calibrated evidence tasks completed, the remaining 86 were
waiting on batch priority while debug had idle capacity. Pending evidence,
join, recruitment and pair-refinement jobs were permitted on either batch or
debug. All 86 remaining evidence tasks were subsequently observed running
concurrently on debug. No inference inputs, parameters or dependencies changed,
and no running task was restarted. The update audit is
`graph-batch-dual-partition-update.json`; later completion must still be verified
from the individual output manifests.

### Gene-internal refinement development comparison

The completed recruitment comparison is archived in
`development/completed-graph-recruitment-batch-v1/`. All eight donors completed
all conditions; graph recruitment scored 40–41/47 four-field genotypes versus
41/47 for genomic IPD. This provides no incremental gain on this development set.

The next refinement version retains the same pair scoring and thresholds, but
admits a fragment only when both mates lie wholly inside the annotated exon
envelope (including introns) at every placement at the selected gene. A boundary
placement rejects the whole fragment rather than deleting only some candidate
placements. Unfiltered competing-locus evidence still determines locus
assignment. This deliberately excludes terminal UTRs as well as flanks; extending
the interval would require explicit annotation. It is a development restriction,
not a calibrated diagnostic-variant likelihood or a validation candidate freeze.

The same eight donors and four conditions run under a new code snapshot and
output directory. Version-one outputs remain intact. Tests cover a flanking mate,
alternative boundary placement, and non-destructive filtering; all 112 tests pass.

All 32 gene-internal comparisons completed. The unsampled additive graph scores
43/47 four-field genotypes, correcting DQA1 in HG03516 and NA19185 without losing
another correct genotype. HPRC remains 41/47; sampled additive remains 40/47.
The unsampled HLA-A error is removed, while the sampled DPB1 error persists.
`development/completed-graph-internal-refinement-v1/` archives the checked
results. This is a development gain, not independent proof. The global
informative-fragment counter still needs a pair-specific discrimination check
before freezing the validation candidate. No reserved outcome has been inspected.

### Candidate selection and pair-specific support

The pair-specific development check requires at least 20 fragments whose pair
emissions differ between the proposed and native genotype, and between the
proposed genotype and its closest competitor. When closest competitors tie, use
the smallest such count. This prevents unrelated nuisance-candidate variation
from inflating the count; it does not turn emissions into independent diagnostic
variant observations. A regression test reproduces that inflation, and all 113
tests pass.

All 32 comparisons completed with unchanged genotype outcomes. The two
unsampled-additive DQA1 rescues retain 370 and 385 fragments distinguishing their
closest competitor. Results are archived in
`development/completed-graph-pairwise-refinement-v1/`.
`validation/PRIMARY_CANDIDATE.json` selects this inference version, records its
hashes and parameters, and retains the original evaluation gates. All 3,515
observed-reference source occurrences were checked against the immutable
400-haplotype exclusion list, with zero overlap. The complete validation execution
driver still needs to be frozen before any outcome inspection.

The reserved genomic-IPD control is submitted behind a pilot gate with pinned
read, cohort, reservation, reference and wrapper hashes. It decodes predictions
on IBEX for subsequent inference but does not load truth or score outcomes.
Monitoring observes its job states only; predictions remain remote. This is a
control and prerequisite, not the final graph validation run.

### Frozen independent graph execution

`validation/FROZEN_VALIDATION.json` fixes the inference workflow before outcome
inspection: 84 reserved donors, original T1K and genomic-IPD controls, unsampled
HPRC and HPRC-plus-Asian graphs, gene-internal paired evidence, pair-specific
20-fragment support and the existing score thresholds. It pins 18 code files,
input metadata and 20 reference/index manifests. Each stage checks the freeze
digest, code, assets and the donor's actual read hashes. The reserved genomic-IPD
control is a prerequisite for native evidence and final refinement.

The validation projection adapter calls the same explicit-path implementation
as development; no mapping or inference parameters change. Stage-coordinate
tests cover all 84 donors, eight loci and both panels exactly once. All 114 tests
pass. A complete one-donor pilot gates the remaining cohort, with arrays capped
at 300 simultaneous tasks and both batch/debug partitions available. Outcome
predictions stay remote and unscored; the monitor reads scheduler state only for
the reserved jobs. The primary four-field and literal two-field gates are
unchanged, and graph-versus-IPD and additive-versus-HPRC contrasts stay separate.

IBEX rejected advance submission of the 1,328-task cohort mapping array with
`QOSMaxSubmitJobPerUserLimit`. No mapping job was created by that request. The
pilot and cohort preparation jobs were retained. `advance_validation_graph.py`
now submits mapping, evidence, joins and refinement only after every predecessor
task has completed successfully, keeping one large remaining stage queued at a
time. The persistent monitor calls this dispatcher; incomplete or failed task
sets do not release their successors. An uncertain submission is recorded for
reconciliation rather than blindly retried. This scheduling adjustment changes
neither the frozen inference nor the endpoint. The suite now has 115 passing tests.

### Frozen endpoint implementation

`validation/EVALUATION_FREEZE.json` pins the endpoint code, truth catalogue and
nomenclature inputs before outcome inspection. `evaluate_reserved.py` accepts an
explicit, complete four-method result grid; importing it does not read reserved
calls or truth. It preserves failed donors as zero-credit eligible results,
rejects missing/pending results, reports ineligible truth separately, and runs
the existing family-paired bootstrap with 10,000 replicates and seed 20260918.
Only the original-T1K contrast carries the primary protocol decision; genomic-IPD
and HPRC contrasts are explicitly exploratory. Tests exercise failed-donor
denominators and incomplete-grid rejection. All 117 tests pass. Provenance
verification and the final snapshot handoff remain necessary before unblinding.

## References

- Song et al. (2023), https://doi.org/10.1101/gr.277585.122
- Siren et al. (2024), https://doi.org/10.1038/s41592-024-02407-2
- Source/paper review: ../research/t1k-review/REVIEW.md
