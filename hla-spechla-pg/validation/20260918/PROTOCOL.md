# DōgoHLA 0.1.0: prospective 32-donor extension

Frozen before new inference on 18 September 2026. Deadline: **09:02:07 UTC**
(12:02:07 Riyadh), five hours after the request. Cohort and method digests are
in `frozen-method.json`; `cohort.tsv` lists every planned donor.

## Design

Evaluate all 32 remaining fold 1–4 donors: 16 East Asian and 16 South Asian,
32 distinct families, none overlapping the eight-donor method-development
pilot. These are an extension of the existing development cohort, **not** the
untouched 228/946 locked benchmark. Use existing complete 1000 Genomes MHC
read recruitments, identical input FASTQs for each paired comparison.
Exclude the entire held-out fold and related family haplotypes from panel
references and graphs. Verify graph path sequences and exclusions before use.

Run native SpecHLA 1.0.12 with its default small-variant workflow, DōgoHLA
0.1.0, and DōgoHLA without graph reconstruction (same repaired phasing,
IPD 3.65 phase database and DRB1 naming correction). The last comparison
isolates the structural graph contribution. Test native SpecHLA's optional
long-indel mode on a pilot donor first; if executable within the time budget,
run it as an additional comparator. Report availability/failures explicitly.

## Outcomes and analysis

Primary: exact whole-gene haplotypes against assembly-derived truth across
8 genes × 2 haplotypes × 32 donors (512 planned sequences), using global
edit distance and minimum-distance unordered diploid assignment. Report
truth availability explicitly. Secondary: total global edits, masked Ns,
call/completion rates, assembly-derived two-field/G-group calls, and the
available independent experimental Gourraud genotypes. Infix scores are
secondary and must never be mixed with global scores.

Report every donor, gene and ancestry stratum, not only improvements. Missing
outputs and no-calls receive no exact/correct credit in planned denominators;
edit-distance comparisons require paired valid outputs and are labelled
conditional on completion. Missing truth is separately accounted for and
excluded identically for every method. Report paired donor differences and
95% percentile bootstrap intervals (10,000 resamples, seed 20260918), sampling
whole donors/families rather than treating correlated alleles as independent.
Show the graph ablation even if it is neutral or worse.

## Execution and stopping rules

Freeze inference code and thresholds before new results. Isolate outputs and
record code, reference and input hashes, Slurm IDs, timing, exit status and
validated completion markers. Do not silently reuse partial runs. Execution
repairs must be logged; changing inference after viewing results creates a
new exploratory version and cannot retain the frozen primary label.

Use Slurm on DDBJ, four CPUs per donor, and reuse existing read/reference
assets. Build the four missing normalised DRB1 graphs from cached raw graphs.
At the deadline, report all completed, failed and pending donors. Do not claim
32-donor validation if fewer completed. No outcome-dependent cohort extension
or early stopping. Native long-indel feasibility is an execution gate, not an
accuracy-based selection. This short validation can support further study;
it cannot establish general superiority or clinical validity.
