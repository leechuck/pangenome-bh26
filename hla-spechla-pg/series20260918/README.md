# DogoHLA: matched ancestry / reference-panel experiment

Work directory: `/ibex/scratch/projects/c2014/rob/dogohla-benchmark/`.
Slurm partition `batch`, account `pi-hohndor` (verified association; `c2014`
is the storage project). Current job IDs are in `ibex-launch.json` and
`ibex-recruit-job.txt`. `monitor_ibex.py` records scheduler state every minute;
the local service is `dogohla-ibex-monitor-20260918`.

## Matched 64-donor cohort

Retain the prospectively selected 32 Asian donors (16 EAS, 16 SAS), and select
16 EUR and 16 AFR donors using metadata and truth availability only. Selection
uses unique known families and population round robin; no prediction outcomes.
`selection.json` pins the cohort hash. The new non-Asian donors are now reserved
for this experiment and cannot be called an untouched test in later development.

Exclude all 64 donors and their known family groups from every panel-derived
reference, including alternate assembly aliases. There are 132 excluded
haplotypes. `select_cohort.py` checks canonical donor/source metadata explicitly.
The same exclusions apply to binning, phasing panel candidates and graph paths.

Three panel arms share the same reads, software, IPD release and calling settings:

- `full`: the expanded AsianHLA panel, including its HPRC component; 622 retained
  haplotypes including the two reference haplotypes.
- `hprc`: HPRC plus the two references; 340 retained haplotypes.
- `asian_matched`: the non-HPRC Asian-enriched component plus a deterministic HPRC
  subset, with exactly the same 340-haplotype size as `hprc`.

The DRB1 graphs add the same SpecHLA coordinate-reference path to each panel:
623, 341 and 341 paths respectively. This prototype uses graph structural calls
at DRB1; panel read recruitment affects all eight genes. Each DogoHLA arm also
retains its no-graph output to measure the incremental structural-graph effect.
The panel comparison changes recruitment references as well as the DRB1 graph;
it must not be presented as a graph-mapping-only ablation.

## Matched database and comparators

The separate shadow SpecHLA environment uses IPD 3.65.0, the version already
used by T1K. Genomic/CDS/ARS candidate and binning databases are regenerated;
original coordinate references are retained. Both updated SpecHLA and DogoHLA
use `nonuse`, disabling the legacy frequency filter so newly added alleles are
not excluded by the older frequency table. Label this baseline
`SpecHLA-IPD365-noFreq`; the frozen v0.1.0 baseline used IPD 3.38 and `Unknown`.
Do not attribute differences between those configurations solely to the graph.

T1K is rerun with `--alleleDigitUnits 4 --alleleDelimiter :`; score its actual
`*_genotype.tsv`, not representative alleles from `*_allele.tsv`. The pinned
reference contains four-field names. Unresolved calls remain unresolved.

Assembly sequence and strict four-field truth remain separate from Gourraud's
experimental, ambiguity-aware exon typing. `abi_rached2018/` supplies a secondary
computational concordance reference, not an additional independent truth cohort.
The expanded Gourraud cohort contains 946 donors in 660 known family groups:
230 EAS, 355 EUR, 188 AFR and 173 AMR. Canonical graph-donor and family overlap
are both zero. Recruitment and four-field T1K arrays are submitted; the persistent
`dogohla-gourraud-scheduler-20260918` service submits SpecHLA and three DogoHLA
panel arms in packed arrays after the corresponding matched-cohort check passes. It keeps
total submitted/running user tasks below 1,900 (the verified QoS limit is 2,000).
Job IDs appear in the ledger as each array is submitted.

The wall-clock deadline was withdrawn on 18 September. The deadline publisher
timer is stopped; monitoring continues until completion. IBEX enforces 1,300
CPUs per user through its partition QoS, and backfill examines at most 64 jobs
per user per cycle. Pending work is therefore regrouped into eight independent
four-thread samples per 32-CPU/64-GB allocation. This retains per-sample methods
and thread counts while allowing about 41 allocations to fill the CPU allowance.
Observed on 18 September at 10:18 Riyadh: 40 packed T1K allocations plus
four method checks used 1,296 CPUs; further jobs waited on the CPU QoS limit.
Running tasks are preserved; only held pending tasks are replaced. `packed/`
records exact donor/index membership and cohort hashes. Experimental truth is
two-field; four-field predictions do not create four-field experimental truth.
Family groups must be resampled together in uncertainty estimates.
HLA-HD still requires an authorized installation/download package.

## Execution and recovery

The original 16 GB workstation-mediated transfer was stopped after measuring
only tens of KB/s. Pinned Conda packages and public reference data are downloaded
directly on IBEX. The 103 MB custom bundle was transferred in checked chunks;
`custom-transfer-complete.json` records its hash locally. Ubuntu 24 runs in
Singularity with an explicit bind preserving the original installation prefixes.
The runtime now excludes the host home and inherited environment and uses the
pinned Python explicitly. This avoids the host user's pyenv Python 3.7 shim.

The first installer attempt lacked CA certificates. The retry binds the host
certificate bundle and preserves TLS verification. Database preparation adapts
IPD's legacy six-field EMBL ID headers to Biopython's seven-field layout; only
headers change, not sequences/features. Allele names come from feature
annotations (or the record description), not the source feature.
The initial PGGB package export described a different environment; the corrected
explicit lock was exported from the actual DDBJ `pggb053/conda-meta` records.
Graph generation uses the original wfmash 0.10.5 command semantics.

Sequence-free deleted-allele records in IPD are excluded from binning FASTA files:
they otherwise generate invalid zero-length SAM headers. `repair_empty_refs.py`
preserves the original FASTA files and hashes, rebuilds all four binning indexes,
and updates preparation manifests. No sequence-bearing candidate is removed.
Dependent inference waits for successful repair and a complete sample from its method.
A subsequent smoke run exposed a read-only container working directory; the
wrapper now creates a writable, separate directory per Slurm allocation. Failed
attempts remain under `attempts/readonly-cwd/`, and only pending tasks are regrouped.

Three execution fixes are also isolated for DDBJ recovery:

1. SpecHap's legacy `-N` insertion command crashes because its current parser
   requires `--protocols`. Translate it to `--protocols nanopore --weights 1`,
   retaining its intended protocol. The failed HG00658 input produced a valid
   445-record VCF in Slurm test 20667505 after this change.
2. Native annotation sometimes needs read depth to distinguish candidate alleles.
   Structural-overlay outputs now link the source BAM and index, checking their
   hashes against the source manifest before annotation.
3. Local assembly tests file non-emptiness for read and selected-contig lists.
   Empty lists now use the existing upstream fallback instead of attempting to
   index a missing assembly or reuse a stale BWA index. The same guard applies
   to every updated IBEX method through the shared shadow environment.

DDBJ failed attempts are retained under `attempts/20260918-execution-fixes/`.
Already verified successful outputs are preserved. `ddbj-recovery-jobs.json`
records the recovery arrays; the frozen primary inference thresholds and truth
remain unchanged. Initial failed IBEX job IDs remain in the launch ledger with
`_initial` suffixes. Scheduler completion alone is insufficient: result manifests
and expected output files must also validate.

## Missing-make recovery (18 September, 11:30 Riyadh)

The first complete-path execution checks failed because the minimal container
lacked GNU make, used by FermiKit local assembly. Strict nested-log checks
rejected these outputs; they never entered accuracy estimates. The exact DDBJ
`/usr/bin/make` binary was copied with its SHA-256 recorded in
`portable-binaries.json`; Slurm preflight 52039301 verified it and the assembly
tools inside the container. `recover_missing_make.py` preserves the failed
attempts, reuses only completed DogoHLA read bins, and restarts each method gate.
Dependent arrays are retargeted to the new gates.

The 946-sample T1K run produced 945 validated completions and one analyzer
segmentation fault (NA10846); the failed attempt is preserved and the same
four-thread, four-field configuration is retried separately. No failed sample
is removed from scoring denominators.
