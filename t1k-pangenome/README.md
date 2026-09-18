# T1K with personalized Asian-enriched graph references

Development preparation; no demonstrated improvement over T1K yet.
The evaluation conditions and controls are in [PROTOCOL.md](PROTOCOL.md).

`build_panel.py` builds observed HLA sequences for HPRC-only and HPRC-plus-Asian
references. It verifies the validation reservation, removes excluded donor/family
paths, validates sequences against their gene/CDS hashes and deduplicates identical
sequences. JSON sidecars retain the allowed source paths, coding coordinates,
unknown-base counts and separate exact-genomic versus exact-CDS labels.
CDS matches are never promoted into four-field labels. No frequency prior is
inferred from duplicated path counts.

```sh
python3 t1k-pangenome/build_panel.py \
  --source hla-analysis/source \
  --metadata hla-typer/source/haplotype_donors.tsv \
  --catalogue hla-analysis/results/sequence_catalogue.tsv \
  --exons hla-spechla-pg/source/exons.tsv \
  --reservation t1k-pangenome/validation \
  --output t1k-pangenome/references/NEW_BUILD
python3 -m unittest discover -s t1k-pangenome -p 'test_*.py'
```

The builder refuses an existing output directory. `COMPLETE.json` records builder,
input and output hashes only after successful preparation. Downloaded/generated
reference FASTAs are local artifacts, excluded from Git; `reference-build.json`
records the latest verified build.

`audit_reads.py` resolves reserved donors through the public 1000 Genomes CRAM
manifest and checks both CRAM and CRAI objects with HTTP HEAD requests. It verifies
the reserved cohort hash and rejects conflicting donor/source mappings. The
2026-09-18 audit in `validation/read-availability.json` found both objects
accessible for all 84 donors. This establishes availability, not successful
decoding; recruitment must still validate the index and extracted paired reads.
The original reservation is unchanged and no outcome files enter this check.

```sh
python3 t1k-pangenome/audit_reads.py \
  --reservation t1k-pangenome/validation \
  --public-manifest hla-targeted/source/C4Investigator/resources/1000Genomes_resources/tgp_full_30x.tsv \
  --output t1k-pangenome/validation/read-availability.json
```

An audit refuses to overwrite an existing report; use a new output name to
record a later availability check.

Current references preserve 86 HPRC/reference source haplotypes per locus and
351–354 source haplotypes in the additive set (some assemblies lack a locus).
These counts are smaller than the old benchmark panels because the new reserved
test families and previously exposed development families are all excluded.

All 16 panel/locus graphs completed construction and sampling-index preparation;
four needed an explicit backbone when building the distance index. Remaining implementation:
competing paralog/decoy references, real-sample k-mer selection, graph-based
candidate-pair refinement, and frozen independent evaluation.
The existing eight target loci alone are not a complete mapping decoy set.

The upstream [vg haplotype sampling documentation](https://github.com/vgteam/vg/wiki/Haplotype-Sampling)
requires a compatible chain structure; do not silently replace a failed sampling
step with the full graph and label it personalized.

`build_graph.py` prepares one locus graph under Slurm using the installed
wfmash/seqwish/smoothxg/gfaffix/odgi/vg tools. Its pilot parameters use 500 bp
segments and a 1 kb minimum alignment block; the older DRB1 builder's 25 kb
cutoff is longer than the new HLA-A sequences. These are development parameters,
not an accuracy-validated choice. The builder rejects changed reference hashes,
preserves failed attempts, records executed commands and binary hashes, and
requires exact preservation of every named input sequence through graph
normalization. It then builds distance, r-index and haplotype-sampling indexes.
A failed topology/index step is a failed build, with no full-graph fallback.
Successful preprocessing does not yet prove sample-specific selection or typing
accuracy. The graph remains a locus preparation artifact without a complete
paralog/decoy reference.

`smoke_personalization.py` tests the integrated vg k-mer counting, haplotype
sampling and paired-read mapping workflow using deterministic synthetic reads
from two training paths. It copies graph artifacts to an isolated directory,
records source/tool/read hashes, and requires one alignment record per read and
at least 95% mapped reads. This gate checks execution and output completeness;
it does not test correct placements, HLA genotype accuracy or held-out
generalization. Validation donors and their truth are not consumed.

The installed vg 1.76.1 integrated workflow still executes external `kmc`.
KMC 3.2.4 is staged separately under `tools/kmc-3.2.4`; its upstream archive and
binary hashes are in `kmc-provenance.json`. The smoke driver sets
`OMP_NUM_THREADS` to its allocated thread count, overriding the benchmark
container's one-thread default for this subprocess only.

`repair_graph_index.py` handles the diagnosed top-level-loop failure by rebuilding
the distance index with `vg index -P` and the original deterministic indexing
backbone. It copies the graph unchanged, verifies the vg binary, preserves the
failed attempt, and records the old manifest hash. Successful repairs do not
replace failed files in place.

The synthetic smoke supports explicit `--kmer-coverage`. For the generated
40-fold coverage per haplotype and 150 bp reads, the expected shared 29-mer
coverage is approximately `2 * 40 * (150 - 29 + 1) / 150 = 65.07`. The v3
experiment uses 65 to test the explicit-coverage workflow. This number is known
from simulation and must not become a fixed assumption for real data. Real-data
coverage estimation and selection calibration remain unresolved.

`t1k_reference.py` prepares the linear sequence-information control. The complete
IPD reference is copied byte-for-byte before adding observed full genomic
contexts. Each added sequence has a short internal `GENE*PG...` identifier and
a sidecar retaining genomic/CDS label alternatives, original graph path and
source metadata. These identifiers are not official HLA allele names. The
resolver never promotes a CDS-only match to four fields and never pads short
labels. Missing coding annotations, unknown bases and exact duplicates already
in IPD are recorded as skipped additions; the IPD fallback remains present.

The adapter uses the available coding intervals as T1K coverage annotations;
these are CDS intervals, not complete transcript exon/UTR annotations. T1K
treats different contexts as separate candidates, and reference multiplicity
can affect its weighting and tie-breaking. Consequently this is an explicit
linear-reference control, not the final allele-level graph inference model.
Synthetic and development testing must precede any held-out evaluation.

`decode_t1k.py` translates genotype tables into explicit two- and four-field
alternative sets using the context sidecar. It retains an unresolved flag when
any candidate lacks the requested resolution, preserves homozygous copy
assignments, and rejects cross-gene labels. These JSON calls require the same
conservative ambiguity handling during evaluation; a matching member of an
unresolved set is not a resolved correct call.

`path_support.py` connects graph alignments to observed path candidates. It
projects each contiguous alignment walk in either orientation onto all
compatible named paths, retaining repeated placements, offsets, mapping scores
and mate metadata. Shared sequence does not force a single allele assignment;
walks absent from the observed paths stay explicitly unsupported. The adapter
checks that alignment node IDs belong to the input graph. Fragment likelihoods,
paralog competition and allele-level inference still need integration.

The corrected synthetic T1K controls completed successfully: baseline, HPRC
contexts and HPRC-plus-Asian contexts all decoded to the same expected HLA-A
four-field pair. The first adapter version incorrectly mixed `A` and `HLA-A`
gene namespaces; v2 derives its namespace from the unchanged IPD reference and
the decoder rejects duplicate normalized gene rows. These training-path tests
establish integration, not improvement. The first graph-evidence projection
also retained compatible observed-path placements for all 4,008 reads.

`fragment_support.py` joins mate-linked graph placements into one evidence record
per fragment. It retains all concordant shared paths, keeps only the best score
for duplicate placements, and explicitly records missing mates and unsupported
pairs. The insert-size limit is configurable (1,000 bp for the synthetic smoke),
not a learned library distribution. All 2,004 synthetic fragments had concordant
placements; this remains an integration check rather than genotype validation.

`prepare_validation_reads.py` prepares the 84 reserved donors in an isolated
`validation-reads/` directory. The input list is bound to the immutable cohort
and public-source audit. Recruitment uses the existing benchmark's region and
mate policy, followed by CRAM sample-ID, complete FASTQ syntax/pair-name/count,
and output checksum verification. It does not load genotype truth or run a
typer. A successful pilot gates the remaining array, capped at eight concurrent
four-CPU tasks. The existing monitor records preparation failures.
