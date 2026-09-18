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

Remaining implementation: graph construction and topology validation, competing
paralog/decoy references, sample-specific k-mer selection, T1K reference adapter,
graph-based candidate-pair refinement, and frozen independent evaluation.
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
