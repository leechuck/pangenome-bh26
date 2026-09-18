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
