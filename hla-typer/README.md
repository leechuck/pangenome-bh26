# hla-typer: short-read HLA typing with the Asian MHC pangenome

Plan: `~/.claude/plans/ultrathink-about-this-make-melodic-moth.md` (panel supplies the backbone,
IPD-IMGT/HLA supplies the exons). Benchmark code lives in `../hla-bench/`. All compute runs on DDBJ
under `/home/leechuck/hla/hla-typer` (partition `asianhla-c32`); `driver/remote.py` uploads and runs.
Status: development (40 hla-asian50 pilot donors only). Nothing here is a test-set result.

## Stages and verified facts

| Step | Code | Output (cluster) | Verified |
|---|---|---|---|
| Graph scan | `loci/scan_graph.py` | `results/graph_scan/ref_nodes.tsv.gz` | 754 haplotypes, reference path 5,169,788 bp, 1.5 min |
| Locus boundaries | `loci/define_loci.py` | `source/loci_spec.tsv` | boundary nodes forward, traversed once by >= 99.6% of haplotypes |
| Locus cut | `loci/extract_loci.py` | `results/loci/<LOCUS>.fa.gz` | 751-754 haplotypes per locus, all ACGT-only; every annotated gene copy inside its locus (DRB1: 139 haplotypes lose <= 20 bp of 3' UTR at the shared DRB345/DRB1 boundary) |
| Labels | `panel/make_labels.py` | `source/haplotype_labels.tsv.gz` | from `hla-analysis/results/sequence_catalogue.tsv` (IPD 3.65.0) |
| Panel arms | `panel/build_panel.py` | `panels/<arm>/` | full 754; hprc, asian_matched, random_matched 466 haplotypes each (units = donor x assembly source) |
| Weights, CDS | `panel/make_weights.py` | `weights/*.bed`, `weights/cds_intervals.tsv` | CDS cut from locus coordinates equals catalogue `cds_sha256` for all 2,398 A/DQ/DRB345 copies checked |
| Recruitment | `stageA/recruit.sbatch` | `reads/<donor>/` | chr6 MHC +/-100 kb, chr6 alts, 525 `HLA-*` contigs, unmapped pairs; no mate fetching (1.5% of region reads have mates elsewhere on NA18620) |
| Stage B | `stageB/setup_db.sbatch`, `stageB/genotype.sbatch`, `stageB/score_calls.py` | `db/<arm>`, `calls/<arm>/<donor>` | `augment` takes seconds to minutes per locus (HLA-A 11 s, DQ 174 s); leave-one-out via `locityper genotype --leave-out` |
| Stage C1 | `stageC1/ipd_cds.py`, `graft.py`, `build_candidates.py`, `rescore.sbatch` | `c1/<arm>-<variant>/<donor>` | graft round-trip tests; candidate CDS self-check in the builder |

Loci (ctg1 = chr6 - 28,410,700): HLA-A 18.9 kb, HLA-C 20.6 kb, HLA-B 23.6 kb, HLA-DRB345 structural
block 93.2 kb in GRCh38 (6.7-176 kb across haplotypes), HLA-DRB1 21.2 kb, HLA-DQ 47.9 kb, HLA-DP 43.0 kb.
The DR block between DRA and the DRB1 3' end has no boundary node shared by all haplotypes, so
DRB3/4/5 are typed jointly as one structural locus.

## Development results (40 hla-asian50 pilot donors, leave-one-donor-out, own-assembly truth)

Stage B without exon weights, `results/stageB/summary.tsv`:

| arm | two-field | G-group | exact CDS | exact full gene |
|---|---|---|---|---|
| full (754) | 406/432 = 94.0% | 94.2% | 93.1% | 70.7% |
| asian_matched (466) | 94.7% | 95.4% | 93.3% | 67.3% |
| hprc (466) | 91.4% | 91.9% | 90.3% | 65.7% |

Classical 8 genes, two-field: full 295/319 (92.5%), asian_matched 92.8%, hprc 88.4%; T1K 1.0.6 on the
same reads with the pinned IPD 3.65.0 database 306/311 (98.4%). DRB3/4/5 with copy number: panel 111/113,
T1K 48/110 (T1K reports one allele per paralog, no copy number; both scoring conventions must be reported).
Disagreements on the classical genes: both right 283, only T1K right 23, only panel right 5. Of the panel's
30 CDS-level errors, 18 are wrong picks with the truth in the panel (8 at DRB1) and 12 have no panel
haplotype carrying the truth CDS; in most misses the correct pair is not among Locityper's retained options.
The recurring DRB1 pick `HIFI032302D.1` is a protein-exact DRB1*08:03 haplotype chosen instead of a second
copy in homozygous donors.

Variants queued from these findings: `calls/<arm>-w` (Stage B with `--reg-weights`; the donor's own rows
are filtered from the BED files because Locityper maps left-out names and their identical substitutes to one
haplotype), `c1/full-ens` (grafting round whose backbones also include leave-one-out panel haplotypes labelled
with T1K's alleles), `c1/full-ens-p` (plus partial/exon-level IPD grafts, for DRB1 where IPD has a complete
CDS for 896 of 4,023 alleles). Candidate sets are capped at 500 per locus (priority: backbones, T1K grafts,
IPD grafts by edit distance, at most 2 per backbone and two-field type); uncapped they reached 16,370 at HLA-A.
Of the 12 panel-absent truth alleles, 11 exist as complete IPD CDS, 9 within 2 edits of a panel backbone and
two DPB1 alleles at 4-5 edits (`--max-edits 6` needed for the no-T1K ablation); one (NA18620 HLA-C) is novel.

## Known issues
- Pair fetch (`samtools view -P`) over S3 read > 5.8 GB and ran > 30 min per donor, so it was dropped.
- IPD 3.65.0 has a complete CDS for only 896 of 4,023 DRB1 alleles; the rest cannot be grafted whole.
- `locityper target` collapses identical haplotypes; checked on HG02155 that `--leave-out` keeps a shared
  sequence under another donor's name (log: `replaced 1 haplotype(s) with identical [HG02155.1->HIFI032097D.2]`).

## Tests
`python3 -m unittest` in `panel/` and `stageC1/` (the latter needs `edlib`; cluster env `hlatyper`).
