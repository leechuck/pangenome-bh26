# pangenome-bh26

Work from the BioHackathon 2026 (Japan) on human pangenomes and the HLA region,
run on the NIG supercomputer BioHackathon node.

- `hla/` CWL workflow: extract the MHC from haplotype assemblies (pgr-tk),
  call HLA/KIR/C4 alleles (Immuannot), fetch every HLA gene across haplotypes
  (pgr-tk) and visualise them (pgr-tk principal bundles, pggb + odgi).
  See `hla/README.md`.
- `data/` provenance notes for the pangenome datasets uploaded to
  `/home/asianhla/data/upload/` on the NIG node (Arab Pangenome Reference,
  Korean pangenome K-PanRef).
- `1000g_ground_truth/` experimental HLA gold-standard labels for validating
  `hla/` pipeline calls: Gourraud et al. 2014 Sanger typing of 1,267 1000
  Genomes samples, and Lai et al. 2024 4-field HPRC labels for 44 samples.
  See `1000g_ground_truth/README.md`.

- **[DōgoHLA](hla-spechla-pg/README.md)**: an AsianHLA pangenome extension of
  SpecHLA, named for Matsuyama’s Dōgo Onsen. Repaired phasing, guarded graph
  reconstruction, reproducible development results and 1000 Genomes validation.

## Latest results: six-step route (2026-11-17)

Plain-language summary with re-runnable checks:
[`scientist/evidence/route-20261114/README.md`](scientist/evidence/route-20261114/README.md)
(made with [showboat](https://github.com/simonw/showboat); `uvx showboat verify` re-runs every block).

- On 66 real Asian genomes, a small Asian haplotype panel, our gene builder and re-fetching reads the standard
  alignment lost together rebuild exactly 45 of 134 HLA gene copies that are missing from the public allele
  database (standard T1K with that database alone: 0), without making ordinary typing worse.
- The re-fetch fixes a hole in DRB1 (regions with no reads 103 -> 0) and passed its frozen test (12 -> 18 of 134).
- More haplotypes from the same population keep helping (the curve is still rising at 40).
- Picking paths from a pangenome graph does not work as a typer, and nothing we tried helps at very low depth.
