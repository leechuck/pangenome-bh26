# BioHackathon 2026 wrap-up: three summary slides

`slides.tex` (Beamer, 16:9), `slides.pdf`, `render/slide-*.png`. Figures in `figures/` are built by
`make_figures.py` from result tables in this repository; `fig_cohorts`, `fig_graph_density`,
`fig_variant_classes` and `fig_hla_dev` are new, the workflow and method schematics are TikZ inside the deck.

```sh
python3 make_figures.py          # needs hla-bench/ and hla-typer/ (scorer, labels) in the repo
pdflatex slides.tex; pdflatex slides.tex
pdftoppm -png -r 110 slides.pdf render/slide
```

## Narrative

**Slide 1, workflow and aims.** Two kinds of input (assemblies from APR, HPRC release 2 and JaSaPaGe;
haplotype walks from the K-PanRef and CPC graphs) go through one MHC extraction step, then split into
Immuannot allele calling and Minigraph-Cactus graph construction. The allele labels feed the graph uses.
Four outputs: allele landscape per cohort, novel alleles with read support, SNV/indel/SV genotyping
(PanGenie), short-read HLA typing (Locityper). The cohort bar chart shows that 390 of 754 haplotypes are
Asian or Arab, which is the reference-diversity claim the other two slides test. The three aims are the
three questions the hackathon set out to answer.

**Slide 2, graph and variant validation.** Left figure: node density and haplotype sharing along the GRCh38
path of the 754-haplotype graph (from `hla-typer/results/graph_scan/ref_nodes.tsv.gz`). The graph is
near-linear except at the HLA genes, C4 and the DR-DQ block, which is where a pangenome can matter. Right
and bottom: PanGenie genotypes of 40 held-out 1000 Genomes donors against their own assemblies, full panel
versus HPRC-only panel on the same sites (`hla-asian50/refined/results/variant_summary.tsv`, universe
`frozen_HPRC_sites`), per variant class: SNV, small complex, SV site, SV-bearing truth genotype. Gain grows
with variant complexity (+0.1 to +3.9 points). Third panel: SNVs shared with the published NYGC GATK callset,
where both graph panels beat the linear callset. Caveats on the slide: assembly-derived truth, panel size
and ancestry confounded.

**Slide 3, HLA typing.** Schematic of the typer: recruited reads, Locityper on seven anchor-bounded loci cut
from the graph with a leave-one-donor-out panel, open-set candidates (IPD-IMGT/HLA exons grafted onto
panel backbones plus T1K calls), naming ladder, outputs. Left figure: per-gene two-field accuracy on the 40
development donors for three panel arms and T1K (`hla-typer/results/stageB/scores.tsv`,
`hla-typer/results/comparators/t1k/`). T1K leads on the classical genes (98.4 vs 92.5 percent); the panel
wins on DRB3/4/5 with copy number (111/113 vs 48/110) and the Asian size-matched panel beats HPRC-only
(+4.4 points). Right figure: SpecHLA with the pangenome in its read-binning step reconstructs more gene
sequences exactly (51 vs 41 of 128 haplotypes; `hla-spechla-pg/results/sequence_scores.tsv`). Status line:
the locked test on 228 leave-one-out donors and 946 Gourraud-typed donors has not run.

## Figures not used and why

- `hla/results/figures/fig5_population_hla.png` (allele landscape): good but a full slide on its own; slide 1
  cites the result in the workflow box instead.
- `hla-asian50/refined/results/comparison.png`: superseded by the per-variant-class figure; its HLA panels
  use the older Locityper setup with the recruitment defect.
- odgi `MHC.viz` rendering of the graph: 7,900 px tall, unreadable on a slide; replaced by the density plot.
- `hla-structural/results/benchmark.png` (RCCX/DRB structural typing): separate line of work, no room.
