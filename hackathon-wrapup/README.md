# BioHackathon wrap-up

Three slides for a short spoken wrap-up: **the build → the payoff → the surprise**.

- `slides.pdf`: presentation, 16:9, vector graphics and selectable text.
- `slides.tex`: editable Beamer/TikZ source; numbers explicitly set from the sources below.
- `render/slide-*.png`: slide previews.
- `speaker-notes.md`: short talk track and evidence details.
- `previous/`: deck and README before the audience-focused redesign.
- `make_figures.py`: detailed supporting plots, retained for deeper discussion. Its T1K scoring now uses the same JaSaPaGe fallback as Stage B, restoring NA19088 to the comparison. These supporting PNGs are not embedded in the new slides.

Build from the repository root (two LaTeX passes are required for page-relative coordinates):

```sh
pdflatex -interaction=nonstopmode -halt-on-error -output-directory hackathon-wrapup hackathon-wrapup/slides.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory hackathon-wrapup hackathon-wrapup/slides.tex
pdftoppm -png -r 140 hackathon-wrapup/slides.pdf hackathon-wrapup/render/slide
```

## Evidence

1. Cohorts: `hla-typer/source/haplotype_donors.tsv`; 754 total, 390 in the specified Asian/Arab cohorts. The graph drawing is a conceptual schematic, explicitly labelled, not a locus topology plot.
2. PanGenie: `hla-asian50/refined/results/variant_summary.tsv` and `variant_paired.tsv`, universe `all_truth`, class `truth_SV_length`, endpoint `all`. EAS: HPRC 1060/1465, full 1117/1465; SAS: HPRC 1188/1577, full 1219/1577. EAS delta 3.89 pp (95% interval 2.40–5.29); SAS 1.97 pp (0.32–3.47). This differs intentionally from the previous deck's frozen-HPRC-site plot, to align the displayed rates and intervals to one endpoint. Donors: 20 EAS, 20 SAS. Donor bootstrap intervals; assembly-derived graph truth; test families excluded from panels but test assemblies retained in graph topology. Panel size and ancestry are confounded.
3. SpecHLA: `hla-spechla-pg/results/sequence_scores.tsv` and `REPORT.md`. Native vs mode A (pangenome added to read-extraction/binning database): 41 vs 51 exact sequences of 128, in eight fold-0 development donors. Native designation gives 61/63 two-field genotypes in both arms. Excluded panel donors include fold relatives. Exact sequence matching uses edlib infix alignment, with N counted as a mismatch; see report for its scope.
4. Separate Locityper Stage B experiment: `hla-typer/results/stageB/scores.tsv`, 40 development donors. Full panel classical genes 295/319 (92.5%), DRB3/4/5 111/113 (98.2%). Corrected T1K comparison: 314/319 (98.4%), DRB3/4/5 48/113 (42.5%). T1K does not report the copy-number endpoint, so its DRB3/4/5 score is omitted from the main slide. These are Stage B scores before grafting or T1K candidate integration.

The 228 leave-one-out / 946 experimentally typed donor validation is future work. Newly completed C1 jobs have not been incorporated as scored results.
