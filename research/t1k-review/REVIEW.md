# T1K and AsianHLA: review, not implementation

Reviewed 2026-09-18. Paper: Song et al., *Efficient and accurate KIR and HLA
genotyping with massively parallel sequencing data*, Genome Research (2023),
doi:10.1101/gr.277585.122. Full text is saved in paper.html and paper.txt;
Figure 1 is figure1.jpg. PDF downloads were blocked; no placeholder PDF retained.
T1K/ contains upstream HEAD; T1K-v1.0.6/ is the tag matching the benchmark's
Conda version. Commits and file hashes are in provenance.json.

## Interpretation of the current gap

Do not conclude that SpecHLA is inherently inferior from this modified pipeline.
The 64-donor two-field results show a severe HLA-A-specific deficit: T1K 64/64,
HPRC DogoHLA 8/64, Asian DogoHLA 7/64, updated SpecHLA-long 5/64. HLA-A
accounts for 56 of the net 88-genotype gap between T1K and HPRC DogoHLA.
HLA-B and HLA-C tie T1K in that DogoHLA arm. See matched-per-gene.tsv.

The deployed SpecHLA whole-gene naming routine crops to fixed coordinates,
aggregates local BLAST HSP lengths/errors and ranks percentage identity without
an explicit full-query coverage penalty. Our updated database configuration uses
nonuse, disabling the legacy frequency filter. A larger candidate set plus this
scoring is a plausible contributor, not a demonstrated sole cause.
For HG00658, an existing HPRC-arm BLAST output gives A*11:335:01 a perfect
3071-base local match, versus 98.533% over 3203 alignment bases for truth-compatible
A*11:01:01:01. The second reconstruction also favours a wrong allele. Thus
reconstruction errors must be investigated as well; correcting naming alone
has not been shown sufficient. Original outputs remain unchanged.

## T1K mechanisms worth retaining

Candidate extraction and paired-end allele alignment feed a joint abundance
model across genes. Equivalence classes and weighted EM preserve competing
read assignments; SQUAREM accelerates convergence. Pair selection evaluates
read support and the code includes missing-exonic-coverage penalties and
cross-gene noise-based confidence. It need not reconstruct two complete genes
before calling known alleles. ReadAssignmentWeight uses discrete similarity
weights; do not describe it as a full base-quality-aware Bayesian model.

The paper used v1.0.2; our baseline uses v1.0.6. In v1.0.6, ParseDatFile.pl
rescues suitable incomplete alleles by filling missing introns with modal gene
sequences, unlike the paper's original exclusion policy. The default DNA
reference keeps 200 bases of intron on each exon flank and 50 bases of UTR.
The source also supports genome mode, so simply using longer references is
not itself a novel graph method. Database multiplicity influences sequence
weights/initialization; these are not measured population frequencies.

## Proposed priorities, requiring a later implementation decision

1. Diagnose the HLA-A regression on development samples: separate read
   assignment, reconstructed sequence, database preparation and naming effects.
   Preserve the frozen benchmark and label any corrective evaluation separately.
2. Use the Asian pangenome's observed intron/exon/flank linkage to distinguish
   full genomic variants within candidate allele families. Keep imputed sequence
   visibly separate from observed sequence. This targets four-field resolution.
3. Represent related HLA genes and pseudogenes jointly, retaining uncertain
   read and mate placements. The paper's HLA-U/A*02:783 example motivates this;
   an exon-only whitelist cannot establish the true locus of a read.
4. Add a graph-derived likelihood stage for candidate diplotypes, preserving
   T1K's joint competition. Map to phased paths, then marginalize compatible
   placements/paths; do not convert the graph mapping into hard per-gene bins.
   Retain all-IPD alternatives and an unknown/novel route when no panel path fits.
5. Explore weak, regularized local haplotype priors, with global fallback and
   family-balanced counts. Do not count duplicated panel paths as independent
   read evidence or exclude an allele because it is absent from Asian donors.

The present DogoHLA graph intervention is chiefly conservative homozygous
noncoding long-indel repair at DRB1; its panel also changes upstream read
collection. It is not a general graph-based joint HLA allele inference model.

## Evaluation proposal

First preserve ordinary T1K and repair our diagnostic baseline. Compare T1K
against graph-derived reference sequences with its existing aligner, then graph
alignment with the joint likelihood, then optional population priors. The
linear graph-derived-reference control is essential to distinguish additional
sequence information from graph alignment. Compare size-controlled HPRC and
Asian panels, using the same IPD, reads, donor/family exclusions and endpoints.
Retain absent-panel alleles, no-calls, ancestry strata, novel sequence and
four-field ambiguity. No short-read evidence means unresolved, not an invented
four-field assignment. Current benchmark outcomes are now observed and cannot
serve as a fresh untouched test set for choices made from this review.

No typer code was changed and no experiment was launched for this review.
Subsequent implementation authorization and the validation reservation are in
../../t1k-pangenome/PROTOCOL.md.

## Personalized pangenome paper and additive panel audit

Siren et al., *Personalized pangenome references*, Nature Methods (2024),
doi:10.1038/s41592-024-02407-2, is saved as paper-PMC12643174.html/.txt.
It motivates selecting a sample-specific subset from a broad source panel using
read k-mer support before alignment. It is a whole-genome variant-calling study,
not evidence of HLA typing improvement over T1K. Its graph-topology, coverage and
blockwise phasing assumptions must be checked for our locus graphs.

The live membership audit confirms that full contains all 340 HPRC/reference
members plus 282 non-HPRC Asian-source haplotypes; it removes no HPRC members.
The replacement panel removes 282 HPRC haplotypes and substitutes those 282.
On the completed matched cohort, full DogoHLA scores 404/509, HPRC 407/509,
and replacement 400/509. Full gains 3 and loses 6 two-field genotypes relative
to HPRC. The no-graph full and HPRC arms both score 402/509. These are exploratory
results; no graph improvement is established. Full is already the additive
comparison requested by the user, and no duplicate experiment is needed.
