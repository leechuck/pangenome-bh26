# Revised scope: separately analysed Asian strata

**Current scope:** Robert has deferred Saudi/Arabian validation and requested 1000 Genomes only. The running cohort is 20 East Asian plus 20 South Asian donors, analysed separately. The three-stratum proposal below is retained as the future extension; Arabian inputs are no longer a blocker. This changes validation donors, not the composition of the full reference graph.

Robert requested a mixed East, West Asian/Arabian and South Asian pilot, with separate analysis because graph representation differs across populations. This supersedes the East-Asian-only accuracy-benchmark proposal in README.md. The existing 50-person selection and representation results remain a historical feasibility experiment; they are not results for the revised cohort.

## Available candidates

Intersecting the existing pedigree metadata with two eligible RCCX assembly annotations gives:

| Stratum | Eligible donors | Population breakdown | Donors in current experimental HLA truth resource |
|---|---:|---|---:|
| East Asian | 55 | JPT 21, KHV 13, CHS 11, CDX 5, CHB 5 | 23 |
| South Asian | 36 | PJL 10, BEB 10, GIH 8, STU 5, ITU 3 | 0 |

Candidate IDs are in `source/mixed_asian_candidates.tsv`. This is an eligibility inventory, not a final selection, confirmation of read availability, or a truth-quality guarantee. No prediction outcomes were used. There are 55 distinct recorded families among East Asian candidates and 36 among South Asian candidates.

The [official 1000 Genomes population metadata](https://github.com/igsr/1000Genomes_data_indexes/blob/master/README_populations.md) has no dedicated Central Asian or Arabian/West Asian stratum. Saudi genomes represent West Asia, not Central Asia. Pakistani PJL remains classified as South Asian; it must not be relabelled Central Asian to fill a quota. Korean reference haplotypes likewise do not create a Korean 1000 Genomes validation population.

Robert confirmed **West Asia / Arabian** as the third stratum. Central Asia is outside this pilot. West Asian/Arabian validation requires another source of donor-matched short reads and reliable truth. Existing paired RCCX annotations include five JaSaPaGe-Saudi donors and 53 APR donors, but the latter are not automatically Saudi or Central Asian, and these counts do not establish independent, usable short-read validation samples. Do not infer ancestry from sample names or assembly-provider labels alone.

## Analysis changes

- Use approximately balanced broad strata within a roughly 50-donor total, subject to confirmed third-stratum eligibility. Do not freeze numerical quotas or replacement samples until matched reads, donor metadata and validation data are verified. Within strata, balance named populations where feasible and report population counts.
- Treat East Asian, South Asian and West Asian/Arabian strata as separate primary summaries. Show each method's accuracy, call coverage, uncertainty, truth denominator and full-panel versus HPRC-only paired difference within each stratum. A pooled figure may accompany these results but cannot replace them. Population-level results with few samples are descriptive.
- Measure training-panel representation after family exclusions: distinct donors/haplotypes, locus completeness, allele/structure availability and sequence distance to available haplotypes. Reference-donor counts alone do not establish how well a test allele is represented.
- Test the hypothesis that poorer represented South Asian haplotypes have different error rates or gains. Do not assume South Asian performance is worse; HLA alleles are shared across populations and HPRC already includes South Asian donors.
- Preserve strict family/alias exclusion and full-panel versus HPRC-only inference with identical inputs. Consider size-matched panel subsampling to distinguish more haplotypes from their population composition. Do not tune either arm on new test predictions.
- Classical experimental HLA accuracy is currently unavailable for South Asian candidates in the saved truth resource. Obtain independent labels or report assembly-derived agreement as a distinct tier. Never compare experimental-label accuracy in one stratum with assembly-label agreement in another as if the endpoints were identical.
- The same frozen caller settings, read recruitment, callable masks and endpoint definitions must apply across strata. Any external third cohort requires explicit sequencing/batch/coverage reporting; population differences may otherwise be technical differences.

**Execution has begun:** see [the initial run](run/README.md) for submitted jobs and remaining validation work. No completed mixed-cohort accuracy result is available. The scripts `select_panel.py` and `check_setup.py` currently reproduce the historical East Asian feasibility experiment only; they must not be used to claim a completed mixed-cohort selection.
