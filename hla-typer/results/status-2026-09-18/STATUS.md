# HLA typing status — 2026-09-18

Live DDBJ check at 11:50 JST: no running jobs for leechuck. Weighted Stage B array
20641990 tasks 18–39 and SpecHLA recruitment 20641142 tasks 8–39 remain user-held.
The previous check found all five 20659273 prior tasks failed on a missing
`discarded_haplotypes.txt`; the normalised fold-0 DRB1 graph build completed.

## Newly scored completed development runs

Downloaded compact C1 result/candidate snapshots from DDBJ to `c1-snapshot.tar.gz`
and `c1/`, preserving older local snapshots. Scored with:

```
python3 hla-typer/stageB/score_calls.py --c1 \
  --calls hla-typer/results/status-2026-09-18/c1 \
  --out hla-typer/results/status-2026-09-18/scores
```

Only donors with COMPLETE markers enter the comparisons below. All 40 full-ens
donors completed. This is development data, not the locked test.

| Endpoint | Stage B full | C1 full-ens |
|---|---:|---:|
| Classical eight genes, two-field | 295/319 (92.5%) | 297/319 (93.1%) |
| Classical eight genes, G group | 296/319 (92.8%) | 303/319 (95.0%) |
| Classical eight genes, exact CDS | 292/320 (91.3%) | 294/320 (91.9%) |
| DRB3/4/5, type plus copy number | 111/113 (98.2%) | 109/113 (96.5%) |
| All genes, two-field | 406/432 (94.0%) | 406/432 (94.0%) |

T1K classical two-field, corrected to the same donor/assembly selection:
314/319 (98.4%). See `hackathon-wrapup/make_figures.py` and its generated
`figures/hla_dev_numbers.txt`. The earlier T1K denominator excluded NA19088
because it lacked Stage B's JaSaPaGe fallback.

The four completed full-ens-p donors score 27/31 classical two-field, versus
29/31 for Stage B on those same donors. This small selected subset gives no
reason to prioritise indiscriminate partial-exon expansion.

Do not use C1 `gene_exact` scores as sequence-reconstruction evidence: grafted
candidate rows have an empty gene_sha256 in candidate_rows(), so that endpoint
needs explicit scoring of the generated sequences.

## Run provenance defect

Despite its name, full-ens is not a T1K ensemble. The live and local
`stageC1/rescore.sbatch` dispatch only `t1k*` to `--t1k`; `ens` and `ens-p`
fall through to `--no-t1k`. Remote job logs show zero t1k_backbones.
All full-ens candidate labels contain only 3,011 backbone rows and 111,777
ipd_near rows; no T1K-labelled candidate sources. The partial arm also includes
ipd_partial. These are no-T1K IPD-grafting runs. Directory names must not be
used as evidence for the experimental condition.

## Interpretation and next development step

We have a working graph-derived haplotype-panel typer: graph-cut loci feed
Locityper, predictions map to allele labels and sequences, and the DRB345
structural block supports presence/copy-number-aware typing. This is distinct
from mapping reads through graph topology and calling recombined/novel paths.

The full development run supports a small classical typing improvement and a
stronger G-group improvement, with regressions elsewhere. It does not yet
justify replacing T1K or claiming a validated new HLA typer.

Priority next experiment: correct and record the actual candidate-source
configuration, verify source counts before genotyping, compare a real T1K-assisted
candidate arm with a pangenome-only arm on the development set, and report paired
rescues/regressions. Keep genotype ranking and candidate coverage separate.
Repair and test prior generation before evaluating it; missing duplicate metadata
must be distinguished from an incomplete database. Priors need review of panel
membership, sequence-collapse handling and candidate-multiplicity effects.

Then freeze the method and execute the locked independent validation. Direct
graph-alignment experiments in SpecHLA remain a separate route; the normalised
DRB1 graph build alone establishes neither mapping performance nor typing gain.
