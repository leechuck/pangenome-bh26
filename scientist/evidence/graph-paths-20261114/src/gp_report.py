"""Summary tables from row-level outputs (no sequences). Runs on the workstation on copied TSVs.
  python gp_report.py RESULTS_DIR"""
import collections
import csv
import sys
from pathlib import Path

R = Path(sys.argv[1])


def load(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def tb(x):
    return x == 'True'


def pct(a, b):
    return f'{a} of {b} ({100 * a / b:.0f}%)' if b else f'{a} of 0'


out = []
P = out.append
place = load(R / 'placement_rows.tsv')
samp = load(R / 'sampling_rows.tsv')
arms = load(R / 'arms_rows.tsv')

# ---------------------------------------------------------------- placement
P('## Placement (GRAPH_TEST.md section 4 measure)\n')
P('Reads that contain the true local version at a site where the true copy differs from the nearest arm-A call; '
  '"exact" = Giraffe placed the read on a route that spells the whole local version. Testable sites only (local version <= 120 bp).\n')
for variant in ('ipd', 'ipd_cohort'):
    P(f'\n**Graph `{variant}`**\n')
    P('| gene | read set | sites | sites with true version in a stored path | reads exact / reads carrying true version | intron-1 run (col 5785) | intron-2 repeat | other length | substitution | sites >= 80% exact | window in one stored path: exact | window only as new route: exact |')
    P('|---|---|---:|---:|---|---|---|---|---|---|---|---|')
    for g in ('DRB1', 'DQA1', 'B', 'all'):
        for rs in ('D-sim', 'D-real'):
            S = [r for r in place if r['variant'] == variant and r['readset'] == rs and (g == 'all' or r['gene'] == g) and tb(r['testable'])]
            if not S:
                continue

            def agg(sel):
                X = [r for r in S if sel(r)]
                return sum(int(r['reads_true_exact']) for r in X), sum(int(r['reads_true']) for r in X)
            a = agg(lambda r: True)
            hp = agg(lambda r: r['named'] == 'homopolymer_5785')
            i2 = agg(lambda r: r['named'] == 'intron2_repeat')
            ol = agg(lambda r: r['named'] in ('', '-') and r['kind'] == 'length')
            su = agg(lambda r: r['named'] in ('', '-') and r['kind'] == 'substitution')
            w1 = agg(lambda r: tb(r['window_in_stored_path']))
            w0 = agg(lambda r: not tb(r['window_in_stored_path']) and tb(r['window_is_walk']))
            s80 = sum(1 for r in S if int(r['reads_true']) and int(r['reads_true_exact']) >= 0.8 * int(r['reads_true']))
            s_cov = sum(1 for r in S if int(r['reads_true']))
            P(f"| {g} | {rs} | {len(S)} | {sum(tb(r['true_in_stored_path']) for r in S)} | {pct(*a)} | {pct(*hp)} | {pct(*i2)} | {pct(*ol)} | {pct(*su)} | {s80} of {s_cov} | {pct(*w1)} | {pct(*w0)} |")
P('\nOld PGGB `hprc_asian` graph (GRAPH_TEST.md, DRB1): simulated 588 of 1,422 (41%), real 845 of 1,625 (52%); '
  'intron-1 run 64/64 and 35/35; intron-2 repeat 402 of 898 (45%) and 688 of 1,271 (54%).\n')

# ---------------------------------------------------------------- sampling
def summarize(rows):
    ab = [r for r in rows if tb(r['ipd_absent'])]
    pr = [r for r in rows if not tb(r['ipd_absent'])]
    donors = collections.defaultdict(list)
    for r in rows:
        donors[r['donor']].append(tb(r['exact']))
    return dict(n=len(rows), ab=len(ab), pr=len(pr),
                ab_exact=sum(tb(r['exact']) for r in ab), pr_exact=sum(tb(r['exact']) for r in pr),
                both=sum(all(v) for v in donors.values()), nd=len(donors),
                wrong_ab=sum(int(r['dist']) for r in ab), wrong_pr=sum(int(r['dist']) for r in pr),
                ab_any=sum(tb(r['exact_any']) for r in ab), pr_any=sum(tb(r['exact_any']) for r in pr),
                halluc=sum(not tb(r['exact']) for r in pr),
                novel_pick=sum(not tb(r['picked_stored']) for r in rows),
                errors=sum(bool(r['error']) for r in rows))


P('\n## Path picking with vg haplotypes\n')
P('Per copy: gene body (first exon start to last exon end) of the paired sampled haplotype vs the true copy, '
  'edlib global distance. dip2 = `--diploid-sampling` (best pair of 32 candidates); hap4/hap8 = 4 or 8 sampled '
  'haplotypes, best one per copy (candidate generator; columns "any"). Input i = all T1K candidate reads, '
  'ii = reads assigned to the gene by bwa with no mate on another gene.\n')
keys = sorted({(r['gene'], r['readset'], r['variant'], r['input'], int(r['subchain']), r['mode']) for r in samp},
              key=lambda k: (('DRB1', 'DQA1', 'B').index(k[0]), k[1], k[2], k[3], -k[4], k[5]))
for g in ('DRB1', 'DQA1', 'B'):
    for rs in ('D-sim', 'D-real'):
        P(f'\n**{g}, {rs}**\n')
        P('| graph | input | subchain | mode | IPD-absent body exact | IPD-present body exact (non-exact = hallucination) | donors both exact | wrong bases absent / present | picked body not a stored path | absent among candidates | present among candidates | errors |')
        P('|---|---|---:|---|---|---|---|---|---:|---|---|---:|')
        for k in keys:
            if k[0] != g or k[1] != rs:
                continue
            X = [r for r in samp if (r['gene'], r['readset'], r['variant'], r['input'], int(r['subchain']), r['mode']) == k]
            s = summarize(X)
            P(f"| {k[2]} | {k[3]} | {k[4]} | {k[5]} | {s['ab_exact']}/{s['ab']} | {s['pr_exact']}/{s['pr']} | {s['both']}/{s['nd']} | "
              f"{s['wrong_ab']} / {s['wrong_pr']} | {s['novel_pick']} | {s['ab_any']}/{s['ab']} | {s['pr_any']}/{s['pr']} | {s['errors']} |")

# ---------------------------------------------------------------- arms
P('\n## Existing arms on the same donors (gene-body distance of the first label of each slot)\n')
P('| gene | read set | arm | IPD-absent body exact | IPD-present body exact | wrong bases absent / present | donors both exact |')
P('|---|---|---|---|---|---|---|')
for g in ('DRB1', 'DQA1', 'B'):
    for rs in ('D-sim', 'D-real'):
        for arm in ('A', 'B', 'C4', 'D4'):
            X = [r for r in arms if r['gene'] == g and r['readset'] == rs and r['arm'] == arm]
            if not X:
                continue
            ab = [r for r in X if tb(r['ipd_absent'])]
            pr = [r for r in X if not tb(r['ipd_absent'])]
            donors = collections.defaultdict(list)
            for r in X:
                donors[r['donor']].append(tb(r['exact']))
            P(f"| {g} | {rs} | {arm} | {sum(tb(r['exact']) for r in ab)}/{len(ab)} | {sum(tb(r['exact']) for r in pr)}/{len(pr)} | "
              f"{sum(int(r['dist']) for r in ab)} / {sum(int(r['dist']) for r in pr)} | {sum(all(v) for v in donors.values())}/{len(donors)} |")
(R / 'TABLES.txt').write_text('\n'.join(out) + '\n')
print('\n'.join(out))

# ---------------------------------------------------------------- alternative picker
altp = R / 'alt_rows.tsv'
if altp.exists():
    alt = load(altp)
    P('\n## Alternative path picking: edge-coverage DP on the clean graph (Giraffe, input ii)\n')
    P('| gene | read set | graph | IPD-absent body exact | IPD-present body exact | donors both exact | wrong bases absent / present | route not a stored path | errors |')
    P('|---|---|---|---|---|---|---|---:|---:|')
    for g in ('DRB1', 'DQA1', 'B'):
        for rs in ('D-sim', 'D-real'):
            for variant in ('ipd', 'ipd_cohort'):
                X = [r for r in alt if r['gene'] == g and r['readset'] == rs and r['variant'] == variant]
                if not X:
                    continue
                ab = [r for r in X if tb(r['ipd_absent'])]
                pr = [r for r in X if not tb(r['ipd_absent'])]
                donors = collections.defaultdict(list)
                for r in X:
                    donors[r['donor']].append(tb(r['exact']))
                P(f"| {g} | {rs} | {variant} | {sum(tb(r['exact']) for r in ab)}/{len(ab)} | {sum(tb(r['exact']) for r in pr)}/{len(pr)} | "
                  f"{sum(all(v) for v in donors.values())}/{len(donors)} | {sum(int(r['dist']) for r in ab)} / {sum(int(r['dist']) for r in pr)} | "
                  f"{sum(not tb(r['picked_stored']) for r in X)} | {sum(bool(r['error']) for r in X)} |")
    (R / 'TABLES.txt').write_text('\n'.join(out) + '\n')
