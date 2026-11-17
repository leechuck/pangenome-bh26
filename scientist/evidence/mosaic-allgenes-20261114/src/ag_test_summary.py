"""Pre-registered endpoints (DESIGN.md section 4) from the row-level outputs of ag_score.py and
ag_proposals_table.py. Numbers only. Family-clustered bootstrap and sign-flip test: frozen score.uncertainty.
  python ag_test_summary.py TAG READSET[,READSET]
Writes results/TEST_SUMMARY_<TAG>.md and results/endpoints_<TAG>.json."""
import collections
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L          # noqa: E402
C = L.C
import score as SC            # noqa: E402

S2 = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114')
GENES = C.GENES
NON_DRB1 = tuple(g for g in GENES if g != 'DRB1')
PRIMARY = ('C4R4n', 'D4R4n')
ARMS = ('A', 'B', 'C4', 'C4R4n', 'C4R1n', 'D4', 'D4R4n', 'D4R1n')
FALLBACK_KINDS = ('timeout', 'bound', 'error')
T = lambda x: x == 'True'
tag, readsets = sys.argv[1], sys.argv[2].split(',')
res = S2 / 'results'
rows = list(csv.DictReader(open(res / f'score_rows_{tag}.tsv'), delimiter='\t'))
tasks = list(csv.DictReader(open(res / f'tasks_{tag}.tsv'), delimiter='\t'))
haps = list(csv.DictReader(open(res / f'haplotypes_{tag}.tsv'), delimiter='\t'))
group_of = {d: tuple(sorted(v)) for d, v in C.exclusion_map().items()}
by = collections.defaultdict(dict)
for r in rows:
    by[(r['readset'], r['arm'])][(r['donor'], r['gene'], r['haplotype'])] = r


def base(arm):
    return 'B' if arm[0] in ('B', 'D') else 'A'


def ed(r):
    return None if r['ed_final'] in ('None', '') else int(r['ed_final'])


def right(r):
    return T(r['recovered_final']) or T(r['final_body_ok'])


def body(r):
    return T(r['body_recovered']) or T(r['recovered'])


def paired(rs, arm, val, sel, scale=1.0):
    b = by.get((rs, base(arm)), {})
    pairs = []
    for k, r in by.get((rs, arm), {}).items():
        if k in b and sel(r):
            x, y = val(r), val(b[k])
            if x is None or y is None:
                continue
            pairs.append(dict(donor=r['donor'], difference=int(x) - int(y)))
    if not pairs:
        return None
    u = SC.uncertainty(pairs, group_of)
    return dict(n=len(pairs), sum=sum(p['difference'] for p in pairs), diff=u['difference_pp'] / 100 * scale,
                ci=[u['bootstrap95_pp'][0] / 100 * scale, u['bootstrap95_pp'][1] / 100 * scale],
                p=u['sign_flip_p'], clusters=u['family_clusters'])


def fmt(u, pp=True):
    if u is None:
        return 'n/a'
    m = 100 if pp else 1
    unit = ' pp' if pp else ''
    return (f"{u['sum']:+d} of {u['n']}; {u['diff'] * m:+.{1 if pp else 2}f}{unit} [{u['ci'][0] * m:+.{1 if pp else 2}f}, "
            f"{u['ci'][1] * m:+.{1 if pp else 2}f}], p={u['p']:.3g}, {u['clusters']} clusters")


def counts(rs, arm, genes):
    R = [r for k, r in by.get((rs, arm), {}).items() if k[1] in genes]
    ab = [r for r in R if T(r['ipd_absent'])]
    pr = [r for r in R if not T(r['ipd_absent']) and r['truth4']]
    eds = [ed(r) for r in R if ed(r) is not None]
    return dict(absent=len(ab), body=sum(body(r) for r in ab), exact_final=sum(T(r['recovered_final']) for r in ab),
                exact_any=sum(T(r['recovered']) for r in ab), named4=sum(T(r['named4']) for r in pr), n4=len(pr),
                novel_present=sum(T(r['final_mos']) for r in R if not T(r['ipd_absent'])),
                right=sum(T(r['final_mos']) and right(r) for r in ab),
                wrong=sum(T(r['final_mos']) and not right(r) for r in ab),
                ed_mean=sum(eds) / len(eds) if eds else None, ed_missing=len(R) - len(eds), copies=len(R))


def cfs(rs, arm, genes):
    b = by.get((rs, base(arm)), {})
    c = collections.Counter()
    for k, r in by.get((rs, arm), {}).items():
        if k[1] in genes and T(r['final_mos']) and k in b and ed(r) is not None and ed(b[k]) is not None:
            c['closer' if ed(r) < ed(b[k]) else 'farther' if ed(r) > ed(b[k]) else 'equal'] += 1
    return c


out, E = [], {}
P = out.append
P(f'# Endpoints, tag {tag}\n')
# f: fallbacks
fb = collections.defaultdict(list)
for t in tasks:
    if t['readset'] in readsets and t['status'] != 'complete':
        fb[(t['readset'], t['donor'], t['arm'])].append('task_' + t['status'])
for h in haps:
    if h['readset'] not in readsets:
        continue
    if h['gene_status'] != 'complete':
        key = (h['readset'], h['donor'], h['arm'])
        if h['failure_kind'] in FALLBACK_KINDS and f"{h['gene']}:{h['failure_kind']}" not in fb[key]:
            fb[key].append(f"{h['gene']}:{h['failure_kind']}")
gene_fail = collections.Counter()
seen = set()
for h in haps:
    k = (h['readset'], h['donor'], h['arm'], h['gene'])
    if h['readset'] in readsets and h['gene_status'] != 'complete' and k not in seen:
        seen.add(k)
        gene_fail[h['failure_kind']] += 1
avail_f = res / f'score_avail_{tag}.tsv'
if avail_f.exists():
    for line in open(avail_f):
        x = line.rstrip('\n').split('\t')
        if len(x) == 4 and x[0] in readsets and x[2][:1] in 'CD' and x[3] != 'complete':
            key = (x[0], x[1], x[2][0])
            if 'scored_' + x[3] not in fb[key]:
                fb[key].append('scored_' + x[3])
n_tasks = sum(t['readset'] in readsets for t in tasks)
fbt = {k: v for k, v in fb.items() if v}
E['f'] = dict(tasks=n_tasks, fallback_tasks=len(fbt), detail={'/'.join(k): v for k, v in fbt.items()},
              gene_failures=dict(gene_fail), pass_=len(fbt) <= 2)
P(f"**f. fallbacks:** {len(fbt)} of {n_tasks} tasks ({'PASS' if len(fbt) <= 2 else 'FAIL'}; limit 2). "
  f"Gene-level non-complete by kind: {dict(gene_fail) or 'none'}. Detail: {E['f']['detail'] or 'none'}\n")
for rs in readsets:
    E[rs] = {}
    P(f'\n## {rs}\n')
    have = [a for a in ARMS if (rs, a) in by]
    # a
    a_ = {arm: counts(rs, arm, GENES)['novel_present'] for arm in have if arm[0] in 'CD'}
    ok_a = all(a_[x] <= 1 for x in PRIMARY if x in a_)
    E[rs]['a'] = dict(novel_present=a_, pass_=ok_a)
    P(f"**a. novel final calls on IPD-present copies (all genes):** {a_} -> {'PASS' if ok_a else 'FAIL'} (primary arms <= 1)\n")
    # b
    bl, ok_b = [], True
    for arm in [x for x in PRIMARY if x in have]:
        ca, cb = counts(rs, arm, GENES), counts(rs, base(arm), GENES)
        tot_ok = ca['named4'] >= cb['named4']
        per = {g: counts(rs, arm, (g,))['named4'] - counts(rs, base(arm), (g,))['named4'] for g in GENES}
        gene_ok = all(v >= -1 for v in per.values())
        ok_b &= tot_ok and gene_ok
        u = paired(rs, arm, lambda r: T(r['named4']), lambda r: not T(r['ipd_absent']) and bool(r['truth4']))
        bl.append(f"{arm} vs {base(arm)}: total {cb['named4']} -> {ca['named4']} of {ca['n4']} ({fmt(u)}); per-gene change {per}")
        E[rs].setdefault('b', {})[arm] = dict(base=cb['named4'], arm=ca['named4'], n=ca['n4'], per_gene=per, total_ok=tot_ok,
                                              gene_ok=gene_ok, interval=u)
    E[rs].setdefault('b', {})['pass_'] = ok_b
    P(f"**b. four-field named accuracy:** {'PASS' if ok_b else 'FAIL'}\n")
    for x in bl:
        P(f'- {x}')
    # c, d (+ secondary non-DRB1)
    for lab, arm in (('c', 'C4R4n'), ('d', 'D4R4n')):
        if arm not in have:
            continue
        for gl, genes in (('all genes', GENES), ('non-DRB1 genes (secondary)', NON_DRB1)):
            u = paired(rs, arm, body, lambda r, genes=genes: T(r['ipd_absent']) and r['gene'] in genes)
            ca, cb = counts(rs, arm, genes), counts(rs, base(arm), genes)
            E[rs].setdefault(lab, {})[gl] = dict(base=cb['body'], arm=ca['body'], absent=ca['absent'], interval=u)
            verdict = ''
            if lab == 'c' and gl == 'all genes':
                ok = u is not None and u['ci'][0] > 0
                E[rs][lab]['pass_'] = ok if rs == 'T-sim' else 'reported'
                verdict = (' -> ' + ('PASS' if ok else 'FAIL')) if rs == 'T-sim' else ' (reported)'
            P(f"\n**{lab}. gene-body recovery of IPD-absent copies, {arm} minus {base(arm)}, {gl}:** {cb['body']} -> {ca['body']} "
              f"of {ca['absent']}; {fmt(u)}{verdict}")
    # e
    P('\n**e. novel final calls on IPD-absent copies: right gene body / not**')
    E[rs]['e'] = {}
    for full, ab_ in (('C4', 'C4R4n'), ('D4', 'D4R4n')):
        if full not in have:
            continue
        c0, c1 = counts(rs, full, GENES), counts(rs, ab_, GENES)
        keep = c1['right'] / c0['right'] if c0['right'] else None
        ok = keep is None or keep >= 0.75
        E[rs]['e'][ab_] = dict(without=(c0['right'], c0['wrong']), with_=(c1['right'], c1['wrong']), keep=keep, pass_=ok)
        sh = lambda c: f"{c['right']}/{c['right'] + c['wrong']}" + (f" = {100 * c['right'] / (c['right'] + c['wrong']):.0f}%" if c['right'] + c['wrong'] else '')
        P(f"- {full}: {sh(c0)} right; {ab_}: {sh(c1)} right; right calls kept {c1['right']}/{c0['right']}"
          f" -> {'PASS' if ok else 'FAIL'}{' (no right calls without abstaining)' if keep is None else ''}")
    # secondary per-gene table
    P(f'\n**Secondary: per gene ({rs})**\n')
    P('| gene | arm | IPD-absent | body | byte-exact final | 4-field named | novel on present | novel absent right/wrong | mean wrong bases per copy | closer/equal/farther |')
    P('|---|---|---|---|---|---|---|---|---|---|')
    for g in list(GENES) + ['all', 'non-DRB1']:
        genes = GENES if g == 'all' else NON_DRB1 if g == 'non-DRB1' else (g,)
        for arm in have:
            c = counts(rs, arm, genes)
            cc = cfs(rs, arm, genes) if arm not in ('A', 'B') else None
            P(f"| {g} | {arm} | {c['absent']} | {c['body']} | {c['exact_final']} | {c['named4']}/{c['n4']} | {c['novel_present']} | "
              f"{c['right']}/{c['wrong']} | {'%.2f' % c['ed_mean'] if c['ed_mean'] is not None else 'n/a'}"
              f"{'' if not c['ed_missing'] else ' (' + str(c['ed_missing']) + ' missing)'} | "
              f"{'' if cc is None else '%d/%d/%d' % (cc['closer'], cc['equal'], cc['farther'])} |")
    P(f'\n**Secondary: mean wrong bases per copy, arm minus base, paired, family-clustered ({rs})**\n')
    for arm in [a for a in have if a not in ('A', 'B')]:
        for gl, genes in (('all genes', GENES), ('non-DRB1', NON_DRB1)):
            for cl, sel in (('all copies', lambda r: True), ('IPD-absent copies', lambda r: T(r['ipd_absent']))):
                u = paired(rs, arm, ed, lambda r, genes=genes, sel=sel: r['gene'] in genes and sel(r))
                P(f'- {arm} minus {base(arm)}, {gl}, {cl}: {fmt(u, pp=False)} bases per copy')
                E[rs].setdefault('ed', {})[f'{arm}|{gl}|{cl}'] = u
    P(f'\n**Secondary: byte-exact recovery (final call) of IPD-absent copies, all genes ({rs})**\n')
    for arm in [a for a in have if a not in ('A', 'B')]:
        u = paired(rs, arm, lambda r: T(r['recovered_final']), lambda r: T(r['ipd_absent']))
        P(f'- {arm} minus {base(arm)}: {fmt(u)}')
(res / f'TEST_SUMMARY_{tag}.md').write_text('\n'.join(out) + '\n')
(res / f'endpoints_{tag}.json').write_text(json.dumps(E, indent=1, default=str) + '\n')
print('\n'.join(out))
