"""Step 2b pre-registered endpoints P1-P6 (DESIGN_TEST.txt): re-pulled reads (NEWRS) versus the old reads (OLDRS),
paired within donor x gene x haplotype, same arm; family-clustered bootstrap / sign-flip = frozen score.uncertainty.
Numbers only.
  python end3.py OUTTAG NEW_ROWS OLD_ROWS NEWRS OLDRS NEW_HAPS OLD_HAPS NEW_TASKS T1K_RUN_ROOT
Writes work/ag/results/SUMMARY_OUTTAG.txt and endpoints_OUTTAG.json."""
import collections
import csv
import glob
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import frozen4               # noqa: E402
frozen4.use()
import common3 as K          # noqa: E402
import score as SC           # noqa: E402
C = K.C

outtag, nrows, orows, NRS, ORS, nhaps, ohaps, ntasks, t1kroot = sys.argv[1:10]
GENES = C.GENES
T = lambda x: x == 'True'
group_of = {d: tuple(sorted(v)) for d, v in C.exclusion_map().items()}


def load(p, rs):
    by = collections.defaultdict(dict)
    for r in csv.DictReader(open(p), delimiter='\t'):
        if r['readset'] == rs:
            by[r['arm']][(r['donor'], r['gene'], r['haplotype'])] = r
    return by


N, O = load(nrows, NRS), load(orows, ORS)
body = lambda r: T(r['body_recovered']) or T(r['recovered'])
right = lambda r: T(r['recovered_final']) or T(r['final_body_ok'])
absent = lambda r: T(r['ipd_absent'])
named_sel = lambda r: not T(r['ipd_absent']) and r['truth4'] != ''


def ed(r):
    return None if r['ed_final'] in ('None', '') else int(r['ed_final'])


def paired(X, Y, val, sel):
    """X minus Y over keys present in both"""
    pairs, sx, sy, n = [], 0, 0, 0
    for k, r in X.items():
        if k in Y and sel(r):
            x, y = val(r), val(Y[k])
            if x is None or y is None:
                continue
            pairs.append(dict(donor=r['donor'], difference=int(x) - int(y)))
            sx += int(x); sy += int(y); n += 1
    if not pairs:
        return dict(n=0)
    u = SC.uncertainty(pairs, group_of)
    return dict(n=n, new=sx, old=sy, diff=sx - sy, pp=round(u['difference_pp'], 2),
                ci95_pp=[round(x, 2) for x in u['bootstrap95_pp']], p=u['sign_flip_p'], clusters=u.get('clusters'))


def fmt(u, a='new', b='old'):
    if not u.get('n'):
        return 'n/a'
    return (f"{u[b]} -> {u[a]} of {u['n']} ({u['diff']:+d}; {u['pp']:+.1f} pp [{u['ci95_pp'][0]:+.1f}, {u['ci95_pp'][1]:+.1f}], "
            f"p={u['p']:.3g})")


E, out = {}, []
P = out.append
P(f'# Step 2b endpoints, tag {outtag}: {NRS} (re-pulled reads) versus {ORS} (old reads)\n')
# P1 / P3
for key, arm in (('P1', 'C4R4n'), ('P3', 'D4R4n')):
    u = paired(N[arm], O[arm], body, absent)
    E[key] = dict(arm=arm, **u, passed=(u.get('n', 0) > 0 and u['ci95_pp'][0] > 0) if key == 'P1' else None)
    P(f"**{key}** gene bodies recovered among IPD-absent copies, all genes, {arm} {NRS} minus {ORS}: {fmt(u)}"
      + (f" -> {'PASS' if E[key]['passed'] else 'FAIL'}" if key == 'P1' else ''))
# P2 (i)
P('\n**P2 safety**')
p2 = {}
for arm in ('A', 'B'):
    u = paired(N[arm], O[arm], lambda r: T(r['named4']), named_sel)
    per = {}
    for g in GENES:
        ug = paired(N[arm], O[arm], lambda r: T(r['named4']), lambda r, g=g: named_sel(r) and r['gene'] == g)
        per[g] = ug.get('diff', 0)
    ok = u.get('n', 0) > 0 and u['diff'] >= 0 and min(per.values()) >= -1
    p2[f'i_{arm}'] = dict(total=u, per_gene_change=per, passed=ok)
    P(f"- (i) arm {arm} four-field named, {NRS} vs {ORS}: {fmt(u)}; per gene change {per} -> {'PASS' if ok else 'FAIL'}")
nov = {}
for arm in ('C4R4n', 'D4R4n'):
    nov[arm] = sum(T(r['final_mos']) for r in N[arm].values() if not absent(r))
ok = all(v <= 1 for v in nov.values())
p2['ii'] = dict(novel_final_calls_on_ipd_present=nov, passed=ok)
P(f"- (ii) novel final calls on IPD-present copies ({NRS}): {nov} -> {'PASS' if ok else 'FAIL'}")
# (iii) fallbacks: T1K arm A/B runs (no COMPLETE.json) + builder tasks per the frozen criterion f
fails = []
donors = sorted({k[0] for k in N['A']})
for d in donors:
    for a in ('A', 'B'):
        if not (Path(t1kroot) / d / a / 'COMPLETE.json').exists():
            fails.append(f'{d}:{a}:t1k_not_complete')
for t in csv.DictReader(open(ntasks), delimiter='\t'):
    if t['readset'] == NRS and t['status'] != 'complete':
        fails.append(f"{t['donor']}:{t['arm']}:task_{t['status']}")
for h in csv.DictReader(open(nhaps), delimiter='\t'):
    if h['readset'] == NRS and h['gene_status'] != 'complete' and h['failure_kind'] in ('timeout', 'bound', 'error'):
        fails.append(f"{h['donor']}:{h['arm']}:{h['gene']}:{h['failure_kind']}")
fails = sorted(set(fails))
n_tasks = len(donors) * 4
ok = len({':'.join(f.split(':')[:2]) for f in fails}) <= 2
p2['iii'] = dict(tasks=n_tasks, fallbacks=fails, passed=ok)
P(f"- (iii) tasks in fallback: {len({':'.join(f.split(':')[:2]) for f in fails})} of {n_tasks} {fails[:10]} -> {'PASS' if ok else 'FAIL'}")
p2['passed'] = all(v['passed'] for v in p2.values() if isinstance(v, dict) and 'passed' in v)
E['P2'] = p2
# P4 within the new read set
P('\n**P4** within ' + NRS)
E['P4'] = {}
for arm, b in (('C4R4n', 'A'), ('D4R4n', 'B')):
    ub = paired(N[arm], N[b], body, absent)
    un = paired(N[arm], N[b], lambda r: T(r['named4']), named_sel)
    E['P4'][f'{arm}-{b}'] = dict(body=ub, named4=un)
    P(f"- {arm} minus {b}: gene body {fmt(ub)}; four-field named {fmt(un)}")
    ub0 = paired(O[arm], O[b], body, absent)
    E['P4'][f'{arm}-{b}_old_reads'] = dict(body=ub0)
    P(f"  (same on {ORS}: gene body {fmt(ub0)})")
# P5 novel calls on IPD-absent copies right / wrong per gene
P('\n**P5** novel final calls on IPD-absent copies, right gene body / not, per gene (old -> new)')
P('| gene | ' + ' | '.join(['C4', 'C4R4n', 'D4', 'D4R4n']) + ' |')
P('|---|---|---|---|---|')
E['P5'] = {}
for g in GENES + ('all',):
    cells = []
    for arm in ('C4', 'C4R4n', 'D4', 'D4R4n'):
        def rw(X):
            rr = [r for r in X[arm].values() if absent(r) and T(r['final_mos']) and (g == 'all' or r['gene'] == g)]
            return sum(right(r) for r in rr), sum(not right(r) for r in rr)
        o, n = rw(O), rw(N)
        E['P5'][f'{g}:{arm}'] = dict(old=o, new=n)
        cells.append(f'{o[0]}/{o[1]} -> {n[0]}/{n[1]}')
    P(f'| {g} | ' + ' | '.join(cells) + ' |')
# P6 DRB1 and DR53-group carriers
P('\n**P6** DRB1 only; DR53 group = donors whose old arm-A DRB1 call has *04, *07 or *09 in a slot')
dr53 = set()
for (d, g, h), r in O['A'].items():
    if g == 'DRB1' and any(x.startswith(('HLA-DRB1*04', 'HLA-DRB1*07', 'HLA-DRB1*09', 'DRB1*04', 'DRB1*07', 'DRB1*09'))
                           for x in r['call'].replace('|', ',').split(',')):
        dr53.add(d)
E['P6'] = dict(dr53_donors=len(dr53))
P(f'DR53-group donors: {len(dr53)}')
P('| subset | arm | IPD-absent | body old -> new | byte-exact final old -> new | mean wrong bases old -> new (copies) | four-field named old -> new |')
P('|---|---|---|---|---|---|---|')
for sub, sel in (('DRB1 all', lambda d: True), ('DRB1 DR53', lambda d: d in dr53), ('DRB1 non-DR53', lambda d: d not in dr53)):
    for arm in ('A', 'B', 'C4', 'C4R4n', 'D4', 'D4R4n'):
        def st(X):
            rr = [r for (d, g, h), r in X[arm].items() if g == 'DRB1' and sel(d)]
            ab = [r for r in rr if absent(r)]
            eds = [ed(r) for r in rr if ed(r) is not None]
            return (len(ab), sum(body(r) for r in ab), sum(T(r['recovered_final']) for r in ab),
                    (sum(eds) / len(eds) if eds else float('nan')), len(eds),
                    sum(T(r['named4']) for r in rr if named_sel(r)), sum(1 for r in rr if named_sel(r)))
        o, n = st(O), st(N)
        E['P6'][f'{sub}:{arm}'] = dict(old=o, new=n)
        P(f'| {sub} | {arm} | {n[0]} | {o[1]} -> {n[1]} | {o[2]} -> {n[2]} | {o[3]:.2f} -> {n[3]:.2f} ({n[4]}) | {o[5]} -> {n[5]} of {n[6]} |')
# blind regions (truth-free haplotype tables)
bl = {}
for lab, p, rs in (('old', ohaps, ORS), ('new', nhaps, NRS)):
    c = collections.Counter()
    for h in csv.DictReader(open(p), delimiter='\t'):
        if h['readset'] == rs and h['n_blind'] not in ('', 'None'):
            c[(h['gene'], h['arm'])] += int(h['n_blind'])
    bl[lab] = c
P('\nBlind regions (sum over haplotypes), old -> new: ' + ', '.join(
    f"{g}/{a}: {bl['old'].get((g, a), 0)} -> {bl['new'].get((g, a), 0)}" for g in GENES for a in ('C', 'D')))
E['P6']['blind_regions'] = {f'{g}/{a}': [bl['old'].get((g, a), 0), bl['new'].get((g, a), 0)] for g in GENES for a in ('C', 'D')}
# pairs added per donor
added = {}
for p in sorted(glob.glob(str(K.READS[NRS] / '*/PAIRS.json'))):
    j = json.load(open(p))
    added[j['donor']] = [j['old_regional_pairs'], j['added_pairs']]
E['P6']['pairs_old_added'] = added
if added:
    v = sorted(a for _, a in added.values())
    P(f'Pairs added per donor: median {v[len(v) // 2]}, range {v[0]}-{v[-1]}, total {sum(v)} over {len(v)} donors')
# 2-/3-field arm A per gene
P('\n**Arm A 2-/3-field correct per gene (all copies), old -> new; decoded (strict)**')
P('| gene | copies | 3-field | 2-field |')
P('|---|---|---|---|')
E['field23_A'] = {}
for g in GENES:
    o = [r for (d, gg, h), r in O['A'].items() if gg == g]
    nn = {k: r for k, r in N['A'].items() if k[1] == g}
    f = lambda rr, c: sum(T(r[c]) for r in rr)
    n_ = list(nn.values())
    E['field23_A'][g] = {c: [f(o, c), f(n_, c)] for c in ('named3_decoded', 'named3_strict', 'named2_decoded', 'named2_strict')}
    P(f"| {g} | {len(n_)} | {f(o, 'named3_decoded')} -> {f(n_, 'named3_decoded')} ({f(o, 'named3_strict')} -> {f(n_, 'named3_strict')}) | "
      f"{f(o, 'named2_decoded')} -> {f(n_, 'named2_decoded')} ({f(o, 'named2_strict')} -> {f(n_, 'named2_strict')}) |")
res = K.ST / 'ag/results'
(res / f'SUMMARY_{outtag}.txt').write_text('\n'.join(out) + '\n')
json.dump(E, open(res / f'endpoints_{outtag}.json', 'w'), indent=1, default=str)
print('\n'.join(out))
