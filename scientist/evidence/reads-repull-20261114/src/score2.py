"""Part B step 7: score old (D-real) vs re-pulled (D-real2) reads on the 20 development donors with the
logic of mosaic-drb1-20261113/src/drb1_score.py (rows and DRB1 table copied unchanged; only the
places where tables are found differ). Truth copies of D-real2 = those of D-real (same donors).
Adds a per-gene all-gene four-field table for arm A old vs new. Numbers and names only.

  python score2.py TAG
"""
from __future__ import annotations

import collections
import os
import sys
from pathlib import Path

sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113/src')
import drb1_score as DS       # noqa: E402
C, SC = DS.C, DS.SC

ST = Path('/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114')
S2 = DS.S2
G = 'DRB1'
NEW = DS.NEW
RS = os.environ.get('RS', 'D-real2')


def arm_table(rs, d, arm):
    if rs == RS:
        if arm == 'A':
            return C.verified_table(ST / 'work' / RS / d / 'A'), []
        DS.S2 = ST
        try:
            return DS.new_arm_table(rs, d, arm)
        finally:
            DS.S2 = S2
    return DS.new_arm_table(rs, d, arm) if arm in NEW else SC.arm_table(rs, d, arm)


def new_post(rs, d, arm, g, calls):
    DS.S2 = ST if rs == RS else S2
    try:
        return DS.new_post(rs, d, arm, g, calls)
    finally:
        DS.S2 = S2


def main():
    tag = sys.argv[1]
    assert os.environ.get('SLURM_JOB_ID')
    readsets, arms = ['D-real', RS], ['A', 'C3', 'C4']
    T = SC.Truth(['D-real'])
    for (rs, d, g), cps in list(T.copies.items()):
        T.copies[(RS, d, g)] = [dict(c, readset=RS) for c in cps]
    body_of = {}
    for r in T.records:
        ex = r['core_exons']
        body_of[(r['donor'], r['gene'], str(r['haplotype']))] = C.sha256_str(r['core_sequence'].upper()[ex[0][0]:ex[-1][1]])
    nidx = DS.Index(T.maps)
    rows, avail = [], []
    for (rs, d, g), copies in sorted(T.copies.items()):
        copies = sorted(copies, key=lambda c: c['haplotype'])
        for arm in arms:
            key = (rs, d, arm)
            if key not in nidx.cache:
                try:
                    t, fas = arm_table(rs, d, arm)
                except Exception as e:
                    t, fas = None, repr(e)
                if t is None:
                    nidx.cache[key] = None
                    avail.append((rs, d, arm, str(fas)))
                else:
                    for fa in fas:
                        if Path(fa).exists():
                            nidx.add_fasta(fa)
                    nidx.cache[key] = C.read_t1k(t)
                    avail.append((rs, d, arm, 'complete'))
            calls = nidx.cache[key]
            complete = calls is not None
            if complete and arm in NEW:
                calls = new_post(rs, d, arm, g, calls)
            slots = calls.get(g, [[], []]) if complete else [[], []]
            gene_bodies = {nidx.body.get(l) for s in slots for l in s if SC.is_novel(l)} - {None}
            for c, r in zip(copies, SC.score_gene(copies, slots, nidx)):
                rows.append(dict(readset=rs, donor=d, gene=g, haplotype=c['haplotype'], arm=arm,
                                 ipd_absent=c['ipd_absent'], truth4=c['truth4'] or '', class_full=c['class_full'],
                                 class_ipd=c['class_ipd'], arm_complete=complete,
                                 body_recovered=body_of.get((d, g, str(c['haplotype']))) in gene_bodies, **r))
    res = ST / 'results'
    res.mkdir(exist_ok=True)
    cols = list(rows[0])
    with open(res / f'score_rows_{tag}.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    out = []
    P = out.append
    P(f'Re-pull scores, tag {tag} (20 development donors; D-real = old reads, {RS} = re-pulled reads)\n')
    P('incomplete arms: ' + (', '.join(f'{a[0]}/{a[1]}/{a[2]}:{a[3]}' for a in avail if a[3] != 'complete') or 'none') + '\n')
    for rs in readsets:
        P(f'\n**Read set {rs}**\n')
        P('| arm | DRB1 IPD-absent copies | exact (any slot) | exact (final call) | gene body exact | mosaic-class: exact / body / n | other absent: exact / body / n | DRB1 IPD-present 4-field correct | novel final calls on IPD-present DRB1 copies | novel final calls on absent copies: exact / body only / incomplete | DRB1 3-field / 2-field correct (all copies) | all genes 4-field correct |')
        P('|---|---:|---:|---:|---:|---|---|---|---:|---|---|---|')
        for arm in arms:
            R = [r for r in rows if r['readset'] == rs and r['arm'] == arm]
            Dr = [r for r in R if r['gene'] == G]
            ab = [r for r in Dr if r['ipd_absent']]
            mo = [r for r in ab if r['class_ipd'] == 'mosaic']
            ot = [r for r in ab if r['class_ipd'] != 'mosaic']
            pr = [r for r in Dr if not r['ipd_absent'] and r['truth4']]
            al = [r for r in R if not r['ipd_absent'] and r['truth4']]
            P(f"| {arm} | {len(ab)} | {sum(r['recovered'] for r in ab)} | {sum(r['recovered_final'] for r in ab)} | "
              f"{sum(r['body_recovered'] or r['recovered'] for r in ab)} | "
              f"{sum(r['recovered'] for r in mo)} / {sum(r['body_recovered'] or r['recovered'] for r in mo)} / {len(mo)} | "
              f"{sum(r['recovered'] for r in ot)} / {sum(r['body_recovered'] or r['recovered'] for r in ot)} / {len(ot)} | "
              f"{sum(r['named4'] for r in pr)}/{len(pr)} | {sum(r['final_mos'] for r in Dr if not r['ipd_absent'])} | "
              f"{sum(r['final_mos'] and r['recovered_final'] for r in ab)} / "
              f"{sum(r['final_mos'] and not r['recovered_final'] and r['body_recovered'] for r in ab)} / "
              f"{sum(r['final_mos'] and not r['recovered_final'] and not r['body_recovered'] for r in ab)} | "
              f"{sum(r['named3_decoded'] for r in Dr)} / {sum(r['named2_decoded'] for r in Dr)} of {len(Dr)} | "
              f"{sum(r['named4'] for r in al)}/{len(al)} |")
        for arm in arms:
            if arm == 'A':
                continue
            base = {(r['donor'], r['gene'], r['haplotype']): r for r in rows if r['readset'] == rs and r['arm'] == 'A'}
            R = {(r['donor'], r['gene'], r['haplotype']): r for r in rows if r['readset'] == rs and r['arm'] == arm}
            _diffs(P, T, R, base, f'{arm} minus A ({rs})')
    # new reads vs old reads, same arm
    P(f'\n**{RS} minus D-real, same arm**\n')
    for arm in arms:
        base = {(r['donor'], r['gene'], r['haplotype']): r for r in rows if r['readset'] == 'D-real' and r['arm'] == arm}
        R = {(r['donor'], r['gene'], r['haplotype']): r for r in rows if r['readset'] == RS and r['arm'] == arm}
        _diffs(P, T, R, base, f'{arm}: {RS} minus D-real')
    P('\n**Arm A four-field named correct per gene (IPD-present copies with a four-field truth name)**\n')
    P(f'| gene | copies | D-real (old reads) | {RS} (new reads) | gained | lost |')
    P('|---|---:|---:|---:|---:|---:|')
    for g in C.GENES:
        o = {(r['donor'], r['haplotype']): r for r in rows if r['readset'] == 'D-real' and r['arm'] == 'A' and r['gene'] == g and not r['ipd_absent'] and r['truth4']}
        n = {(r['donor'], r['haplotype']): r for r in rows if r['readset'] == RS and r['arm'] == 'A' and r['gene'] == g and not r['ipd_absent'] and r['truth4']}
        P(f"| {g} | {len(o)} | {sum(r['named4'] for r in o.values())} | {sum(r['named4'] for r in n.values())} | "
          f"{sum(n[k]['named4'] and not o[k]['named4'] for k in o if k in n)} | {sum(o[k]['named4'] and not n[k]['named4'] for k in o if k in n)} |")
    (res / f'SCORES_{tag}.md').write_text('\n'.join(out) + '\n')
    print('\n'.join(out))


def _diffs(P, T, R, base, label0):
    for field, sel, label in (('recovered', lambda r: r['gene'] == G and r['ipd_absent'], 'DRB1 exact recovery'),
                              ('body', lambda r: r['gene'] == G and r['ipd_absent'], 'DRB1 gene-body recovery'),
                              ('named4', lambda r: r['gene'] == G and not r['ipd_absent'] and r['truth4'], 'DRB1 4-field named'),
                              ('named4', lambda r: not r['ipd_absent'] and r['truth4'], 'all genes 4-field named')):
        pairs = []
        for k, r in R.items():
            if k not in base or not sel(r):
                continue
            val = (lambda x: (x['body_recovered'] or x['recovered'])) if field == 'body' else (lambda x: x[field])
            pairs.append(dict(donor=r['donor'], difference=int(val(r)) - int(val(base[k]))))
        if not pairs:
            continue
        u = SC.uncertainty(pairs, T.group_of)
        P(f"- {label0}, {label}: {sum(p['difference'] for p in pairs):+d} of {len(pairs)} copies; "
          f"{u['difference_pp']:+.1f} pp, 95% CI [{u['bootstrap95_pp'][0]:+.1f}, {u['bootstrap95_pp'][1]:+.1f}], "
          f"sign-flip p = {u['sign_flip_p']:.3g}")


if __name__ == '__main__':
    main()
