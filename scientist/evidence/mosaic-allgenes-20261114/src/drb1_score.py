"""Scorer for the DRB1 stage: the frozen stage's scorer (score.py) with one more arm (C3) and one
more column: gene-body recovery (the candidate equals the true copy from first-exon start to
last-exon end; allele ends are ignored). Numbers and names only.

  python drb1_score.py TAG READSETS ARMS      e.g.  dev D-sim,D-real A,C,E,C3
"""
from __future__ import annotations

import collections
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L          # noqa: E402
C = L.C
import score as SC            # noqa: E402

S2 = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113')
G = 'DRB1'
NEW = ('C3', 'C4', 'D3', 'D4')   # x3: T1K pass 2 decides; x4: the patched allele is the call
DIR = {'C3': 'C3', 'C4': 'C3', 'D3': 'D3', 'D4': 'D3'}   # C: IPD panel; D: IPD + cohort panel (on top of arm B)


class Index(SC.NovelIndex):
    def __init__(self, maps):
        super().__init__(maps)
        self.body = {}

    def add_fasta(self, path):
        if str(path) in self.cache:
            return
        self.cache[str(path)] = True
        for name, ints, seq in C.read_fasta_headers_and_seqs(path):
            if not SC.is_novel(name):
                continue
            g = C.gene_of_name(name)
            if g not in C.GENES:
                continue
            ex = C.header_exons(ints)
            cds = C.cds_of(seq, ex)
            cn, pn = self.maps[g]
            self.cache[name] = (C.sha256_str(seq), set(cn.get(cds, ())), set(pn.get(C.translate(cds), ())))
            self.body[name] = C.sha256_str(seq[ex[0][0]:ex[-1][1]])


def new_arm_table(rs, d, arm):
    w = S2 / 'work' / rs / d / DIR[arm]
    if not (w / 'COMPLETE.json').exists():
        return None, 'failed' if (w / 'FAILURE.json').exists() else 'missing'
    rec = C.load_json(w / 'COMPLETE.json')
    fas = [w / 'additions.fa'] + ([Path(rec['pass1_reference'])] if arm[0] == 'D' else [])
    if arm.endswith('4'):
        return C.verified_table(Path(rec['pass1_run'])), fas
    t = Path(rec['calls_table'])
    assert C.sha256_file(t) == rec['calls_table_sha256']
    return t, fas


def new_post(rs, d, arm, g, calls):
    """Acceptance rule (as the frozen scorer) and renames for the new arm."""
    w = S2 / 'work' / rs / d / DIR[arm]
    calls = dict(calls)
    a = w / 'ACCEPT.json'
    if arm.endswith('3') and a.exists():
        rej = {x['slot']: x for x in C.load_json(a)['genes'].get(g, []) if not x['accepted']}
        if rej:
            slots = [list(s) for s in calls.get(g, [[], []])]
            for i, x in rej.items():
                if i < len(slots):
                    slots[i] = [x['displaced']] + slots[i]
            calls[g] = slots
    o = w / 'OVERRIDES.json'
    if g == G and o.exists():
        slots = [list(s) for s in calls.get(g, [[], []])]
        for ov in C.load_json(o)['overrides']:
            if arm.endswith('3') and ov.get('kind') != 'rename':
                continue
            done = False
            for i, s in enumerate(slots):
                if s and ov['pass1_label'] in s and not any('*MOS_' in l for l in s):
                    if done and not ov.get('also_other_slot'):
                        continue
                    slots[i] = list(ov['new_labels'])
                    done = True
        calls[g] = slots
    return calls


def main():
    tag, readsets, arms = sys.argv[1], sys.argv[2].split(','), sys.argv[3].split(',')
    T = SC.Truth(readsets)
    body_of = {}
    for r in T.records:
        ex = r['core_exons']
        body_of[(r['donor'], r['gene'], str(r['haplotype']))] = C.sha256_str(r['core_sequence'].upper()[ex[0][0]:ex[-1][1]])
    nidx = Index(T.maps)
    rows, avail = [], []
    for (rs, d, g), copies in sorted(T.copies.items()):
        copies = sorted(copies, key=lambda c: c['haplotype'])
        for arm in arms:
            if arm in ('B', 'D', 'D3', 'D4') and rs == 'D-sim':
                continue
            key = (rs, d, arm)
            if key not in nidx.cache:
                try:
                    t, fas = new_arm_table(rs, d, arm) if arm in NEW else SC.arm_table(rs, d, arm)
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
            if complete and arm in ('C', 'D'):
                calls = SC.apply_acceptance(rs, d, arm, g, calls)
            if complete and arm in NEW:
                calls = new_post(rs, d, arm, g, calls)
            slots = calls.get(g, [[], []]) if complete else [[], []]
            gene_bodies = {nidx.body.get(l) for s in slots for l in s if SC.is_novel(l)} - {None}
            for c, r in zip(copies, SC.score_gene(copies, slots, nidx)):
                rows.append(dict(readset=rs, donor=d, gene=g, haplotype=c['haplotype'], arm=arm,
                                 ipd_absent=c['ipd_absent'], truth4=c['truth4'] or '', class_full=c['class_full'],
                                 class_ipd=c['class_ipd'], arm_complete=complete,
                                 body_recovered=body_of.get((d, g, str(c['haplotype']))) in gene_bodies, **r))
    res = S2 / 'results'
    res.mkdir(exist_ok=True)
    cols = list(rows[0])
    with open(res / f'score_rows_{tag}.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    out = []
    P = out.append
    P(f'DRB1 stage scores, tag {tag}\n')
    P('incomplete arms: ' + (', '.join(f'{a[0]}/{a[1]}/{a[2]}:{a[3]}' for a in avail if a[3] != 'complete') or 'none') + '\n')
    for rs in readsets:
        P(f'\n**Read set {rs}**\n')
        P('| arm | DRB1 IPD-absent copies | exact (any slot) | exact (final call) | gene body exact | mosaic-class: exact / body / n | other absent: exact / body / n | DRB1 IPD-present 4-field correct | novel final calls on IPD-present DRB1 copies | novel final calls on absent copies: exact / body only / incomplete | DRB1 3-field / 2-field correct (all copies) | all genes 4-field correct |')
        P('|---|---:|---:|---:|---:|---|---|---|---:|---|---|---|')
        for arm in arms:
            R = [r for r in rows if r['readset'] == rs and r['arm'] == arm]
            if not R:
                continue
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
            ref_arm = 'B' if (arm in ('D', 'D3', 'D4') and 'B' in arms) else 'A'
            base = {(r['donor'], r['gene'], r['haplotype']): r for r in rows if r['readset'] == rs and r['arm'] == ref_arm}
            R = {(r['donor'], r['gene'], r['haplotype']): r for r in rows if r['readset'] == rs and r['arm'] == arm}
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
                P(f"- {arm} minus {ref_arm}, {label}: {sum(p['difference'] for p in pairs):+d} of {len(pairs)} copies; "
                  f"{u['difference_pp']:+.1f} pp, 95% CI [{u['bootstrap95_pp'][0]:+.1f}, {u['bootstrap95_pp'][1]:+.1f}], "
                  f"sign-flip p = {u['sign_flip_p']:.3g}")
    (res / f'SCORES_{tag}.md').write_text('\n'.join(out) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
