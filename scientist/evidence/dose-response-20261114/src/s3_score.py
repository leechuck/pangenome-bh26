"""Scoring for the dose-response stage, using the frozen scorer's rules (mosaic-typer score.py) unchanged:
slot assignment by best permutation, exact = byte equality of the final call with the truth core,
named4 for IPD-present copies with one four-field truth label, family-clustered bootstrap (4000, seed 20261020).
Additions: gene-body exactness (first exon start to last exon end), supply by the panel, conversion.

  python s3_score.py --tag pilot --manifest work/manifest_pilot.json
  python s3_score.py --tag full  --manifest work/manifest_full.json
Writes results/dose_rows_<tag>.tsv, dose_summary_<tag>.tsv (long), dose_contrasts_<tag>.tsv, dose_by_gene_<tag>.tsv.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s3_lib as L  # noqa: E402

C = L.C
sys.path.insert(0, str(L.MT / 'src'))
import score as SC  # noqa: E402  frozen scorer (imported read-only)
import itertools  # noqa: E402

R = L.STAGE / 'results'


class Index(SC.NovelIndex):
    """NovelIndex plus gene-body hash of each novel candidate."""

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


def score_gene2(copies, slots, nidx):
    """Frozen SC.score_gene logic (same permutation search and weights) plus gene-body fields."""
    ev = [SC.slot_eval(s, nidx) for s in slots]
    gene_shas = set().union(*(e['all_shas'] for e in ev)) if ev else set()
    gene_bodies = {nidx.body[l] for s in slots for l in s if SC.is_novel(l) and l in nidx.body}
    best = None
    for perm in itertools.permutations(range(len(ev)), len(copies)):
        ms = [SC.copy_eval(c, ev[p]) for c, p in zip(copies, perm)]
        sc = sum(1000 * m['exact'] + 100 * m['named4'] + 10 * m['named3_decoded'] + m['named2_decoded']
                 for m in ms)
        if best is None or sc > best[0]:
            best = (sc, perm, ms)
    rows = []
    for c, p, m in zip(copies, best[1], best[2]):
        s = ev[p]
        fin_bodies = {nidx.body[l] for l in s['typed'] if SC.is_novel(l) and l in nidx.body}
        rows.append(dict(m, recovered=c['sha'] in gene_shas, body_final=c['body'] in fin_bodies,
                         body_any=c['body'] in gene_bodies, no_call=not s['labels'],
                         call='|'.join(s['typed']) if s['typed'] else '.'))
    return rows


def boot_share(clusters, reps=4000, seed=20261020):
    """clusters: list of (k, n). Family-clustered percentile bootstrap of sum k / sum n (percent)."""
    clusters = [c for c in clusters if c[1] > 0]
    n = sum(c[1] for c in clusters)
    if not n:
        return None, None, None
    k = sum(c[0] for c in clusters)
    rng = random.Random(seed)
    b = []
    for _ in range(reps):
        d = [rng.choice(clusters) for _ in clusters]
        dn = sum(x[1] for x in d)
        if dn:
            b.append(100 * sum(x[0] for x in d) / dn)
    b.sort()
    return 100 * k / n, b[int(0.025 * len(b))], b[int(0.975 * len(b)) - 1]


def load_tsv(p):
    return list(csv.DictReader(open(p), delimiter='\t'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', required=True)
    ap.add_argument('--manifest', required=True)
    a = ap.parse_args()
    man = C.load_json(a.manifest)
    donors = {r['donor']: r for r in load_tsv(R / 'donors.tsv')}
    cells = {(r['donor'], r['cell']): r for r in load_tsv(R / 'cells.tsv')}
    supply = {(r['donor'], r['cell'], r['gene'], r['copy']): r for r in load_tsv(R / 'supply_rows.tsv')}
    cellsets = {t['donor']: [c for c in t['cells'] if not c.startswith('REPRO')] for t in man['tasks']}
    T = SC.Truth(['T-sim', 'T-real', 'D-real'])
    nidx = Index(T.maps)
    body_of = {}
    for r in T.records:
        ex = r['core_exons']
        body_of[(r['donor'], r['gene'], str(r['haplotype']))] = C.sha256_str(
            r['core_sequence'].upper()[ex[0][0]:ex[-1][1]])
    runs = L.donor_runs()
    rows, avail = [], []
    for donor in sorted(cellsets):
        d = donors[donor]
        rs = d['readset']
        for cell in ['n0'] + cellsets[donor] + ['B_existing']:
            key = (donor, cell)
            try:
                if cell == 'n0':
                    table, fas = C.verified_table(Path(runs[donor]['armA_run'])), []
                    meta = dict(source='IPD', n_requested=0, n_drawn=0)
                elif cell == 'B_existing':
                    table, fas = SC.arm_table(rs, donor, 'B')
                    meta = dict(source='B_existing', n_requested='', n_drawn='')
                else:
                    w = L.STAGE / 'work' / donor / cell
                    rec = C.load_json(w / 'COMPLETE.json')
                    assert C.sha256_file(w / 'calls_genotype.tsv') == rec['table_sha256']
                    table, fas = w / 'calls_genotype.tsv', [w / 'additions.fa']
                    meta = dict(source=rec['source'], n_requested=rec['n_requested'], n_drawn=rec['n_drawn'])
                for fa in fas:
                    if Path(fa).exists():
                        nidx.add_fasta(fa)
                calls = C.read_t1k(table)
                state = 'complete'
            except Exception as e:  # noqa: BLE001
                calls, state = None, 'missing_or_failed: ' + repr(e)[:200]
                meta = dict(source=cells[key]['source'] if key in cells else '', n_requested='', n_drawn='')
            avail.append((donor, cell, state))
            for g in C.GENES:
                copies = sorted(T.copies[(rs, donor, g)], key=lambda c: c['haplotype'])
                for c in copies:
                    c['body'] = body_of[(donor, g, str(c['haplotype']))]
                slots = calls.get(g, [[], []]) if calls is not None else [[], []]
                for c, r in zip(copies, score_gene2(copies, slots, nidx)):
                    sp = supply.get((donor, cell, g, str(c['haplotype'])), {})
                    rows.append(dict(donor=donor, group=d['group'], readtype=d['readtype'], readset=rs, cell=cell,
                                     source=meta['source'], n_requested=meta['n_requested'], n_drawn=meta['n_drawn'],
                                     gene=g, copy=c['haplotype'], ipd_absent=int(c['ipd_absent']),
                                     has_truth4=int(bool(c['truth4'])),
                                     supplied_byte=int(sp['supplied_byte']) if sp else (0 if cell == 'n0' else -1),
                                     supplied_body=int(sp['supplied_body']) if sp else (0 if cell == 'n0' else -1),
                                     call=r['call'], exact=int(r['exact']), recovered_any=int(r['recovered']),
                                     body_exact=int(r['body_final']), body_any=int(r['body_any']),
                                     named4=int(r['named4']), no_call=int(r['no_call']),
                                     complete=int(calls is not None)))
    cols = list(rows[0])
    with open(R / f'dose_rows_{a.tag}.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    C.save_json(R / f'dose_avail_{a.tag}.json', [dict(donor=x, cell=y, state=z) for x, y, z in avail
                                                 if z != 'complete'])
    print('rows', len(rows), 'incomplete', sum(1 for x in avail if x[2] != 'complete'))
    summarize(a.tag, rows, T.group_of)


# metric: name -> (row filter, row value)
METRICS = collections.OrderedDict([
    ('absent_exact_any', (lambda r: r['ipd_absent'], lambda r: r['recovered_any'])),
    ('absent_exact_final', (lambda r: r['ipd_absent'], lambda r: r['exact'])),
    ('absent_body_any', (lambda r: r['ipd_absent'], lambda r: r['body_any'])),
    ('absent_body_final', (lambda r: r['ipd_absent'], lambda r: r['body_exact'])),
    ('present_named4', (lambda r: (not r['ipd_absent']) and r['has_truth4'], lambda r: r['named4'])),
    ('all_correct', (lambda r: r['ipd_absent'] or r['has_truth4'],
                     lambda r: r['exact'] if r['ipd_absent'] else r['named4'])),
    ('conv_exact_any', (lambda r: r['ipd_absent'] and r['supplied_byte'] == 1, lambda r: r['recovered_any'])),
    ('conv_exact_final', (lambda r: r['ipd_absent'] and r['supplied_byte'] == 1, lambda r: r['exact'])),
    ('conv_body_any', (lambda r: r['ipd_absent'] and r['supplied_body'] == 1, lambda r: r['body_any'])),
    ('conv_body_final', (lambda r: r['ipd_absent'] and r['supplied_body'] == 1, lambda r: r['body_exact'])),
    ('supplied_byte', (lambda r: r['ipd_absent'] and r['cell'] != 'B_existing', lambda r: r['supplied_byte'])),
    ('supplied_body', (lambda r: r['ipd_absent'] and r['cell'] != 'B_existing', lambda r: r['supplied_body'])),
    ('absent_exact_any_unsupplied', (lambda r: r['ipd_absent'] and r['supplied_byte'] != 1 and r['cell'] != 'B_existing',
                                     lambda r: r['recovered_any'])),
])


def cast(rows):
    for r in rows:
        for k in ('ipd_absent', 'has_truth4', 'exact', 'recovered_any', 'body_exact', 'body_any', 'named4', 'no_call'):
            r[k] = int(r[k])
        for k in ('supplied_byte', 'supplied_body'):
            r[k] = int(r[k]) if r[k] not in ('', None) else -1
    return rows


def summarize(tag, rows, group_of):
    comp = lambda d: group_of.get(d, (d,))  # noqa: E731
    out, bygene = [], []
    keys = sorted({(r['group'], r['readtype'], r['cell']) for r in rows})
    for (grp, rt, cell) in keys:
        sel = [r for r in rows if r['group'] == grp and r['cell'] == cell]
        nr = [r['n_requested'] for r in sel if r['n_requested'] != '']
        nd = [int(r['n_drawn']) for r in sel if r['n_drawn'] != '']
        for m, (flt, val) in METRICS.items():
            cl = collections.defaultdict(lambda: [0, 0])
            for r in sel:
                if flt(r):
                    x = cl[comp(r['donor'])]
                    x[0] += int(val(r))
                    x[1] += 1
            pct, lo, hi = boot_share([tuple(v) for v in cl.values()])
            k = sum(v[0] for v in cl.values())
            n = sum(v[1] for v in cl.values())
            out.append(dict(group=grp, readtype=rt, cell=cell, source=sel[0]['source'],
                            n_requested=nr[0] if nr else '', mean_n_drawn=round(sum(nd) / len(nd), 1) if nd else '',
                            test_donors=len({r['donor'] for r in sel}), metric=m, k=k, n=n,
                            pct='' if pct is None else round(pct, 2), lo='' if lo is None else round(lo, 2),
                            hi='' if hi is None else round(hi, 2)))
        for g in C.GENES:
            for m in ('absent_exact_any', 'absent_body_any', 'present_named4', 'supplied_byte'):
                flt, val = METRICS[m]
                s2 = [r for r in sel if r['gene'] == g and flt(r)]
                if not s2:
                    continue
                k = sum(int(val(r)) for r in s2)
                bygene.append(dict(group=grp, readtype=rt, cell=cell, gene=g, metric=m, k=k, n=len(s2),
                                   pct=round(100 * k / len(s2), 2)))
    cols = list(out[0])
    for name, data in (('dose_summary', out), ('dose_by_gene', bygene)):
        c = list(data[0])
        with open(R / f'{name}_{tag}.tsv', 'w') as f:
            f.write('\t'.join(c) + '\n')
            for r in data:
                f.write('\t'.join(str(r[x]) for x in c) + '\n')
    # ---- paired contrasts (same copies in both cells), cluster bootstrap + sign flip from the frozen scorer
    idx = {(r['donor'], r['gene'], r['copy'], r['cell']): r for r in rows}
    cellset = collections.defaultdict(set)
    for r in rows:
        cellset[r['group']].add(r['cell'])
    contr = []
    for grp in sorted(cellset):
        cs = cellset[grp]
        pairs = []
        for n in (5, 10, 20, 40):
            if f'match_n{n:02d}' in cs:
                pairs.append((f'match_n{n:02d}', 'n0', 'matched vs IPD only'))
        for a_, b_ in ((10, 5), (20, 10), (40, 20), (40, 5), (40, 10)):
            if f'match_n{a_:02d}' in cs and f'match_n{b_:02d}' in cs:
                pairs.append((f'match_n{a_:02d}', f'match_n{b_:02d}', 'trend within matched'))
        for c_ in sorted(cs):
            if c_.startswith('x_'):
                n = int(c_[-2:])
                pairs.append((f'match_n{n:02d}', c_, 'matched vs cross-population, same n'))
        for (ca, cb, kind) in pairs:
            for m in ('absent_exact_any', 'absent_body_any', 'present_named4', 'all_correct'):
                flt, val = METRICS[m]
                prs = []
                for r in rows:
                    if r['group'] != grp or r['cell'] != ca or not flt(r):
                        continue
                    b = idx.get((r['donor'], r['gene'], r['copy'], cb))
                    if b is None:
                        continue
                    prs.append(dict(donor=r['donor'], difference=int(val(r)) - int(val(b))))
                if not prs:
                    continue
                u = SC.uncertainty(prs, group_of)
                contr.append(dict(group=grp, kind=kind, cell_a=ca, cell_b=cb, metric=m, n_pairs=len(prs),
                                  gains=sum(p['difference'] > 0 for p in prs), losses=sum(p['difference'] < 0 for p in prs),
                                  diff_pp=round(u['difference_pp'], 2), lo=round(u['bootstrap95_pp'][0], 2),
                                  hi=round(u['bootstrap95_pp'][1], 2), sign_flip_p=round(u['sign_flip_p'], 4)))
    if contr:
        c = list(contr[0])
        with open(R / f'dose_contrasts_{tag}.tsv', 'w') as f:
            f.write('\t'.join(c) + '\n')
            for r in contr:
                f.write('\t'.join(str(r[x]) for x in c) + '\n')


if __name__ == '__main__':
    main()
