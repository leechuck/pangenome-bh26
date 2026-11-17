"""DEVELOPMENT DONORS ONLY: per-copy gene-body distance of the existing arms (A, B, C4, D4) for the 20
development donors, genes DRB1, DQA1, B. Output results/arms_rows.tsv (names, distances, flags only).
  python gp_arms.py"""
import sys
from pathlib import Path

import edlib

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gp_lib as G   # noqa: E402
C = G.C
import drb1_score as DS   # noqa: E402
SC = DS.SC


def bodies_from(fa, out):
    for name, ints, seq in C.read_fasta_headers_and_seqs(fa):
        ex = C.header_exons(ints)
        if ex:
            out.setdefault(name, seq[ex[0][0]:ex[-1][1]])


def ed(a, b):
    if a == b:
        return 0
    if not a or not b:
        return max(len(a), len(b))
    return edlib.align(a, b, mode='NW', task='distance')['editDistance']


def main():
    records = C.load_records()
    classes = C.load_classes()
    dev = [tuple(l.split()) for l in open(G.E / 'mosaic-drb1-20261113/work/dev_tasks.txt') if l.strip()]
    body = {}
    bodies_from(C.GENOMIC_FA, body)
    rows, missing = [], []
    for rs, d in dev:
        for arm in ('A', 'B', 'C4', 'D4'):
            if arm in ('B', 'D4') and rs == 'D-sim':
                continue
            try:
                t, fas = DS.new_arm_table(rs, d, arm) if arm in ('C4', 'D4') else SC.arm_table(rs, d, arm)
            except Exception as e:
                t, fas = None, repr(e)
            if t is None:
                missing.append((rs, d, arm, str(fas)))
                continue
            for fa in fas:
                if Path(fa).exists():
                    bodies_from(fa, body)
            calls = C.read_t1k(t)
            if arm in ('C4', 'D4'):
                calls = DS.new_post(rs, d, arm, 'DRB1', calls)
            for g in G.GENES3:
                if arm in ('C4', 'D4') and g != 'DRB1':
                    continue
                truth = []
                for r in sorted((r for r in records if r['donor'] == d and r['gene'] == g), key=lambda r: str(r['haplotype'])):
                    core, ex = r['core_sequence'].upper(), r['core_exons']
                    truth.append((r, core[ex[0][0]:ex[-1][1]]))
                slots = calls.get(g, [[], []])
                cand = [s[0] if s else None for s in slots] + [None] * (2 - len(slots))
                cb = [body.get(x, '') if x else '' for x in cand[:2]]
                D = [[ed(tb, b) for b in cb] for _, tb in truth]
                pick = (0, 1) if D[0][0] + D[1][1] <= D[0][1] + D[1][0] else (1, 0)
                for i, (r, tb) in enumerate(truth):
                    lab = cand[pick[i]] or ''
                    rows.append(dict(readset=rs, donor=d, gene=g, arm=arm, copy=r['haplotype'],
                                     ipd_absent=not r['genomic_labels'], cls=classes[r['source_id']]['ipd_label_L150'],
                                     body_len=len(tb), label=lab, novel_label=SC.is_novel(lab) if lab else False,
                                     label_has_body=bool(lab) and lab in body, dist=D[i][pick[i]],
                                     exact=D[i][pick[i]] == 0))
    res = G.S / 'results'
    res.mkdir(exist_ok=True)
    cols = list(rows[0])
    with open(res / 'arms_rows.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    C.save_json(res / 'arms_COMPLETE.json', dict(status='complete', rows=len(rows), missing=missing))
    print(len(rows), 'rows; missing', missing)


if __name__ == '__main__':
    main()
