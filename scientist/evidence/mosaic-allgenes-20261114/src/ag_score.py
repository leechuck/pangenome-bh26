"""DEVELOPMENT scorer for the all-gene builder (stage mosaic-allgenes-20261114). Generalises
drb1_score.py (frozen scorer score.py underneath) to every gene and adds, per copy, the edit distance
(edlib, global) between the final call's gene body and the true gene body. Numbers/names only are
written; sequences stay in memory on IBEX.

  python ag_score.py TAG READSETS ARMS        e.g. main D-sim,D-real A,B,C4,D4,C4R1,D4R1
ARMS: A, B (existing T1K runs); C4[Rk[n]] / D4[Rk[n]]: pass-1 table of arm A / B plus the builder's
overrides of every gene (OVERRIDES.json, or OVERRIDES_Rk[n].json for an abstain rule).
"""
from __future__ import annotations

import collections
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L          # noqa: E402
C = L.C
import score as SC            # noqa: E402

S2 = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114')
TEST = ('T-sim', 'T-real')
FALLBACK = set()


def is_new(arm):
    return arm[:2] in ('C4', 'D4')


def ov_file(arm):
    return 'OVERRIDES.json' if arm in ('C4', 'D4') else f'OVERRIDES_{arm[2:]}.json'


class Index(SC.NovelIndex):
    def __init__(self, maps):
        super().__init__(maps)
        self.body = {}
        self.body_seq = {}

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
            self.body_seq[name] = seq[ex[0][0]:ex[-1][1]]


def task_dir(tag, rs, d, arm):
    return S2 / 'work' / tag / rs / d / arm[0]


def new_arm_table(tag, rs, d, arm):
    w = task_dir(tag, rs, d, arm)
    if not (w / 'COMPLETE.json').exists():
        if (w / 'FAILURE.json').exists():
            f = C.load_json(w / 'FAILURE.json')
            if f.get('fallback_calls_table'):      # explicit fallback: pass-1 calls, no overrides
                return Path(f['fallback_calls_table']), [Path(f['fallback_pass1_reference'])] if arm[0] == 'D' else []
        # no record at all, or a failure without a fallback table (e.g. killed by Slurm after the retry):
        # the pass-1 calls stand (explicit fallback, counted as a failure by ag_test_summary.py)
        import run_donor as RD
        spec = RD.pass1_spec(rs, d, arm[0])
        FALLBACK.add((rs, d, arm[0]))
        return C.verified_table(spec['run']), [Path(spec['ref'])] if arm[0] == 'D' else []
    fas, run = [], None
    for g in C.GENES:
        gc = w / g / 'COMPLETE.json'
        if gc.exists():
            rec = C.load_json(gc)
            run = run or Path(rec['pass1_run'])
            if arm[0] == 'D':
                fas.append(Path(rec['pass1_reference']))
            fas.append(w / g / 'additions.fa')
    if run is None:
        return None, 'no_gene_complete'
    return C.verified_table(run), sorted(set(fas), key=str)


def new_post(tag, rs, d, arm, g, calls):
    o = task_dir(tag, rs, d, arm) / g / ov_file(arm)
    if not o.exists():
        return calls
    calls = dict(calls)
    slots = [list(s) for s in calls.get(g, [[], []])]
    for ov in C.load_json(o)['overrides']:
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
    import edlib
    tag, readsets, arms = sys.argv[1], sys.argv[2].split(','), sys.argv[3].split(',')
    assert not set(readsets) & set(TEST) or (S2 / 'FREEZE.json').exists(), 'test read sets only after FREEZE.json'
    T = SC.Truth(readsets)
    body_of, body_seq = {}, {}
    for r in T.records:
        if r['donor'] not in {d for (rs, d, g) in T.copies}:
            continue
        ex = r['core_exons']
        b = r['core_sequence'].upper()[ex[0][0]:ex[-1][1]]
        body_of[(r['donor'], r['gene'], str(r['haplotype']))] = C.sha256_str(b)
        body_seq[(r['donor'], r['gene'], str(r['haplotype']))] = b
    # named allele -> gene body
    named_body = {}
    for g, msa in T.msas.items():
        for i, names in enumerate(msa.ipd_names):
            if names and msa.ipd_exons[i]:
                ex = msa.ipd_exons[i]
                for n in names:
                    named_body[n] = msa.seqs[i][ex[0][0]:ex[-1][1]]
    nidx = Index(T.maps)
    ed_cache = {}

    def ed(label, tkey):
        s = nidx.body_seq.get(label) if SC.is_novel(label) else named_body.get(label)
        if s is None:
            return None
        k = (label, tkey)
        if k not in ed_cache:
            ed_cache[k] = edlib.align(s, body_seq[tkey], mode='NW', task='distance')['editDistance']
        return ed_cache[k]
    rows, avail = [], []
    for (rs, d, g), copies in sorted(T.copies.items()):
        copies = sorted(copies, key=lambda c: c['haplotype'])
        for arm in arms:
            if arm[0] in ('B', 'D') and rs == 'D-sim':
                continue
            key = (rs, d, arm)
            if key not in nidx.cache:
                try:
                    t, fas = new_arm_table(tag, rs, d, arm) if is_new(arm) else SC.arm_table(rs, d, arm)
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
                    avail.append((rs, d, arm, 'fallback_pass1' if is_new(arm) and (rs, d, arm[0]) in FALLBACK else 'complete'))
            calls = nidx.cache[key]
            complete = calls is not None
            gene_status = ''
            if complete and is_new(arm):
                w = task_dir(tag, rs, d, arm) / g
                gene_status = 'complete' if (w / 'COMPLETE.json').exists() else \
                    'failed' if (w / 'FAILURE.json').exists() else 'missing'
                calls = new_post(tag, rs, d, arm, g, calls)
            slots = calls.get(g, [[], []]) if complete else [[], []]
            gene_bodies = {nidx.body.get(l) for s in slots for l in s if SC.is_novel(l)} - {None}
            # edit distance: own copy-to-slot assignment minimising the summed distance (the frozen
            # scorer's assignment maximises exact/named matches and is arbitrary when nothing matches)
            firsts = [(SC.typed(s_)[0] if s_ else None) for s_ in slots]
            ed_of = {}
            best = None
            for perm in itertools.permutations(range(len(slots)), len(copies)):
                v = [ed(firsts[p], (d, g, str(c['haplotype']))) if firsts[p] else None for c, p in zip(copies, perm)]
                tot = sum(x if x is not None else 10 ** 6 for x in v)
                if best is None or tot < best[0]:
                    best = (tot, v)
            if best:
                ed_of = {str(c['haplotype']): x for c, x in zip(copies, best[1])}
            for c, r in zip(copies, SC.score_gene(copies, slots, nidx)):
                tkey = (d, g, str(c['haplotype']))
                first = r['call'].split('|')[0] if r['call'] != '.' else None
                rows.append(dict(readset=rs, donor=d, gene=g, haplotype=c['haplotype'], arm=arm,
                                 ipd_absent=c['ipd_absent'], truth4=c['truth4'] or '', class_full=c['class_full'],
                                 class_ipd=c['class_ipd'], arm_complete=complete, gene_status=gene_status,
                                 body_recovered=body_of.get(tkey) in gene_bodies,
                                 final_body_ok=(first is not None and SC.is_novel(first)
                                                and nidx.body.get(first) == body_of.get(tkey)),
                                 ed_final=ed_of.get(str(c['haplotype'])), **r))
    res = S2 / 'results'
    res.mkdir(exist_ok=True)
    cols = list(rows[0])
    with open(res / f'score_rows_{tag}.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    with open(res / f'score_avail_{tag}.tsv', 'w') as f:
        for a in avail:
            f.write('\t'.join(a) + '\n')
    print(len(rows), 'rows;', sum(a[3] != 'complete' for a in avail), 'incomplete arm tables')


if __name__ == '__main__':
    main()
