"""DEVELOPMENT ONLY: why do sampled haplotypes differ from the truth? For every dip2 configuration of one
donor, align picked gene body to the true gene body and say where the wrong bases are: inside read holes
(true-body stretches with no read k-mer at all in that read input) or in covered sequence; deletion vs
insertion vs substitution. Output work/<rs>/<donor>/diag_rows.tsv and holes.tsv (counts only).
  python gp_diag.py READSET DONOR"""
import collections
import re
import sys
from pathlib import Path

import edlib

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gp_lib as G   # noqa: E402
import gp_donor as GD   # noqa: E402
C = G.C


def covered_mask(body, cnt, k=29):
    n = len(body)
    cov = [0] * n
    for i in range(n - k + 1):
        x = body[i:i + k]
        y = G.rc(x)
        if cnt.get(x if x < y else y, 0) > 0:
            for j in range(i, i + k):
                cov[j] = 1
    return cov


def holes(cov, min_len=50):
    out, i = [], 0
    while i < len(cov):
        if not cov[i]:
            j = i
            while j < len(cov) and not cov[j]:
                j += 1
            if j - i >= min_len:
                out.append((i, j))
            i = j
        else:
            i += 1
    return out


def main():
    rs, donor = sys.argv[1], sys.argv[2]
    out = G.S / 'work' / rs / donor
    assert (out / 'COMPLETE.json').exists()
    records = C.load_records()
    msas = C.load_msas(records, genes=G.GENES3)
    fold = G.GT.fold_of(donor)
    spec = G.RD.pass1_spec(rs, donor, 'A')
    fq = [spec['run'] / 'calls_candidate_1.fq', spec['run'] / 'calls_candidate_2.fq']
    rows, hrows = [], []
    for g in G.GENES3:
        msa = msas[g]
        truth = []
        for r in sorted((r for r in records if r['donor'] == donor and r['gene'] == g), key=lambda r: str(r['haplotype'])):
            core, ex = r['core_sequence'].upper(), r['core_exons']
            bc = msa.base_cols(msa.row_of_source[r['source_id']])
            truth.append(dict(hap=r['haplotype'], body=core[ex[0][0]:ex[-1][1]], c0=bc[ex[0][0]], c1=bc[ex[-1][1] - 1],
                              absent=not r['genomic_labels']))
        tk = {i: G.kmers(t['body']) for i, t in enumerate(truth)}
        inputs = {'i': fq, 'ii': [out / g / 'ii_1.fq', out / g / 'ii_2.fq']}
        cov = {}
        for inp, files in inputs.items():
            cnt = GD.count_targets(files, tk)
            for i, t in enumerate(truth):
                cm = covered_mask(t['body'], cnt[i])
                cov[(inp, i)] = cm
                hs = holes(cm)
                hrows.append(dict(readset=rs, donor=donor, gene=g, copy=t['hap'], ipd_absent=t['absent'], input=inp,
                                  body_len=len(t['body']), uncovered_bases=len(cm) - sum(cm), holes_ge50=len(hs),
                                  hole_bases_ge50=sum(b - a for a, b in hs), holes=';'.join(f'{a}-{b}' for a, b in hs)))
        for variant in G.VARIANTS:
            gdir = G.graph_dir(g, variant, fold)
            node_seq, node_cols = G.load_nodes(gdir)
            wd = out / g / variant
            for gaf in sorted(wd.glob('*_dip2.gaf')):
                inp, Lc, _ = gaf.stem.split('_')
                haps = []
                for line in open(gaf):
                    if line.startswith('@') or not line.strip():
                        continue
                    f = line.rstrip('\n').split('\t')
                    haps.append(GD.hap_bases(G.parse_walk(f[5]), node_seq, node_cols))
                if not haps:
                    continue
                if len(haps) == 1:
                    haps = haps * 2
                B = [[GD.body_of(h, t['c0'], t['c1']) for h in haps] for t in truth]
                D = [[GD.ed(truth[i]['body'], b) for b in B[i]] for i in range(2)]
                pick = (0, 1) if D[0][0] + D[1][1] <= D[0][1] + D[1][0] else (1, 0)
                for i, t in enumerate(truth):
                    pb = B[i][pick[i]]
                    ctr = collections.Counter()
                    if pb != t['body']:
                        a = edlib.align(pb, t['body'], mode='NW', task='path')
                        tpos = 0
                        for n, op in re.findall(r'(\d+)([=XID])', a['cigar']):
                            n = int(n)
                            # query = picked, target = truth: I = extra base in picked, D = truth base missing in picked
                            if op in '=X':
                                if op == 'X':
                                    for p in range(tpos, tpos + n):
                                        ctr['sub_hole' if not cov[(inp, i)][p] else 'sub_cov'] += 1
                                tpos += n
                            elif op == 'D':
                                h = sum(1 for p in range(tpos, tpos + n) if not cov[(inp, i)][p])
                                ctr['del_hole'] += h
                                ctr['del_cov'] += n - h
                                tpos += n
                            else:
                                p = min(tpos, len(t['body']) - 1)
                                ctr['ins_hole' if not cov[(inp, i)][p] else 'ins_cov'] += n
                    rows.append(dict(readset=rs, donor=donor, gene=g, variant=variant, input=inp, subchain=Lc[1:],
                                     copy=t['hap'], ipd_absent=t['absent'], dist=D[i][pick[i]], body_len=len(t['body']),
                                     picked_len=len(pb), **{k: ctr[k] for k in ('sub_cov', 'sub_hole', 'del_cov', 'del_hole', 'ins_cov', 'ins_hole')}))
    for name, R in (('diag_rows.tsv', rows), ('holes.tsv', hrows)):
        cols = list(R[0])
        with open(out / name, 'w') as f:
            f.write('\t'.join(cols) + '\n')
            for r in R:
                f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    C.save_json(out / 'DIAG_COMPLETE.json', dict(status='complete', rows=len(rows)))
    print('done', rs, donor, len(rows))


if __name__ == '__main__':
    main()
