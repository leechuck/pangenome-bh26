"""Step 3: try the local-version patch builder on the development donors (truth for evaluation
only). Numbers and names only.   python drb1_diag3.py READSET DONOR"""
from __future__ import annotations

import itertools
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L          # noqa: E402
import drb1_patch as PT       # noqa: E402
import drb1_diag2 as D2       # noqa: E402
C, EV, PR = L.C, L.EV, L.PR
import run_donor as RD        # noqa: E402

S2 = D2.S2
G = 'DRB1'


def main():
    readset, donor = sys.argv[1], sys.argv[2]
    assert readset in ('D-sim', 'D-real'), 'development read sets only'
    assert os.environ.get('SLURM_JOB_ID'), 'run under Slurm'
    out = S2 / 'work/diag3' / f'{readset}_{donor}'
    out.mkdir(parents=True, exist_ok=True)
    X = D2.setup(readset, donor, out)
    panel, pu, msa = X['panel'], X['pu'], X['msa']
    sa, sb = X['sa'], X['sb']
    classes = C.load_classes()
    frz = PR.Params(**C.load_json(RD.FROZEN_PARAMS))
    proj = L.project_all(pu, X['rd_u'], 6)
    ev1 = L.make_evidence(panel.C, proj, 1, X['bad_q'], True)
    D = L.tables_quant(ev1)[3]['D']
    res = dict(readset=readset, donor=donor, pass1=[X['names'][sa[0]], X['names'][sb[0]]], D=D,
               same_family=X['same'], family_sizes=[len(X['fam_a']), len(X['fam_b'])], variants={})
    trecs = sorted([r for r in X['records'] if r['donor'] == donor and r['gene'] == G],
                   key=lambda r: str(r['haplotype']))
    tp = PR.Panel([r['source_id'] for r in trecs],
                  [msa.aligned[msa.row_of_source[r['source_id']]] for r in trecs])
    truth = []
    for i, r in enumerate(trecs):
        core, ex = r['core_sequence'].upper(), r['core_exons']
        cl = classes[r['source_id']]
        nz = np.nonzero(tp.P[i] != L.OUT)[0]
        truth.append(dict(codes=tp.P[i], sha=C.sha256_str(core), interior=core[ex[0][0]:ex[-1][1]],
                          label=cl['ipd_label_L150'], full=cl['full_label_L150'], absent=not r['genomic_labels'],
                          hap=str(r['haplotype']), start=int(nz[0]), end=int(nz[-1])))
    if len(truth) == 1:
        truth = truth * 2
    for vname, maxedit, near in (('E3fam', 3, None), ('E6fam', 6, None), ('E3near', 3, 0.004)):
        reads = [x for x, p in zip(X['rd_u'], proj)
                 if p[2] is not None and p[2] <= maxedit and x['qname'] not in X['bad_q']]
        haps = []
        fams = (X['fam_a'], X['fam_b'])
        p1 = (sa[0], sb[0])
        info = []
        for h in (0, 1):
            fam = list(fams[h])
            if near is not None:
                rate = L.subst_rate(panel.P, p1[h])
                fam = [r for r in fam if rate[r] <= near or r == p1[h]]
            fp = L.subpanel(panel, fam)
            ra = fam.index(p1[h])
            blocks = PT.build_blocks(fp)
            PT.count_blocks(blocks, reads)
            partner = panel.P[p1[1 - h]]
            patches, amb = PT.decide(blocks, fp, ra, partner, D, same_family=X['same'])
            path = PT.patched_path(fp, ra, patches, blocks)
            codes = path.codes(fp.P)
            exr = {i: msa.exons(X['rows'][r]) for i, r in enumerate(fam)}
            seq = PR.path_sequence(fp, path)
            exo = PR.path_exons(fp, path, exr)
            sp = PT.spell(fp, codes, frz) if patches else (0, 0)
            haps.append(dict(codes=codes, seq_sha=C.sha256_str(seq),
                             interior=(seq[exo[0][0]:exo[-1][1]] if exo else None), p1codes=fp.P[ra]))
            nz = np.nonzero(codes != L.OUT)[0]
            info.append(dict(family_rows=len(fam), blocks=len(blocks),
                             resolvable=sum(b['resolvable'] for b in blocks),
                             too_long=sum(b['longest'] > 110 for b in blocks),
                             substring_conflict=sum(b['substring_conflict'] for b in blocks),
                             n_patches=len(patches), ambiguous=amb[:10], n_ambiguous=len(amb),
                             patches=[{k: v for k, v in p.items() if k not in ('new',)} for p in patches][:30],
                             spell_mismatch=sp[0], spell_switches=sp[1], exons_ok=exo is not None,
                             start=int(nz[0]), end=int(nz[-1])))
        perm = min(itertools.permutations((0, 1)),
                   key=lambda p: sum(L.mism(haps[k]['p1codes'], truth[p[k]]['codes']) for k in (0, 1)))
        for k in (0, 1):
            t = truth[perm[k]]
            info[k].update(truth_hap=t['hap'], truth_label=t['label'], truth_full=t['full'], truth_absent=t['absent'],
                           mismatch_before=L.mism(haps[k]['p1codes'], t['codes']),
                           mismatch_after=L.mism(haps[k]['codes'], t['codes']),
                           byte_exact=haps[k]['seq_sha'] == t['sha'],
                           interior_exact=haps[k]['interior'] is not None and haps[k]['interior'] == t['interior'],
                           start_end_offset=[info[k]['start'] - t['start'], info[k]['end'] - t['end']])
        res['variants'][vname] = dict(reads=len(reads), haps=info)
        C.save_json(out / 'DIAG3.json', res)
    res['complete'] = True
    C.save_json(out / 'DIAG3.json', res)


if __name__ == '__main__':
    main()
