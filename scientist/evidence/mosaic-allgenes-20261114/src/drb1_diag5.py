"""Step 5: evaluate patch builder v3.1 on the development donors (truth for evaluation only).
  python drb1_diag5.py READSET DONOR"""
from __future__ import annotations

import itertools
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L          # noqa: E402
import drb1_patch2 as PT      # noqa: E402
import drb1_diag2 as D2       # noqa: E402
C, EV, PR = L.C, L.EV, L.PR
import run_donor as RD        # noqa: E402

S2 = D2.S2
G = 'DRB1'


def main():
    readset, donor = sys.argv[1], sys.argv[2]
    assert readset in ('D-sim', 'D-real'), 'development read sets only'
    assert os.environ.get('SLURM_JOB_ID'), 'run under Slurm'
    out = S2 / 'work/diag5' / f'{readset}_{donor}'
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
    reads = [x for x in X['rd'] if x['qname'] not in X['bad_q']]
    p1 = (sa[0], sb[0])
    fams = (X['fam_a'], X['fam_b'])
    for vname, kw in (('v31_e2', dict(max_edit=2)), ('v31_e1', dict(max_edit=1)),
                      ('v31_e2_strict', dict(max_edit=2, sup_share=0.35, min_sup=4))):
        haps, info = [], []
        for h in (0, 1):
            Hh = PT.build_haplotype(panel, fams[h], p1[h], p1[1 - h], reads, D, **dict(kw))
            fp, path = Hh['fp'], Hh['path']
            codes = path.codes(fp.P)
            exr = {i: msa.exons(X['rows'][r]) for i, r in enumerate(Hh['fam'])}
            seq = PR.path_sequence(fp, path)
            exo = PR.path_exons(fp, path, exr)
            sp = PT.spell(fp, codes, frz) if Hh['patches'] else (0, 0)
            haps.append(dict(codes=codes, seq_sha=C.sha256_str(seq), events=Hh['events'],
                             interior=(seq[exo[0][0]:exo[-1][1]] if exo else None), p1codes=fp.P[Hh['ra']]))
            nz = np.nonzero(codes != L.OUT)[0]
            ev = Hh['events']
            info.append(dict(family_rows=len(Hh['fam']), events=len(ev),
                             resolvable=sum(e['resolvable'] for e in ev),
                             too_long=sum(e['longest'] > 120 for e in ev),
                             substring_conflict=sum(e['substring_conflict'] for e in ev),
                             n_patches=len(Hh['patches']),
                             patches=[{k: v for k, v in p.items() if k not in ('new',)} for p in Hh['patches']][:30],
                             notes=Hh['notes'][:20], n_notes=len(Hh['notes']),
                             in_reference=seq in X['ref_by_gene'].get(G, set()),
                             spell_mismatch=sp[0], spell_switches=sp[1], exons_ok=exo is not None,
                             start=int(nz[0]), end=int(nz[-1])))
        perm = min(itertools.permutations((0, 1)),
                   key=lambda p: sum(L.mism(haps[k]['p1codes'], truth[p[k]]['codes']) for k in (0, 1)))
        for k in (0, 1):
            t = truth[perm[k]]
            # where do the remaining differences to the truth lie, relative to the events?
            diff = np.nonzero((haps[k]['codes'] != L.OUT) & (t['codes'] != L.OUT) & (haps[k]['codes'] != t['codes']))[0]
            dbg = []
            for c in diff[:12].tolist():
                hit = [e for e in haps[k]['events'] if e['cl'] <= c <= e['cr']]
                near = [e for e in hit if e['first'] <= c <= e['last']] or hit
                if not near:
                    dbg.append(dict(col=c, event=None))
                    continue
                e = near[0]
                ts = PT.row_string(t['codes'], e['cl'], e['cr'])
                dbg.append(dict(col=c, first=e['first'], last=e['last'], cl=e['cl'], cr=e['cr'], anchor=e['anchor'],
                                resolvable=e['resolvable'], conflict=e['substring_conflict'], longest=e['longest'],
                                n_versions=len(e['alleles']), truth_version_known=ts in e['alleles'],
                                truth_is_cur=ts == e['cur'], raw={str(a): b for a, b in e['raw'].items()},
                                cur_total=e['cur_total'],
                                sup_truth=e['sup'].get(ts), sup_lens={len(a): b for a, b in e['sup'].items()}))
            info[k]['debug_remaining'] = dbg
            info[k].update(truth_hap=t['hap'], truth_label=t['label'], truth_full=t['full'], truth_absent=t['absent'],
                           mismatch_before=L.mism(haps[k]['p1codes'], t['codes']),
                           mismatch_after=L.mism(haps[k]['codes'], t['codes']),
                           remaining_diff_cols=diff[:12].tolist(),
                           byte_exact=haps[k]['seq_sha'] == t['sha'],
                           interior_exact=haps[k]['interior'] is not None and haps[k]['interior'] == t['interior'],
                           start_end_offset=[info[k]['start'] - t['start'], info[k]['end'] - t['end']])
        res['variants'][vname] = dict(reads=len(reads), haps=info)
        C.save_json(out / 'DIAG5.json', res)
    res['complete'] = True
    C.save_json(out / 'DIAG5.json', res)


if __name__ == '__main__':
    main()
