"""Step 1: measure why DRB1 fails and try read/evidence/search variants. DEVELOPMENT DONORS ONLY
(D-sim, D-real). Truth is used for evaluation only, never by the search. Output: numbers and
allele names, no sequence.
  python drb1_diag.py READSET DONOR
"""
from __future__ import annotations

import itertools
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L          # noqa: E402
C, EV, PR = L.C, L.EV, L.PR
import panel_select as PS     # noqa: E402
import run_donor as RD        # noqa: E402

S2 = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113')
G = 'DRB1'
THETA = 0.01


def main():
    readset, donor = sys.argv[1], sys.argv[2]
    assert readset in ('D-sim', 'D-real'), 'development read sets only'
    assert os.environ.get('SLURM_JOB_ID'), 'run under Slurm'
    out = S2 / 'work/diag1' / f'{readset}_{donor}'
    out.mkdir(parents=True, exist_ok=True)
    t00 = time.time()
    C.verify_inputs()
    records = C.load_records()
    msa = C.load_msas(records, genes=(G,))[G]
    classes = C.load_classes()
    spec = RD.pass1_spec(readset, donor, 'C')
    table = C.verified_table(spec['run'])
    geno = C.read_t1k_raw(table)
    wanted, slots = set(), []
    for g, row in geno.items():
        for sl in (row['slots'][:row['copies']] if row['copies'] else []):
            if sl['abundance'] <= 0:
                continue
            wanted |= set(sl['labels'][:3])
            if g == G:
                slots.append(sl['labels'][:3])
    seqs, _ = RD.read_reference(spec['ref'], wanted, {G})
    target = {n: seqs[n] for n in sorted(wanted) if n in seqs}
    reads, mstats, gene_of_q = L.map_reads2(target, spec['run'] / 'calls_candidate_1.fq',
                                            spec['run'] / 'calls_candidate_2.fq', out / 'map', 4, RD.ENV)
    shutil.rmtree(out / 'map', ignore_errors=True)
    shas = [C.sha256_str(s) for s in msa.seqs]
    row_is_ipd = [bool(n) for n in msa.ipd_names]
    rows = PS.select_rows('C', row_is_ipd, shas, None)
    assert all(row_is_ipd[r] for r in rows)
    names = [msa.ipd_names[r][0] for r in rows]
    panel = PR.Panel(names, [msa.aligned[r] for r in rows])
    local = {r: i for i, r in enumerate(rows)}

    def label_row(lab):
        r = msa.row_of_allele.get(lab)
        return local.get(r) if r is not None else None
    p1 = [[r for r in (label_row(l) for l in labs) if r is not None] for labs in slots]
    p1 = [x for x in p1 if x]
    assert p1, 'no pass-1 DRB1 call'
    if len(p1) == 1:
        p1 = p1 * 2
    sa, sb = p1[0], p1[1]
    rd, unplaced = L.windows(panel, reads.get(G, []), label_row)
    bad_q = {x['qname'] for x in reads.get(G, []) if gene_of_q[x['qname']] - {G}}
    other = {}
    for x in reads.get(G, []):
        for g in gene_of_q[x['qname']] - {G}:
            other[g] = other.get(g, 0) + 1
    fam_a = sorted(set(L.family_rows(panel, sa[0], THETA)) | set(sa))
    fam_b = sorted(set(L.family_rows(panel, sb[0], THETA)) | set(sb))
    same = sb[0] in fam_a or sa[0] in fam_b
    if same:
        fam_a = fam_b = sorted(set(fam_a) | set(fam_b))
    U = sorted(set(fam_a) | set(fam_b))
    pu = L.subpanel(panel, U)
    uidx = {r: i for i, r in enumerate(U)}
    rd_u = [dict(x, prefer_row=uidx.get(x['prefer_row'])) for x in rd]
    t0 = time.time()
    proj = L.project_all(pu, rd_u, 6)
    t_proj = time.time() - t0
    t0 = time.time()
    ev0 = EV.build_evidence(panel, rd, 6)
    t_ev0 = time.time() - t0
    pa, pb = L.subpanel(panel, fam_a), L.subpanel(panel, fam_b)
    ia, ib = {r: i for i, r in enumerate(fam_a)}, {r: i for i, r in enumerate(fam_b)}
    st_a, st_b = [ia[r] for r in sa], [ib[r] for r in sb]
    ex_a = {i: msa.exons(rows[r]) for i, r in enumerate(fam_a)}
    ex_b = {i: msa.exons(rows[r]) for i, r in enumerate(fam_b)}

    # ---- truth, evaluation only
    trecs = sorted([r for r in records if r['donor'] == donor and r['gene'] == G],
                   key=lambda r: str(r['haplotype']))
    tp = PR.Panel([r['source_id'] for r in trecs], [msa.aligned[msa.row_of_source[r['source_id']]] for r in trecs])
    truth = []
    frz = PR.Params(**C.load_json(RD.FROZEN_PARAMS))
    for i, r in enumerate(trecs):
        tc = tp.P[i]
        cl = classes[r['source_id']]
        core, ex = r['core_sequence'].upper(), r['core_exons']
        info = dict(haplotype=str(r['haplotype']), ipd_absent=not r['genomic_labels'],
                    ipd_label=cl['ipd_label_L150'], full_label=cl['full_label_L150'],
                    truth4=(r['genomic_labels'] or [''])[0] if r['genomic_labels'] else '',
                    length=len(core))
        nz = np.nonzero(tc != L.OUT)[0]
        info['start_col'], info['end_col'] = int(nz[0]), int(nz[-1])
        # fewest-switch spelling of the truth from allowed rows
        T = np.where((np.arange(PR.NCODE)[None, :] != tc[:, None]), 1000.0, 0.0).astype(np.float32)
        T[:, L.OUT] = 0
        T[tc == L.OUT, :] = 0
        dec = {}
        best_fam = None
        for tag, pn in (('all', panel), ('fam_a', pa), ('fam_b', pb)):
            c, path = PR.viterbi(pn, T, frz, set(), True)
            dec[tag] = dict(mismatch=int(round(c // 1000)), switches=path.switches,
                            sources=[pn.names[r_] for r_, _, _ in path.segs][:8])
            c0, p0 = PR.viterbi(pn, T, frz, set(), False)
            dec[tag]['nearest_row'] = pn.names[p0.segs[0][0]]
            dec[tag]['nearest_mismatch'] = int(round(c0 // 1000))
            nzr = np.nonzero(pn.P[p0.segs[0][0]] != L.OUT)[0]
            dec[tag]['nearest_start_end_offset'] = [int(nzr[0]) - info['start_col'], int(nzr[-1]) - info['end_col']]
        info['spelling'] = dec
        info['pass1_mismatch'] = [L.mism(panel.P[sa[0]], tc), L.mism(panel.P[sb[0]], tc)]
        info['pass1_start_end_offset'] = []
        for r_ in (sa[0], sb[0]):
            nzr = np.nonzero(panel.P[r_] != L.OUT)[0]
            info['pass1_start_end_offset'].append([int(nzr[0]) - info['start_col'], int(nzr[-1]) - info['end_col']])
        info['_codes'] = tc
        info['_sha'] = C.sha256_str(core)
        info['_interior'] = core[ex[0][0]:ex[-1][1]]
        truth.append(info)
    if len(truth) == 1:
        truth = truth * 2

    evsets = {'R0': ev0}
    for E in (1, 2, 3):
        evsets[f'E{E}'] = L.make_evidence(panel.C, proj, E, bad_q, True)
    evsets['E2np'] = L.make_evidence(panel.C, proj, 2, frozenset(), False)
    res = dict(readset=readset, donor=donor, pass1=[names[sa[0]], names[sb[0]]],
               reads_assigned=len(reads.get(G, [])), unplaced=unplaced, mate_on_other_gene=len(bad_q),
               mate_other_genes=other, mapping={k: v for k, v in mstats.items() if 'DRB' in k or k in ('reads', 'unmapped')},
               family_sizes=[len(fam_a), len(fam_b)], same_family=same, panel_rows=panel.R,
               family_first_fields=[sorted({L.first_field(names[r]) for r in fam_a}),
                                    sorted({L.first_field(names[r]) for r in fam_b})],
               t_proj=round(t_proj, 1), t_ev0=round(t_ev0, 1),
               truth=[{k: v for k, v in t.items() if not k.startswith('_')} for t in truth], configs={})
    C.save_json(out / 'DIAG.json', res)
    plan = {'R0': ('F', 'Q'), 'E1': ('Q', 'G0', 'G3'), 'E2': ('Q', 'G0', 'G3'), 'E3': ('Q', 'G0', 'G3'),
            'E2np': ('G0',)}
    for evn, ev in evsets.items():
        for rule in plan[evn]:
            if rule == 'F':
                sup, contra, depth, info = L.tables_frozen(ev, frz)
            elif rule == 'Q':
                sup, contra, depth, info = L.tables_quant(ev)
            else:
                sup, contra, depth, info = L.tables_nogap(ev, pu.P, int(rule[1:]))
            pos = depth[depth > 0]
            N = PR.pair_evidence_cost(sup, contra, truth[0]['_codes'], truth[1]['_codes'])
            uc = PR.unexplained_columns(sup, contra, truth[0]['_codes'], truth[1]['_codes'])
            ngap = int(sum(1 for c in uc if sup[c, 4] or truth[0]['_codes'][c] == 4 or truth[1]['_codes'][c] == 4))
            Z1 = PR.pair_evidence_cost(sup, contra, panel.P[sa[0]], panel.P[sb[0]])
            base = dict(reads=ev.stats, median_depth=float(np.median(pos)) if len(pos) else 0, info=info,
                        N_truth_pair=N, N_cols=int(len(uc)), N_gap_cols=ngap, Z_pass1=Z1)
            for lam in (0.1, 0.5, 1.0):
                t0 = time.time()
                prm = PR.Params(**dict(frz.as_dict(), lam=lam, w_link=0.0, max_iter=3))
                best0, best = L.search(pa, pb, sup, contra, st_a, st_b, prm)
                cfg = dict(base, lam=lam, Z_best0=best0['cost'], best_cost=best['cost'],
                           gain=best0['cost'] - best['cost'],
                           best0_pair=[pa.names[best0['ha'].segs[0][0]], pb.names[best0['hb'].segs[0][0]]],
                           switches=[best['ha'].switches, best['hb'].switches], haps=[])
                hc = [(best['ha'], pa, ex_a), (best['hb'], pb, ex_b)]
                h0 = [best0['ha'].codes(pa.P), best0['hb'].codes(pb.P)]
                codes = [h.codes(pn.P) for h, pn, _ in hc]
                perm = min(itertools.permutations((0, 1)),
                           key=lambda p: sum(L.mism(codes[k], truth[p[k]]['_codes']) for k in (0, 1)))
                perm0 = min(itertools.permutations((0, 1)),
                            key=lambda p: sum(L.mism(h0[k], truth[p[k]]['_codes']) for k in (0, 1)))
                for k in (0, 1):
                    h, pn, exr = hc[k]
                    t = truth[perm[k]]
                    seq = PR.path_sequence(pn, h)
                    exo = PR.path_exons(pn, h, exr)
                    nz = np.nonzero(codes[k] != L.OUT)[0]
                    cfg['haps'].append(dict(
                        truth_hap=t['haplotype'], truth_label=t['ipd_label'], truth_ipd_absent=t['ipd_absent'],
                        switches=h.switches, sources=[pn.names[r_] for r_, _, _ in h.segs][:8],
                        mismatch_to_truth=L.mism(codes[k], t['_codes']),
                        best0_mismatch_to_truth=L.mism(h0[k], truth[perm0[k]]['_codes']) if perm0 == perm else None,
                        byte_exact=C.sha256_str(seq) == t['_sha'],
                        interior_exact=(exo is not None and seq[exo[0][0]:exo[-1][1]] == t['_interior']),
                        start_end_offset=[int(nz[0]) - t['start_col'], int(nz[-1]) - t['end_col']]))
                cfg['seconds'] = round(time.time() - t0, 1)
                res['configs'][f'{evn}|{rule}|{lam}'] = cfg
                C.save_json(out / 'DIAG.json', res)
    res['seconds_total'] = round(time.time() - t00, 1)
    res['complete'] = True
    C.save_json(out / 'DIAG.json', res)


if __name__ == '__main__':
    main()
