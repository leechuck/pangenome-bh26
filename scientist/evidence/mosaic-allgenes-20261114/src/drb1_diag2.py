"""Step 2: look at the evidence at the positions that distinguish each true DRB1 copy from its
nearest family allele (development donors only; truth used for evaluation only). Counts and
positions only, no sequence.   python drb1_diag2.py READSET DONOR"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L          # noqa: E402
C, EV, PR = L.C, L.EV, L.PR
import panel_select as PS     # noqa: E402
import run_donor as RD        # noqa: E402

S2 = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114')
G = 'DRB1'


def load_shared(readset, donor, out, arm='C'):
    """Gene-independent part of setup(): truth-free pass-1 calls of all genes and the reads mapped
    once to the pass-1 alleles of all genes (the mapping target never depended on the gene)."""
    C.verify_inputs()
    records = C.load_records()
    spec = RD.pass1_spec(readset, donor, 'C' if arm == 'C' else 'D')
    table = C.verified_table(spec['run'])
    geno = C.read_t1k_raw(table)
    aliases = C.load_json(spec['aliases']) if spec.get('aliases') else {}
    wanted, slots = set(), {}
    for g, row in geno.items():
        for sl in (row['slots'][:row['copies']] if row['copies'] else []):
            if sl['abundance'] <= 0:
                continue
            wanted |= set(sl['labels'][:3])
            slots.setdefault(g, []).append(sl['labels'][:3])
    seqs, ref_by_gene = RD.read_reference(spec['ref'], wanted, set(C.GENES))
    target = {n: seqs[n] for n in sorted(wanted) if n in seqs}
    reads, mstats, gene_of_q = L.map_reads2(target, spec['run'] / 'calls_candidate_1.fq',
                                            spec['run'] / 'calls_candidate_2.fq', out / 'map', 4, RD.ENV)
    shutil.rmtree(out / 'map', ignore_errors=True)
    return dict(records=records, spec=spec, aliases=aliases, slots=slots, reads=reads, gene_of_q=gene_of_q,
                ref_by_gene=ref_by_gene, mstats=mstats)


def setup(readset, donor, out, arm='C', gene=G, shared=None):
    """Shared loading: returns a dict with panel, families, reads, pass-1 rows (no truth).
    arm 'C': panel = IPD alleles, pass 1 = the genomic T1K run.
    arm 'D': panel = IPD + the cohort alleles of the donor's fold reference (arm B), pass 1 = arm B.
    gene: any of C.GENES (stage mosaic-allgenes: was the constant DRB1)."""
    G = gene
    if shared is None:
        shared = load_shared(readset, donor, out, arm)
    records, spec, aliases = shared['records'], shared['spec'], shared['aliases']
    msa = C.load_msas(records, genes=(G,))[G]
    slots = shared['slots'].get(G, [])
    reads, gene_of_q, ref_by_gene = shared['reads'], shared['gene_of_q'], shared['ref_by_gene']
    shas = [C.sha256_str(s) for s in msa.seqs]
    row_is_ipd = [bool(n) for n in msa.ipd_names]
    pg_shas = {v['sequence_sha256'] for v in aliases.values()} if arm == 'D' else None
    rows = PS.select_rows(arm, row_is_ipd, shas, pg_shas)
    excluded = RD.excluded_for(readset, donor, C.exclusion_map())
    comp_absent = {msa.row_of_source[r['source_id']] for r in records
                   if r['gene'] == G and r['donor'] in excluded and not r['genomic_labels']}
    PS.check_no_leakage(rows, row_is_ipd, msa.cohort_donors, excluded, donor, comp_absent, arm)
    label_of_sha = {}
    for lab, v in aliases.items():
        label_of_sha.setdefault(v['sequence_sha256'], lab)
    names, all_names = [], []
    for r in rows:
        if msa.ipd_names[r]:
            names.append(msa.ipd_names[r][0])
            all_names.append(list(msa.ipd_names[r]))
        else:
            n = label_of_sha.get(shas[r], 'COH_' + shas[r][:12])
            names.append(n)
            all_names.append([n])
    panel = PR.Panel(names, [msa.aligned[r] for r in rows])
    local = {r: i for i, r in enumerate(rows)}
    row_by_sha = {shas[r]: local[r] for r in rows}

    def label_row(lab):
        r = msa.row_of_allele.get(lab)
        if r is not None and r in local:
            return local[r]
        if lab in aliases and aliases[lab]['sequence_sha256'] in row_by_sha:
            return row_by_sha[aliases[lab]['sequence_sha256']]
        return None
    p1 = [[r for r in (label_row(l) for l in labs) if r is not None] for labs in slots]
    p1 = [x for x in p1 if x]
    assert p1, f'no pass-1 {G} call'
    if len(p1) == 1:
        p1 = p1 * 2
    sa, sb = p1[0], p1[1]
    rd, unplaced = L.windows(panel, reads.get(G, []), label_row)
    bad_q = {x['qname'] for x in reads.get(G, []) if gene_of_q[x['qname']] - {G}}
    fam_a = sorted(set(L.family_rows(panel, sa[0], 0.01)) | set(sa))
    fam_b = sorted(set(L.family_rows(panel, sb[0], 0.01)) | set(sb))
    same = sb[0] in fam_a or sa[0] in fam_b
    if same:
        fam_a = fam_b = sorted(set(fam_a) | set(fam_b))
    U = sorted(set(fam_a) | set(fam_b))
    pu = L.subpanel(panel, U)
    uidx = {r: i for i, r in enumerate(U)}
    rd_u = [dict(x, prefer_row=uidx.get(x['prefer_row'])) for x in rd]
    return dict(records=records, msa=msa, spec=spec, rows=rows, names=names, all_names=all_names, panel=panel,
                sa=sa, sb=sb, rd=rd, rd_u=rd_u, bad_q=bad_q, fam_a=fam_a, fam_b=fam_b, same=same, pu=pu, U=U,
                ref_by_gene=ref_by_gene, n_assigned=len(reads.get(G, [])), unplaced=unplaced, arm=arm, gene=G,
                label_row=label_row, rows_msa=rows)


def main():
    readset, donor = sys.argv[1], sys.argv[2]
    assert readset in ('D-sim', 'D-real'), 'development read sets only'
    assert os.environ.get('SLURM_JOB_ID'), 'run under Slurm'
    out = S2 / 'work/diag2' / f'{readset}_{donor}'
    out.mkdir(parents=True, exist_ok=True)
    X = setup(readset, donor, out)
    panel, pu, msa = X['panel'], X['pu'], X['msa']
    sa, sb = X['sa'], X['sb']
    classes = C.load_classes()
    proj = L.project_all(pu, X['rd_u'], 6)
    # reads projected through the pass-1 rows only (edit distance to the nearer pass-1 allele)
    p1rows = sorted(set(sa) | set(sb))
    pp = L.subpanel(panel, p1rows)
    pidx = {r: i for i, r in enumerate(p1rows)}
    rd_p = [dict(x, prefer_row=pidx.get(x['prefer_row'])) for x in X['rd']]
    proj_p = L.project_all(pp, rd_p, 6)
    evs = {'E1': L.make_evidence(panel.C, proj, 1, X['bad_q'], True),
           'E3': L.make_evidence(panel.C, proj, 3, X['bad_q'], True),
           'P2': L.make_evidence(panel.C, proj_p, 2, X['bad_q'], True),
           'P4': L.make_evidence(panel.C, proj_p, 4, X['bad_q'], True)}
    both_gap = np.ones(panel.C, dtype=bool)
    for r in p1rows:
        both_gap &= (panel.P[r] >= 4)
    trecs = sorted([r for r in X['records'] if r['donor'] == donor and r['gene'] == G],
                   key=lambda r: str(r['haplotype']))
    tp = PR.Panel([r['source_id'] for r in trecs],
                  [msa.aligned[msa.row_of_source[r['source_id']]] for r in trecs])
    tcs = [tp.P[i] for i in range(len(trecs))]
    if len(tcs) == 1:
        tcs = tcs * 2
        trecs = trecs * 2
    frz = PR.Params(**C.load_json(RD.FROZEN_PARAMS))
    res = dict(readset=readset, donor=donor, pass1=[X['names'][sa[0]], X['names'][sb[0]]],
               both_pass1_gap_columns=int(both_gap.sum()), evidence={}, truth=[])
    tabs = {}
    for k, ev in evs.items():
        sup, contra, depth, info = L.tables_quant(ev)
        sup2, contra2 = sup.copy(), contra.copy()
        sup2[:, 4] = False
        contra2[:, 4] = False
        sup2[both_gap] = False
        contra2[both_gap] = False
        tabs[k] = (ev, sup, contra, depth, info, sup2, contra2)
        N = PR.pair_evidence_cost(sup, contra, tcs[0], tcs[1])
        N2 = PR.pair_evidence_cost(sup2, contra2, tcs[0], tcs[1])
        Z = PR.pair_evidence_cost(sup, contra, panel.P[sa[0]], panel.P[sb[0]])
        Z2 = PR.pair_evidence_cost(sup2, contra2, panel.P[sa[0]], panel.P[sb[0]])
        uc = PR.unexplained_columns(sup, contra, tcs[0], tcs[1])
        runs = []
        if len(uc):
            st = pv = int(uc[0])
            for c in uc[1:].tolist():
                if c > pv + 3:
                    runs.append([st, pv])
                    st = c
                pv = c
            runs.append([st, pv])
        ex = []
        for c in uc[:12].tolist():
            ex.append(dict(col=c, counts=ev.counts[c].sum(axis=1).tolist(), truth=[int(tcs[0][c]), int(tcs[1][c])],
                           pass1=[int(panel.P[sa[0]][c]), int(panel.P[sb[0]][c])]))
        res['evidence'][k] = dict(reads=ev.stats, info=info, N=N, Z=Z, N_nogap=N2, Z_nogap=Z2,
                                  N_runs=[[a, b, b - a + 1] for a, b in runs][:15], N_examples=ex)
    for i, r in enumerate(trecs[:2]):
        tc = tcs[i]
        cl = classes[r['source_id']]
        other = tcs[1 - i]
        # nearest family row
        best = None
        for tag, fam in (('a', X['fam_a']), ('b', X['fam_b'])):
            for row in fam:
                m = L.mism(panel.P[row], tc)
                if best is None or m < best[0]:
                    best = (m, row, tag)
        m, row, tag = best
        rc = panel.P[row]
        v = (rc != L.OUT) & (tc != L.OUT) & (rc != tc)
        sites = []
        for c in np.nonzero(v)[0].tolist():
            s = dict(col=c, truth=int(tc[c]), nearest=int(rc[c]), partner=int(other[c]))
            for k, (ev, sup, contra, depth, info, sup2, contra2) in tabs.items():
                s[k] = dict(counts=ev.counts[c].sum(axis=1).tolist()[:6], D=info['D'],
                            truth_sup=bool(sup[c, tc[c]]), nearest_sup=bool(sup[c, rc[c]]),
                            nearest_contra=bool(contra[c, rc[c]]),
                            overdeep=bool(depth[c] > 1.6 * info['D']))
            sites.append(s)
        # can the search find it when the partner is right? (single Viterbi, partner fixed to truth)
        fam = X['fam_a'] if tag == 'a' else X['fam_b']
        pf = L.subpanel(panel, fam)
        orc = {}
        for k in ('E1', 'E3', 'P2', 'P4'):
            ev, sup, contra, depth, info, sup2, contra2 = tabs[k]
            for rule, (s_, c_) in (('Q', (sup, contra)), ('Qng', (sup2, contra2))):
                for lam in (0.1, 0.5):
                    prm = PR.Params(**dict(frz.as_dict(), lam=lam, w_link=0.0))
                    T = PR.pair_cost_table(s_, c_, other)
                    cost, path = PR.viterbi(pf, T, prm, {fam.index(row)}, True)
                    hc = path.codes(pf.P)
                    orc[f'{k}|{rule}|{lam}'] = dict(switches=path.switches, mismatch_to_truth=L.mism(hc, tc),
                                                   pair_cost=PR.pair_evidence_cost(s_, c_, hc, other),
                                                   zero_cost=PR.pair_evidence_cost(s_, c_, rc, other))
        res['truth'].append(dict(haplotype=str(r['haplotype']), ipd_label=cl['ipd_label_L150'],
                                 ipd_absent=not r['genomic_labels'], nearest_row=X['names'][row],
                                 nearest_mismatch=m, sites=sites[:40], n_sites=len(sites),
                                 partner_fixed_search=orc))
    C.save_json(out / 'DIAG2.json', res)


if __name__ == '__main__':
    main()
