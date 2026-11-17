"""Step 4: where do the reads for each local version come from? For blocks of interest (given as
MSA column ranges), count reads containing each known local version among ALL T1K candidate reads,
split by the gene bwa placed them on and by filter stage. Development donors only. Counts only.
  python drb1_diag4.py READSET DONOR"""
from __future__ import annotations

import collections
import os
import shutil
import subprocess
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
COMP = str.maketrans('ACGTN', 'TGCAN')
TARGETS = [(4497, 4527), (4621, 4661), (5757, 5804), (4914, 4944), (11858, 12023), (17176, 17241)]


def main():
    readset, donor = sys.argv[1], sys.argv[2]
    assert readset in ('D-sim', 'D-real')
    assert os.environ.get('SLURM_JOB_ID')
    out = S2 / 'work/diag4' / f'{readset}_{donor}'
    out.mkdir(parents=True, exist_ok=True)
    X = D2.setup(readset, donor, out)
    panel = X['panel']
    sa, sb = X['sa'], X['sb']
    proj = L.project_all(X['pu'], X['rd_u'], 6)
    edit_of = {}
    for x, p in zip(X['rd_u'], proj):
        edit_of[(x['qname'], x['reverse'], x['seq'][:20])] = p[2]
    # all candidate reads with their placement gene (second mapping pass, cheap)
    spec = X['spec']
    geno = C.read_t1k_raw(C.verified_table(spec['run']))
    wanted = set()
    for g, row in geno.items():
        for sl in (row['slots'][:row['copies']] if row['copies'] else []):
            if sl['abundance'] > 0:
                wanted |= set(sl['labels'][:3])
    seqs, _ = RD.read_reference(spec['ref'], wanted, {G})
    work = out / 'map'
    work.mkdir(exist_ok=True)
    fa = work / 't.fa'
    fa.write_text(''.join(f'>{n}\n{seqs[n]}\n' for n in sorted(wanted) if n in seqs))
    subprocess.run([str(RD.ENV / 'bwa'), 'index', str(fa)], check=True, capture_output=True)
    sam = work / 'a.sam'
    with open(sam, 'w') as f:
        subprocess.run([str(RD.ENV / 'bwa'), 'mem', '-t', '4', '-v', '1', str(fa),
                        str(spec['run'] / 'calls_candidate_1.fq'), str(spec['run'] / 'calls_candidate_2.fq')],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)
    import pysam
    allreads = []
    mate_gene = collections.defaultdict(set)
    with pysam.AlignmentFile(str(sam)) as af:
        for rec in af.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary:
                continue
            g = 'unmapped' if rec.is_unmapped else C.gene_of_name(rec.reference_name)
            mate_gene[rec.query_name].add(g)
            allreads.append((rec.query_name, rec.query_sequence, g, rec.mapping_quality,
                             rec.get_tag('NM') if rec.has_tag('NM') else -1))
    shutil.rmtree(work, ignore_errors=True)
    res = dict(readset=readset, donor=donor, pass1=[X['names'][sa[0]], X['names'][sb[0]]], blocks=[])
    for h, (fam, p1) in enumerate(((X['fam_a'], sa[0]), (X['fam_b'], sb[0]))):
        rate = L.subst_rate(panel.P, p1)
        near = [r for r in fam if rate[r] <= 0.004 or r == p1]
        fp = L.subpanel(panel, near)
        blocks = PT.build_blocks(fp)
        for b in blocks:
            if not any(b['cl'] <= t1 and b['cr'] >= t0 for t0, t1 in TARGETS):
                continue
            cur = PT.row_string(fp.P[near.index(p1)], b['cl'], b['cr'])
            vers = sorted(b['alleles'], key=lambda s: (s != cur, s))
            rows = []
            for vi, s in enumerate(vers):
                rc = s[::-1].translate(COMP)
                cnt = collections.Counter()
                for q, seq, g, mq, nm in allreads:
                    if s in seq or rc in seq:
                        cnt['all'] += 1
                        cnt['gene:' + g] += 1
                        if g == G:
                            cnt['drb1_mate_other_gene'] += bool(mate_gene[q] - {G, 'unmapped'})
                            cnt['drb1_nm0'] += nm == 0
                            cnt['drb1_nm_le1'] += nm <= 1
                            cnt['drb1_nm_ge3'] += nm >= 3
                rows.append(dict(version=vi, is_pass1_version=s == cur, length=len(s),
                                 n_rows=len(b['alleles'][s]),
                                 example_rows=[fp.names[r] for r in b['alleles'][s][:3]], counts=dict(cnt)))
            res['blocks'].append(dict(hap=h, pass1=X['names'][p1], cl=b['cl'], cr=b['cr'], resolvable=b['resolvable'],
                                      versions=rows))
    C.save_json(out / 'DIAG4.json', res)


if __name__ == '__main__':
    main()
