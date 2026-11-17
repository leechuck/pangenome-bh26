"""Part A step 2 (and step 8): own-copy read coverage along each true DRB1 copy of a development donor.

A read counts for true copy h ("own_unique") when it aligns to h over its full length with <= 1 edit
(edlib, infix) and aligns worse to the donor's other true DRB1 copy; "own_strict" additionally requires
it to align worse to every pass-1 DRB1 allele (only informative where h differs from the pass-1 allele);
"any_le1" = <= 1 edit to h regardless of the other copy.
Coverage is counted per 100-base window of copy h (reads overlapping the window).
Read sources: the full recruited read set of the read set (r1/r2 of the arm-A T1K run) and
T1K's own candidate reads of that run (calls_candidate_*.fq).
Also: for reads whose best target is copy h with 2..10 edits, the positions on h where they differ
(to see whether the real genome differs from the truth sequence).
Development donors only. Output: counts and positions only.

  python cov.py READSET DONOR [RUN_DIR]     READSET in D-sim, D-real, D-real2
"""
from __future__ import annotations

import collections
import gzip
import json
import os
import sys
import time
from pathlib import Path

import edlib
import numpy as np

sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src')
import common as C          # noqa: E402
import run_donor as RD      # noqa: E402

ST = Path('/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114')
G = 'DRB1'
K = 25
WIN = 100
COMP = str.maketrans('ACGTN', 'TGCAN')
DEV = set(l.split()[1] for l in open('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113/work/dev_tasks.txt'))


def rc(s):
    return s[::-1].translate(COMP)


def fastq_seqs(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                yield line.strip().upper()


def run_spec(readset, donor, run_dir=None):
    if run_dir:
        run = Path(run_dir)
        rec = C.load_json(run / 'COMPLETE.json')
        cmd = RD.parse_cmd(rec['command'])
        return dict(run=run, ref=RD.to_host(cmd['-f']) if cmd['-f'].startswith('/benchmark') else Path(cmd['-f']),
                    r1=cmd['-1'], r2=cmd['-2'])
    return RD.pass1_spec(readset, donor, 'A')


def host(p):
    p = str(p)
    for k in RD.BIND:
        if p == k or p.startswith(k + '/'):
            return RD.to_host(p)
    return Path(p)


def main():
    readset, donor = sys.argv[1], sys.argv[2]
    run_dir = sys.argv[3] if len(sys.argv) > 3 else None
    assert donor in DEV, 'development donors only'
    assert os.environ.get('SLURM_JOB_ID'), 'run under Slurm'
    out = ST / 'work/cov' / f'{readset}_{donor}'
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    records = C.load_records()
    msa = C.load_msas(records, genes=(G,))[G]
    trecs = sorted([r for r in records if r['donor'] == donor and r['gene'] == G], key=lambda r: str(r['haplotype']))
    classes = C.load_classes()
    copies = []
    for r in trecs:
        seq = r['core_sequence'].upper()
        row = msa.aligned[msa.row_of_source[r['source_id']]]
        col_of = [c for c, ch in enumerate(row) if ch != '-']
        assert len(col_of) == len(seq)
        cl = classes.get(r['source_id'], {})
        copies.append(dict(hap=str(r['haplotype']), seq=seq, col_of=col_of, absent=not r['genomic_labels'],
                           label=cl.get('ipd_label_L150'), labels=list(r['genomic_labels'])[:3]))
    spec = run_spec(readset, donor, run_dir)
    geno = C.read_t1k_raw(C.verified_table(spec['run']))
    p1 = []
    row = geno.get(G)
    if row and row['copies']:
        for sl in row['slots'][:row['copies']]:
            if sl['abundance'] > 0:
                p1 += sl['labels'][:3]
    p1 = sorted(set(p1))
    seqs, _ = RD.read_reference(spec['ref'], set(p1), {G})
    others = [seqs[n] for n in p1 if n in seqs]
    # k-mer index over both copies (forward); reads tested in both orientations
    kidx = set()
    for c in copies:
        s = c['seq']
        for i in range(len(s) - K + 1):
            kidx.add(s[i:i + K])
    sources = {'recruited': [host(spec['r1']), host(spec['r2'])],
               't1k_candidates': [spec['run'] / 'calls_candidate_1.fq', spec['run'] / 'calls_candidate_2.fq']}
    res = dict(readset=readset, donor=donor, run=str(spec['run']), pass1=p1,
               copies=[dict(hap=c['hap'], length=len(c['seq']), ipd_absent=c['absent'], class_label=c['label'],
                            genomic_labels=c['labels']) for c in copies], sources={})
    for sname, files in sources.items():
        n_reads = n_cand = 0
        cov_u = [np.zeros(len(c['seq']), dtype=np.int32) for c in copies]
        cov_any = [np.zeros(len(c['seq']), dtype=np.int32) for c in copies]
        cov_s = [np.zeros(len(c['seq']), dtype=np.int32) for c in copies]
        diffs = [collections.Counter() for _ in copies]
        near_cov = [np.zeros(len(c['seq']), dtype=np.int32) for c in copies]
        n_own = [0, 0]
        for fpath in files:
            for s in fastq_seqs(fpath):
                n_reads += 1
                L_ = len(s)
                if L_ < 60:
                    continue
                hit = False
                for q in (s, rc(s)):
                    for i in range(0, L_ - K + 1, 8):
                        if q[i:i + K] in kidx:
                            hit = True
                            break
                    if hit:
                        break
                if not hit:
                    continue
                n_cand += 1
                best = []
                for ci, c in enumerate(copies):
                    b = None
                    for q in (s, rc(s)):
                        a = edlib.align(q, c['seq'], mode='HW', task='locations', k=12)
                        if a['editDistance'] >= 0 and (b is None or a['editDistance'] < b[0]):
                            b = (a['editDistance'], a['locations'][0], q)
                    best.append(b)
                oth = 99
                for t in others:
                    for q in (s, rc(s)):
                        a = edlib.align(q, t, mode='HW', task='distance', k=12)
                        if a['editDistance'] >= 0:
                            oth = min(oth, a['editDistance'])
                for ci, b in enumerate(best):
                    if b is None:
                        continue
                    ed, (st, en), q = b
                    ed_other_copy = min([x[0] for j, x in enumerate(best) if j != ci and x is not None] or [99])
                    if len(copies) == 2 and copies[0]['seq'] == copies[1]['seq']:
                        ed_other_copy = 99
                    if ed <= 1:
                        cov_any[ci][st:en + 1] += 1
                        if ed < ed_other_copy:
                            cov_u[ci][st:en + 1] += 1
                            n_own[ci] += 1
                            if ed < oth:
                                cov_s[ci][st:en + 1] += 1
                    elif ed <= 10 and ed <= ed_other_copy and ed <= oth:
                        near_cov[ci][st:en + 1] += 1
                        a = edlib.align(q, copies[ci]['seq'][st:en + 1], mode='NW', task='path')
                        pos = st
                        for n, op in _cigar(a['cigar']):
                            if op in '=M':
                                pos += n
                            elif op == 'X':
                                for j in range(n):
                                    diffs[ci][(pos + j, 'X')] += 1
                                pos += n
                            elif op == 'D':   # base in target missing from read
                                diffs[ci][(pos, f'D{n}')] += 1
                                pos += n
                            elif op == 'I':   # extra bases in read
                                diffs[ci][(pos, f'I{n}')] += 1
        per_copy = []
        for ci, c in enumerate(copies):
            n = len(c['seq'])
            wins = []
            for w0 in range(0, n, WIN):
                w1 = min(n, w0 + WIN)
                wins.append(dict(start=w0, end=w1, col0=c['col_of'][w0], col1=c['col_of'][w1 - 1],
                                 own_unique=int(cov_u[ci][w0:w1].max()), own_unique_min=int(cov_u[ci][w0:w1].min()),
                                 own_strict=int(cov_s[ci][w0:w1].max()), any_le1=int(cov_any[ci][w0:w1].max()), near_2_10=int(near_cov[ci][w0:w1].max())))
            top = sorted(((v, p, k) for (p, k), v in diffs[ci].items() if v >= 3), reverse=True)[:60]
            per_copy.append(dict(hap=c['hap'], own_unique_reads=n_own[ci],
                                 median_window_own_unique=float(np.median([w['own_unique'] for w in wins])),
                                 windows=wins,
                                 recurrent_diffs=[dict(pos=p, col=c['col_of'][min(p, n - 1)], kind=k, reads=v,
                                                       own_unique_depth_here=int(cov_u[ci][min(p, n - 1)]))
                                                  for v, p, k in top]))
        res['sources'][sname] = dict(files=[str(x) for x in files], reads=n_reads, kmer_candidates=n_cand,
                                     copies=per_copy)
        C.save_json(out / 'COV.json', res)
    res.update(complete=True, seconds=round(time.time() - t0))
    C.save_json(out / 'COV.json', res)
    C.save_json(out / 'COMPLETE.json', dict(status='complete', job=os.environ['SLURM_JOB_ID']))


def _cigar(c):
    n = ''
    for ch in c:
        if ch.isdigit():
            n += ch
        else:
            yield int(n), ch
            n = ''


if __name__ == '__main__':
    try:
        main()
    except BaseException as e:
        import traceback
        out = ST / 'work/cov' / f'{sys.argv[1]}_{sys.argv[2]}'
        out.mkdir(parents=True, exist_ok=True)
        C.save_json(out / 'FAILURE.json', dict(status='failed', error=repr(e), traceback=traceback.format_exc()))
        raise
