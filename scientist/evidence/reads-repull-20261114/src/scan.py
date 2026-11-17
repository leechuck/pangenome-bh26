"""Part A step 3: look for the donor's own DRB1 intron-1 reads anywhere in the original WGS alignment.

  python scan.py kmers DONOR    -> work/scan/DONOR/patterns.txt, KMERS.json
  python scan.py analyze DONOR  -> work/scan/DONOR/SCAN.json  (reads work/scan/DONOR/hits.sam)

k-mers (31): from each true DRB1 copy of the donor, the stretch aligned to MSA columns COL0..COL1
(the intron-1 hole seen on real reads; widened by the hole runs measured in step 2), keeping k-mers
absent from the donor's other true copy. Both orientations are searched in every read of the CRAM
(mapped anywhere, secondary/supplementary records included, unmapped included) with grep -F.
A hit read is "own" for copy h when it aligns to copy h over its full length with <= 2 edits and
better than to the other true copy (whether it also fits better than the pass-1 alleles of the D-real
arm-A run is reported separately).
Development donors only (truth used for diagnosis). Output: counts, contig names, positions, flags.
"""
from __future__ import annotations

import collections
import gzip
import json
import os
import sys
from pathlib import Path

import edlib

sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src')
import common as C          # noqa: E402
import run_donor as RD      # noqa: E402

ST = Path('/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114')
SUB = os.environ.get('SCAN_SUB', 'scan')     # 'scan' = hole k-mers; 'scan_diag' = diagnostic k-mers (whole copy)
DIAG = SUB == 'scan_diag'
G = 'DRB1'
K = 31
COL0, COL1 = 3300, 5600
COMP = str.maketrans('ACGTN', 'TGCAN')
DEV = set(l.split()[1] for l in open('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113/work/dev_tasks.txt'))
MHC = ('chr6', 28500000, 33500000)


def rc(s):
    return s[::-1].translate(COMP)


def truth(donor):
    records = C.load_records()
    msa = C.load_msas(records, genes=(G,))[G]
    out = []
    for r in sorted([r for r in records if r['donor'] == donor and r['gene'] == G], key=lambda r: str(r['haplotype'])):
        seq = r['core_sequence'].upper()
        row = msa.aligned[msa.row_of_source[r['source_id']]]
        col_of = [c for c, ch in enumerate(row) if ch != '-']
        out.append(dict(hap=str(r['haplotype']), seq=seq, col_of=col_of, absent=not r['genomic_labels']))
    return out


def pass1_kmers(donor):
    spec = RD.pass1_spec('D-real', donor, 'A')
    geno = C.read_t1k_raw(C.verified_table(spec['run']))
    p1 = sorted({l for sl in geno[G]['slots'][:geno[G]['copies']] if sl['abundance'] > 0 for l in sl['labels'][:3]})
    seqs, _ = RD.read_reference(spec['ref'], set(p1), {G})
    out = set()
    for n in p1:
        q = seqs[n]
        out |= {q[i:i + K] for i in range(len(q) - K + 1)}
    return out


def hole_cols(donor):
    """MSA column range: fixed 3300-5600, widened by runs of 100-base windows where D-real recruited
    own-unique coverage is <= 1 while D-sim has >= 5 (step 2 output), restricted to 2000-8000."""
    lo, hi = COL0, COL1
    try:
        real = C.load_json(ST / f'work/cov/D-real_{donor}/COV.json')['sources']['recruited']['copies']
        sim = C.load_json(ST / f'work/cov/D-sim_{donor}/COV.json')['sources']['recruited']['copies']
        for cr, cs in zip(real, sim):
            for wr, ws in zip(cr['windows'], cs['windows']):
                if wr['own_unique'] <= 1 and ws['own_unique'] >= 5 and 2000 <= wr['col0'] <= 8000:
                    lo, hi = min(lo, wr['col0']), max(hi, wr['col1'])
    except FileNotFoundError:
        pass
    return lo, hi


def kmers(donor):
    out = ST / 'work' / SUB / donor
    out.mkdir(parents=True, exist_ok=True)
    T = truth(donor)
    lo, hi = (0, 10 ** 9) if DIAG else hole_cols(donor)
    allk = [set(c['seq'][i:i + K] for i in range(len(c['seq']) - K + 1)) for c in T]
    pats, info = set(), []
    for h, c in enumerate(T):
        idx = [i for i, col in enumerate(c['col_of']) if lo <= col <= hi]
        if not idx:
            info.append(dict(hap=c['hap'], kmers=0))
            continue
        a, b = idx[0], idx[-1]
        other = set().union(*[allk[j] for j in range(len(T)) if j != h and T[j]['seq'] != c['seq']]) if len(T) > 1 else set()
        ks = {c['seq'][i:i + K] for i in range(a, b - K + 2)} - other
        ks = {k for k in ks if 'N' not in k}
        if DIAG:   # only k-mers that separate the true copy from the pass-1 alleles of the D-real arm-A run
            ks -= pass1_kmers(donor)
        pats |= ks | {rc(k) for k in ks}
        info.append(dict(hap=c['hap'], copy_start=a, copy_end=b, kmers=len(ks)))
    (out / 'patterns.txt').write_text(''.join(p + '\n' for p in sorted(pats)))
    C.save_json(out / 'KMERS.json', dict(donor=donor, col_range=[lo, hi], copies=info, patterns=len(pats)))


def filter_(donor):
    out = ST / 'work' / SUB / donor
    raw = [l.strip() for l in open(out / 'patterns.raw.txt') if l.strip()]
    rep = {l.strip() for l in open(out / 'repeat_kmers.txt') if l.strip()}
    rep |= {rc(k) for k in rep}
    keep = [k for k in raw if k not in rep]
    (out / 'patterns.txt').write_text(''.join(k + '\n' for k in keep))
    info = C.load_json(out / 'KMERS.json')
    info.update(patterns_before_repeat_filter=len(raw), repeat_patterns_removed=len(raw) - len(keep), patterns=len(keep),
                repeat_filter='k-mers (either orientation) seen in reads aligned to chr1:1-60,000,000 of the same CRAM')
    C.save_json(out / 'KMERS.json', info)


def in_recruited_region(rname, pos):
    if rname == '*':
        return False
    if rname.startswith('HLA-') or rname.startswith('chr6_'):
        return True
    return rname == MHC[0] and MHC[1] <= pos - 1 < MHC[2]


def analyze(donor):
    out = ST / 'work' / SUB / donor
    T = truth(donor)
    lo, hi = C.load_json(out / 'KMERS.json')['col_range']
    pats = set(l.strip() for l in open(out / 'patterns.txt') if l.strip())
    spec = RD.pass1_spec('D-real', donor, 'A')
    geno = C.read_t1k_raw(C.verified_table(spec['run']))
    p1 = sorted({l for sl in geno[G]['slots'][:geno[G]['copies']] if sl['abundance'] > 0 for l in sl['labels'][:3]})
    seqs, _ = RD.read_reference(spec['ref'], set(p1), {G})
    others = [seqs[n] for n in p1 if n in seqs]
    names = set()
    for f in (spec['r1'], spec['r2']):
        with gzip.open(RD.to_host(f), 'rt') as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    names.add(line[1:].split()[0].removesuffix('/1').removesuffix('/2'))
    recs = []
    with open(out / 'hits.sam') as f:
        for line in f:
            v = line.rstrip('\n').split('\t')
            if len(v) < 11:
                continue
            seq = v[9].upper()
            nk = sum(1 for i in range(len(seq) - K + 1) if seq[i:i + K] in pats)
            if nk == 0:
                continue
            recs.append(dict(q=v[0], flag=int(v[1]), rname=v[2], pos=int(v[3]), mapq=int(v[4]), rnext=v[6],
                             pnext=int(v[7]), seq=seq, nk=nk))
    res = dict(donor=donor, col_range=[lo, hi], pass1=p1, hit_records=len(recs), copies=[])
    per_copy = [collections.defaultdict(list) for _ in T]
    k25 = set()
    for c in T:
        k25 |= {c['seq'][i:i + 25] for i in range(len(c['seq']) - 24)}
    res['primary_hit_records'] = sum(1 for r in recs if not r['flag'] & 0x900)
    res['prefilter'] = 'edlib only for reads with >= 50% of their 25-mers (either orientation) in the true copies'
    n_pass = 0
    for r in recs:
        if r['flag'] & 0x900:
            continue
        sq = r['seq']
        n25 = max(1, len(sq) - 24)
        f1 = sum(1 for i in range(len(sq) - 24) if sq[i:i + 25] in k25) / n25
        sr = rc(sq)
        f2 = sum(1 for i in range(len(sr) - 24) if sr[i:i + 25] in k25) / n25
        if max(f1, f2) < 0.5:
            continue
        n_pass += 1
        best = []
        for c in T:
            b = 99, None
            for q in (r['seq'], rc(r['seq'])):
                a = edlib.align(q, c['seq'], mode='HW', task='locations', k=10)
                if 0 <= a['editDistance'] < b[0]:
                    b = (a['editDistance'], a['locations'][0])
            best.append(b)
        oth = 99
        for t in others:
            for q in (r['seq'], rc(r['seq'])):
                a = edlib.align(q, t, mode='HW', task='distance', k=10)
                if a['editDistance'] >= 0:
                    oth = min(oth, a['editDistance'])
        for h, (ed, loc) in enumerate(best):
            od = min([best[j][0] for j in range(len(T)) if j != h and T[j]['seq'] != T[h]['seq']] or [99])
            if ed <= 2 and ed < od and (not DIAG or ed < oth):
                mate_rname = r['rname'] if r['rnext'] == '=' else r['rnext']
                unm, munm = bool(r['flag'] & 4), bool(r['flag'] & 8)
                pulled = (not unm and in_recruited_region(r['rname'], r['pos'])) or \
                         (not munm and in_recruited_region(mate_rname, r['pnext']))
                col = T[h]['col_of'][loc[0]]
                per_copy[h][r['q']].append(dict(flag=r['flag'], rname=r['rname'] if not unm else '*unmapped*',
                                                pos=r['pos'], mapq=r['mapq'], mate=mate_rname if not munm else '*unmapped*',
                                                mate_pos=r['pnext'], proper=bool(r['flag'] & 2), ed=ed, col=col,
                                                in_hole=lo <= col <= hi, pulled_by_rule=pulled, better_than_pass1=ed < oth,
                                                in_recruited_fastq=r['q'] in names))
    res['prefilter_passed'] = n_pass
    for h, c in enumerate(T):
        reads = [x for v in per_copy[h].values() for x in v]
        hole = [x for x in reads if x['in_hole']]
        def tab(rs, key):
            return dict(collections.Counter(str(x[key]) for x in rs).most_common(30))
        contig_class = lambda x: ('unmapped' if x['rname'] == '*unmapped*' else 'MHC chr6 28.5-33.5M' if x['rname'] == 'chr6' and in_recruited_region('chr6', x['pos'])
                                  else 'chr6 outside MHC' if x['rname'] == 'chr6' else 'HLA contig' if x['rname'].startswith('HLA-')
                                  else 'chr6 alt' if x['rname'].startswith('chr6_') else 'decoy' if x['rname'].endswith('_decoy') or x['rname'].startswith('chrUn_JTFH')
                                  else 'other:' + x['rname'].split('_')[0])
        for x in hole:
            x['contig_class'] = contig_class(x)
        mq = collections.Counter(('0' if x['mapq'] == 0 else '1-19' if x['mapq'] < 20 else '20-59' if x['mapq'] < 60 else '60') for x in hole)
        res['copies'].append(dict(
            hap=c['hap'], ipd_absent=c['absent'], own_reads_all=len(reads), own_reads_in_hole=len(hole),
            own_pairs_in_hole=len({q for q, v in per_copy[h].items() if any(x['in_hole'] for x in v)}),
            in_hole_in_recruited_fastq=sum(x['in_recruited_fastq'] for x in hole),
            in_hole_pulled_by_rule=sum(x['pulled_by_rule'] for x in hole),
            in_hole_contig_class=tab(hole, 'contig_class'), in_hole_contig=tab(hole, 'rname'),
            in_hole_mate_contig=tab(hole, 'mate'), in_hole_mapq=dict(mq),
            in_hole_proper_pair=sum(x['proper'] for x in hole),
            in_hole_better_than_pass1=sum(x['better_than_pass1'] for x in hole),
            in_hole_not_recruited=sum(not x['in_recruited_fastq'] for x in hole),
            in_hole_not_recruited_contig_class=tab([x for x in hole if not x['in_recruited_fastq']], 'contig_class'),
            in_hole_not_recruited_positions=sorted([(x['rname'], x['pos'], x['mapq'], x['mate'], x['mate_pos'], x['flag'], x['col'])
                                                    for x in hole if not x['in_recruited_fastq']])[:300],
            in_hole_recruited_positions_sample=sorted([(x['rname'], x['pos'], x['mapq'], x['col'])
                                                       for x in hole if x['in_recruited_fastq']])[:100]))
    if DIAG:
        for h, c in enumerate(T):
            reads = [x for v in per_copy[h].values() for x in v]
            bins = collections.defaultdict(lambda: [0, 0])
            for x in reads:
                bins[x['col'] // 500 * 500][0 if x['in_recruited_fastq'] else 1] += 1
            res['copies'][h]['diag_reads_by_col500_recruited_vs_not'] = {str(k): v for k, v in sorted(bins.items())}
            res['copies'][h]['diag_reads_recruited'] = sum(x['in_recruited_fastq'] for x in reads)
            res['copies'][h]['diag_reads_not_recruited'] = sum(not x['in_recruited_fastq'] for x in reads)
            res['copies'][h]['diag_not_recruited_contig_class'] = dict(collections.Counter(
                ('unmapped' if x['rname'] == '*unmapped*' else 'chr6 MHC' if in_recruited_region(x['rname'], x['pos']) and x['rname'] == 'chr6'
                 else x['rname'] if x['rname'].startswith(('HLA-', 'chr6_')) else x['rname'].split('_')[0] + ('_decoy' if 'decoy' in x['rname'] else ''))
                for x in reads if not x['in_recruited_fastq']).most_common(30))
    C.save_json(out / 'SCAN.json', res)


if __name__ == '__main__':
    mode, donor = sys.argv[1], sys.argv[2]
    assert donor in DEV
    assert os.environ.get('SLURM_JOB_ID')
    {'kmers': kmers, 'filter': filter_, 'analyze': analyze}[mode](donor)
