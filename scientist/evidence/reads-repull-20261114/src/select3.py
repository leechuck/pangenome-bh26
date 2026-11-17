"""Part B, rule 2 (truth-free), read set D-real3 for one development donor.

D-real3 = old regional pull (D-real: chr6:28.5-33.5 Mb + chr6_* alt + HLA-* contigs, with mates)
          UNION every read pair of the whole-genome data in which either mate carries >= 3 bait 31-mers
          among the 31-mers starting at every 4th position (either strand, canonical k-mers).
Bait 31-mers (work/bait/bait_kmers.txt) = 31-mers of the IPD genomic alleles of all genes in genomic.fa
that do not occur anywhere in GRCh38 outside the MHC (chr6:28.5-33.5 Mb masked, alt/HLA contigs left out;
decoys, unplaced contigs and EBV included). Whole-genome read pairs are scanned through T1K's own
full-WGS candidate set (baselines/results/wgs.<donor>/t1k/t1k_genomic_ipd/calls_candidate_*.fq,
fastq-extractor -s 0.97 against genomic.fa), which contains every read with substantial IPD sequence.
No truth is read.

  python select3.py DONOR    -> work/reads3/DONOR/{r1.fq.gz,r2.fq.gz,PAIRS.json}
"""
import gzip
import os
import sys
from pathlib import Path

sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src')
import common as C          # noqa: E402
import run_donor as RD      # noqa: E402

ST = Path('/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114')
K, STRIDE, MIN_HITS = 31, 4, 3
COMP = str.maketrans('ACGTN', 'TGCAN')


def fq(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as f:
        while True:
            h = f.readline()
            if not h:
                return
            s, p, q = f.readline(), f.readline(), f.readline()
            yield h, s, q


def qname(h):
    n = h[1:].split()[0]
    return n[:-2] if n.endswith(('/1', '/2')) else n


def hits(seq, bait):
    s = seq.strip().upper()
    r = s[::-1].translate(COMP)
    L = len(s)
    n = 0
    for i in range(0, L - K + 1, STRIDE):
        a = s[i:i + K]
        b = r[L - i - K:L - i]
        if (a if a < b else b) in bait:
            n += 1
            if n >= MIN_HITS:
                return True
    return False


def main():
    donor = sys.argv[1]
    assert os.environ.get('SLURM_JOB_ID')
    bait = set(l.strip() for l in open(ST / 'work/bait/bait_kmers.txt'))
    old = RD.pass1_spec('D-real', donor, 'A')
    o1, o2 = RD.to_host(old['r1']), RD.to_host(old['r2'])
    cand = C.ZEN / 'baselines/results' / f'wgs.{donor}' / 't1k/t1k_genomic_ipd'
    out = ST / 'work/reads3' / donor
    out.mkdir(parents=True, exist_ok=True)
    names = set()
    n_old = n_add = n_cand = 0
    with gzip.open(out / 'r1.fq.gz', 'wt', compresslevel=4) as w1, gzip.open(out / 'r2.fq.gz', 'wt', compresslevel=4) as w2:
        for (h1, s1, q1), (h2, s2, q2) in zip(fq(o1), fq(o2)):
            n1 = qname(h1)
            assert n1 == qname(h2)
            names.add(n1)
            w1.write(f'{h1}{s1}+\n{q1}')
            w2.write(f'{h2}{s2}+\n{q2}')
            n_old += 1
        for (h1, s1, q1), (h2, s2, q2) in zip(fq(cand / 'calls_candidate_1.fq'), fq(cand / 'calls_candidate_2.fq')):
            n_cand += 1
            n1 = qname(h1)
            assert n1 == qname(h2)
            if n1 in names:
                continue
            if hits(s1, bait) or hits(s2, bait):
                names.add(n1)
                w1.write(f'{h1}{s1}+\n{q1}')
                w2.write(f'{h2}{s2}+\n{q2}')
                n_add += 1
    C.save_json(out / 'PAIRS.json', dict(
        status='paired_fastq_verified', donor=donor, readset='D-real3', pairs=n_old + n_add,
        old_regional_pairs=n_old, added_pairs=n_add, wgs_candidate_pairs_scanned=n_cand,
        rule=__doc__.split('\n\n')[1], bait=C.load_json(ST / 'work/bait/BAIT.json'),
        r1=str(out / 'r1.fq.gz'), r2=str(out / 'r2.fq.gz'),
        r1_sha256=C.sha256_file(out / 'r1.fq.gz'), r2_sha256=C.sha256_file(out / 'r2.fq.gz'),
        truth_loaded=False, job=os.environ['SLURM_JOB_ID']))


if __name__ == '__main__':
    main()
