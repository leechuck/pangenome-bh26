"""Step 2b: rule-2 read set for one donor (rule exactly as select3.py, step 2; same bait file, same
hits() test, same union with the old regional pull; imported from select3.py).
Whole-genome candidate pairs: T1K full-WGS candidate files of the baseline run
(baselines/results/wgs.<donor>/t1k/t1k_genomic_ipd) when they exist; otherwise this stage's own
extraction with the same fastq-extractor call over all primary pairs of the CRAM (work/wgsx/<donor>, wgsx.sbatch).

  python select_rule2.py RS DONOR     RS: D-real3 | T-real3   -> READS[RS]/DONOR/{r1,r2}.fq.gz, PAIRS.json
"""
import gzip
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common3 as K          # noqa: E402
import select3 as S          # noqa: E402
C, RD = K.C, K.RD


def cand_dir(donor):
    base = C.ZEN / 'baselines/results' / f'wgs.{donor}' / 't1k'
    if (base / 't1k_genomic_ipd/calls_candidate_1.fq').exists():
        rec = C.load_json(base / 'COMPLETE.json')
        assert rec['status'] == 'complete' and rec['donor'] == donor
        log = (base / 't1k_genomic_ipd/run.log').read_text()
        assert 'fastq-extractor' in log and '-s 0.97' in log and log.rstrip().endswith('Finish.')
        return base / 't1k_genomic_ipd', 'baseline full-WGS T1K run'
    own = K.ST / 'work/wgsx' / donor
    rec = C.load_json(own / 'COMPLETE.json')
    assert rec['status'] == 'complete' and rec['donor'] == donor
    return own, 'this stage, wgsx.sbatch (same fastq-extractor call)'


def main():
    rs, donor = sys.argv[1], sys.argv[2]
    assert os.environ.get('SLURM_JOB_ID') and rs in K.OLD
    K.require_freeze(rs)
    bait = set(l.strip() for l in open(K.ST / 'work/bait/bait_kmers.txt'))
    old = K._orig_pass1(K.OLD[rs], donor, 'A')
    o1, o2 = RD.to_host(old['r1']), RD.to_host(old['r2'])
    cand, cand_src = cand_dir(donor)
    out = K.READS[rs] / donor
    out.mkdir(parents=True, exist_ok=True)
    names = set()
    n_old = n_add = n_cand = 0
    with gzip.open(out / 'r1.fq.gz', 'wt', compresslevel=4) as w1, gzip.open(out / 'r2.fq.gz', 'wt', compresslevel=4) as w2:
        for (h1, s1, q1), (h2, s2, q2) in zip(S.fq(o1), S.fq(o2)):
            n1 = S.qname(h1)
            assert n1 == S.qname(h2)
            names.add(n1)
            w1.write(f'{h1}{s1}+\n{q1}')
            w2.write(f'{h2}{s2}+\n{q2}')
            n_old += 1
        for (h1, s1, q1), (h2, s2, q2) in zip(S.fq(cand / 'calls_candidate_1.fq'), S.fq(cand / 'calls_candidate_2.fq')):
            n_cand += 1
            n1 = S.qname(h1)
            assert n1 == S.qname(h2)
            if n1 in names:
                continue
            if S.hits(s1, bait) or S.hits(s2, bait):
                names.add(n1)
                w1.write(f'{h1}{s1}+\n{q1}')
                w2.write(f'{h2}{s2}+\n{q2}')
                n_add += 1
    C.save_json(out / 'PAIRS.json', dict(
        status='paired_fastq_verified', donor=donor, readset=rs, pairs=n_old + n_add, old_readset=K.OLD[rs],
        old_regional_pairs=n_old, added_pairs=n_add, wgs_candidate_pairs_scanned=n_cand, candidates=str(cand),
        candidates_source=cand_src, rule=S.__doc__.split('\n\n')[1], bait_sha256=C.sha256_file(K.ST / 'work/bait/bait_kmers.txt'),
        r1=str(out / 'r1.fq.gz'), r2=str(out / 'r2.fq.gz'),
        r1_sha256=C.sha256_file(out / 'r1.fq.gz'), r2_sha256=C.sha256_file(out / 'r2.fq.gz'),
        truth_loaded=False, job=os.environ['SLURM_JOB_ID']))


if __name__ == '__main__':
    main()
