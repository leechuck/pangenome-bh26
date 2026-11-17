"""Step 2b wgsx check (development donor): compare this stage's whole-genome T1K candidate extraction
(work/wgsx_check/DONOR) with the baseline run's candidate files, as sets of (name, read 1, read 2). Counts only.
  python wgsx_cmp.py DONOR   -> work/wgsx_check/DONOR/CMP.json"""
import hashlib
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common3 as K          # noqa: E402
import select3 as S          # noqa: E402
C = K.C


def pairs(d):
    out = set()
    for (h1, s1, q1), (h2, s2, q2) in zip(S.fq(d / 'calls_candidate_1.fq'), S.fq(d / 'calls_candidate_2.fq')):
        out.add(hashlib.sha1(f'{S.qname(h1)}\t{s1.strip()}\t{s2.strip()}'.encode()).digest())
    return out


donor = sys.argv[1]
assert os.environ.get('SLURM_JOB_ID') and donor in {l.split()[1] for l in open(K.ST / 'src/tasks_cov.txt')}
a = pairs(K.ST / 'work/wgsx_check' / donor)
b = pairs(C.ZEN / 'baselines/results' / f'wgs.{donor}' / 't1k/t1k_genomic_ipd')
C.save_json(K.ST / 'work/wgsx_check' / donor / 'CMP.json',
            dict(donor=donor, new_pairs=len(a), baseline_pairs=len(b), only_new=len(a - b), only_baseline=len(b - a),
                 identical_sets=a == b))
print(len(a), len(b), len(a - b), len(b - a))
