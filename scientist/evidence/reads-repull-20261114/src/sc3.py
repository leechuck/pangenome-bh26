"""Step 2b: frozen step-4 scorer (mosaic-allgenes-20261114/src/ag_score.py, verified against its FREEZE.json,
imported unchanged) for re-pulled read sets, plus a reproduction mode for the old outputs.
Redirections only (no scoring logic changed):
  - truth copies of D-real3 / T-real3 = those of D-real / T-real (same donors; frozen score.Truth);
  - arm A/B tables of D-real3 / T-real3 = this stage's T1K runs (common3.py);
  - builder outputs: re-pulled read sets from this stage (work/ag/work/TAG/...); old read sets from the
    frozen step-4 stage (work/OLDTAG/...), read-only;
  - output: work/ag/results/score_rows_OUTTAG.tsv, score_avail_OUTTAG.tsv.

  python sc3.py OUTTAG TAG OLDTAG READSETS ARMS
  e.g. python sc3.py repro_old - test T-real A,B,C4,C4R4n,C4R1n,D4,D4R4n,D4R1n
"""
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import frozen4               # noqa: E402
frozen4.use()
import common3 as K          # noqa: E402
import ag_score as AS        # noqa: E402  (frozen step-4 file)
assert Path(AS.__file__).resolve().parent == frozen4.S4 / 'src'
SC, C = AS.SC, AS.C

outtag, tag, oldtag, readsets, arms = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
for rs in readsets.split(','):
    K.require_freeze(rs)
assert os.environ.get('SLURM_JOB_ID')

_Truth = SC.Truth


class Truth3(_Truth):
    def __init__(self, rsl):
        base = sorted({K.OLD.get(r, r) for r in rsl})
        super().__init__(base)
        for (r0, d, g), cps in list(self.copies.items()):
            for r in rsl:
                if K.OLD.get(r) == r0:
                    self.copies[(r, d, g)] = [dict(c, readset=r) for c in cps]
        for k in [k for k in self.copies if k[0] not in rsl]:
            del self.copies[k]


_arm_table = SC.arm_table


def arm_table(rs, d, arm):
    if rs in K.OLD and arm in ('A', 'B'):
        run = K.run_dir(rs, d, arm)
        if arm == 'A':
            return C.verified_table(run), []
        spec = K.pass1_spec(rs, d, 'D')
        return C.verified_table(run), [spec['ref']]
    return _arm_table(rs, d, arm)


def task_dir(tag_, rs, d, arm):
    if rs in K.OLD:
        return K.ST / 'ag/work' / tag / rs / d / arm[0]
    return frozen4.S4 / 'work' / oldtag / rs / d / arm[0]


SC.Truth = Truth3
SC.arm_table = arm_table
AS.task_dir = task_dir
AS.TEST = ()                  # this stage's own gate: require_freeze above
AS.S2 = K.ST / 'ag'
(AS.S2 / 'results').mkdir(parents=True, exist_ok=True)
sys.argv = ['ag_score.py', outtag, readsets, arms]
AS.main()
