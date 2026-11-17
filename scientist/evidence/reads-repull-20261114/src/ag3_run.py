"""Step 2b: frozen step-4 builder (mosaic-allgenes-20261114/src/ag_run.py, verified against its FREEZE.json,
imported unchanged) on re-pulled reads. Only redirections: run_donor.pass1_spec for D-real3/T-real3 returns
this stage's arm-A/B runs (common3.py); output root is this stage (work/ag/work/TAG/RS/DONOR/ARM);
the read set is admitted to ag_run's development list only after this stage's own check
(T-real3 requires FREEZE_TEST.json).

  python ag3_run.py RS DONOR ARM TAG --confounded-dir <step-4 inputs/confounded_v1>
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import frozen4               # noqa: E402
frozen4.use()
import common3 as K          # noqa: E402
import ag_run as AG          # noqa: E402  (frozen step-4 file)

rs = sys.argv[1]
assert rs in K.OLD
K.require_freeze(rs)
assert Path(AG.__file__).resolve().parent == frozen4.S4 / 'src'
AG.S2 = K.ST / 'ag'
AG.DEV = tuple(AG.DEV) + (rs,)
AG.main()
