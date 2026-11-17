"""Item 1: reproduce the T-sim APR051 arm-D hang with the FROZEN DRB1 builder (imported read-only from
the DRB1 stage), truth-free and unscored. Times each step and dumps Python stacks every 120 s.
Output goes to this stage only (work/hang/...). No truth is read beyond what the frozen setup()
loads for its leakage check; nothing is scored.
  python hang_probe.py READSET DONOR ARM [frozen|new]"""
from __future__ import annotations

import faulthandler
import functools
import os
import sys
import time
from pathlib import Path

NEW = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114')
FROZEN_SRC = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113/src')
readset, donor, arm = sys.argv[1], sys.argv[2], sys.argv[3]
which = sys.argv[4] if len(sys.argv) > 4 else 'frozen'
assert os.environ.get('SLURM_JOB_ID'), 'run under Slurm'
sys.path.insert(0, str(FROZEN_SRC if which == 'frozen' else NEW / 'src'))
faulthandler.dump_traceback_later(120, repeat=True, file=sys.stderr)

import drb1_run as R          # noqa: E402
import drb1_lib as L          # noqa: E402
import drb1_patch2 as PT      # noqa: E402
import drb1_diag2 as D2       # noqa: E402

R.S2 = NEW
D2.S2 = NEW
T0 = time.time()
LOG = []


def timed(mod, name):
    f = getattr(mod, name)

    @functools.wraps(f)
    def g(*a, **k):
        t = time.time()
        print(f'[{t - T0:8.1f}s] enter {mod.__name__}.{name}', file=sys.stderr, flush=True)
        try:
            return f(*a, **k)
        finally:
            dt = time.time() - t
            LOG.append((name, round(dt, 2)))
            print(f'[{time.time() - T0:8.1f}s] leave {mod.__name__}.{name} {dt:.1f}s', file=sys.stderr, flush=True)
    setattr(mod, name, g)


for mod, names in ((D2, ('setup',)), (L, ('map_reads2', 'project_all', 'make_evidence', 'family_rows', 'blind_intervals')),
                   (PT, ('build_haplotype', 'build_events', 'count_events', 'decide', 'spell', 'structural_filter'))):
    for n in names:
        if hasattr(mod, n):
            timed(mod, n)
for n in ('path_exons', 'path_sequence', 'viterbi'):
    timed(L.PR, n)

out = NEW / 'work' / 'hang' / f'{readset}_{donor}_{arm}_{which}_{os.environ["SLURM_JOB_ID"]}'
import shutil       # noqa: E402
shutil.rmtree(out, ignore_errors=True)
out.mkdir(parents=True)
# the frozen confounded list (the probe reproduces the frozen builder)
cf = NEW / 'inputs/confounded_events.json'
if which == 'frozen' and not cf.exists():
    shutil.copy(NEW / 'inputs/confounded_DRB1_frozen.json', cf)
kw = {'arm': arm}
if which != 'frozen':
    kw['gene'] = 'DRB1'
res = R.run(readset, donor, out, False, **kw) if which == 'frozen' else R.run(readset, donor, out, False, **kw)
import json   # noqa: E402
(out / 'PROBE.json').write_text(json.dumps(dict(result={k: str(v) for k, v in res.items()}, steps=LOG,
                                                total=time.time() - T0), indent=1))
print('done', time.time() - T0, file=sys.stderr)
