"""Step 2b: access to the FROZEN step-4 stage mosaic-allgenes-20261114 (read-only).
verify() checks every file listed in its FREEZE.json against the recorded sha256 and fails otherwise.
  python frozen4.py   -> work/VERIFY_STEP4.json"""
import hashlib
import json
import os
import sys
from pathlib import Path

S4 = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114')
ST = Path('/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114')


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 23), b''):
            h.update(c)
    return h.hexdigest()


def verify():
    fz = json.loads((S4 / 'FREEZE.json').read_text())
    bad = {k: v for k, v in fz['sha256'].items() if sha(S4 / k) != v}
    assert not bad, f'step-4 frozen files changed: {sorted(bad)}'
    return dict(freeze_json_sha256=sha(S4 / 'FREEZE.json'), frozen=fz['frozen'], files_checked=len(fz['sha256']))


def use():
    """verify, then put the frozen step-4 src first on sys.path"""
    v = verify()
    sys.path.insert(0, str(S4 / 'src'))
    return v


if __name__ == '__main__':
    assert os.environ.get('SLURM_JOB_ID')
    v = verify()
    (ST / 'work/VERIFY_STEP4.json').write_text(json.dumps(v, indent=1))
    print(v)
