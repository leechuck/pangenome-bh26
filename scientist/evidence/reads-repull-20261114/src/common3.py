"""Step 2b shared definitions: re-pulled read sets (rule 2) and where their files live.
D-real3 = 20 development donors (dry run); T-real3 = 66 test donors (pre-registered test)."""
import sys
from pathlib import Path

sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src')
import common as C          # noqa: E402
import run_donor as RD      # noqa: E402

ST = Path('/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114')
OLD = {'D-real3': 'D-real', 'T-real3': 'T-real'}          # old regional read set of the same donors
READS = {'D-real3': ST / 'work/reads3', 'T-real3': ST / 'work/reads_T-real3'}
RUNS = {'D-real3': ST / 'work/D-real3', 'T-real3': ST / 'work/T-real3'}
FREEZE = ST / 'FREEZE_TEST.json'
RD.BIND[str(ST)] = ST       # identity bind so that T1K in the container can read this stage


def require_freeze(rs):
    if rs == 'T-real3':
        assert FREEZE.exists(), 'T-real3 only after FREEZE_TEST.json'


def run_dir(rs, donor, arm):
    """arm A (IPD genomic) or B (IPD + fold cohort panel) T1K run on the re-pulled reads"""
    return RUNS[rs] / donor / arm


_orig_pass1 = RD.pass1_spec


def pass1_spec(readset, donor, arm):
    if readset not in OLD:
        return _orig_pass1(readset, donor, arm)
    a = 'A' if arm in ('A', 'C', 'E', 'A_repro') else 'B' if arm == 'D' else None
    assert a, arm
    run = run_dir(readset, donor, a)
    rec = C.load_json(run / 'COMPLETE.json')
    assert rec['status'] == 'complete' and rec['donor'] == donor and not rec.get('truth_loaded')
    cmd = RD.parse_cmd(rec['command'])
    ref = RD.to_host(cmd['-f'])
    aliases = None
    if a == 'B':
        aliases = ref.parent / 'aliases.json'
        assert C.sha256_file(ref) == rec['reference_sha256'], 'arm-B reference changed'
    return dict(run=run, ref=ref, r1=cmd['-1'], r2=cmd['-2'], threads=int(cmd['-t']), aliases=aliases, receipt=rec)


RD.pass1_spec = pass1_spec
