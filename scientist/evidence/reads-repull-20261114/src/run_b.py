"""Part B step 7, builder arm C3/C4 on re-pulled reads (D-real2), frozen code unchanged:
mosaic-drb1-20261113/src/drb1_run.run() is imported and called as is. Only redirections:
  - run_donor.pass1_spec('D-real2', donor, 'C') returns this stage's arm-A run (work/D-real2/DONOR/A);
  - this stage's folder is bound into the container (identity bind) for the pass-2 T1K run;
  - output folder work/D-real2/DONOR/C3 (same layout as the frozen stage).
The confounded-event list is still read from the frozen stage (drb1_run.S2/inputs).

  python run_b.py DONOR
"""
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113/src')
import drb1_run as DR       # noqa: E402
RD, C = DR.RD, DR.C

ST = Path('/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114')
RD.BIND[str(ST)] = ST
RS = os.environ.get('RS', 'D-real2')        # re-pulled read set: D-real2 (rule 1) or D-real3 (rule 2)
READS = {'D-real2': 'reads2', 'D-real3': 'reads3'}[RS]
_orig = RD.pass1_spec


def pass1_spec(readset, donor, arm):
    if readset != RS:
        return _orig(readset, donor, arm)
    assert arm in ('A', 'C')
    run = ST / 'work' / RS / donor / 'A'
    rec = C.load_json(run / 'COMPLETE.json')
    assert rec['status'] == 'complete' and rec['donor'] == donor and not rec.get('truth_loaded')
    cmd = RD.parse_cmd(rec['command'])
    return dict(run=run, ref=RD.to_host(cmd['-f']), r1=cmd['-1'], r2=cmd['-2'], threads=int(cmd['-t']),
                aliases=None, receipt=rec)


RD.pass1_spec = pass1_spec


def main():
    donor = sys.argv[1]
    assert os.environ.get('SLURM_JOB_ID')
    dev = {l.split()[1] for l in open(DR.S2 / 'work/dev_tasks.txt')}
    assert donor in dev, 'development donors only'
    RD.verify_runtime()
    out = ST / 'work' / RS / donor / 'C3'
    if (out / 'COMPLETE.json').exists():
        return
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rec = dict(status='running', readset=RS, donor=donor, arm='C3', job_id=os.environ['SLURM_JOB_ID'],
               code_sha256={p.name: C.sha256_file(p) for p in sorted((DR.S2 / 'src').glob('*.py'))},
               wrapper_sha256=C.sha256_file(Path(__file__)), started=time.time())
    try:
        rec.update(DR.run(RS, donor, out, True, 'C'), status='complete', finished=time.time())
        C.save_json(out / 'COMPLETE.json', rec)
    except BaseException as e:
        rec.update(status='failed', error=repr(e), traceback=traceback.format_exc(), finished=time.time())
        C.save_json(out / 'FAILURE.json', rec)
        raise


if __name__ == '__main__':
    main()
