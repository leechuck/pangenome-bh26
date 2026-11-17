"""All-gene driver for the patch builder (stage mosaic-allgenes-20261114). One task = one read set x
donor x panel arm; reads are mapped once to the pass-1 alleles of all genes (as before), then the
builder (drb1_run.run, gene as a parameter, no T1K second pass) runs for each gene. No truth is read
beyond the frozen setup's leakage check.

  python ag_run.py READSET DONOR ARM TAG [--genes A,B,...] [--confounded-dir DIR]

ARM: C (IPD panel, on top of arm A) or D (IPD + fold cohort panel, on top of arm B).
TAG: output folder name under work/ (e.g. learn, main).
Output: work/TAG/READSET/DONOR/ARM/<gene>/{proposals.json,OVERRIDES*.json,additions.fa,COMPLETE.json|FAILURE.json}
        work/TAG/READSET/DONOR/ARM/COMPLETE.json | FAILURE.json (task level, per-gene status and seconds)
"""
from __future__ import annotations

import faulthandler
import os
import shutil
import signal
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_run as R          # noqa: E402
import drb1_diag2 as D2       # noqa: E402
C = R.C
S2 = R.S2
DEV = ('D-sim', 'D-real')
GENE_TIMEOUT = 1800     # seconds per gene (wall clock); on expiry the gene falls back to its pass-1 call


def _alarm(signum, frame):
    raise TimeoutError(f'gene exceeded {GENE_TIMEOUT} s')


def main():
    readset, donor, arm, tag = sys.argv[1:5]
    genes = list(C.GENES)
    cdir = None
    a = sys.argv[5:]
    if '--genes' in a:
        genes = a[a.index('--genes') + 1].split(',')
    if '--confounded-dir' in a:
        cdir = Path(a[a.index('--confounded-dir') + 1])
    assert os.environ.get('SLURM_JOB_ID'), 'run under Slurm'
    assert readset in DEV or (S2 / 'FREEZE.json').exists(), 'test read sets only after FREEZE.json'
    assert arm in ('C', 'D') and not (arm == 'D' and readset == 'D-sim'), 'no cohort-panel run for D-sim'
    faulthandler.dump_traceback_later(600, repeat=True, file=sys.stderr)   # stacks in the log if slow
    signal.signal(signal.SIGALRM, _alarm)
    t_v = time.time()
    R.RD.verify_runtime()
    t_v = round(time.time() - t_v, 2)
    out = S2 / 'work' / tag / readset / donor / arm
    if (out / 'COMPLETE.json').exists():
        return
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    code = {p.name: C.sha256_file(p) for p in sorted(Path(__file__).resolve().parent.glob('*.py'))}
    rec = dict(status='running', readset=readset, donor=donor, arm=arm, tag=tag, job_id=os.environ['SLURM_JOB_ID'],
               task_id=os.environ.get('SLURM_ARRAY_TASK_ID'), code_sha256=code, started=time.time(), genes={}, verify_runtime_seconds=t_v)
    try:
        t0 = time.time()
        signal.alarm(GENE_TIMEOUT)
        try:
            shared = D2.load_shared(readset, donor, out, arm)
        finally:
            signal.alarm(0)
        rec['shared_seconds'] = round(time.time() - t0, 2)
        rec['map_stats'] = shared['mstats']
        for g in genes:
            go = out / g
            go.mkdir()
            t1 = time.time()
            grec = dict(gene=g, started=t1)
            signal.alarm(GENE_TIMEOUT)
            try:
                cf = (cdir / f'confounded_{g}.json') if cdir else None
                res = R.run(readset, donor, go, False, arm, gene=g, shared=shared, confounded_path=cf)
                grec.update(res, status='complete', seconds=round(time.time() - t1, 2))
                C.save_json(go / 'COMPLETE.json', grec)
            except Exception as e:      # per-gene failure is recorded, the other genes still run
                kind = ('bound' if isinstance(e, R.PT.BoundExceeded) else 'timeout' if isinstance(e, TimeoutError)
                        else 'no_pass1_call' if 'no pass-1' in str(e) else 'error')
                grec.update(status='failed', failure_kind=kind, fallback='pass-1 call kept (no override written)',
                            error=repr(e), traceback=traceback.format_exc(), seconds=round(time.time() - t1, 2),
                            work=dict(R.PT.WORK))
                for f in go.glob('OVERRIDES*.json'):
                    f.unlink()
                C.save_json(go / 'FAILURE.json', grec)
            finally:
                signal.alarm(0)
            rec['genes'][g] = dict(status=grec['status'], seconds=grec['seconds'], error=grec.get('error'),
                                   failure_kind=grec.get('failure_kind'))
        rec.update(status='complete', finished=time.time())
        C.save_json(out / 'COMPLETE.json', rec)
    except BaseException as e:
        rec.update(status='failed', error=repr(e), traceback=traceback.format_exc(), finished=time.time(),
                   failure_kind='timeout' if isinstance(e, TimeoutError) else 'error')
        try:   # explicit fallback: the pass-1 calls of this arm stand for every gene
            spec = R.RD.pass1_spec(readset, donor, arm)
            rec.update(fallback='pass-1 calls kept for all genes', fallback_calls_table=str(C.verified_table(spec['run'])),
                       fallback_pass1_reference=str(spec['ref']))
        except Exception as e2:
            rec['fallback_error'] = repr(e2)
        C.save_json(out / 'FAILURE.json', rec)
        raise


if __name__ == '__main__':
    main()
