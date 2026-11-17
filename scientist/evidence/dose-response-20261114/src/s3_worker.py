"""Per-donor worker: runs the listed cells of one donor with T1K (reference content is the only thing that
changes).  python s3_worker.py --manifest M --index I      (index defaults to SLURM_ARRAY_TASK_ID)

Manifest: {"tasks": [{"donor": D, "cells": [names]}]}. Cell names: REPRO_A (existing arm-A reference, rerun),
REPRO_B (existing arm-B reference, rerun), or a panel cell from work/DRAWS.json.
Every cell writes work/<donor>/<cell>/COMPLETE.json or FAILURE.json.  Candidate-read and aligned-read files
are deleted after the run (tables, vcf, log and additions.fa are kept).
"""
import argparse
import collections
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s3_lib as L  # noqa: E402

C = L.C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--index', type=int, default=None)
    a = ap.parse_args()
    assert os.environ.get('SLURM_JOB_ID'), 'run under Slurm'
    idx = a.index if a.index is not None else int(os.environ['SLURM_ARRAY_TASK_ID'])
    task = C.load_json(a.manifest)['tasks'][idx]
    donor = task['donor']
    version = L.verify_runtime()
    draws = C.load_json(L.STAGE / 'work/DRAWS.json')[donor]
    runs = L.donor_runs()[donor]
    ctx = {}
    failed = []
    for cell in task['cells']:
        out = L.STAGE / 'work' / donor / cell
        if (out / 'COMPLETE.json').exists():
            continue
        if out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True)
        rec = dict(status='running', donor=donor, cell=cell, readset=runs['readset'],
                   job_id=os.environ.get('SLURM_JOB_ID'), task=idx, t1k_version=version,
                   r1=runs['r1'], r2=runs['r2'], started=time.time())
        try:
            threads = runs['threads']
            if cell == 'REPRO_A':
                ref = L.to_host(runs['armA_ref'])
                want = runs['armA_table_sha256']
                rec.update(reference_sha256=C.sha256_file(ref), n_added=0)
            elif cell == 'REPRO_B':
                b = L.armB_spec(runs['readset'], donor)
                assert (b['r1'], b['r2']) == (runs['r1'], runs['r2']), 'arm A and B used different reads'
                ref, want = Path(b['ref']), b['table_sha256']
                rec.update(reference_sha256=C.sha256_file(ref), n_added=None)
            else:
                d = draws[cell]
                if 'records' not in ctx:
                    ctx['records'] = C.load_records()
                    ctx['by_donor'] = collections.defaultdict(list)
                    for r in ctx['records']:
                        ctx['by_donor'][r['donor']].append(r)
                    ctx['ipd'] = L.ipd_sequence_set()
                    ctx['na'] = None
                if d['source'] == L.NONASIAN:
                    if ctx['na'] is None:
                        ctx['na'] = L.load_nonasian()
                    src = ctx['na']
                else:
                    src = ctx['by_donor']
                assert donor not in d['drawn']
                recs = [r for x in d['drawn'] for r in src[x]]
                adds, st = L.build_additions(recs, ctx['ipd'])
                (out / 'additions.fa').write_text(''.join(x['fa'] for x in adds))
                ref = out / 'reference.fa'
                with open(ref, 'wb') as f:
                    with open(C.GENOMIC_FA, 'rb') as g0:
                        shutil.copyfileobj(g0, f)
                    f.write((out / 'additions.fa').read_bytes())
                want = None
                rec.update(reference_sha256=C.sha256_file(ref), n_added=len(adds), build_stats=st,
                           source=d['source'], n_requested=d['n_requested'], n_drawn=d['n_drawn'],
                           drawn_donors=d['drawn'])
            t0 = time.time()
            cmd = L.t1k(runs['r1'], runs['r2'], ref, out, threads, out / 'run.log')
            rec['t1k_seconds'] = round(time.time() - t0, 1)
            if cell not in ('REPRO_A', 'REPRO_B'):
                ref.unlink()
            for p in list(out.glob('calls_candidate_*.fq')) + list(out.glob('calls_aligned_*.fa')):
                p.unlink()
            got = C.sha256_file(out / 'calls_genotype.tsv')
            rec.update(command=cmd, table_sha256=got, threads=threads)
            if want is not None:
                rec.update(original_table_sha256=want, identical=got == want)
            rec.update(status='complete', finished=time.time())
            C.save_json(out / 'COMPLETE.json', rec)
        except BaseException as e:  # noqa: BLE001
            rec.update(status='failed', finished=time.time(), error=repr(e), traceback=traceback.format_exc())
            C.save_json(out / 'FAILURE.json', rec)
            failed.append(cell)
            print(traceback.format_exc(), file=sys.stderr)
    if failed:
        raise SystemExit('failed cells: ' + ','.join(failed))


if __name__ == '__main__':
    main()
