"""Part B step 7, arm A on re-pulled reads (D-real2): T1K with the IPD full-gene reference, the same
command as mosaic-typer run_donor.t1k (genomic.fa, hla-wgs preset, 4 threads, 4-digit units).
Only change: this stage's folder is bound into the container so that the new reads can be read.

  python t1k_a.py DONOR      -> work/D-real2/DONOR/A/{calls_*,COMPLETE.json}
"""
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src')
import common as C          # noqa: E402
import run_donor as RD      # noqa: E402

ST = Path('/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114')
RD.BIND[str(ST)] = ST       # identity bind for this stage
RS = os.environ.get('RS', 'D-real2')        # re-pulled read set: D-real2 (rule 1) or D-real3 (rule 2)
READS = {'D-real2': 'reads2', 'D-real3': 'reads3'}[RS]


def main():
    donor = sys.argv[1]
    assert os.environ.get('SLURM_JOB_ID')
    out = ST / 'work' / RS / donor / 'A'
    if (out / 'COMPLETE.json').exists():
        return
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rec = dict(status='running', donor=donor, panel='genomic', readset=RS, job_id=os.environ['SLURM_JOB_ID'],
               started=time.time())
    try:
        RD.verify_runtime()
        reads = ST / 'work' / READS / donor
        pr = C.load_json(reads / 'PAIRS.json')
        assert pr['status'] == 'paired_fastq_verified' and pr['donor'] == donor
        r1, r2 = str(reads / 'r1.fq.gz'), str(reads / 'r2.fq.gz')
        assert C.sha256_file(r1) == pr['r1_sha256'] and C.sha256_file(r2) == pr['r2_sha256']
        cmd = RD.t1k(r1, r2, C.GENOMIC_FA, out, 4, out / 'run.log')
        rec.update(status='complete', command=cmd[cmd.index(RD.RUN_T1K):], truth_loaded=False,
                   reference_sha256=C.sha256_file(C.GENOMIC_FA),
                   native_table_sha256=C.sha256_file(out / 'calls_genotype.tsv'), finished=time.time())
        C.save_json(out / 'COMPLETE.json', rec)
    except BaseException as e:
        rec.update(status='failed', error=repr(e), traceback=traceback.format_exc())
        C.save_json(out / 'FAILURE.json', rec)
        raise


if __name__ == '__main__':
    main()
