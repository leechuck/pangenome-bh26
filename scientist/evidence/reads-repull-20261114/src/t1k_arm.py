"""Step 2b: T1K pass 1 on re-pulled reads, arm A (IPD genomic.fa) or arm B (the donor's existing fold
reference genomic_core_hprc_asian: same FASTA file, checked by sha256, as the old arm-B run of the donor).
Same command as mosaic-typer run_donor.t1k (hla-wgs preset, same thread count as the old run of that arm).

  python t1k_arm.py RS DONOR ARM      RS: D-real3 | T-real3, ARM: A | B  -> work/RS/DONOR/ARM/
"""
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common3 as K         # noqa: E402
C, RD = K.C, K.RD


def main():
    rs, donor, arm = sys.argv[1:4]
    assert os.environ.get('SLURM_JOB_ID') and arm in ('A', 'B') and rs in K.OLD
    K.require_freeze(rs)
    out = K.run_dir(rs, donor, arm)
    if (out / 'COMPLETE.json').exists():
        return
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rec = dict(status='running', donor=donor, readset=rs, arm=arm, job_id=os.environ['SLURM_JOB_ID'],
               task_id=os.environ.get('SLURM_ARRAY_TASK_ID'), started=time.time())
    try:
        RD.verify_runtime()
        old = K._orig_pass1(K.OLD[rs], donor, 'A' if arm == 'A' else 'D')
        ref = old['ref']
        if arm == 'B':
            assert C.sha256_file(ref) == old['receipt']['reference_sha256']
            assert (ref.parent / 'aliases.json').exists()
        else:
            assert ref == C.GENOMIC_FA and C.sha256_file(ref) == C.GENOMIC_SHA
        reads = K.READS[rs] / donor
        pr = C.load_json(reads / 'PAIRS.json')
        assert pr['status'] == 'paired_fastq_verified' and pr['donor'] == donor
        r1, r2 = str(reads / 'r1.fq.gz'), str(reads / 'r2.fq.gz')
        assert C.sha256_file(r1) == pr['r1_sha256'] and C.sha256_file(r2) == pr['r2_sha256']
        cmd = RD.t1k(r1, r2, ref, out, old['threads'], out / 'run.log')
        rec.update(status='complete', panel='genomic' if arm == 'A' else 'genomic_core_hprc_asian',
                   command=cmd[cmd.index(RD.RUN_T1K):], truth_loaded=False, reference_sha256=C.sha256_file(ref),
                   old_run=str(old['run']), native_table_sha256=C.sha256_file(out / 'calls_genotype.tsv'),
                   finished=time.time())
        C.save_json(out / 'COMPLETE.json', rec)
    except BaseException as e:
        rec.update(status='failed', error=repr(e), traceback=traceback.format_exc(), finished=time.time())
        C.save_json(out / 'FAILURE.json', rec)
        raise


if __name__ == '__main__':
    main()
