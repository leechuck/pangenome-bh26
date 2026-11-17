"""Worker for stage thin-data-20261114: one donor.

  thin_donor.py <manifest.txt> [--repro] [--keep-reads]     (task index = SLURM_ARRAY_TASK_ID)

manifest line: "<readset> <donor>".
1. nested subsampling of the donor's full recruited read pair files (hash of read name, seed 20261114,
   keep if hash mod 16 < 8/4/2/1 for depth 1/2, 1/4, 1/8, 1/16);
2. T1K (identical command to run_donor.t1k) for arm A, B (and M for Arab/Saudi donors) at every level
   and N at 1/4 and 1/16, each into work/runs/<rs>/<donor>/<arm>/d<N>/;
3. with --repro: also rerun every arm at FULL depth (4 threads) into work/repro/ and compare the genotype table
   with the existing run.
Every run writes COMPLETE.json or FAILURE.json; the task writes work/tasks/<rs>/<donor>/COMPLETE|FAILURE.json.
"""
import gzip, hashlib, json, os, shutil, sys, time, traceback
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib as L
C = L.C
W = L.TS / 'work'


def sha(p):
    return C.sha256_file(p)


def subsample(r1, r2, outdir: Path):
    """Single pass over the pair files; returns {denominator: {pairs, files, sha}}."""
    outdir.mkdir(parents=True, exist_ok=True)
    f1 = {n: gzip.open(outdir / f'd{n}_r1.fq.gz', 'wt', compresslevel=1) for n in L.LEVELS}
    f2 = {n: gzip.open(outdir / f'd{n}_r2.fq.gz', 'wt', compresslevel=1) for n in L.LEVELS}
    cnt = {n: 0 for n in L.LEVELS}
    total = 0
    with gzip.open(r1, 'rt') as a, gzip.open(r2, 'rt') as b:
        while True:
            x = [a.readline() for _ in range(4)]
            y = [b.readline() for _ in range(4)]
            if not x[0]:
                assert not y[0], 'read files differ in length'
                break
            assert y[0], 'read files differ in length'
            n1, n2 = x[0].split()[0].removesuffix('/1'), y[0].split()[0].removesuffix('/2')
            assert n1 == n2, f'pair names differ: {n1} {n2}'
            total += 1
            k = L.read_key(x[0])
            for n, thr in L.LEVELS.items():
                if k < thr:
                    f1[n].write(''.join(x)); f2[n].write(''.join(y)); cnt[n] += 1
    for n in L.LEVELS:
        f1[n].close(); f2[n].close()
    rec = {'full': dict(pairs=total)}
    for n in L.LEVELS:
        rec[str(n)] = dict(pairs=cnt[n], fraction=cnt[n] / total,
                           sha_r1=sha(outdir / f'd{n}_r1.fq.gz'), sha_r2=sha(outdir / f'd{n}_r2.fq.gz'))
    return rec


def run_one(rs, donor, arm, tag, r1, r2, ref, outdir: Path, compare_sha=None):
    if (outdir / 'COMPLETE.json').exists():
        return True
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)
    rec = dict(readset=rs, donor=donor, arm=arm, level=tag, job_id=os.environ.get('SLURM_JOB_ID'),
               r1=str(r1), r2=str(r2), reference=str(ref), threads=L.THREADS, started=time.time())
    try:
        cmd = L.t1k(r1, r2, ref, outdir, L.THREADS, outdir / 'run.log')
        rec.update(command=cmd, status='complete', finished=time.time(),
                   table_sha256=sha(outdir / 'calls_genotype.tsv'), **L.log_stats(outdir / 'run.log'))
        if compare_sha:
            rec.update(original_table_sha256=compare_sha, identical=rec['table_sha256'] == compare_sha)
        for junk in ('calls_aligned_1.fa', 'calls_aligned_2.fa', 'calls_candidate_1.fq',
                     'calls_candidate_2.fq'):
            p = outdir / junk
            if p.exists():
                if junk.startswith('calls_candidate_1'):
                    rec['candidate_pairs'] = sum(1 for _ in open(p)) // 4
                p.unlink()
        C.save_json(outdir / 'COMPLETE.json', rec)
        return True
    except BaseException as e:
        rec.update(status='failed', finished=time.time(), error=repr(e), traceback=traceback.format_exc())
        C.save_json(outdir / 'FAILURE.json', rec)
        print(traceback.format_exc(), file=sys.stderr)
        return False


def main():
    man = sys.argv[1]
    repro = '--repro' in sys.argv
    keep = '--keep-reads' in sys.argv
    assert os.environ.get('SLURM_JOB_ID'), 'run under Slurm'
    idx = int(os.environ['SLURM_ARRAY_TASK_ID'])
    rs, donor = open(man).read().split('\n')[idx].split()
    tdir = W / 'tasks' / rs / donor
    tdir.mkdir(parents=True, exist_ok=True)
    if (tdir / 'COMPLETE.json').exists() and not repro:
        return
    t0 = time.time()
    task = dict(readset=rs, donor=donor, job_id=os.environ.get('SLURM_JOB_ID'), started=t0, runs={})
    ok = True
    try:
        assert rs in ('D-real', 'T-real', 'T-sim'), rs
        info = {a: L.arm_info(rs, donor, a) for a in L.ARMS_ALL}
        arms = [a for a in L.ARMS_ALL if info[a] is not None]
        assert 'A' in arms and 'B' in arms, 'missing arm A or B receipts'
        for a in arms:
            assert (info[a]['r1'], info[a]['r2']) == (info['A']['r1'], info['A']['r2']), 'arms use different reads'
        rdir = W / 'reads' / rs / donor
        sub = rdir / 'SUBSAMPLE.json'
        if not sub.exists():
            ts = time.time()
            rec = subsample(info['A']['r1'], info['A']['r2'], rdir)
            rec['seconds'] = time.time() - ts
            rec['source_r1'] = str(info['A']['r1']); rec['source_r2'] = str(info['A']['r2'])
            C.save_json(sub, rec)
        task['subsample'] = C.load_json(sub)
        for arm in arms:
            for n in L.levels_for_arm(arm):
                od = W / 'runs' / rs / donor / arm / f'd{n}'
                ok &= run_one(rs, donor, arm, f'1/{n}', rdir / f'd{n}_r1.fq.gz', rdir / f'd{n}_r2.fq.gz',
                              info[arm]['ref'], od)
                task['runs'][f'{arm}-{n}'] = str(od)
            if repro:
                od = W / 'repro' / rs / donor / arm
                ok &= run_one(rs, donor, arm, 'full', info[arm]['r1'], info[arm]['r2'], info[arm]['ref'], od,
                              compare_sha=info[arm]['table_sha'])
                task['runs'][f'{arm}-full-repro'] = str(od)
        if not keep and ok:
            shutil.rmtree(rdir / 'x', ignore_errors=True)
            for p in rdir.glob('d*_r*.fq.gz'):
                p.unlink()
    except BaseException as e:
        ok = False
        task.update(error=repr(e), traceback=traceback.format_exc())
        print(traceback.format_exc(), file=sys.stderr)
    task.update(finished=time.time(), seconds=time.time() - t0)
    if ok:
        task['status'] = 'complete'
        C.save_json(tdir / 'COMPLETE.json', task)
        if (tdir / 'FAILURE.json').exists():
            (tdir / 'FAILURE.json').unlink()
    else:
        task['status'] = 'failed'
        C.save_json(tdir / 'FAILURE.json', task)
        raise SystemExit('task failed')


if __name__ == '__main__':
    main()
