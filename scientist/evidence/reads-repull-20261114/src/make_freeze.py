"""Step 2b: write FREEZE_TEST.json (sha256 of src/, bait file, DESIGN_TEST.txt, step-4 FREEZE.json; development
check results; timestamp). Refuses to overwrite.   python make_freeze.py"""
import datetime
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import frozen4               # noqa: E402

ST = frozen4.ST
assert os.environ.get('SLURM_JOB_ID')
out = ST / 'FREEZE_TEST.json'
assert not out.exists(), 'FREEZE_TEST.json exists'
v = frozen4.verify()
files = {f'src/{p.name}': frozen4.sha(p) for p in sorted((ST / 'src').iterdir()) if p.is_file()}
files['work/bait/bait_kmers.txt'] = frozen4.sha(ST / 'work/bait/bait_kmers.txt')
files['DESIGN_TEST.txt'] = frozen4.sha(ST / 'DESIGN_TEST.txt')
checks = {}
for name, p in (('prefilter_check_HG00558', 'work/pref/HG00558/PREF.json'), ('prefilter_check_HG00438', 'work/pref/HG00438/PREF.json'),
                ('wgsx_check_HG00438', 'work/wgsx_check/HG00438/CMP.json'), ('scorer_reproduction_T-real', 'ag/results/CMP_repro_old.json'),
                ('dry_run_old_D-real_rows', 'ag/results/CMP_dry_old_D-real.json')):
    q = ST / p
    checks[name] = json.loads(q.read_text()) if q.exists() else 'missing'
rec = dict(frozen=datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
           note='step 2b: frozen before any T-real3 job; no code change after this file',
           step4_freeze_json_sha256=v['freeze_json_sha256'], step4_files_verified=v['files_checked'],
           sha256=files, development_checks=checks, job=os.environ['SLURM_JOB_ID'])
out.write_text(json.dumps(rec, indent=1))
print(json.dumps(rec, indent=1))
