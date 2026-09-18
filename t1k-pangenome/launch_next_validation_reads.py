"""Stage and submit truth-blind read preparation for the reserved HGSVC cohort."""
import csv
import json
from pathlib import Path
import shlex
import subprocess
import time
from build_graph import sha
from launch_development_evidence import B, HOST, remote, submit

HERE = Path(__file__).resolve().parent
ROOT = '/home/leechuck/hla/t1k-pangenome/hgsvc-validation-v1'


def main():
    source = HERE/'validation-next'
    ledger = source/'READ_PREPARATION_LAUNCH.json'
    if ledger.exists():
        raise FileExistsError('Reconcile the existing submission before any retry')
    lock = json.loads((source/'RESERVATION.json').read_text())
    assert sha(source/'cohort.tsv') == lock['cohort_sha256']
    for name, digest in lock['audit_sha256'].items():
        assert sha(source/name) == digest, name
    cohort = list(csv.DictReader((source/'cohort.tsv').open(), delimiter='\t'))
    audit = json.loads((source/'READ_AVAILABILITY.json').read_text())
    assert [r['donor'] for r in cohort] == [r['donor'] for r in audit['donors']]
    assert all(r[k]['accessible'] for r in audit['donors'] for k in ('cram', 'crai'))
    stage = source/'read-preparation'
    stage.mkdir(exist_ok=False)
    for name in ('RESERVATION.json', 'cohort.tsv'):
        (stage/name).write_bytes((source/name).read_bytes())
    audit['cohort_sha256'] = lock['cohort_sha256']
    audit['source_audit_sha256'] = sha(source/'READ_AVAILABILITY.json')
    (stage/'read-availability.json').write_text(json.dumps(audit, indent=2)+'\n')
    with (stage/'recruit.tsv').open('w') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(['donor', 'cram'])
        writer.writerows((r['donor'], r['cram_url']) for r in audit['donors'])
    script = (HERE/'validation/recruit.sh').read_text()
    old = 'out=/home/leechuck/hla/t1k-pangenome/validation-reads/$donor'
    assert script.count(old) == 1
    script = script.replace(old, 'out='+ROOT+'/validation-reads/$donor')
    # Each task writes its own region list, avoiding a concurrent shared-file writer.
    script = script.replace('bed=reads/regions.bed', 'bed="$out/regions.bed"')
    (stage/'recruit.sh').write_text(script)
    for name in ('prepare_validation_reads.py', 'build_graph.py'):
        (stage/name).write_bytes((HERE/name).read_bytes())
    hashes = {p.name:sha(p) for p in stage.iterdir()}
    destination = B+'/hla/t1k-pangenome/hgsvc-validation-v1/validation'
    remote('mkdir -p '+shlex.quote(destination))
    subprocess.run(['scp', *map(str,stage.iterdir()), HOST+':'+destination+'/'], check=True)
    actual = remote('sha256sum '+' '.join(shlex.quote(destination+'/'+name) for name in hashes))
    assert {Path(line.split()[1]).name:line.split()[0] for line in actual.splitlines()} == hashes
    record = dict(status='submitting', created=time.time(), root=ROOT,
                  scope='Read preparation only; no inference or truth access', donors=len(cohort),
                  code_and_input_sha256=hashes, job_name='t1k-hgsvc28-reads')
    def save():
        ledger.write_text(json.dumps(record, indent=2)+'\n')
    save()
    job = submit(record['job_name'],4,'8G','04:00:00',
                 ['python3', ROOT+'/validation/prepare_validation_reads.py','--root',ROOT,
                  '--index','$SLURM_ARRAY_TASK_ID'], ['--array=0-'+str(len(cohort)-1)+'%28'])
    record.update(status='submitted', job=job)
    save()
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
