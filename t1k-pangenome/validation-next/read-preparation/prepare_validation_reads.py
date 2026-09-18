"""Recruit reserved public reads and verify pairing without accessing truth."""
import argparse
import csv
import gzip
import itertools
import json
import os
from pathlib import Path
import subprocess
import time
from build_graph import sha


def fastq(stream):
    while True:
        name = stream.readline()
        if not name:
            return
        sequence, plus, quality = [stream.readline().rstrip('\r\n') for _ in range(3)]
        if not name.startswith('@') or not plus.startswith('+') or not sequence or len(sequence)!=len(quality):
            raise ValueError('Malformed FASTQ record')
        identifier = name[1:].split()[0]
        if identifier.endswith(('/1','/2')):
            identifier = identifier[:-2]
        yield identifier


def verify_pairs(first, second):
    count = 0
    for a,b in itertools.zip_longest(fastq(first),fastq(second)):
        if a is None or b is None or a != b:
            raise ValueError('FASTQ mates differ in names or counts')
        count += 1
    if not count:
        raise ValueError('No paired reads recruited')
    return count


def prepare(root, index):
    if not os.environ.get('SLURM_CPUS_PER_TASK'):
        raise RuntimeError('Requires Slurm allocation')
    inputs = root/'validation'
    lock = json.loads((inputs/'RESERVATION.json').read_text())
    if sha(inputs/'cohort.tsv') != lock['cohort_sha256']:
        raise ValueError('Reserved cohort changed')
    rows = list(csv.DictReader((inputs/'recruit.tsv').open(),delimiter='\t'))
    cohort = list(csv.DictReader((inputs/'cohort.tsv').open(),delimiter='\t'))
    availability = json.loads((inputs/'read-availability.json').read_text())
    if availability['cohort_sha256'] != lock['cohort_sha256']:
        raise ValueError('Read availability audit belongs to a different cohort')
    sources = {r['donor']:r['cram_url'] for r in availability['donors']}
    if [r['donor'] for r in rows] != [r['donor'] for r in cohort] or any(r['cram']!=sources[r['donor']] for r in rows):
        raise ValueError('Recruitment donor/source list differs from reserved audit')
    if not 0 <= index < len(rows):
        raise ValueError('Invalid cohort index')
    donor = rows[index]['donor']
    out = root/'validation-reads'/donor
    out.mkdir(parents=True,exist_ok=True)
    if (out/'VERIFIED.json').exists():
        raise FileExistsError('Already verified: '+donor)
    record = dict(status='running',started=time.time(),donor=donor,cram=rows[index]['cram'],
                  cohort_sha256=lock['cohort_sha256'],script_sha256=sha(inputs/'recruit.sh'),
                  driver_sha256=sha(Path(__file__)),scope='Read preparation only; no typing or truth access')
    def save():
        (out/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        env=dict(os.environ,DONORS=str(inputs/'recruit.tsv'),SLURM_ARRAY_TASK_ID=str(index))
        subprocess.run(['bash',str(inputs/'recruit.sh')],env=env,check=True)
        header = subprocess.check_output(['samtools','view','-H',str(out/'recruited.cram')],text=True)
        samples = {field[3:] for line in header.splitlines() if line.startswith('@RG\t')
                   for field in line.split('\t') if field.startswith('SM:')}
        if samples != {donor}:
            raise ValueError('CRAM sample identity differs from reserved donor: '+repr(samples))
        with gzip.open(out/'r1.fq.gz','rt') as first, gzip.open(out/'r2.fq.gz','rt') as second:
            pairs = verify_pairs(first,second)
        record.update(status='complete',finished=time.time(),pairs=pairs,sample_ids=sorted(samples),
                      output_sha256={name:sha(out/name) for name in ('r1.fq.gz','r2.fq.gz')})
        save()
        (out/'VERIFIED.json').write_text(json.dumps(record,indent=2)+'\n')
    except Exception as error:
        record.update(status='failed',finished=time.time(),error=repr(error));save()
        raise


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--index',type=int,default=None)
    a=p.parse_args()
    prepare(a.root,a.index if a.index is not None else int(os.environ['SLURM_ARRAY_TASK_ID']))
