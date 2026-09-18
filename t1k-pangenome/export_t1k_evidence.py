"""Isolated T1K 1.0.6 instrumentation: export existing fragment weights only.

The seven columns are fragment, allele, start, end, weight, qual, adjustWeight.
They are pre-coalescing assignment values, not posterior allele probabilities.
No alignment, filtering, quantification or genotype logic is changed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from build_graph import sha

SOURCE_SHA = 'b45577d98f93383ab769cd18a5e6a6d623b63273a67e0841df0c5ca3c9098c85'


def patch(source):
    if hashlib.sha256(source.encode()).hexdigest() != SOURCE_SHA:
        raise ValueError('Expected pinned unmodified T1K 1.0.6 source')
    old_format = '"%s\\t%s\\t%d\\t%d\\n"'
    new_format = '"%s\\t%s\\t%d\\t%d\\t%.9g\\t%.9g\\t%.9g\\n"'
    old_args = 'assignments[j].start, assignments[j].end) ;'
    new_args = ('assignments[j].start, assignments[j].end, assignments[j].weight, '
                'assignments[j].qual, assignments[j].adjustWeight) ;')
    if source.count(old_format)!=2 or source.count(old_args)!=2:
        raise ValueError('Unexpected assignment export sites')
    return source.replace(old_format,new_format).replace(old_args,new_args)


def build(source, output, reads, reference):
    if not os.environ.get('SLURM_CPUS_PER_TASK'):
        raise RuntimeError('Build and smoke test require Slurm')
    if output.exists():
        raise FileExistsError(output)
    original = (source/'Genotyper.cpp').read_text()
    instrumented = patch(original)
    compiler = shutil.which('g++')
    if compiler is None:
        raise FileNotFoundError('g++ compiler unavailable')
    output.mkdir(parents=True)
    record = dict(status='running', source_sha256=SOURCE_SHA,
                  driver_sha256=sha(Path(__file__)), compiler=compiler,
                  compiler_version=subprocess.check_output([compiler,'--version'],text=True),
                  input_sha256={str(p):sha(p) for p in [reference,reads/'r1.fq',reads/'r2.fq']},
                  scope='Export-only instrumentation and synthetic preservation check')
    def save():
        (output/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        for variant,code in (('original',original),('export',instrumented)):
            folder = output/variant
            folder.mkdir()
            for p in source.iterdir():
                if p.suffix in ('.h','.hpp'):
                    shutil.copy2(p,folder/p.name)
            (folder/'Genotyper.cpp').write_text(code)
            with (folder/'build.log').open('w') as log:
                subprocess.run([compiler,'-O3','-o',str(folder/'genotyper'),
                                str(folder/'Genotyper.cpp'),'-lpthread','-lz'],
                               stdout=log,stderr=subprocess.STDOUT,check=True)
            command = [str(folder/'genotyper'),'-f',str(reference),
                       '-1',str(reads/'r1.fq'),'-2',str(reads/'r2.fq'),'-s','0.97',
                       '--alleleDigitUnits','4','--alleleDelimiter',':',
                       '--outputReadAssignment','-t','2','-o',str(folder/'synthetic')]
            with (folder/'run.log').open('w') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        original_dir, exported_dir = output/'original', output/'export'
        genotype = 'synthetic_genotype.tsv'
        if (original_dir/genotype).read_bytes() != (exported_dir/genotype).read_bytes():
            raise ValueError('Instrumentation changed synthetic genotype output')
        native = (original_dir/'synthetic_assign.tsv').read_text().splitlines()
        exported = (exported_dir/'synthetic_assign.tsv').read_text().splitlines()
        if not exported or native != ['\t'.join(line.split('\t')[:4]) for line in exported]:
            raise ValueError('Instrumentation changed assignment rows')
        if any(len(line.split('\t'))!=7 for line in exported):
            raise ValueError('Invalid evidence columns')
        record.update(status='complete',assignment_rows=len(exported),
                      genotype_unchanged=True,assignment_coordinates_unchanged=True,
                      output_sha256={str(p.relative_to(output)):sha(p) for p in
                         [original_dir/'genotyper',exported_dir/'genotyper',
                          original_dir/genotype,exported_dir/genotype,
                          exported_dir/'synthetic_assign.tsv']})
        save()
        (output/'COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')
    except BaseException as error:
        record.update(status='failed',error=repr(error))
        save()
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('source','output','reads','reference'):
        p.add_argument('--'+name,type=Path,required=True)
    a = p.parse_args()
    build(a.source,a.output,a.reads,a.reference)
