"""Export native assignments from an existing development baseline's candidate reads."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
from build_graph import sha
from t1k_evidence import audit


def run(baseline, reads, tool, reference, output, threads):
    if not 1 <= threads <= int(os.environ.get('SLURM_CPUS_PER_TASK','0')):
        raise ValueError('Insufficient Slurm allocation')
    if output.exists():
        raise FileExistsError(output)
    original = json.loads((baseline/'manifest.json').read_text())
    config = original['configuration']
    donor = config['donor']
    marker = json.loads((baseline/'COMPLETE').read_text())
    if original['status']!='complete':
        raise ValueError('Baseline incomplete')
    if marker['configuration_sha256']!=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest():
        raise ValueError('Baseline configuration changed')
    if config['alleleDigitUnits']!=4 or config['alleleDelimiter']!=':' or config['preset']!='hla-wgs':
        raise ValueError('Unsupported baseline configuration')
    for name,digest in marker['outputs'].items():
        if sha(baseline/name)!=digest:
            raise ValueError('Baseline output changed')
    raw = {str(i):sha(reads/('r'+str(i)+'.fq.gz')) for i in (1,2)}
    if raw!=config['reads'] or sha(reference)!=config['database_sha256']:
        raise ValueError('Baseline input mismatch')
    build = json.loads((tool/'COMPLETE.json').read_text())
    exe = tool/'export/genotyper'
    if build['status']!='complete' or sha(exe)!=build['output_sha256']['export/genotyper']:
        raise ValueError('Instrumented executable changed')
    candidates = [baseline/(donor+'_candidate_'+str(i)+'.fq') for i in (1,2)]
    output.mkdir(parents=True)
    record = dict(status='running',started=time.time(),donor=donor,
                  baseline_manifest_sha256=sha(baseline/'manifest.json'),
                  baseline_completion_sha256=sha(baseline/'COMPLETE'),
                  instrumented_build_sha256=sha(tool/'COMPLETE.json'),
                  input_sha256={str(reads/('r'+str(i)+'.fq.gz')):raw[str(i)] for i in (1,2)},
                  candidate_reads_sha256={str(p):sha(p) for p in candidates},
                  driver_sha256=sha(Path(__file__)),reference_sha256=sha(reference),
                  scope='Native assignment export; baseline genotype preservation required')
    command = [str(exe),'-s','0.97','--alleleDigitUnits','4','--alleleDelimiter',':',
               '--outputReadAssignment','-f',str(reference),'-1',str(candidates[0]),
               '-2',str(candidates[1]),'-t',str(threads),'-o',str(output/donor)]
    record['command'] = command
    def save():
        (output/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        with (output/'run.log').open('w') as log:
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        genotype = donor+'_genotype.tsv'
        if (output/genotype).read_bytes()!=(baseline/genotype).read_bytes():
            raise ValueError('Instrumented genotype differs from frozen baseline; inspect before using evidence')
        assignment = donor+'_assign.tsv'
        audit(output/assignment,output/'evidence-audit.json')
        record.update(status='complete',finished=time.time(),genotype_unchanged=True,
                      output_sha256={name:sha(output/name) for name in (genotype,assignment,'evidence-audit.json')})
        save()
        (output/'COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')
    except BaseException as error:
        record.update(status='failed',finished=time.time(),error=repr(error));save()
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline','reads','tool','reference','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--threads',type=int,default=4)
    a = p.parse_args()
    run(a.baseline,a.reads,a.tool,a.reference,a.output,a.threads)
