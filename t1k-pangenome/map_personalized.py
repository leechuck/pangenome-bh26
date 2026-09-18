"""Map unbinned paired reads to a coverage-informed observed HLA graph.

This emits alignment evidence only. It is not a standalone HLA typer; native
all-IPD evidence and joint competing-locus inference remain necessary.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
from build_graph import sha


def sampling_decision(coverage, gene, disabled=False):
    if disabled:
        return None,'selection_disabled_control'
    record = coverage['genes'][gene]
    if not record['usable']:
        return None,record['fallback_reason']
    depth = record['personalization_coverage']
    if not isinstance(depth,int) or depth < 1:
        raise ValueError('Invalid measured personalization coverage')
    return depth,None


def fragment_parameters(calibration, hashes):
    if calibration['status'] != 'complete' or not calibration['estimate']['usable']:
        raise ValueError('Incomplete library calibration')
    if calibration['reads_sha256'] != hashes:
        raise ValueError('Library calibration used different reads')
    fit = calibration['estimate']
    values = [fit['fragment_mean'], fit['fragment_stdev']]
    if any(isinstance(v, bool) or not isinstance(v, (float, int)) or
           not math.isfinite(v) or v <= 0 for v in values):
        raise ValueError('Invalid insert distribution')
    return ['--fragment-mean', str(values[0]), '--fragment-stdev', str(values[1])]


def run(graph, coverage_file, gene, reads, output, threads, disable_selection=False,
        calibration_file=None):
    if not 1 <= threads <= int(os.environ.get('SLURM_CPUS_PER_TASK','0')):
        raise ValueError('Requires sufficient Slurm CPU allocation')
    if output.exists():
        raise FileExistsError(output)
    graph,coverage_file = graph.resolve(),coverage_file.resolve()
    reads = [p.resolve() for p in reads]
    graph_record = json.loads((graph/'COMPLETE.json').read_text())
    coverage = json.loads(coverage_file.read_text())
    if graph_record['status']!='complete' or coverage['status']!='complete':
        raise ValueError('Incomplete inputs')
    if graph_record['gene'] != gene:
        raise ValueError('Graph locus mismatch')
    for name,digest in graph_record['output_sha256'].items():
        if sha(graph/name)!=digest:
            raise ValueError('Graph product changed: '+name)
    hashes = {str(p):sha(p) for p in reads}
    if hashes != coverage['reads_sha256']:
        raise ValueError('Coverage was estimated from different reads')
    insert_args = []
    if calibration_file is not None:
        calibration = json.loads(calibration_file.read_text())
        insert_args = fragment_parameters(calibration, hashes)
    depth,fallback = sampling_decision(coverage,gene,disable_selection)
    vg = Path('/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin/vg')
    kmc = Path('/home/leechuck/hla/t1k-pangenome/tools/kmc-3.2.4/kmc')
    output.mkdir(parents=True)
    output = output.resolve()
    (output/'tmp').mkdir()
    env = dict(os.environ,OMP_NUM_THREADS=str(threads),TMPDIR=str(output/'tmp'))
    record = dict(status='running',started=time.time(),gene=gene,threads=threads,
                  graph_manifest_sha256=sha(graph/'COMPLETE.json'),coverage_sha256=sha(coverage_file),
                  reads_sha256=hashes,driver_sha256=sha(Path(__file__)),vg_sha256=sha(vg),
                  measured_kmer_coverage=depth,selection_disabled=disable_selection,
                  sampling_fallback=fallback,commands=[],
                  scope='Per-locus alignment evidence only; retain all-IPD and competing-locus evidence')
    if calibration_file is not None:
        record['library_calibration_sha256'] = sha(calibration_file)
        record['fragment_parameters'] = insert_args
    def save():
        (output/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    def execute(command, target=None):
        record['commands'].append(command)
        save()
        with (output/'run.log').open('ab') as log:
            if target is None:
                subprocess.run(command,cwd=output,env=env,stdout=log,stderr=log,check=True)
            else:
                with target.open('wb') as stream:
                    subprocess.run(command,cwd=output,env=env,stdout=stream,stderr=log,check=True)
    save()
    try:
        for name in ('graph.gbz','graph.hapl'):
            shutil.copy2(graph/name,output/name)
        mapping_graph = 'graph.gbz'
        if depth is not None:
            record['kmc_sha256'] = sha(kmc)
            (output/'reads.list').write_text('\n'.join(map(str,reads))+'\n')
            execute([str(kmc),'-k29','-m4','-okff','-t'+str(threads),'-hp',
                     '@reads.list','counts','tmp'])
            execute([str(vg),'haplotypes','-v','2','-t',str(threads),'--coverage',str(depth),
                     '--num-haplotypes','8','--include-reference','-i','graph.hapl',
                     '-k','counts.kff','-g','sampled.gbz','graph.gbz'])
            mapping_graph = 'sampled.gbz'
        command = [str(vg),'giraffe','-Z',mapping_graph,'-f',str(reads[0]),'-f',str(reads[1]),'-t',str(threads)]
        command += insert_args
        if mapping_graph=='graph.gbz':
            # The haplotype-preparation index intentionally omits distances;
            # giraffe/minimizer construction requires a full distance index.
            execute([str(vg),'index','-j','mapping.dist','-t',str(threads),
                     '-P',graph_record['indexing_backbone'],'graph.gbz'])
            command += ['-d','mapping.dist']
        execute(command,output/'mapping.gam')
        execute([str(vg),'view','-aj',str(output/'mapping.gam')],output/'alignments.jsonl')
        count = mapped = 0
        with (output/'alignments.jsonl').open() as stream:
            for line in stream:
                alignment = json.loads(line)
                count += 1
                mapped += bool(alignment.get('path',{}).get('mapping'))
        if count != coverage['reads']:
            raise ValueError('Alignment output did not retain every input read')
        record.update(status='complete',finished=time.time(),alignments=count,mapped=mapped,
                      output_sha256={name:sha(output/name) for name in
                                     ('mapping.gam','alignments.jsonl',mapping_graph)})
        save()
        (output/'COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')
    except BaseException as error:
        record.update(status='failed',finished=time.time(),error=repr(error))
        save()
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('graph','coverage','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--gene',required=True)
    p.add_argument('--reads',type=Path,nargs=2,required=True)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--disable-selection',action='store_true')
    p.add_argument('--calibration',type=Path)
    a = p.parse_args()
    run(a.graph,a.coverage,a.gene,a.reads,a.output,a.threads,a.disable_selection,a.calibration)
