#!/usr/bin/env python3
"""Native SpecHLA -v True control, with documented execution compatibility fixes."""
import argparse
import csv
from pathlib import Path
import shutil
import re
import time
from experiment import ENV, command, environment, failed, finish, json_write, sha, start, validate_outputs


def patch(source):
    old='port=$(date +%N|cut -c5-9)'
    new="port=$($python_bin -c 'import socket; s=socket.socket(); s.bind((\"127.0.0.1\",0)); print(s.getsockname()[1])')"
    if source.count(old)!=1: raise ValueError('Unexpected upstream port command')
    source=source.replace(old,new)
    old='hla_ref=$db/ref/HLA_$hla.fa'
    if source.count(old)!=1: raise ValueError('Unexpected upstream phase-reference command')
    return source.replace(old,'hla_ref=$db/HLA/HLA_$hla/HLA_$hla.fa')


def patch_insertion_phaser(source):
    # Upstream SpecHap now dereferences --protocols unconditionally; legacy -N
    # is still accepted by its parser but no longer initializes that argument.
    old = 'SpecHap --ncs --window_size 15000 -N --vcf'
    new = 'SpecHap --ncs --window_size 15000 --protocols nanopore --weights 1 --vcf'
    if source.count(old) != 1:
        raise ValueError('Unexpected upstream insertion-phasing command')
    return source.replace(old, new)


def patch_empty_assembly(source):
    # An existing but empty list means no contig was selected. Upstream's -f
    # check incorrectly enters BWA and can reuse an index from a previous region.
    for name in ('extract.fa', 'id.list'):
        old = 'if [ ! -f "$outdir/' + name + '" ]'
        new = 'if [ ! -s "$outdir/' + name + '" ]'
        if source.count(old)==0 and source.count(new)==1:continue
        if source.count(old) != 1 or source.count(new):
            raise ValueError('Unexpected local-assembly guard: ' + name)
        source = source.replace(old, new)
    return source


def run(root,donor,fold,threads):
    out=root/'runs'/donor/'native-long'
    reads=root/'reads'/donor
    if not (reads/'COMPLETE').exists(): raise ValueError('Incomplete recruitment')
    original=ENV/'share/spechla/script'
    rewritten=patch((original/'whole/SpecHLA.sh').read_text())
    import hashlib
    config=dict(donor=donor,fold=fold,arm='native-long',threads=threads,
                inputs={str(reads/f'r{i}.fq.gz'):sha(reads/f'r{i}.fq.gz') for i in (1,2)},
                source_sha256=sha(original/'whole/SpecHLA.sh'),driver_sha256=sha(__file__),
                patched_pipeline_sha256=hashlib.sha256(rewritten.encode()).hexdigest(),
                compatibility_fixes=['OS-selected valid TCP port at ScanIndel startup','Correct installed path to the same per-gene reference FASTA','Honor allocated threads in ScanIndel BWA','Translate deprecated SpecHap -N to explicit nanopore protocol with unit weight','Use upstream no-contig fallback for empty read/contig lists; prevent stale BWA index reuse'])
    if not start(out,config): return
    t0=time.monotonic()
    try:
        script=out/'vendor/script'; shutil.copytree(original,script)
        (script/'whole/SpecHLA.sh').write_text(rewritten)
        phaser=script/'phase_variants.py'
        phaser.write_text(patch_insertion_phaser(phaser.read_text()))
        assembly=script/'run.assembly.realign.sh'
        assembly.write_text(patch_empty_assembly(assembly.read_text()))
        scan=script/'ScanIndel/ScanIndel.py'
        scan_source=scan.read_text()
        if scan_source.count('bwa mem -M -t8') != 1: raise ValueError('Unexpected ScanIndel thread option')
        scan.write_text(scan_source.replace('bwa mem -M -t8',f'bwa mem -M -t{threads}'))
        env=environment(script); tmp=root/'runs'/donor/'native-long_work'
        command(['timeout','--kill-after=30s','1800','spechla','-n',donor,'-1',reads/'r1.fq.gz',
                 '-2',reads/'r2.fq.gz','-o',tmp,'-j',threads,'-u','0','-p','Unknown','-v','True'],out/'spechla.log',env)
        for p in (tmp/donor).iterdir(): shutil.move(str(p),out/p.name)
        for log in [out/'spechla.log',*out.glob('*.log'),*(out/'Scanindel').glob('*.log')]:
            content=log.read_text(errors='replace')
            if any(s in content for s in ('Traceback (most recent call last)','SyntaxError:','command not found','module not found','Segmentation fault','Execution failed for','ERROR - could not find','No such file or directory')):
                raise ValueError(f'Nested failure: {log}')
        segments=re.split(r'Analyzing sample: '+re.escape(donor)+r'\.',(out/'spechla.log').read_text())[1:]
        if len(segments)!=8 or {s.splitlines()[0].strip() for s in segments} != {'A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1'}:
            raise ValueError('Expected eight ScanIndel gene stages')
        for segment in segments:
            if not any(marker in segment for marker in ('ScanIndel running done:', 'We have not found softclip and unmapped reads.')):
                raise ValueError('Unverified ScanIndel stage: '+segment.splitlines()[0])
        outputs=validate_outputs(out,donor)
        (out/'timing.tsv').write_text(f'stage\tseconds\ntotal\t{time.monotonic()-t0:.3f}\n')
        finish(out,config,outputs)
    except BaseException as exc:
        failed(out,exc); raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True); p.add_argument('--cohort-index',type=int,required=True)
    p.add_argument('--threads',type=int,default=4)
    a=p.parse_args(); root=a.root.resolve()
    rows=list(csv.DictReader(open(root/'cohort.tsv'),delimiter='\t')); r=rows[a.cohort_index]
    run(root,r['donor'],int(r['fold']),a.threads)
