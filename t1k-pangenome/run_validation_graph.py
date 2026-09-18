"""Execute frozen graph validation stages without loading evaluation truth."""
import argparse
import csv
import json
from pathlib import Path
from evidence_io import sha
from run_development_mapping import ROOT,GENES,graph_path
from map_personalized import run as mapping
from run_development_evidence import run_paths as evidence
from kmer_coverage import measure
from library_calibration import run as calibrate
from run_native_evidence import run as native
from join_evidence import build as join
from run_graph_pair_refinement import run as refine

STAGES=('coverage','native','bootstrap','bootstrap_evidence','calibration','mapping','evidence','join','refine')


def coordinates(stage,index,count):
    if stage not in STAGES:raise ValueError('Unknown stage')
    width=16 if stage in ('mapping','evidence') else 8 if stage in ('bootstrap','bootstrap_evidence') else 2 if stage in ('join','refine') else 1
    if not 0<=index<count*width:raise ValueError('Index outside frozen cohort')
    donor,offset=divmod(index,width)
    panel=('hprc','hprc_asian')[offset//8 if width==16 else offset] if width in (16,2) else 'hprc_asian'
    gene=GENES[offset%8] if width in (8,16) else None
    return donor,panel,gene


def run(stage,index,assets,freeze_sha256):
    frozen=assets/'FROZEN_VALIDATION.json'
    if sha(frozen)!=freeze_sha256:raise ValueError('Validation freeze changed')
    plan=json.loads(frozen.read_text())
    for filename,digest in plan['code_sha256'].items():
        if sha(assets/filename)!=digest:raise ValueError('Frozen code changed: '+filename)
    for filename,digest in plan['metadata_sha256'].items():
        if sha(assets/filename)!=digest:raise ValueError('Frozen metadata changed: '+filename)
    for filename,digest in plan['asset_sha256'].items():
        if sha(ROOT/filename)!=digest:raise ValueError('Frozen reference asset changed: '+filename)
    cohort=list(csv.DictReader((assets/'cohort.tsv').open(),delimiter='\t'))
    if len(cohort)!=plan['donors']:raise ValueError('Frozen cohort count changed')
    i,panel,gene=coordinates(stage,index,len(cohort));donor=cohort[i]['donor']
    if not donor or Path(donor).name!=donor or donor in ('.','..'):raise ValueError('Unsafe donor')
    prep=json.loads((assets/'READ_PREPARATION.json').read_text())
    entries=[r for r in prep['donors'] if r['donor']==donor]
    if len(entries)!=1 or prep['cohort_sha256']!=plan['metadata_sha256']['cohort.tsv']:
        raise ValueError('Read preparation/cohort mismatch')
    reads=ROOT/'validation-reads'/donor
    if {n:sha(reads/n) for n in ('r1.fq.gz','r2.fq.gz')}!=entries[0]['reads_sha256']:
        raise ValueError('Reserved reads changed')
    out=ROOT/'validation-runs/graph-v1';coverage=out/'coverage'/(donor+'.json')
    calibration=out/'library-calibration-v1'/donor
    baseline=ROOT/'validation-runs/genomic-ipd-v1'/donor
    if stage=='coverage':
        measure(ROOT/'coverage/markers-genome-v1.json',[reads/'r1.fq.gz',reads/'r2.fq.gz'],coverage)
    elif stage=='native':
        native(baseline,reads,ROOT/'tools/t1k-evidence-v1',ROOT/'t1k-references/genome-v2/reference.fa',
               out/'native-genome-evidence-v1'/donor,4)
    elif stage in ('bootstrap','mapping'):
        version='v2' if stage=='bootstrap' else 'v3'
        mapping(graph_path(panel,gene),coverage,gene,[reads/'r1.fq.gz',reads/'r2.fq.gz'],
                out/('mapping-'+version)/donor/panel/'unsampled'/gene,4,disable_selection=True,
                calibration_file=None if stage=='bootstrap' else calibration/'COMPLETE.json')
    elif stage in ('bootstrap_evidence','evidence'):
        version='v2' if stage=='bootstrap_evidence' else 'v3'
        evidence(out/('mapping-'+version)/donor/panel/'unsampled'/gene,graph_path(panel,gene),
                 ROOT/'references/observed-v2'/panel/('HLA-'+gene+'.fa'),
                 out/('fragments-'+version)/donor/panel/'unsampled'/gene)
    elif stage=='calibration':
        calibrate(out/'fragments-v2'/donor/'hprc_asian/unsampled',
                  out/'mapping-v2'/donor/'hprc_asian/unsampled',calibration,8)
    elif stage=='join':
        source=out/'native-genome-evidence-v1'/donor
        fragments=out/'fragments-v3'/donor/panel/'unsampled'
        mappings=out/'mapping-v3'/donor/panel/'unsampled'
        join(source,source/(donor+'_assign.tsv'),[fragments/g/'fragments' for g in GENES],
             [fragments/g/'support' for g in GENES],[mappings/g for g in GENES],
             out/'join-v3'/donor/panel/'unsampled')
    elif stage=='refine':
        refine(baseline,out/'join-v3'/donor/panel/'unsampled',ROOT/'references/observed-v2'/panel,
               calibration/'COMPLETE.json',out/'calls'/donor/panel/'unsampled',
               internal_pairs=True,pair_specific=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=STAGES,required=True);p.add_argument('--index',type=int,required=True)
    p.add_argument('--assets',type=Path,required=True);p.add_argument('--freeze-sha256',required=True)
    a=p.parse_args();run(a.stage,a.index,a.assets,a.freeze_sha256)
