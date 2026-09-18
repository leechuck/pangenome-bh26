"""Execute one stage of the fixed, balanced graph development expansion."""
import argparse
import json
from pathlib import Path
from evidence_io import sha
from run_development_mapping import ROOT,GENES,graph_path
from map_personalized import run as mapping
from run_development_evidence import run as evidence
from kmer_coverage import measure
from library_calibration import run as calibrate
from run_native_evidence import run as native
from join_evidence import build as join
from graph_recruitment import run as recruit

CONDITIONS=(('hprc','sampled'),('hprc','unsampled'),('hprc_asian','sampled'),('hprc_asian','unsampled'))
STAGES=('coverage','native','bootstrap','bootstrap_evidence','calibration','mapping','evidence','join','recruit')


def task_coordinates(stage,index,count):
    if stage not in STAGES:raise ValueError('Unknown stage')
    per_donor=32 if stage in ('mapping','evidence') else 8 if stage in ('bootstrap','bootstrap_evidence') else 4 if stage in ('join','recruit') else 1
    if not 0<=index<count*per_donor:raise ValueError('Task index outside fixed cohort')
    donor_index,offset=divmod(index,per_donor)
    condition=offset//8 if per_donor==32 else offset if per_donor==4 else 3
    gene_index=offset%8 if per_donor in (8,32) else None
    return donor_index,condition,gene_index


def run(stage,index,config):
    plan=json.loads(config.read_text())
    cohort=Path('/home/leechuck/hla/dogohla-series20260918/cohort.tsv')
    if sha(cohort)!=plan['cohort_sha256']:raise ValueError('Development cohort changed')
    i,condition,gene_index=task_coordinates(stage,index,len(plan['donors']))
    donor=plan['donors'][i]['donor']
    if not donor or Path(donor).name!=donor or donor in ('.','..'):raise ValueError('Unsafe donor identifier')
    panel,selection=CONDITIONS[condition]
    reads=Path('/home/leechuck/hla/hla-typer/reads')/donor
    coverage=ROOT/'coverage/development-genome-screen-v1'/(donor+'.json')
    development=ROOT/'development'
    calibration=development/'library-calibration-v1'/donor
    if stage=='coverage':
        measure(ROOT/'coverage/markers-genome-v1.json',[reads/'r1.fq.gz',reads/'r2.fq.gz'],coverage)
    elif stage=='native':
        native(development/'linear-v1/ipd_genome'/donor,reads,ROOT/'tools/t1k-evidence-v1',
               ROOT/'t1k-references/genome-v2/reference.fa',development/'native-genome-evidence-v1'/donor,4)
    elif stage in ('bootstrap','mapping'):
        gene=GENES[gene_index];version='v2' if stage=='bootstrap' else 'v3'
        mapping(graph_path(panel,gene),coverage,gene,[reads/'r1.fq.gz',reads/'r2.fq.gz'],
                development/('mapping-'+version)/donor/panel/selection/gene,4,
                disable_selection=selection=='unsampled',
                calibration_file=None if stage=='bootstrap' else calibration/'COMPLETE.json')
    elif stage in ('bootstrap_evidence','evidence'):
        evidence(panel,selection,gene_index,'v2' if stage=='bootstrap_evidence' else 'v3',donor)
    elif stage=='calibration':
        calibrate(development/'fragments-v2'/donor/'hprc_asian/unsampled',
                  development/'mapping-v2'/donor/'hprc_asian/unsampled',calibration,8)
    elif stage=='join':
        source=development/'native-genome-evidence-v1'/donor
        fragments=development/'fragments-v3'/donor/panel/selection
        mappings=development/'mapping-v3'/donor/panel/selection
        join(source,source/(donor+'_assign.tsv'),[fragments/g/'fragments' for g in GENES],
             [fragments/g/'support' for g in GENES],[mappings/g for g in GENES],
             development/'join-v3'/donor/panel/selection)
    elif stage=='recruit':
        recruit(development/'linear-v1/ipd_genome'/donor,reads,ROOT/'t1k-references/genome-v2',
                ROOT/'tools/t1k-evidence-v1',development/'graph-recruitment-v1'/donor/panel/selection,
                'graph',development/'join-v3'/donor/panel/selection,4)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=STAGES,required=True)
    p.add_argument('--index',type=int,required=True)
    p.add_argument('--config',type=Path,required=True)
    a=p.parse_args();run(a.stage,a.index,a.config)
