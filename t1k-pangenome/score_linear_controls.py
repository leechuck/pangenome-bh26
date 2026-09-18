"""Score development controls with the existing fixed, ambiguity-aware endpoint."""
import argparse
import csv
import json
from pathlib import Path
import sys
from build_graph import sha
from decode_t1k import decode

HERE = Path(__file__).resolve().parent
SERIES = HERE.parent/'hla-spechla-pg/series20260918'
PANELS = ('hprc','hprc_asian','ipd_genome','genome_hprc','genome_hprc_asian')
sys.path.insert(0,str(SERIES.parent))
from score_comparators import genomic_truth
from score import GENES, REPO, Nomenclature, compare, pred_pair, truth_a, truth_slots
from experiment import completed
from score_validation import write_tsv


def verify_cohort_metadata(planned, execution):
    """Allow relocated CRAM URLs only; every cohort-selection field must agree."""
    if len(planned) != len(execution) or not planned:
        raise ValueError('Execution cohort size changed')
    if len({row['donor'] for row in execution}) != len(execution):
        raise ValueError('Repeated execution donor')
    for left,right in zip(planned,execution):
        if set(left) != set(right) or any(left[key] != right[key] for key in left if key != 'cram'):
            raise ValueError('Execution cohort metadata/order changed')


def prediction(calls, gene, fields):
    slots = calls.get(gene,{}).get('resolutions',{}).get(str(fields))
    if slots is None:
        return [None,None]
    if len(slots) != 2:
        raise ValueError('Expected exactly two allele slots')
    if fields == 2:
        names = [None if slot['unresolved'] or not slot['alleles'] else
                 ','.join(slot['alleles']) for slot in slots]
        return pred_pair(names,'two_field',None,gene)
    return [None if slot['unresolved'] or not slot['alleles'] else
            [frozenset([allele]) for allele in slot['alleles']] for slot in slots]


def main(snapshot, output):
    cohort = SERIES/'cohort.tsv'
    donors = list(csv.DictReader(cohort.open(),delimiter='\t'))
    execution_cohort = HERE/'execution-cohort.tsv'
    execution = list(csv.DictReader(execution_cohort.open(),delimiter='\t'))
    verify_cohort_metadata(donors,execution)
    catalogue = REPO/'hla-analysis/results/sequence_catalogue.tsv'
    truth = truth_a(catalogue,{d['donor'] for d in donors})
    nomenclature = Nomenclature()
    rows, states, hashes = [], [], {}
    for donor_row in donors:
        donor = donor_row['donor']
        baseline = SERIES/'snapshot/t1k4'/donor
        baseline_manifest = json.loads((baseline/'manifest.json').read_text())
        config = baseline_manifest['configuration']
        if baseline_manifest['status'] != 'complete' or not completed(baseline,config):
            raise ValueError('Existing baseline is not verified complete: '+donor)
        table = baseline/(donor+'_genotype.tsv')
        hashes[str(table)] = sha(table)
        baseline_calls = decode(table.read_text(),{})
        methods = {'T1K':('complete',baseline_calls)}
        for panel in PANELS:
            folder = snapshot/panel/donor
            state, calls = 'not_started', {}
            if (folder/'manifest.json').exists():
                record = json.loads((folder/'manifest.json').read_text())
                state = record['status']
                if record['donor'] != donor or record['cohort_sha256'] != sha(execution_cohort):
                    raise ValueError('Control cohort changed')
                if record['reads_sha256'] != config['reads']:
                    raise ValueError('Control and baseline used different reads: '+donor)
                if state == 'complete':
                    if json.loads((folder/'COMPLETE.json').read_text()) != record:
                        raise ValueError('Completion/manifest mismatch')
                    for name,digest in record['output_sha256'].items():
                        if sha(folder/name) != digest:
                            raise ValueError('Changed completed output')
                        hashes[str(folder/name)] = digest
                    calls = json.loads((folder/'calls.json').read_text())
            methods[panel] = state,calls
        for method,(state,calls) in methods.items():
            states.append(dict(donor=donor,method=method,state=state))
            for gene in GENES:
                records = truth.get((donor,gene))
                for fields in (2,4):
                    slots = (genomic_truth(records) if fields==4 else
                             truth_slots(records,'two_field',nomenclature,gene)) if records else None
                    called,correct,matches = compare(prediction(calls,gene,fields),slots) if slots is not None else (0,0,0)
                    rows.append(dict(donor=donor,stratum=donor_row['stratum'],family=donor_row['family'],
                                     method=method,gene=gene,fields=fields,state=state,eligible=int(slots is not None),
                                     called=called,correct=correct,allele_matches=matches))
    reference = {(r['donor'],r['gene'],r['fields']):r for r in rows if r['method']=='T1K'}
    summary = []
    for method in ('T1K',*PANELS):
        for stratum in ('ALL',*sorted({d['stratum'] for d in donors})):
            donor_ids = {d['donor'] for d in donors if stratum=='ALL' or d['stratum']==stratum}
            for fields in (2,4):
                subset = [r for r in rows if r['method']==method and r['fields']==fields and r['donor'] in donor_ids]
                wins = losses = 0
                for row in subset:
                    base = reference[row['donor'],row['gene'],fields]
                    wins += row['correct'] > base['correct']
                    losses += row['correct'] < base['correct']
                summary.append(dict(method=method,stratum=stratum,fields=fields,
                    planned_donors=len(donor_ids),completed_donors=sum(s['state']=='complete' for s in states if s['method']==method and s['donor'] in donor_ids),
                    **{k:sum(r[k] for r in subset) for k in ('eligible','called','correct','allele_matches')},
                    gains_vs_T1K=wins,losses_vs_T1K=losses))
    output.mkdir(parents=True,exist_ok=True)
    for name,data in (('gene_scores',rows),('summary',summary),('status',states)):
        write_tsv(output/(name+'.tsv'),data)
    (output/'provenance.json').write_text(json.dumps(dict(
        scope='Development only; no independent improvement claim',cohort_sha256=sha(cohort),
        execution_cohort_sha256=sha(execution_cohort),
        catalogue_sha256=sha(catalogue),scorer_sha256=sha(Path(__file__)),input_sha256=hashes,
        scoring='All predicted alternatives must agree with eligible truth; unresolved/failed calls get zero credit. Pending results retain denominators; do not interpret incomplete-method contrasts.'),indent=2)+'\n')
    print(json.dumps([r for r in summary if r['stratum']=='ALL'],indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot',type=Path,default=HERE/'work/linear-snapshot')
    p.add_argument('--output',type=Path,default=HERE/'development/linear-analysis')
    a = p.parse_args()
    main(a.snapshot,a.output)
