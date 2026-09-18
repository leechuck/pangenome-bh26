"""Fixed-cohort endpoint calculation; call only after complete provenance review.

This module consumes explicitly supplied predictions and truth. Importing it
does not load reserved predictions or truth. File collection/unblinding is separate.
"""
from score_linear_controls import GENES,prediction,compare
from hgsvc_truth import truth_slots as genomic_slots
from paired_evaluation import evaluate

def truth_slots(records, fields):
    slots = genomic_slots(records, fields)
    if slots is None or fields == 4:
        return slots
    # Existing two-field prediction parser returns numeric names without gene.
    return tuple(frozenset(name.split('*', 1)[1] for name in slot) for slot in slots)


METHODS=('T1K','ipd_genome','graph_hprc','graph_hprc_asian')


def score(cohort,results,truth,repetitions=10000,seed=20260918):
    donors={r['donor'] for r in cohort}
    if not donors or len(donors)!=len(cohort):raise ValueError('Empty or repeated cohort donor')
    expected={(d,m) for d in donors for m in METHODS}
    if set(results)!=expected:raise ValueError('Incomplete or unexpected donor/method grid')
    for result in results.values():
        if result['state'] not in ('complete','failed'):raise ValueError('Nonterminal validation result')
        if result['state']=='failed' and result['calls']:raise ValueError('Failed result must not contribute calls')
    rows=[]
    for donor in cohort:
        if not donor.get('family') or not donor.get('stratum'):raise ValueError('Missing cohort metadata')
        d=donor['donor']
        for gene in GENES:
            records=truth.get((d,gene))
            for fields in (2,4):
                slots=truth_slots(records,fields) if records else None
                for method in METHODS:
                    result=results[d,method]
                    called,correct,matches=compare(prediction(result['calls'],gene,fields),slots) if slots is not None else (0,0,0)
                    rows.append(dict(donor=d,family=donor['family'],stratum=donor['stratum'],gene=gene,
                        method=method,fields=fields,state=result['state'],eligible=int(slots is not None),
                        called=called,correct=correct,allele_matches=matches))
    contrasts={baseline:evaluate(rows,cohort,'graph_hprc_asian',baseline,repetitions,seed)
               for baseline in ('T1K','ipd_genome','graph_hprc')}
    # Only the original-T1K contrast has the prespecified primary hypothesis.
    for baseline in ('ipd_genome','graph_hprc'):
        contrasts[baseline]['interpretation']='Exploratory incremental contrast; primary protocol_pass does not apply to this comparator'
        contrasts[baseline].pop('protocol_pass')
        contrasts[baseline].pop('protocol_checks')
    summary=[]
    for method in METHODS:
        for fields in (2,4):
            subset=[r for r in rows if r['method']==method and r['fields']==fields]
            summary.append(dict(method=method,fields=fields,donors=len(donors),
                failed_donors=len({r['donor'] for r in subset if r['state']=='failed'}),
                ineligible_genotypes=sum(not r['eligible'] for r in subset),
                **{k:sum(r[k] for r in subset) for k in ('eligible','called','correct','allele_matches')}))
    return rows,dict(summary=summary,contrasts=contrasts,
        interpretation='Independent proof additionally requires verified frozen execution, exclusion audit and complete failure accounting')
