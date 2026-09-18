"""Paired family-cluster endpoint analysis of a complete, fixed score table.

This module takes already-scored rows. It never loads truth, chooses a candidate,
or turns development evidence into an independent validation result.
"""
import collections
import random

GENES=('A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1')
TERMINAL={'complete','failed'}


def interval(groups,repetitions=10000,seed=20260918):
    if repetitions<100:raise ValueError('Too few bootstrap replicates')
    values=[groups[key] for key in sorted(groups)]
    denominator=sum(n for _,n in values)
    if not denominator:return None
    rng=random.Random(seed);samples=[]
    for _ in range(repetitions):
        chosen=rng.choices(values,k=len(values));n=sum(x[1] for x in chosen)
        if n:samples.append(sum(x[0] for x in chosen)/n)
    if not samples:raise ValueError('No eligible bootstrap samples')
    samples.sort()
    return dict(difference=sum(d for d,_ in values)/denominator,
                lower_95=samples[int(.025*len(samples))],
                upper_95=samples[min(len(samples)-1,int(.975*len(samples)))],
                eligible=denominator,families=len(values),replicates=repetitions,seed=seed)


def evaluate(rows,cohort,candidate,baseline='T1K',repetitions=10000,seed=20260918):
    if candidate==baseline:raise ValueError('Candidate equals baseline')
    donors={row['donor']:row for row in cohort}
    if not donors or len(donors)!=len(cohort):raise ValueError('Empty or duplicate cohort')
    if any(not r.get('family') or not r.get('stratum') for r in cohort):
        raise ValueError('Every donor requires family and stratum metadata')
    indexed={}
    for row in rows:
        if row['method'] not in (candidate,baseline):continue
        if row['donor'] not in donors:raise ValueError('Score donor outside fixed cohort')
        key=(row['method'],row['donor'],row['gene'],int(row['fields']))
        if key in indexed:raise ValueError('Duplicate score row')
        if row['state'] not in TERMINAL:raise ValueError('Nonterminal result; wait for complete comparison')
        value={k:int(row[k]) for k in ('eligible','called','correct')}
        if any(x not in (0,1) for x in value.values()):raise ValueError('Scores must be binary per genotype')
        if value['correct']>value['called'] or value['correct']>value['eligible']:
            raise ValueError('Impossible genotype score')
        if row['state']=='failed' and (value['called'] or value['correct']):
            raise ValueError('Failed result cannot receive credit')
        metadata=donors[row['donor']]
        if row.get('family',metadata['family'])!=metadata['family'] or row.get('stratum',metadata['stratum'])!=metadata['stratum']:
            raise ValueError('Score metadata differs from fixed cohort')
        indexed[key]=dict(row,**value)
    expected={(m,d,g,f) for m in (candidate,baseline) for d in donors for g in GENES for f in (2,4)}
    if set(indexed)!=expected:raise ValueError('Missing or unexpected genotype rows')
    contrasts={}
    for stratum in ('ALL',*sorted({r['stratum'] for r in cohort})):
        for fields in (2,4):
            groups=collections.defaultdict(lambda:[0,0]);wins=losses=correct_candidate=correct_baseline=0
            failures={candidate:set(),baseline:set()}
            for donor,metadata in donors.items():
                if stratum!='ALL' and metadata['stratum']!=stratum:continue
                for gene in GENES:
                    a=indexed[candidate,donor,gene,fields];b=indexed[baseline,donor,gene,fields]
                    if a['eligible']!=b['eligible']:raise ValueError('Unequal truth eligibility')
                    for method,r in ((candidate,a),(baseline,b)):
                        if r['state']=='failed':failures[method].add(donor)
                    # Retain all families, including zero-eligible families, in resampling.
                    group=groups[metadata['family']]
                    if not a['eligible']:continue
                    delta=a['correct']-b['correct'];group[0]+=delta;group[1]+=1
                    wins+=delta>0;losses+=delta<0
                    correct_candidate+=a['correct'];correct_baseline+=b['correct']
            contrasts[stratum+'/'+str(fields)]=dict(interval=interval(groups,repetitions,seed),
                gains=wins,losses=losses,candidate_correct=correct_candidate,baseline_correct=correct_baseline,
                failed_donors={m:len(s) for m,s in failures.items()})
    four=contrasts['ALL/4']['interval'];two=contrasts['ALL/2']['interval']
    checks=dict(four_field_lower_bound_positive=four is not None and four['lower_95']>0,
                two_field_point_not_decreased=two is not None and two['difference']>=0,
                two_field_excludes_loss_over_one_point=two is not None and two['lower_95']>-.01)
    return dict(candidate=candidate,baseline=baseline,contrasts=contrasts,protocol_checks=checks,
                protocol_pass=all(checks.values()),
                interpretation='Endpoint calculations only. Independent improvement also requires a pre-outcome candidate freeze, untouched cohort, verified run provenance and complete failure accounting.')
