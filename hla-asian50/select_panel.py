#!/usr/bin/env python3
"""Deterministic feasibility cohort, selected without read prediction outcomes."""
import csv, collections, hashlib, json
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SEED = 'asian50-20260917-v1'
def rows(path, delimiter='\t'):
    return list(csv.DictReader(open(path), delimiter=delimiter))
def rank(d):
    return hashlib.sha256((SEED+':'+d).encode()).hexdigest()
def write(name, records):
    with open(ROOT/'source'/name, 'w') as f:
        w=csv.DictWriter(f, fieldnames=list(records[0]), delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(records)
def main():
    ped={r['SampleID']:r for r in rows(REPO/'hla-pilot/source/1000G.ped',' ')}
    truth={r['id']:r for r in rows(REPO/'hla-pilot/source/1000G_HLA.txt',' ')}
    cat=rows(REPO/'hla-structural/source/catalogue.tsv')
    grouped=collections.defaultdict(list)
    for r in cat:
        if r['locus']=='RCCX' and r['eligible']=='1':grouped[r['donor']].append(r)
    eligible=[d for d,rs in grouped.items() if len(rs)==2 and ped.get(d,{}).get('Superpopulation')=='EAS']
    # Retain every donor with independent classical labels, then population-balanced fill.
    chosen=sorted([d for d in eligible if d in truth],key=rank)
    families={grouped[d][0]['family'] for d in chosen}
    assert len(families)==len(chosen)
    while len(chosen)<50:
        counts=collections.Counter(ped[d]['Population'] for d in chosen)
        remaining=[d for d in eligible if d not in chosen and grouped[d][0]['family'] not in families]
        assert remaining, 'Fewer than 50 eligible distinct families'
        d=min(remaining,key=lambda d:(counts[ped[d]['Population']],rank(d)))
        chosen.append(d);families.add(grouped[d][0]['family'])
    # Distribute population groups across five folds, balancing total fold size.
    fold={};sizes=collections.Counter();popsize=collections.Counter()
    for pop in sorted({ped[d]['Population'] for d in chosen}):
        for d in sorted([d for d in chosen if ped[d]['Population']==pop],key=rank):
            k=min(range(5),key=lambda k:(popsize[pop,k],sizes[k],k))
            fold[d]=k;sizes[k]+=1;popsize[pop,k]+=1
    old=set((REPO/'hla-structural/source/evaluation_donors.txt').read_text().split())
    old.update(r['donor'] for r in rows(REPO/'hla-targeted/source/validation_donors.tsv'))
    records=[dict(donor=d,population=ped[d]['Population'],superpopulation='EAS',family=grouped[d][0]['family'],fold=fold[d],experimental_HLA_labels=int(d in truth),paired_RCCX_annotations=1,previous_structural_read_evaluation=int(d in old),selection_sha256=rank(d)) for d in sorted(chosen,key=lambda d:(fold[d],ped[d]['Population'],d))]
    write('donors.tsv',records)
    (ROOT/'source/donors.txt').write_text('\n'.join(r['donor'] for r in records)+'\n')
    hla=[]
    for r in records:
        if r['donor'] not in truth:continue
        for g in ['A','B','C','DRB1','DQB1']:
            t=truth[r['donor']]
            hla.append(dict(donor=r['donor'],gene='HLA-'+g,allele1=t[g],allele2=t[g+'.1'],source='Gourraud2014_experimental_typing',resolution='twofield_with_reported_ambiguities'))
    write('experimental_HLA_truth.tsv',hla)
    manifest=rows(REPO/'hla-audit/2026-09-16/panel_manifest.tsv')
    family_by_donor={r['donor']:r['family'] for r in cat}
    def fam(d):return family_by_donor.get(d,ped.get(d,{}).get('FamilyID',d))
    exclusions=[]
    for k in range(5):
        excluded={r['family'] for r in records if r['fold']==k}
        for p in manifest:
            if fam(p['donor_id']) in excluded:
                exclusions.append(dict(fold=k,sample=p['sample'],hap_id=p['hap_id'],donor=p['donor_id'],family=fam(p['donor_id'])))
    write('excluded_paths.tsv',exclusions)
    # Catalogue-level oracle only, not short-read accuracy or confirmed biological truth.
    oracle=[]
    for r in records:
        excluded={x['family'] for x in records if x['fold']==r['fold']}
        for locus in ['DRB','RCCX']:
            training=[x for x in cat if x['locus']==locus and x['eligible']=='1' and x['family'] not in excluded]
            for p in [x for x in cat if x['donor']==r['donor'] and x['locus']==locus and x['eligible']=='1']:
                for arm in ['full_panel','HPRC_only']:
                    armrows=[x for x in training if arm=='full_panel' or x['cohort'].startswith('HPRC')]
                    oracle.append(dict(donor=r['donor'],fold=r['fold'],locus=locus,hap_id=p['hap_id'],arm=arm,training_haplotypes=len(armrows),signature=p['structural_signature'],signature_represented=int(any(x['structural_signature']==p['structural_signature'] for x in armrows))))
    write('structural_representation_oracle.tsv',oracle)
    summary=dict(seed=SEED,eligible_EAS_donors=len(eligible),selected=len(records),populations=dict(collections.Counter(r['population'] for r in records)),fold_sizes=dict(sizes),experimental_HLA_donors=sum(r['experimental_HLA_labels'] for r in records),experimental_HLA_donor_gene_rows=len(hla),previously_evaluated_donors=sum(r['previous_structural_read_evaluation'] for r in records),oracle={})
    for locus in ['DRB','RCCX']:
        for arm in ['full_panel','HPRC_only']:
            rr=[r for r in oracle if r['locus']==locus and r['arm']==arm]
            summary['oracle'][locus+'_'+arm]=dict(represented=sum(r['signature_represented'] for r in rr),total=len(rr))
    assert len(records)==50 and len({r['family'] for r in records})==50
    assert all(sum(r['fold']==k for r in records)==10 for k in range(5))
    (ROOT/'source/cohort_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
