#!/usr/bin/env python3
"""Independent two-field concordance; resample known families together."""
import collections,csv,json,random,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from score import Nomenclature,compare,pred_pair,read_result
from score_comparators import t1k_calls
from score_validation import write_tsv
from experiment import completed,sha
from truth import gourraud,GOURRAUD_GENES
import score_series as matched

def family_interval(groups,reps=10000):
    """Paired accuracy difference: cluster-sampled numerator / denominator."""
    values=list(groups.values())
    if not values or not sum(n for _,n in values):return None
    rng=random.Random(20260918);estimates=[]
    for _ in range(reps):
        sample=rng.choices(values,k=len(values));den=sum(n for _,n in sample)
        if den:estimates.append(sum(g for g,_ in sample)/den)
    estimates.sort()
    return dict(difference=sum(g for g,_ in values)/sum(n for _,n in values),lower_95=estimates[int(.025*len(estimates))],upper_95=estimates[min(len(estimates)-1,int(.975*len(estimates)))],families=len(values),replicates=reps)

def main():
    matched.ROOT=HERE/'snapshot/gourraud';out=HERE/'analysis/gourraud';out.mkdir(parents=True,exist_ok=True)
    donors=list(csv.DictReader((HERE/'cohort.gourraud.tsv').open(),delimiter='\t'))
    methods=[m for m in matched.METHODS if '-long-' not in m]
    nom=Nomenclature();truth,reasons=gourraud('two_field',nom)
    rows=[];states=[];hashes={}
    for d in donors:
        for method in methods:
            p=matched.locate(d['donor'],method);state='not_started';names={}
            if (p/'manifest.json').exists():
                try:
                    manifest=json.loads((p/'manifest.json').read_text());state=manifest['status']
                    if completed(p,manifest['configuration']):
                        state='complete';f=p/(d['donor']+'_genotype.tsv' if method=='T1K-four-field' else 'hla.result.txt')
                        names=t1k_calls(f) if method=='T1K-four-field' else read_result(f)
                        hashes[str(f.relative_to(matched.ROOT))]=sha(f)
                    elif state=='complete':state='invalid'
                except Exception:state='invalid'
            states.append(dict(donor=d['donor'],method=method,state=state))
            for gene in GOURRAUD_GENES:
                slots=truth.get((d['donor'],gene));pair=names.get(gene,[None,None])
                called,correct,matches=compare(pred_pair(pair,'two_field',nom,gene),slots) if slots else (0,0,0)
                rows.append(dict(donor=d['donor'],family=d['family'],stratum=d['stratum'],method=method,gene=gene,state=state,eligible=int(slots is not None),called=called,correct=correct,allele_matches=matches,allele1=pair[0] or 'NO_CALL',allele2=pair[1] or 'NO_CALL',exclusion_reason=reasons.get((d['donor'],gene),'')))
    summary=[]
    for method in methods:
        for stratum in ('ALL','EAS','EUR','AFR','AMR'):
            ds={d['donor'] for d in donors if stratum=='ALL' or d['stratum']==stratum}
            rs=[r for r in rows if r['method']==method and r['donor'] in ds]
            summary.append(dict(method=method,stratum=stratum,planned=len(ds),complete=sum(s['state']=='complete' for s in states if s['method']==method and s['donor'] in ds),**{k:sum(r[k] for r in rs) for k in ('eligible','called','correct','allele_matches')}))
    write_tsv(out/'gene_scores.tsv',rows);write_tsv(out/'summary.tsv',summary);write_tsv(out/'status.tsv',states)
    # Until both methods complete the whole planned cohort, withhold intervals:
    # scheduling speed must not determine which families enter a comparison.
    finished={s['method'] for s in summary if s['stratum']=='ALL' and s['complete']==len(donors)}
    index={(r['method'],r['donor'],r['gene']):r for r in rows};intervals={}
    target='full:DogoHLA'
    if target in finished:
        for baseline in sorted(finished-{target}):
            groups=collections.defaultdict(lambda:[0,0])
            for r in rows:
                if r['method']!=target or not r['eligible']:continue
                b=index[(baseline,r['donor'],r['gene'])]
                groups[r['family']][0]+=r['correct']-b['correct'];groups[r['family']][1]+=1
            intervals[baseline]=family_interval(groups)
    (out/'bootstrap.json').write_text(json.dumps(intervals,indent=2)+'\n')
    (out/'provenance.json').write_text(json.dumps(dict(cohort_sha256=sha(HERE/'cohort.gourraud.tsv'),output_sha256=hashes,truth='Gourraud 2014 experimental exon typing, ambiguity-aware two-field; no four-field truth',bootstrap='Paired family clusters; fixed seed; 10000 replicates; exploratory, no multiplicity adjustment'),indent=2)+'\n')
    report=['# Gourraud experimental-truth benchmark','','946 donors; 660 family groups; graph-donor and known-family overlap zero. Four-field calls are retained, but this truth supports two-field scoring only.','','| Method | Complete | Two-field genotype pairs |','|---|---:|---:|']
    for s in summary:
        if s['stratum']=='ALL':report.append(f"| {s['method']} | {s['complete']}/946 | "+(f"{s['correct']}/{s['eligible']}" if s['complete']==946 else 'pending')+' |')
    report+=['','Pending/failed calls retain planned eligible denominators. Summary TSV includes ancestry strata. Intervals appear only after both compared methods complete the entire cohort.','']
    (out/'REPORT.md').write_text('\n'.join(report));print('\n'.join(report))

if __name__=='__main__':main()
