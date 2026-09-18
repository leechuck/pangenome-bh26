#!/usr/bin/env python3
"""Score the matched cohort with fixed truth denominators and explicit resolution."""
import collections,csv,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from score_comparators import genomic_truth,four_pair,t1k_calls
from score import GENES,REPO,Nomenclature,compare,pred_pair,read_result,truth_a,truth_slots
from experiment import completed,sha
from score_validation import write_tsv
ROOT=HERE/'snapshot';OUT=HERE/'analysis'
METHODS=['SpecHLA-IPD365-noFreq','SpecHLA-IPD365-long-noFreq','T1K-four-field']+[f'{a}:{m}' for a in ('full','hprc','asian_matched') for m in ('DogoHLA','DogoHLA-no-graph')]

def terminal_failure(path):
 marker=path/'TERMINAL_FAILURE.json'
 if not marker.exists():return False
 record=json.loads(marker.read_text())
 return record['manifest_sha256']==sha(path/'manifest.json') and json.loads((path/'manifest.json').read_text())['status']=='failed'

def locate(d,m):
 if m=='T1K-four-field':return ROOT/'t1k4'/d
 if m.startswith('SpecHLA'):return ROOT/'arms/full/runs'/d/('native-long' if '-long-' in m else 'native')
 a,n=m.split(':');return ROOT/'arms'/a/'final/runs'/d/n

def main():
 donors=list(csv.DictReader((HERE/'cohort.tsv').open(),delimiter='\t'))
 truth=truth_a(REPO/'hla-analysis/results/sequence_catalogue.tsv',{d['donor'] for d in donors});nom=Nomenclature()
 rows=[];states=[];hashes={}
 for d in donors:
  donor=d['donor']
  for method in METHODS:
   p=locate(donor,method);names={};state='not_started';error='NONE'
   if (p/'manifest.json').exists():
    try:
     manifest=json.loads((p/'manifest.json').read_text());state=manifest['status']
     if (p/'COMPLETE').exists() and completed(p,manifest['configuration']):
      state='complete';f=p/(donor+'_genotype.tsv' if method=='T1K-four-field' else 'hla.result.txt')
      names=t1k_calls(f) if method=='T1K-four-field' else read_result(f)
      hashes[str(f.relative_to(ROOT))]=sha(f)
     elif state=='complete':state='incomplete'
    except Exception as e:state='invalid';error=repr(e)
   states.append(dict(donor=donor,method=method,state=state,error=error))
   for gene in GENES:
    recs=truth.get((donor,gene))
    for level in ('two_field','four_field'):
     slots=(genomic_truth(recs) if level=='four_field' else truth_slots(recs,level,nom,gene)) if recs else None
     pair=names.get(gene,[None,None]);pred=four_pair(pair) if level=='four_field' else pred_pair(pair,level,nom,gene)
     called,correct,matches=compare(pred,slots) if slots is not None else (0,0,0)
     rows.append(dict(donor=donor,stratum=d['stratum'],family=d['family'],method=method,gene=gene,level=level,state=state,eligible=int(slots is not None),called=called,correct=correct,allele_matches=matches,allele1=pair[0] or 'NO_CALL',allele2=pair[1] or 'NO_CALL'))
 summary=[]
 for method in METHODS:
  for stratum in ('ALL','EAS','SAS','EUR','AFR'):
   ds={d['donor'] for d in donors if stratum=='ALL' or d['stratum']==stratum}
   for level in ('two_field','four_field'):
    rs=[r for r in rows if r['method']==method and r['donor'] in ds and r['level']==level]
    summary.append(dict(method=method,stratum=stratum,level=level,planned_donors=len(ds),completed_donors=sum(s['state']=='complete' for s in states if s['method']==method and s['donor'] in ds),**{k:sum(r[k] for r in rs) for k in ('eligible','called','correct','allele_matches')}))
 OUT.mkdir(exist_ok=True);write_tsv(OUT/'gene_scores.tsv',rows);write_tsv(OUT/'summary.tsv',summary);write_tsv(OUT/'status.tsv',states)
 (OUT/'provenance.json').write_text(json.dumps(dict(cohort_sha256=sha(HERE/'cohort.tsv'),output_sha256=hashes,truth='Exact CDS matches for two fields; exact genomic matches for four fields. No padding of unresolved labels. Ambiguities preserved; expression suffix ignored at the four-numeric-field endpoint.'),indent=2)+'\n')
 report=['# Matched 64-donor benchmark status','','All methods use IPD 3.65. Four-field truth requires exact genomic identity for both haplotypes; two-field truth requires exact CDS identity. These are assembly-derived labels. Independent Gourraud concordance is reported separately.','','| Method | Complete | Two-field genotype pairs | Four-field genotype pairs |','|---|---:|---:|---:|']
 for m in METHODS:
  rs=[s for s in summary if s['method']==m and s['stratum']=='ALL'];n=rs[0]['completed_donors']
  scores=[f"{s['correct']}/{s['eligible']}" if n==64 else 'pending' for s in rs]
  report.append(f"| {m} | {n}/64 | {' | '.join(scores)} |")
 report += ['','`summary.tsv` also separates EAS, SAS, EUR and AFR. Pending methods retain all planned truth denominators in TSV output; headline accuracy is withheld until the method finishes all 64 donors. Calls on eligible truth that remain unresolved are failures, not removed from the denominator.','']
 (OUT/'REPORT.md').write_text('\n'.join(report));print('\n'.join(report))

if __name__=='__main__':main()
