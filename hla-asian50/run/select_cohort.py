#!/usr/bin/env python3
"""Freeze available strata without using prediction outcomes; reserve Arabian slots."""
import csv,hashlib,json,datetime,collections
from pathlib import Path
R=Path(__file__).resolve().parent;P=R.parent.parent
seed='mixed-asian-pilot-20260917-v1'
def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def rank(d):return hashlib.sha256((seed+':'+d).encode()).hexdigest()
def write(p,rs):
 with open(R/'source'/p,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
candidates=read(R.parent/'source/mixed_asian_candidates.tsv')
urls=dict(x.strip().split('\t') for x in open(P/'hla-targeted/source/C4Investigator/resources/1000Genomes_resources/tgp_full_30x.tsv'))
selected=[];families=set()
for stratum in ['EAS','SAS']:
 ss=[]
 while len(ss)<20:
  counts=collections.Counter(r['population'] for r in ss)
  rr=[r for r in candidates if r['stratum']==stratum and r['family'] not in families and r['donor'] in urls]
  assert rr
  r=min(rr,key=lambda r:(counts[r['population']],-int(r['experimental_HLA_labels']),rank(r['donor'])))
  ss.append(dict(r));families.add(r['family'])
 for i,r in enumerate(sorted(ss,key=lambda r:(r['population'],rank(r['donor'])))):
  r['fold']=i%5;r['public_cram_url']='https://s3.amazonaws.com/1000genomes/'+urls[r['donor']]
  r['selection_sha256']=rank(r['donor']);selected.append(r)
selected.sort(key=lambda r:(r['fold'],r['stratum'],r['population'],r['donor']))
write('donors.tsv',selected)
(R/'source/donors.txt').write_text('\n'.join(r['donor'] for r in selected)+'\n')
cat=read(P/'hla-structural/source/catalogue.tsv');fam={r['donor']:r['family'] for r in cat}
manifest=read(P/'hla-audit/2026-09-16/panel_manifest.tsv');ex=[]
for k in range(5):
 excluded={r['family'] for r in selected if r['fold']==k}
 for r in manifest:
  if fam.get(r['donor_id'],r['donor_id']) in excluded:ex.append(dict(fold=k,sample=r['sample'],hap_id=r['hap_id'],donor=r['donor_id']))
write('excluded_paths.tsv',ex)
meta=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),seed=seed,status='initial_40_available_strata_Arabian_pending',reserved_Arabian_slots=10,strata=dict(collections.Counter(r['stratum'] for r in selected)),populations=dict(collections.Counter(r['population'] for r in selected)),experimental_HLA_donors=sum(int(r['experimental_HLA_labels']) for r in selected),fold_sizes=dict(collections.Counter(r['fold'] for r in selected)),selection='Population-balanced within stratum, experimental-label availability then fixed hash; no prediction outcomes',files={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (R/'source').glob('*.tsv')})
(R/'source/cohort_freeze.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
