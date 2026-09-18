#!/usr/bin/env python3
"""Select 32 non-Asian donors using metadata/truth availability only, never predictions."""
import collections,csv,hashlib,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1]
sys.path.insert(0,str(HERE.parent))
from score import truth_a,GENES
from experiment import sha
seed='DogoHLA-ancestry-20260918'

def rank(d):return hashlib.sha256((seed+d).encode()).hexdigest()

def write(path,rows):
 with open(path,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

asia=list(csv.DictReader(open(HERE.parent/'validation/20260918/cohort.tsv'),delimiter='\t'))
source=list(csv.DictReader(open(REPO/'hla-bench/source/cohorts.tsv'),delimiter='\t'))
truth=truth_a(REPO/'hla-analysis/results/sequence_catalogue.tsv',{r['donor'] for r in source})
used={r['family'] for r in asia};selected=[]
for ancestry in ('EUR','AFR'):
 pools=collections.defaultdict(list)
 for r in source:
  if r['cohort']=='test_A_loo' and r['superpopulation']==ancestry and r['dev']=='0' and all((r['donor'],g) in truth for g in GENES):pools[r['population']].append(r)
 for pool in pools.values():pool.sort(key=lambda r:rank(r['donor']))
 chosen=[]
 while len(chosen)<16:
  before=len(chosen)
  for pop in sorted(pools,key=rank):
   while pools[pop] and pools[pop][0]['family'] in used:pools[pop].pop(0)
   if pools[pop] and len(chosen)<16:
    r=pools[pop].pop(0);chosen.append(r);used.add(r['family'])
  if len(chosen)==before:raise ValueError('Insufficient unrelated donors')
 selected+=chosen
rows=[dict(donor=r['donor'],stratum=r['stratum'],population=r['population'],family=r['family'],fold=0,cohort='Asian',cram=r['public_cram_url'],selection_sha256=rank(r['donor'])) for r in asia]
rows += [dict(donor=r['donor'],stratum=r['superpopulation'],population=r['population'],family=r['family'],fold=0,cohort='nonAsian',cram=r['cram'],selection_sha256=rank(r['donor'])) for r in selected]
assert len(rows)==64 and len({r['family'] for r in rows})==64
write(HERE/'cohort.tsv',rows)
write(HERE/'donors_folds.tsv',[{k:r[k] for k in ('donor','fold','cram')} for r in rows])
write(HERE/'donors_nonasian.tsv',[{k:r[k] for k in ('donor','cram')} for r in rows[32:]])
groups=list(csv.DictReader(open(HERE.parent/'source/validation_groups.tsv'),delimiter='\t'))
donors={r['donor'] for r in rows}
metadata=list(csv.DictReader(open(REPO/'hla-typer/source/haplotype_donors.tsv'),delimiter='\t'))
known={r['hap_id'] for r in groups}
pedigree={r['SampleID']:r for r in csv.DictReader(open(REPO/'hla-pilot/source/1000G.ped'),delimiter=' ')}
for r in metadata:
 if r['hap_id'] not in known:
  groups.append(dict(hap_id=r['hap_id'],donor=r['donor_id'],group=pedigree.get(r['donor_id'],{}).get('FamilyID',r['donor_id']),pedigree_status='canonical_donor_metadata',cohort=r['cohort']))
write(HERE/'validation_groups.tsv',groups)
own_haps={r['hap_id'] for r in metadata if r['donor_id'] in donors}
hit={r['group'] for r in groups if r['donor'] in donors or r['hap_id'] in own_haps or r['group'] in used}
excluded=[dict(fold=0,sample=r['hap_id'].split('#')[0],hap_id=r['hap_id'],donor=r['donor']) for r in groups if r['group'] in hit]
assert not own_haps-{r['hap_id'] for r in excluded}, 'Missing alternate-assembly donor exclusions'
write(HERE/'excluded_paths.tsv',excluded)
(HERE/'selection.json').write_text(json.dumps(dict(seed=seed,cohort_sha256=sha(HERE/'cohort.tsv'),selection='Metadata-only, unique families, population round robin within 16 EUR and 16 AFR; complete assembly truth across eight genes',excluded_haplotypes=len(excluded),strata=dict(collections.Counter(r['stratum'] for r in rows)),scope='New non-Asian donors come from test_A_loo; now reserved for this prospective experiment, not available for a later untouched test.'),indent=2)+'\n')
print((HERE/'selection.json').read_text())
