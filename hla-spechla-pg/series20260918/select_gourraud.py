#!/usr/bin/env python3
"""Stage the existing experimental-truth cohort, with independent leakage checks."""
import collections,csv,hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1]
def read(p,delim='\t'):
 with p.open() as f:return list(csv.DictReader(f,delimiter=delim))
rows=[r for r in read(REPO/'hla-bench/source/cohorts.tsv') if r['cohort']=='test_B_gourraud']
assert len(rows)==946 and len({r['donor'] for r in rows})==946
ped={r['SampleID']:r['FamilyID'] for r in read(REPO/'hla-pilot/source/1000G.ped',' ')}
panel=read(REPO/'hla-typer/source/haplotype_donors.tsv')
donors={r['donor_id'] for r in panel if r['cohort']!='REF'}
families={ped.get(d,d) for d in donors}
urls=dict(line.split('\t')[:2] for line in (REPO/'hla-targeted/source/C4Investigator/resources/1000Genomes_resources/tgp_full_30x.tsv').read_text().splitlines())
out=[]
for r in rows:
 d=r['donor'];assert d in ped and d not in donors and ped[d] not in families
 out.append(dict(donor=d,stratum=r['superpopulation'],population=r['population'],family=ped[d],fold=0,cohort='Gourraud2014',cram='https://s3.amazonaws.com/1000genomes/'+urls[d]))
for name,keys in [('cohort.gourraud.tsv',list(out[0])),('recruit.gourraud.tsv',['donor','cram'])]:
 with (HERE/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(out)
record=dict(donors=len(out),families=len({r['family'] for r in out}),ancestry=dict(collections.Counter(r['stratum'] for r in out)),graph_donor_overlap=0,graph_family_overlap=0,cohort_sha256=hashlib.sha256((HERE/'cohort.gourraud.tsv').read_bytes()).hexdigest(),truth='Experimental Gourraud 2014, ambiguity-aware exon typing; no four-field truth')
(HERE/'gourraud-selection.json').write_text(json.dumps(record,indent=2)+'\n');print(record)
