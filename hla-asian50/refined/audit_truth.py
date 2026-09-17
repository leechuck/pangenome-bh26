#!/usr/bin/env python3
"""Describe overlap between frozen truth sources; never change scoring labels."""
from pathlib import Path
from score_hla import table,compare,write
R=Path(__file__).resolve().parent/'results'
rows=[x for x in table(R/'hla_scores.tsv') if x['method']=='T1K']
assembly={(x['donor'],x['gene']):x for x in rows if x['endpoint']=='assembly_exact_CDS_twofield'}
out=[]
for x in rows:
 if x['endpoint']!='experimental_twofield':continue
 a=assembly.get((x['donor'],x['gene']))
 if a is None:continue
 pair=tuple(a['truth'].split(';'));assert all('/' not in v for v in pair)
 _,correct,matches=compare(pair,[set(v.split('/')) for v in x['truth'].split(';')])
 out.append(dict(donor=x['donor'],gene=x['gene'],experimental_truth=x['truth'],assembly_truth=a['truth'],compatible=correct,allele_matches=matches))
write(R/'truth_source_consistency.tsv',out)
print('Compatible truth labels:',sum(x['compatible'] for x in out),'/',len(out))
for x in out:
 if not x['compatible']:print(x)
