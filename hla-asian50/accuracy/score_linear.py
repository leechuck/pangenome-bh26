#!/usr/bin/env python3
"""Three-way comparison restricted to sequence-matched single-base VCF sites.

NYGC calls are a baseline, never truth. Missing records are never filled as 0/0.
No SNV decomposition inside complex graph bubbles is attempted.
"""
import argparse,csv,collections,json
from pathlib import Path
from score_panels import records,gt,predictions,write,LOW,HIGH,OFFSET

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    ds=list(csv.DictReader(open(a.run/'source/donors.tsv'),delimiter='\t'));truth={};audit=collections.Counter()
    for v,ss in records(a.run/'graph/full.top.vcf.gz'):
        pos=int(v[1])+OFFSET;aa=[v[3]]+v[4].split(',')
        if not LOW<=pos<=HIGH or not all(len(x)==1 and x in 'ACGT' for x in aa):continue
        conflicts=set()
        for tag in v[7].split(';'):
            if tag.startswith('CONFLICT='):conflicts.update(tag.split('=',1)[1].split(','))
        gi=v[8].split(':').index('GT');idx={s:i+9 for i,s in enumerate(ss)}
        truth[pos]=(v[3],{d['donor']:None if d['donor'] in conflicts else gt(v[idx[d['donor']]].split(':')[gi],aa) for d in ds})
    linear={};seen=set();duplicates=set()
    for v,ss in records(a.run/'accuracy/source/linear.vcf.gz'):
        pos=int(v[1]);aa=[v[3]]+v[4].split(',')
        if pos in seen:duplicates.add(pos)
        seen.add(pos)
        if pos not in truth:continue
        if v[6]!='PASS':audit['not_PASS']+=1;continue
        if not all(len(x)==1 and x in 'ACGT' for x in aa):audit['not_pure_SNV']+=1;continue
        assert v[3]==truth[pos][0],('reference mismatch',pos)
        gi=v[8].split(':').index('GT');idx={s:i+9 for i,s in enumerate(ss)}
        linear[pos]={d['donor']:gt(v[idx[d['donor']]].split(':')[gi],aa) for d in ds}
    for pos in duplicates:linear.pop(pos,None)
    audit['duplicate_positions_excluded']=len(duplicates);audit['graph_SNV_sites']=len(truth);audit['matched_linear_PASS_SNV_sites']=len(linear)
    rows=[]
    for d in ds:
        s=d['donor'];pr={arm:predictions(a.run/f'pangenie/{s}/{arm}/calls_genotyping.vcf') for arm in ['full','hprc']};counts=collections.defaultdict(collections.Counter)
        for pos,(ref,ts) in truth.items():
            t=ts[s]
            if t is None:continue
            key=(pos-OFFSET,ref)
            if not all(key in pr[x] for x in pr):continue
            gs={arm:pr[arm][key][0] for arm in pr};gs['linear']=linear[pos][s] if pos in linear else None
            for universe in ['graph_shared_SNV']+(['three_way_shared_PASS_SNV'] if pos in linear else []):
                for arm,g in gs.items():
                    c=counts[universe,arm];c['n']+=1;c['called']+=g is not None;c['correct']+=g==t
                    if t!=(ref,ref):c['nonref_n']+=1;c['nonref_correct']+=g==t
        for (universe,arm),c in sorted(counts.items()):
            rows.append(dict(donor=s,stratum=d['stratum'],population=d['population'],universe=universe,variant_class='SNV',arm=arm,**{k:c[k] for k in ['n','called','correct','nonref_n','nonref_correct']}))
        print(s,'linear scored',flush=True)
    write(a.out/'linear_per_donor.tsv',rows);(a.out/'linear_audit.json').write_text(json.dumps(dict(audit),indent=2)+'\n')
if __name__=='__main__':main()
