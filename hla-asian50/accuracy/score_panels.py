#!/usr/bin/env python3
"""Sequence-aware, unphased held-out assembly concordance. Standard library only."""
import csv,gzip,json,re,argparse,collections
from pathlib import Path
OFFSET=28410700 # extracted reference is chr6:28410701-33580488, inclusive
LOW,HIGH=28510121,33480577 # intersection with primary read recruitment BED

def gt(value,alleles):
    bits=re.split(r'[/|]',value)
    if len(bits)!=2 or any(not x.isdigit() for x in bits): return None
    return tuple(sorted(alleles[int(x)] for x in bits))

def records(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt') as f:
        for line in f:
            if line.startswith('##'):continue
            v=line.rstrip('\n').split('\t')
            if line.startswith('#CHROM'): samples=v[9:];continue
            yield v,samples

def predictions(path):
    out={}
    for v,ss in records(path):
        key=(int(v[1]),v[3]); aa=[v[3]]+v[4].split(',')
        assert key not in out,('duplicate prediction',key)
        gi=v[8].split(':').index('GT')
        out[key]=(gt(v[9].split(':')[gi],aa),set(aa))
    return out

def write(path,rows):
    with open(path,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    donors=list(csv.DictReader(open(a.run/'source/donors.tsv'),delimiter='\t'))
    truth={};audit=collections.Counter()
    for v,ss in records(a.run/'graph/full.top.vcf.gz'):
        pos=int(v[1]);ref=v[3];aa=[ref]+v[4].split(',')
        if pos+OFFSET<LOW or pos+OFFSET+len(ref)-1>HIGH: audit['outside_recruited_primary_interval']+=1;continue
        key=(pos,ref);assert key not in truth
        kind='SNV' if all(len(x)==1 and x in 'ACGT' for x in aa) else 'SV_length_site' if any(abs(len(x)-len(ref))>=50 for x in aa) else 'small_complex'
        conflicts=set()
        for tag in v[7].split(';'):
            if tag.startswith('CONFLICT='):conflicts.update(tag.split('=',1)[1].split(','))
        gi=v[8].split(':').index('GT');indices={s:i+9 for i,s in enumerate(ss)}
        calls={}
        for d in donors:
            s=d['donor'];g=gt(v[indices[s]].split(':')[gi],aa)
            if s in conflicts: g=None;audit['donor_conflict_masked']+=1
            if g and any(set(x)-set('ACGT') for x in g):g=None
            calls[s]=g
        truth[key]=(kind,calls)
    write(a.out/'reference_interval.tsv',[dict(contig='chr6',start_1based=LOW,end_inclusive=HIGH,graph_offset=OFFSET)])
    allrows=[];discord=[]
    for d in donors:
        s=d['donor'];pred={arm:predictions(a.run/f'pangenie/{s}/{arm}/calls_genotyping.vcf') for arm in ['full','hprc']}
        counts=collections.defaultdict(collections.Counter)
        for key,(kind,tcalls) in truth.items():
            t=tcalls[s]
            if t is None:continue
            ref=key[1];nonref=t!=(ref,ref);shared=all(key in pred[arm] for arm in pred)
            cats=[kind,'ALL']
            if any(abs(len(x)-len(ref))>=50 for x in t):cats.append('truth_SV_length')
            for universe in ['all_truth']+(['shared_sites'] if shared else []):
                for arm in pred:
                    item=pred[arm].get(key);g=item[0] if item else None;represented=bool(item and all(x in item[1] for x in t))
                    for cat in cats:
                        c=counts[universe,cat,arm];c['n']+=1;c['site_present']+=item is not None;c['represented']+=represented;c['called']+=g is not None;c['correct']+=g==t
                        if nonref:c['nonref_n']+=1;c['nonref_correct']+=g==t
            if shared and kind=='SV_length_site':
                f=pred['full'][key][0];h=pred['hprc'][key][0]
                if (f==t)!=(h==t):
                    discord.append(dict(donor=s,stratum=d['stratum'],chr6_pos=key[0]+OFFSET,ref_length=len(ref),truth_lengths=','.join(map(str,map(len,t))),full_correct=int(f==t),hprc_correct=int(h==t),truth_nonref=int(nonref)))
        for (universe,kind,arm),c in sorted(counts.items()):
            allrows.append(dict(donor=s,stratum=d['stratum'],population=d['population'],fold=d['fold'],universe=universe,variant_class=kind,arm=arm,**{k:c[k] for k in ['n','site_present','represented','called','correct','nonref_n','nonref_correct']}))
        print(s,'scored',flush=True)
    write(a.out/'panel_per_donor.tsv',allrows)
    if discord:write(a.out/'sv_differential.tsv',discord)
    (a.out/'truth_audit.json').write_text(json.dumps(dict(audit),indent=2)+'\n')
if __name__=='__main__':main()
