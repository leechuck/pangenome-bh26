#!/usr/bin/env python3
"""Frozen-universe before/after comparison using SHA256 allele-sequence identity.
Never use allele numbers across VCFs. No-call and absent-site failures remain counted.
"""
import csv,gzip,hashlib,re,collections,argparse,json
from pathlib import Path
LOW,HIGH,OFFSET=28510121,33480577,28410700

def digest(s):return hashlib.sha256(s.encode()).digest()
def records(p):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt') as f:
        for line in f:
            if line.startswith('##'):continue
            v=line.rstrip('\n').split('\t')
            if line.startswith('#CHROM'):ss=v[9:];continue
            yield v,ss

def indices(s):
    x=re.split(r'[/|]',s)
    return tuple(map(int,x)) if len(x)==2 and all(i.isdigit() for i in x) else None

def prediction(path):
    out={}
    for v,ss in records(path):
        aa=[v[3]]+v[4].split(',');hs=[digest(x) for x in aa];key=(int(v[1]),hs[0]);assert key not in out
        ix=indices(v[9].split(':')[v[8].split(':').index('GT')]);gt=tuple(sorted(hs[i] for i in ix)) if ix else None
        out[key]=(gt,set(hs))
    return out

def write(p,rows):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--original',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    with open(a.original/'source/donors.tsv') as f:ds=list(csv.DictReader(f,delimiter='\t'))
    truth={};audit=collections.Counter()
    for v,ss in records(a.original/'graph/full.top.vcf.gz'):
        pos=int(v[1]);ref=v[3];aa=[ref]+v[4].split(',');hs=[digest(x) for x in aa];key=(pos,hs[0])
        if pos+OFFSET<LOW or pos+OFFSET+len(ref)-1>HIGH:continue
        assert key not in truth
        kind='SNV' if all(len(x)==1 and x in 'ACGT' for x in aa) else 'SV_length_site' if any(abs(len(x)-len(ref))>=50 for x in aa) else 'small_complex'
        conflicts=set()
        for tag in v[7].split(';'):
            if tag.startswith('CONFLICT='):conflicts.update(tag.split('=',1)[1].split(','))
        idx={s:i+9 for i,s in enumerate(ss)};gi=v[8].split(':').index('GT');ts={}
        for d in ds:
            s=d['donor'];ii=indices(v[idx[s]].split(':')[gi]);g=None;sv=False
            if ii and s not in conflicts and all(set(aa[i])<=set('ACGT') for i in ii):
                g=tuple(sorted(hs[i] for i in ii));sv=any(abs(len(aa[i])-len(ref))>=50 for i in ii)
            ts[s]=(g,sv)
        truth[key]=(kind,ts)
    linear={};seen=set();dups=set()
    for v,ss in records(a.original/'accuracy/source/linear.vcf.gz'):
        pos=int(v[1])-OFFSET
        if pos in seen:dups.add(pos)
        seen.add(pos);aa=[v[3]]+v[4].split(',');key=(pos,digest(v[3]))
        if key not in truth or truth[key][0]!='SNV' or v[6]!='PASS' or not all(len(x)==1 and x in 'ACGT' for x in aa):continue
        hs=[digest(x) for x in aa];gi=v[8].split(':').index('GT');idx={s:i+9 for i,s in enumerate(ss)};ts={}
        for d in ds:
            ii=indices(v[idx[d['donor']]].split(':')[gi]);ts[d['donor']]=tuple(sorted(hs[i] for i in ii)) if ii else None
        linear[key]=ts
    linear={k:v for k,v in linear.items() if k[0] not in dups}
    rows=[]
    for d in ds:
        s=d['donor'];pred={}
        for prefix,root in [('old_',a.original),('',a.run)]:
            for arm in ['full','hprc']:pred[prefix+arm]=prediction(root/f'pangenie/{s}/{arm}/calls_genotyping.vcf')
        # Native missing support must preserve every previously emitted site in the analysed interval.
        for arm in ['full','hprc']:
            lost=(set(pred['old_'+arm])&set(truth))-set(pred[arm]);assert not lost,(s,arm,len(lost))
        assert (set(pred['hprc'])&set(truth))<=set(pred['full'])
        counts=collections.defaultdict(collections.Counter)
        for key,(kind,ts) in truth.items():
            t,sv=ts[s]
            if t is None:continue
            nonref=t!=(key[1],key[1]);cats=[kind,'ALL']+(['truth_SV_length'] if sv else [])
            universes=['all_truth']
            if key in pred['old_hprc']:universes.append('frozen_HPRC_sites')
            if key in pred['full'] and key in pred['hprc']:universes.append('shared_new_sites')
            if kind=='SNV' and key in pred['old_full'] and key in pred['old_hprc'] and key in linear:universes.append('frozen_threeway_SNV')
            for u in universes:
                arms=list(pred)+(['linear'] if u=='frozen_threeway_SNV' else [])
                for arm in arms:
                    item=(linear[key][s],set()) if arm=='linear' else pred[arm].get(key);g=item[0] if item else None
                    for cat in cats:
                        c=counts[u,cat,arm];c['n']+=1;c['site_present']+=item is not None;c['called']+=g is not None;c['correct']+=g==t
                        if nonref:c['nonref_n']+=1;c['nonref_correct']+=g==t
        for (u,kind,arm),c in sorted(counts.items()):rows.append(dict(donor=s,stratum=d['stratum'],population=d['population'],fold=d['fold'],universe=u,variant_class=kind,arm=arm,**{k:c[k] for k in ['n','site_present','called','correct','nonref_n','nonref_correct']}))
        print(s,'frozen variants scored',flush=True)
    write(a.out/'variant_per_donor.tsv',rows)
    # Independent implementation must reproduce every original all-truth count exactly.
    with open(a.original/'accuracy/results/panel_per_donor.tsv') as f:old=list(csv.DictReader(f,delimiter='\t'))
    rr={(r['donor'],r['variant_class'],r['arm']):r for r in rows if r['universe']=='all_truth'}
    for x in old:
        if x['universe']!='all_truth':continue
        got=rr[x['donor'],x['variant_class'],'old_'+x['arm']]
        for k in ['n','site_present','called','correct','nonref_n','nonref_correct']:assert int(x[k])==got[k],(x['donor'],x['variant_class'],k)
    (a.out/'variant_checks.json').write_text(json.dumps(dict(donors=len(ds),old_metrics_exactly_reproduced=True,all_previous_output_sites_preserved=True,full_output_contains_HPRC_output_sites=True),indent=2)+'\n')
if __name__=='__main__':main()
