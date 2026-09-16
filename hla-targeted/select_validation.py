#!/usr/bin/env python3
"""Select fresh donors by availability and a fixed hash, before read inspection."""
import csv,hashlib,json,datetime
from pathlib import Path
from collections import defaultdict
ROOT=Path(__file__).resolve().parent
def main():
    rows=list(csv.DictReader(open(ROOT.parent/'hla-structural/source/catalogue.tsv'),delimiter='\t'))
    seen=set((ROOT.parent/'hla-structural/source/evaluation_donors.txt').read_text().splitlines())
    used_families={r['family'] for r in rows if r['donor'] in seen}
    manifest=dict(line.strip().split('\t') for line in open(ROOT/'source/C4Investigator/resources/1000Genomes_resources/tgp_full_30x.tsv'))
    groups=defaultdict(list)
    for r in rows:
        if r['locus']=='RCCX':groups[r['donor']].append(r)
    eligible=[d for d,rs in groups.items() if d not in seen and d in manifest and len(rs)==2 and all(r['family'] not in used_families for r in rs)]
    ranked=sorted(eligible,key=lambda d:hashlib.sha256(('targeted-c4-20260916:'+d).encode()).hexdigest())
    selected=[];families=set()
    for d in (['HG02735'] if 'HG02735' in eligible else [])+ranked:
        fam=groups[d][0]['family']
        if fam in families:continue
        selected.append(d);families.add(fam)
        if len(selected)==24:break
    out=[dict(donor=d,family=groups[d][0]['family'],cohort=groups[d][0]['cohort'],url='https://s3.amazonaws.com/1000genomes/'+manifest[d],selection='prespecified_matched_HLA_case' if d=='HG02735' else 'fixed_hash_no_read_outcomes') for d in selected]
    with open(ROOT/'source/validation_donors.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
    (ROOT/'source/validation_selection.json').write_text(json.dumps(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),eligible=len(eligible),selected=len(out),excluded_previous_donors=sorted(seen),excluded_previous_families=sorted(used_families),selection='HG02735 if eligible, then SHA256 ranking, one donor per new family; freeze all selected families out of marker discovery and path references'),indent=2)+'\n')
    print('Eligible fresh donors',len(eligible),'selected',len(out));print(' '.join(selected))
if __name__=='__main__':main()
