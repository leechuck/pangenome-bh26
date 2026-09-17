#!/usr/bin/env python3
"""Verify cohort independence, paired ablations and archive provenance."""
import csv,hashlib,json,collections,xml.etree.ElementTree as ET
from pathlib import Path
R=Path(__file__).resolve().parent
def read(name):return list(csv.DictReader(open(R/'source'/name),delimiter='\t'))
ds=read('donors.tsv');ex=read('excluded_paths.tsv');oracle=read('structural_representation_oracle.tsv')
assert len(ds)==len({r['donor'] for r in ds})==len({r['family'] for r in ds})==50
assert collections.Counter(r['fold'] for r in ds)=={str(k):10 for k in range(5)}
for d in ds:
    matching=[e for e in ex if e['fold']==d['fold'] and e['donor']==d['donor']]
    assert len(matching)>=2
    assert all(e['family']==d['family'] for e in matching)
for locus in ['DRB','RCCX']:
    for arm in ['full_panel','HPRC_only']:
        rr=[r for r in oracle if r['locus']==locus and r['arm']==arm]
        assert len(rr)==100 and all(n==2 for n in collections.Counter(r['donor'] for r in rr).values())
for name in ['public_reads.tsv','remote_availability.tsv']:
    assert {r['donor'] for r in read(name)}=={r['donor'] for r in ds}
assert len({r['donor'] for r in read('experimental_HLA_truth.tsv')})==23
for m in json.loads((R/'literature/download_manifest.json').read_text()):
    p=R/'literature'/(m['name']+'.xml')
    assert hashlib.sha256(p.read_bytes()).hexdigest()==m['sha256']
    root=ET.parse(p).getroot()
    dois={''.join(e.itertext()) for e in root.findall('./front/article-meta/article-id') if e.get('pub-id-type')=='doi'}
    assert m['doi'] in dois,(m['name'],dois)
print('PASS: 50 donors, five family-separated folds, both ablations, truth tiers, read manifests, five full-text articles.')
