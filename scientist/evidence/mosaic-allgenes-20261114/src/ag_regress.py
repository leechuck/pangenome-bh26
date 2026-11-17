"""Regression check (item 5): the generalised builder, DRB1 only, frozen confounded list, abstain off,
against the frozen DRB1 stage outputs. Compares per haplotype decisions, patched-sequence hashes, patch
and note lists, and the OVERRIDES / additions files. Hashes and names only.
  python ag_regress.py TAG"""
import hashlib
import json
import sys
from pathlib import Path

E = Path('/encrypted/e3008/hohndor-dogohla-20261020')
S2 = E / 'mosaic-allgenes-20261114'
OLD = E / 'mosaic-drb1-20261113' / 'work'
tag = sys.argv[1]
NEWKEYS = {'pass1_msa_row', 'unresolved', 'unresolved_rules'}
out = ['readset\tdonor\tarm\thap\told_decision\tnew_decision\told_sha16\tnew_sha16\tsame_haplotype_record\tsame_overrides\tsame_additions']
n = dict(tasks=0, hap=0, same_hap=0, same_sha=0, same_ov=0, same_add=0, missing=0)
for rs in ('D-sim', 'D-real'):
    for d in sorted((OLD / rs).iterdir()):
        for arm, sub in (('C', 'C3'), ('D', 'D3')):
            o = d / sub / 'proposals.json'
            p = S2 / 'work' / tag / rs / d.name / arm / 'DRB1' / 'proposals.json'
            if not o.exists():
                continue
            n['tasks'] += 1
            if not p.exists():
                n['missing'] += 1
                out.append(f'{rs}\t{d.name}\t{arm}\tNA\tmissing')
                continue
            jo, jn = json.loads(o.read_text()), json.loads(p.read_text())
            ovo = json.loads((d / sub / 'OVERRIDES.json').read_text())
            ovn = json.loads((p.parent / 'OVERRIDES.json').read_text())
            same_ov = ovo['overrides'] == ovn['overrides']
            h1 = hashlib.sha256((d / sub / 'additions.fa').read_bytes()).hexdigest()
            h2 = hashlib.sha256((p.parent / 'additions.fa').read_bytes()).hexdigest()
            n['same_ov'] += same_ov
            n['same_add'] += h1 == h2
            for h, (a, b) in enumerate(zip(jo['haplotypes'], jn['haplotypes'])):
                b2 = {k: v for k, v in b.items() if k not in NEWKEYS}
                same = a == b2
                n['hap'] += 1
                n['same_hap'] += same
                n['same_sha'] += a.get('sha256') == b.get('sha256')
                out.append('\t'.join(map(str, (rs, d.name, arm, h, a.get('decision', a.get('skipped')), b.get('decision', b.get('skipped')),
                                               (a.get('sha256') or '')[:16], (b.get('sha256') or '')[:16], same, same_ov, h1 == h2))))
            if len(jo['haplotypes']) != len(jn['haplotypes']):
                out.append(f'{rs}\t{d.name}\t{arm}\tNA\thaplotype count differs')
            for k in ('pass1', 'same_family', 'reads_assigned', 'reads_used', 'expected_depth', 'confounded_windows_loaded'):
                if jo.get(k) != jn.get(k):
                    out.append(f'{rs}\t{d.name}\t{arm}\tNA\tfield {k} differs: {jo.get(k)} vs {jn.get(k)}')
(S2 / 'results' / f'regression_{tag}.tsv').write_text('\n'.join(out) + '\n')
print(json.dumps(n))
