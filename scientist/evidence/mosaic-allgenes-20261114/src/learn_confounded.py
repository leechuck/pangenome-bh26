"""Truth-free list of 'confounded' positions: events at which, in a donor whose two DRB1 alleles
belong to different families, reads consistent with ONE haplotype support two versions at once
(own version and another), or pile up. One haplotype cannot carry two versions, so such reads come
from a second copy elsewhere. Learned from read evidence only (proposals.json), no truth.
  python learn_confounded.py READSET[,READSET]"""
import json
import sys
from pathlib import Path

S2 = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113')
SUB = sys.argv[2] if len(sys.argv) > 2 else 'C3'
wins, donors = [], 0
for rs in sys.argv[1].split(','):
    for d in sorted((S2 / 'work' / rs).iterdir()):
        p = d / SUB / 'proposals.json'
        if not p.exists():
            continue
        j = json.loads(p.read_text())
        if j['same_family']:
            continue
        donors += 1
        for h in j['haplotypes']:
            for n in h.get('notes', []):
                if n['kind'] == 'pile_up' or (n['kind'] == 'own_version_still_supported' and not n.get('shared')
                                              and n['cur_total'] >= n.get('thr', 3)):
                    wins.append((n['cl'], n['cr']))
wins = sorted(set(wins))
merged = []
for a, b in wins:
    if merged and a <= merged[-1][1]:
        merged[-1][1] = max(merged[-1][1], b)
    else:
        merged.append([a, b])
out = S2 / 'inputs/confounded_events.json'
out.parent.mkdir(exist_ok=True)
out.write_text(json.dumps(dict(source_readsets=sys.argv[1], donors_used=donors, n_raw=len(wins), windows=merged), indent=1) + '\n')
print(len(wins), 'raw windows ->', len(merged), 'merged; donors', donors)
print(merged)
