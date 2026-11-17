"""Count builder patches blocked by the blind-region guard (notes kind near_blind_region), blind regions,
and kept patches per donor from proposals.json files.  python guard_table.py DIR_A LABEL_A [DIR_B LABEL_B ...]
DIR/<donor>/C3/proposals.json"""
import json
import sys
from pathlib import Path
args = sys.argv[1:]
sets = [(Path(args[i]), args[i + 1]) for i in range(0, len(args), 2)]
donors = sorted({p.parent.parent.name for d, _ in sets for p in d.glob('*/C3/proposals.json')})
print('| donor | ' + ' | '.join(f'{l}: blind regions / blocked by guard / patches kept / decisions' for _, l in sets) + ' |')
print('|---|' + '---|' * len(sets))
tot = {l: [0, 0, 0] for _, l in sets}
for d in donors:
    cells = []
    for base, l in sets:
        p = base / d / 'C3/proposals.json'
        if not p.exists():
            cells.append('missing')
            continue
        a = json.loads(p.read_text())
        hs = [h for h in a['haplotypes'] if 'skipped' not in h]
        nb = sum(h['n_blind'] for h in hs)
        bl = sum(1 for h in hs for n in h['notes'] if n.get('kind') == 'near_blind_region')
        kp = sum(h['n_patches'] for h in hs)
        tot[l][0] += nb; tot[l][1] += bl; tot[l][2] += kp
        cells.append(f"{nb} / {bl} / {kp} / {','.join(h['decision'] for h in hs)}")
    print(f'| {d} | ' + ' | '.join(cells) + ' |')
print('| total | ' + ' | '.join(f'{v[0]} / {v[1]} / {v[2]} / ' for v in tot.values()) + ' |')
