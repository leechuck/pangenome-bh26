"""Truth-free per-gene list of 'confounded' positions (same rule as learn_confounded.py, gene as a
parameter): events at which, in a donor whose two alleles of the gene belong to different families,
reads consistent with ONE haplotype support two versions at once (own version and another), or pile
up. Learned from read evidence only (proposals.json of the builder), no truth.
  python ag_learn_confounded.py TAG READSET[,READSET] ARM OUTDIR [GENES]"""
import json
import sys
from pathlib import Path

S2 = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114')
GENES = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')
tag, readsets, arm, outdir = sys.argv[1], sys.argv[2].split(','), sys.argv[3], Path(sys.argv[4])
genes = sys.argv[5].split(',') if len(sys.argv) > 5 else GENES
assert not {'T-sim', 'T-real'} & set(readsets)
outdir.mkdir(parents=True, exist_ok=True)
for g in genes:
    wins, donors, missing = [], 0, 0
    for rs in readsets:
        for d in sorted((S2 / 'work' / tag / rs).iterdir()):
            p = d / arm / g / 'proposals.json'
            if not p.exists():
                missing += 1
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
    (outdir / f'confounded_{g}.json').write_text(json.dumps(dict(gene=g, source_tag=tag, source_readsets=sys.argv[2],
                                                                 arm=arm, donors_used=donors, missing=missing,
                                                                 n_raw=len(wins), windows=merged), indent=1) + '\n')
    print(g, len(wins), 'raw windows ->', len(merged), 'merged; donors', donors, 'missing', missing)
