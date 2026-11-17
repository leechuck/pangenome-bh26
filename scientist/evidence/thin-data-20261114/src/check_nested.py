"""Pilot check: downsampled files are nested, complete pairs, expected counts. usage: check_nested.py rs donor ..."""
import gzip, json, sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import lib as L
C = L.C
out = []
for arg in sys.argv[1:]:
    rs, d = arg.split(':')
    rd = L.TS / 'work/reads' / rs / d
    sub = C.load_json(rd / 'SUBSAMPLE.json')
    names = {}
    for n in L.LEVELS:
        for m in (1, 2):
            with gzip.open(rd / f'd{n}_r{m}.fq.gz', 'rt') as f:
                names[(n, m)] = [l.split()[0].removesuffix(f'/{m}') for i, l in enumerate(f) if i % 4 == 0]
    full = []
    with gzip.open(L.arm_info(rs, d, 'A')['r1'], 'rt') as f:
        full = [l.split()[0].removesuffix('/1') for i, l in enumerate(f) if i % 4 == 0]
    rec = dict(readset=rs, donor=d, full_pairs=len(full), levels={})
    ok = True
    for n, thr in L.LEVELS.items():
        r1, r2 = names[(n, 1)], names[(n, 2)]
        exp = [x for x in full if L.read_key(x) < thr]
        rec['levels'][n] = dict(pairs=len(r1), expected_by_rule=len(exp), fraction=len(r1) / len(full),
                                mates_match=r1 == r2, equals_rule=r1 == exp,
                                manifest_pairs=sub[str(n)]['pairs'])
        ok &= r1 == r2 and r1 == exp and sub[str(n)]['pairs'] == len(r1)
    for a, b in ((16, 8), (8, 4), (4, 2)):
        sa, sb = set(names[(a, 1)]), set(names[(b, 1)])
        rec['levels'][a]['subset_of_%d' % b] = sa <= sb
        ok &= sa <= sb
    rec['all_ok'] = ok
    out.append(rec)
    print(json.dumps(rec))
C.save_json(L.TS / 'results/pilot_nested_check.json', out)
