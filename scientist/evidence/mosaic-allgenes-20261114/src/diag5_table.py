import json, glob, collections, sys
d = sys.argv[1]; show = sys.argv[2] if len(sys.argv) > 2 else None
vs = sorted({v for f in glob.glob(d + '/*.json') for v in json.load(open(f))['variants']})
for rs in ['D-sim', 'D-real']:
    for v in vs:
        agg = collections.defaultdict(collections.Counter)
        for f in sorted(glob.glob(f'{d}/{rs}_*.json')):
            j = json.load(open(f))
            if v not in j['variants']: continue
            for h in j['variants'][v]['haps']:
                k = 'MOSAIC' if (h['truth_label'] == 'mosaic' and h['truth_absent']) else ('EXACT' if not h['truth_absent'] else 'OTHER')
                a = agg[k]; a['n'] += 1; a['patched'] += h['n_patches'] > 0; a['int_exact'] += h['interior_exact']; a['byte'] += h['byte_exact']
                a['mm0'] += h['mismatch_after'] == 0; a['better'] += h['mismatch_after'] < h['mismatch_before']; a['worse'] += h['mismatch_after'] > h['mismatch_before']
                a['mm0_before'] += h['mismatch_before'] == 0
                a['patched_not_truth'] += (h['n_patches'] > 0 and h['mismatch_after'] != 0)
        print(rs, v, ' | '.join(f"{k} " + ' '.join(f'{x}={a[x]}' for x in ['n', 'patched', 'better', 'worse', 'mm0', 'mm0_before', 'int_exact', 'byte', 'patched_not_truth']) for k, a in sorted(agg.items())))
if show:
    for rs in ['D-real', 'D-sim']:
        print('==', rs, show)
        for f in sorted(glob.glob(f'{d}/{rs}_*.json')):
            j = json.load(open(f))
            for h in j['variants'][show]['haps']:
                ps = ' '.join(f"[{p['cl']}-{p['cr']} {p['len_cur']}>{p['len_new']} cur{p['cur_total']} new{p['n_new']}{p['fwd_rev']} exp{p['exp_h']} {p['kind'][:6]}]" for p in h['patches'][:7])
                nt = ' '.join(f"({a['kind'][:9]} {a['cl']} cur{a['cur_total']} new{a['n_new']} exp{a['exp_h']})" for a in h['notes'][:4])
                print(j['donor'], j['pass1'][0][9:] + '|' + j['pass1'][1][9:], 'S' if j['same_family'] else '', h['truth_label'][:5] + ('*' if h['truth_absent'] else ''), 'fam', h['family_rows'], 'ev', h['events'], 'res', h['resolvable'], '| mm', h['mismatch_before'], '>', h['mismatch_after'], 'int', int(h['interior_exact']), 'byte', int(h['byte_exact']), 'inref', int(h['in_reference']), 'spell', h['spell_mismatch'], h['spell_switches'], 'rem', h['remaining_diff_cols'][:5], '|', ps, nt)
