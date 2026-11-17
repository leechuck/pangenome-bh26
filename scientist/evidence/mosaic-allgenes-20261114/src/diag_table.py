"""Local summary of DIAG.json files (numbers and names only). usage: diag_table.py DIR"""
import collections
import glob
import json
import sys

d = sys.argv[1]
J = [json.load(open(f)) for f in sorted(glob.glob(d + '/*.json'))]
J = [j for j in J if j.get('configs')]
print('donor-readsets:', len(J), 'complete:', sum(bool(j.get('complete')) for j in J))
print('\n== reads and truth')
for j in J:
    t = j['truth']
    c0 = j['configs'].get('R0|F|0.1', {})
    print(j['readset'], j['donor'], '|'.join(n.replace('HLA-DRB1*', '') for n in j['pass1']),
          'fam', j['family_sizes'], 'same' if j['same_family'] else '', 'assigned', j['reads_assigned'],
          'mateOther', j['mate_on_other_gene'], j['mate_other_genes'],
          '| truth:', ' ; '.join(
              f"{x['ipd_label']}{'*' if x['ipd_absent'] else ''} sw(all/famA/famB)="
              f"{x['spelling']['all']['switches']}/{x['spelling']['fam_a']['switches']}m{x['spelling']['fam_a']['mismatch']}"
              f"/{x['spelling']['fam_b']['switches']}m{x['spelling']['fam_b']['mismatch']} p1mm={x['pass1_mismatch']}"
              f" endoff_nearest={x['spelling']['all']['nearest_start_end_offset']} endoff_p1={x['pass1_start_end_offset']}"
              for x in t))
keys = sorted({k for j in J for k in j['configs']})
print('\n== per config, per read set: noise floor under the truth pair (N), signal (Z_pass1 - N), outcomes')
for rs in ('D-sim', 'D-real'):
    print('--', rs)
    for k in keys:
        rows = [j['configs'][k] for j in J if j['readset'] == rs and k in j['configs']]
        if not rows:
            continue
        N = sorted(r['N_truth_pair'] for r in rows)
        S = sorted(r['Z_pass1'] - r['N_truth_pair'] for r in rows)
        dep = sorted(r['median_depth'] for r in rows)
        hm = [h for r in rows for h in r['haps']]
        mos = [h for h in hm if h['truth_label'] == 'mosaic' and h['truth_ipd_absent']]
        exa = [h for h in hm if not h['truth_ipd_absent']]
        oth = [h for h in hm if h['truth_ipd_absent'] and h['truth_label'] != 'mosaic']
        def s(hs):
            return (f"n={len(hs)} int_exact={sum(h['interior_exact'] for h in hs)} byte={sum(h['byte_exact'] for h in hs)} "
                    f"mm0={sum(h['mismatch_to_truth'] == 0 for h in hs)} switched={sum(h['switches'] > 0 for h in hs)} "
                    f">3sw={sum(h['switches'] > 3 for h in hs)}")
        print(f"{k:14s} n={len(rows)} depth med {dep[len(dep)//2]:.0f} | N med {N[len(N)//2]:.0f} max {N[-1]:.0f} | "
              f"S med {S[len(S)//2]:.0f} | MOSAIC {s(mos)} | EXACT {s(exa)} | OTHER-ABSENT {s(oth)}")
