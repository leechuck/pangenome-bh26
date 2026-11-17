"""DEVELOPMENT ONLY: is each patch (and each rejected patch) right? Compares the local version
chosen by the builder with the true copy's local version (hashes only).
  python drb1_eval_patches.py READSET[,READSET]"""
import collections
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L
import drb1_patch2 as PT
C, PR = L.C, L.PR
import panel_select as PS

S2 = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113')
G = 'DRB1'
SUB = sys.argv[2] if len(sys.argv) > 2 else 'C3'
records = C.load_records()
msa = C.load_msas(records, genes=(G,))[G]
shas = [C.sha256_str(s) for s in msa.seqs]
rows = PS.select_rows('C', [bool(n) for n in msa.ipd_names], shas, None)
names = [msa.ipd_names[r][0] for r in rows]
panel = PR.Panel(names, [msa.aligned[r] for r in rows])
rown = {n: i for i, n in enumerate(names)}
for r in rows:
    for n in msa.ipd_names[r]:
        rown.setdefault(n, rows.index(r))
out = ['readset\tdonor\tpass1\ttruth_class\ttruth_absent\tstatus\tkind\tcl\tcr\tlen_cur\tlen_new\tcur_total\tn_new\tfwd\trev\texp_h\tn_rows_new\tn_rows_cur\tverdict']
tot = collections.Counter()
for rs in sys.argv[1].split(','):
    assert rs in ('D-sim', 'D-real')
    for d in sorted((S2 / 'work' / rs).iterdir()):
        p = d / SUB / 'proposals.json'
        if not p.exists():
            continue
        j = C.load_json(p)
        tr = [r for r in records if r['donor'] == d.name and r['gene'] == G]
        tp = PR.Panel([r['source_id'] for r in tr], [msa.aligned[msa.row_of_source[r['source_id']]] for r in tr])
        cl_ = C.load_classes()
        for h in j['haplotypes']:
            if 'skipped' in h:
                continue
            a = panel.P[rown[h['pass1_allele']]]
            k = min(range(len(tr)), key=lambda i: L.mism(a, tp.P[i]))
            tc = tp.P[k]
            for status, items in (('applied', h['patches']), ('rejected', h['notes'])):
                for x in items:
                    ts = PT.row_string(tc, x['cl'], x['cr'])
                    tsha = C.sha256_str(ts)[:16] if ts is not None else None
                    v = 'new_is_truth' if tsha == x['new_sha'] else 'cur_is_truth' if tsha == x['cur_sha'] else 'neither'
                    tot[(rs, status, x['kind'], v)] += 1
                    out.append('\t'.join(map(str, (rs, d.name, h['pass1_allele'], cl_[tr[k]['source_id']]['ipd_label_L150'],
                                                   not tr[k]['genomic_labels'], status, x['kind'], x['cl'], x['cr'], x['len_cur'],
                                                   x['len_new'], x['cur_total'], x['n_new'], x['fwd_rev'][0], x['fwd_rev'][1],
                                                   x['exp_h'], x.get('n_rows_new'), x.get('n_rows_cur'), v))))
(S2 / 'results/dev_patch_verdicts.tsv').write_text('\n'.join(out) + '\n')
for k in sorted(tot):
    print(*k, tot[k], sep='\t')
