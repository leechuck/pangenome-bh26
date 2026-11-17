"""DEVELOPMENT ONLY: is each patch (and each rejected patch) right? Generalisation of
drb1_eval_patches.py to every gene and both panel arms: compares the local version chosen by the
builder with the true copy's local version (hashes only), the true copy being the donor's copy
nearest to the pass-1 allele (MSA row stored in proposals.json).
  python ag_eval_patches.py TAG READSET[,READSET] ARM"""
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L
import drb1_patch2 as PT
C, PR = L.C, L.PR

S2 = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114')
tag, readsets, arm = sys.argv[1], sys.argv[2].split(','), sys.argv[3]
assert not {'T-sim', 'T-real'} & set(readsets)
records = C.load_records()
msas = C.load_msas(records)
cl_ = C.load_classes()
out = ['readset\tdonor\tgene\tpass1\ttruth_class\ttruth_absent\tstatus\tkind\tcl\tcr\tlen_cur\tlen_new\tcur_total\tn_new\tfwd\trev\texp_h\tn_rows_new\tn_rows_cur\tverdict']
tot = collections.Counter()
for rs in readsets:
    for d in sorted((S2 / 'work' / tag / rs).iterdir()):
        for g in C.GENES:
            p = d / arm / g / 'proposals.json'
            if not p.exists():
                continue
            j = C.load_json(p)
            msa = msas[g]
            tr = [r for r in records if r['donor'] == d.name and r['gene'] == g]
            if not tr:
                continue
            for h in j['haplotypes']:
                if 'skipped' in h:
                    continue
                prow = h['pass1_msa_row']
                tp = PR.Panel(['p1'] + [r['source_id'] for r in tr],
                              [msa.aligned[prow]] + [msa.aligned[msa.row_of_source[r['source_id']]] for r in tr])
                k = min(range(len(tr)), key=lambda i: L.mism(tp.P[0], tp.P[i + 1]))
                tc = tp.P[k + 1]
                for status, items in (('applied', h['patches']), ('rejected', h['notes'])):
                    for x in items:
                        ts = PT.row_string(tc, x['cl'], x['cr'])
                        tsha = C.sha256_str(ts)[:16] if ts is not None else None
                        v = 'new_is_truth' if tsha == x['new_sha'] else 'cur_is_truth' if tsha == x['cur_sha'] else 'neither'
                        tot[(rs, g, status, x['kind'], v)] += 1
                        out.append('\t'.join(map(str, (rs, d.name, g, h['pass1_allele'], cl_[tr[k]['source_id']]['ipd_label_L150'],
                                                       not tr[k]['genomic_labels'], status, x['kind'], x['cl'], x['cr'], x['len_cur'],
                                                       x['len_new'], x['cur_total'], x['n_new'], x['fwd_rev'][0], x['fwd_rev'][1],
                                                       x['exp_h'], x.get('n_rows_new'), x.get('n_rows_cur'), v))))
(S2 / 'results' / f'patch_verdicts_{tag}_{arm}.tsv').write_text('\n'.join(out) + '\n')
for k in sorted(tot):
    print(*k, tot[k], sep='\t')
