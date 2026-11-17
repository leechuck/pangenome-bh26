"""Diagnostic for one disagreement (names, lengths, distances only)."""
import sys
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rs_common as R
sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src')
import common as C, score as SC
T = SC.Truth(['T-sim'])
for r in T.records:
    if r['donor'] == 'APR049' and r['gene'] == 'DRB1':
        t = R.Seq('t', r['core_sequence'].upper(), [tuple(e) for e in r['core_exons']])
        print('copy hap', r['haplotype'], 'labels', r['genomic_labels'], 'core_len', len(t.seq), 'gb_len', len(t.body), 'cod_len', len(t.cod), 'nexons', len(t.exons))
        m = T.msas['DRB1']
        for lab in r['genomic_labels'] + ['HLA-DRB1*16:09:01', 'HLA-DRB1*16:01:01:02']:
            i = m.row_of_allele.get(lab)
            if i is None:
                print(' label', lab, 'not in msa'); continue
            a = R.Seq(lab, m.seqs[i], m.exons(i))
            print(' label', lab, 'len', len(a.seq), 'gb_len', len(a.body), 'cod_len', len(a.cod), 'nexons', len(a.exons), 'same_row_as_truth', m.row_of_seq.get(t.seq) == i,
                  'd_gb', R.ed(a.body, t.body), 'd_cod', R.ed(a.cod, t.cod), 'd_full', R.ed(a.seq, t.seq), 'names_for_row', m.ipd_names[i][:6])
