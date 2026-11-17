"""Spot-check 5 confirmatory copies: lengths and distances only. Recomputes distances from the sequences with an
edlib alignment path (CIGAR mismatch/indel count) and compares with the row-level file. Prints no sequence."""
import json, random, sys
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import edlib, rs_common as R, rs_confirm as RC
import rs_tables as T
res = Path(sys.argv[1])
rows = T.load_rows(res / 'rows_confirm.tsv')
rng = random.Random(1114)
def pick(f):
    c = [r for r in rows if f(r)]
    return rng.choice(c)
chosen = [pick(lambda r: r['arm'] == 'ipd' and r['missing'] == 0 and r['n_tied'] == 1 and r['cod_dist'] and r['cod_dist'] > 0),
          pick(lambda r: r['arm'] == 'genomic_core_hprc_asian' and r['ipd_absent'] == 1 and r['gb_exact'] == 1 and r['n_tied'] == 1),
          pick(lambda r: r['arm'] == 'genomic' and r['ipd_absent'] == 1 and r['missing'] == 0 and r['gb_dist'] and r['gb_dist'] > 0 and r['n_tied'] == 1),
          pick(lambda r: r['arm'] == 'genomic_core_nonasian' and r['ipd_absent'] == 0 and r['gb_dist'] is not None and r['gb_dist'] > 0 and r['n_tied'] == 1),
          pick(lambda r: r['arm'] == 'genomic' and r['ipd_absent'] == 0 and r['gb_exact'] == 1 and r['n_tied'] == 1)]
plan = json.loads((RC.C / 'code/PLAN.json').read_text())
fold = {x['donor']: x['panel_fold'] for x in plan['donors']}
def cig(a, b):
    r = edlib.align(a, b, mode='NW', task='path')
    c = r['cigar'] or ''
    n = 0; num = ''
    for ch in c:
        if ch.isdigit(): num += ch
        else:
            if ch in 'XID': n += int(num)
            num = ''
    return r['editDistance'], n
for r in chosen:
    d, g, arm = r['donor'], r['gene'], r['arm']
    T_ = json.loads((RC.C / 'truth' / d / 'TRUTH.json').read_text())
    genes, feats = RC.gtf_features(RC.C / 'truth' / d / 'immuannot.gtf.gz')
    rec = [x for x in T_['records'] if x['gene'] == g][r['copy'] - 1]
    t = R.Seq('t', rec['sequence'].upper(), RC.truth_exons(rec, genes, feats, g))
    ref = RC.C / ('refs/ipd.original.fa' if arm == 'ipd' else 'refs/genomic.fa' if arm == 'genomic' else f'refs/fold{fold[d]}/{arm}/reference.fa')
    name = r['call']
    al = R.read_fasta_needed(ref, {name})[name]
    a = R.Seq(name, al[0], al[1])
    out = dict(donor=d, gene=g, arm=arm, copy=r['copy'], call=name, ipd_absent=r['ipd_absent'],
               true_core_len=len(t.seq), true_gb_len=len(t.body), true_cod_len=len(t.cod), call_seq_len=len(a.seq),
               call_gb_len=len(a.body) if a.has_introns else 'NA(exon-only)', call_cod_len=len(a.cod),
               file_gb=r['gb_dist'], file_cod=r['cod_dist'])
    out['recomputed_cod'] = cig(a.cod, t.cod)
    out['recomputed_gb'] = cig(a.body, t.body) if a.has_introns and len(a.exons) == len(t.exons) else 'NA'
    print(out)
