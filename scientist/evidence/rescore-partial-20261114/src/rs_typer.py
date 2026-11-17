"""Partial-credit rows for the mosaic-typer stage (all genes) and the DRB1 patch-builder stage (DRB1).
Usage: rs_typer.py MODE READSET ARMS OUTDIR      MODE = typer | drb1
Reads existing T1K runs/additions through the frozen stages' own loaders (read-only); no typing, no truth
is used except for scoring. Frozen per-copy flags are carried along for validation (columns fz_*)."""
import collections
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rs_common as R

OLD = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src')
NEWD = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113/src')
sys.path.insert(0, str(OLD))
sys.path.insert(0, str(NEWD))
import common as C   # noqa: E402
import score as SC   # noqa: E402

EXTRA = ['arm_complete', 'body_final_any', 'class_full', 'class_ipd', 'truth4', 'fz_recovered', 'fz_recovered_final',
         'fz_named4', 'fz_body_recovered', 'fz_final_mos']


class Idx(SC.NovelIndex):
    def __init__(self, maps):
        super().__init__(maps)
        self.seq = {}
        self.collisions = 0

    def add_fasta(self, path):
        if str(path) in self.cache:
            return
        self.cache[str(path)] = True
        for name, ints, seq in C.read_fasta_headers_and_seqs(path):
            if not SC.is_novel(name):
                continue
            g = C.gene_of_name(name)
            if g not in C.GENES:
                continue
            ex = C.header_exons(ints)
            cds = C.cds_of(seq, ex)
            cn, pn = self.maps[g]
            if name in self.seq and self.seq[name].seq != seq:
                self.collisions += 1
            self.cache[name] = (C.sha256_str(seq), set(cn.get(cds, ())), set(pn.get(C.translate(cds), ())))
            self.seq[name] = R.Seq(name, seq, ex)


def main(mode, readset, arms_s, outdir):
    arms = arms_s.split(',')
    out = Path(outdir)
    C.verify_inputs()
    if mode == 'drb1':
        import drb1_score as D
        genes = ('DRB1',)
        NEW = D.NEW
    else:
        D = None
        genes = C.GENES
        NEW = ()
    T = SC.Truth([readset])
    seq_of = {}
    for r in T.records:
        seq_of[(r['donor'], r['gene'], str(r['haplotype']))] = R.Seq(
            r['source_id'], r['core_sequence'].upper(), [tuple(e) for e in r['core_exons']])
    body_hash = {k: C.sha256_str(v.body) for k, v in seq_of.items()}
    nidx = Idx(T.maps)
    rows, avail = [], []

    def ipd_seq(label):
        g = C.gene_of_name(label)
        if g not in T.msas:
            return None
        m = T.msas[g]
        i = m.row_of_allele.get(label)
        if i is None:
            i = m.row_of_allele.get(label if label.startswith('HLA-') else 'HLA-' + label)
        if i is None:
            return None
        return R.Seq(label, m.seqs[i], m.exons(i))

    def labseq(label):
        if SC.is_novel(label):
            return nidx.seq.get(label)
        return ipd_seq(label)

    for (rs, d, g), copies in sorted(T.copies.items()):
        if g not in genes:
            continue
        copies = sorted(copies, key=lambda c: c['haplotype'])
        for arm in arms:
            if arm in ('B', 'D', 'D3', 'D4') and rs == 'D-sim':
                continue
            if arm in ('B', 'D') and rs in C.SYN_SETS:
                continue
            key = (rs, d, arm)
            if key not in nidx.cache:
                try:
                    t, fas = D.new_arm_table(rs, d, arm) if arm in NEW else SC.arm_table(rs, d, arm)
                except Exception as e:
                    t, fas = None, repr(e)
                if t is None:
                    nidx.cache[key] = None
                    avail.append((rs, d, arm, str(fas)))
                else:
                    for fa in fas:
                        if Path(fa).exists():
                            nidx.add_fasta(fa)
                    nidx.cache[key] = C.read_t1k(t)
                    avail.append((rs, d, arm, 'complete'))
            calls = nidx.cache[key]
            complete = calls is not None
            if complete and arm in ('C', 'D') and SC.ACCEPTANCE:
                calls = SC.apply_acceptance(rs, d, arm, g, calls)
            if complete and arm in NEW:
                calls = D.new_post(rs, d, arm, g, calls)
            slots = calls.get(g, [[], []]) if complete else [[], []]
            typed = [SC.typed(s) for s in slots]
            alleles, unres = [], []
            for tl in typed:
                al = [labseq(l) for l in tl]
                unres.append(sum(a is None for a in al))
                alleles.append([a for a in al if a is not None])
            truth = [seq_of[(d, g, str(c['haplotype']))] for c in copies]
            res = R.assign(alleles, truth)
            fz = SC.score_gene(copies, slots, nidx)
            gene_bodies = set()
            final_bodies = set()
            if mode == 'drb1':
                for s_ in typed:
                    for l in s_:
                        if SC.is_novel(l) and l in nidx.seq and nidx.seq[l].body is not None:
                            final_bodies.add(C.sha256_str(nidx.seq[l].body))
                for s in slots:
                    for l in s:
                        if SC.is_novel(l) and l in nidx.seq and nidx.seq[l].body is not None:
                            gene_bodies.add(C.sha256_str(nidx.seq[l].body))
            for ci, (c, t, (p, rr), f) in enumerate(zip(copies, truth, res, fz)):
                base = dict(set=mode, readset=rs, donor=d, gene=g, copy=ci + 1, haplotype=c['haplotype'],
                            ipd_absent=int(c['ipd_absent']))
                names = typed[p] if p is not None else []
                al = alleles[p] if p is not None else []
                nov = bool(names) and all(SC.is_novel(l) for l in names)
                row = R.make_row(base, arm, t, rr, names, al, unres[p] if p is not None else 0, novel_final=nov)
                row.update(arm_complete=int(complete), body_final_any=int(body_hash[(d, g, str(c['haplotype']))] in final_bodies), class_full=c['class_full'], class_ipd=c['class_ipd'],
                           truth4=c['truth4'] or '', fz_recovered=int(f['recovered']),
                           fz_recovered_final=int(f['recovered_final']), fz_named4=int(f['named4']),
                           fz_final_mos=int(f['final_mos']),
                           fz_body_recovered=int(f['recovered'] or body_hash[(d, g, str(c['haplotype']))] in gene_bodies))
                rows.append(row)
    R.write_rows(out / f'rows_{mode}_{readset}.tsv', rows, R.ROW_COLS + EXTRA)
    inc = collections.Counter((a[2]) for a in avail if a[3] != 'complete')
    print('rows', len(rows), 'incomplete arm-donors', dict(inc), 'novel name collisions', nidx.collisions)
    print('unresolved labels rows', sum(1 for r in rows if r['unresolved'] and int(r['unresolved']) > 0))


if __name__ == '__main__':
    main(*sys.argv[1:5])
