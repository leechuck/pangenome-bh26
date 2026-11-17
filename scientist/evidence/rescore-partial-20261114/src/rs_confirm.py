"""Partial-credit rows for the confirmatory test (33 donors, real reads). Read-only on the stage.
Usage: rs_confirm.py OUTDIR
Truth gene body / exons: sequence = TRUTH.json core; coding exons = Immuannot CDS + stop_codon features of the
copy's gene (mapped into core coordinates). Validated against IPD header exons for copies identical to an
IPD genomic allele (printed)."""
import collections
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rs_common as R

B = Path('/ibex/scratch/projects/c2014/rob/dogohla-benchmark')
C = B / 'confirmatory-test-20261113'
sys.path.insert(0, str(C / 'code'))
from endpoints import named_pair, unmatched_cores   # noqa: E402
from score_native import decode_slots               # noqa: E402

GENES = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
load = lambda p: json.loads(Path(p).read_text())
ALIAS = {}


def gtf_features(gtf):
    genes = {}
    feats = collections.defaultdict(list)
    with gzip.open(gtf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            x = line.rstrip('\n').split('\t')
            if len(x) < 9 or x[2] not in ('gene', 'CDS', 'stop_codon'):
                continue
            gid = re.search(r'gene_id "([^"]+)"', x[8])[1]
            if x[2] == 'gene':
                nm = re.search(r'gene_name "([^"]+)"', x[8])[1]
                genes[(x[0], int(x[3]), int(x[4]), nm)] = gid
            else:
                feats[gid].append((int(x[3]), int(x[4])))
    return genes, feats


def truth_exons(rec, genes, feats, gname):
    gid = genes.get((rec['contig'], rec['start_1based'], rec['end_inclusive'], 'HLA-' + gname))
    if gid is None or not feats.get(gid):
        return []
    s0, e0 = rec['start_1based'], rec['end_inclusive']
    iv = []
    for s, e in feats[gid]:
        if rec['strand'] == '+':
            iv.append((s - s0, e - s0 + 1))
        else:
            iv.append((e0 - e, e0 - s + 1))
    iv.sort()
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def decode_line(line):
    v = line.rstrip('\n').split('\t')
    copies = int(v[1])
    tokens = [v[2] if float(v[4]) > 0 else '.'] if copies in (1, 2) else ['.']
    tokens = [tokens[0], tokens[0] if copies == 1 else v[5] if copies == 2 and float(v[7]) > 0 else '.']
    return [[] if t == '.' else t.split(',') for t in tokens]


def main(outdir):
    out = Path(outdir)
    plan = load(C / 'code' / 'PLAN.json')
    arms = plan['arms']
    donors = [r['donor'] for r in plan['donors']]
    fold = {r['donor']: r['panel_fold'] for r in plan['donors']}
    pop = {r['donor']: r['population'] for r in plan['donors']}
    screen = load(C / 'receipts' / 'RELATEDNESS_SCREEN.json')
    excluded = sorted({f['a'] for f in screen['flags'] if f['b_role'] != 'test'})
    parent = {d: d for d in donors}

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for g in plan['clusters']:
        for d in g[1:]:
            parent[find(d)] = find(g[0])
    for f in screen['flags']:
        if f['b_role'] == 'test':
            parent[find(f['a'])] = find(f['b'])
    comp = collections.defaultdict(list)
    for d in donors:
        comp[find(d)].append(d)
    group_of = {d: sorted(v) for v in comp.values() for d in v}
    (out / 'groups_confirm.json').write_text(json.dumps(dict(groups=group_of, excluded_setP=excluded)))

    def refpath(arm, d):
        if arm == 'ipd':
            return C / 'refs/ipd.original.fa'
        if arm == 'genomic':
            return C / 'refs/genomic.fa'
        return C / 'refs' / f'fold{fold[d]}' / arm / 'reference.fa'

    # ---- pass 1: calls and names needed
    calls = {}
    need = collections.defaultdict(set)
    for d in donors:
        for arm in arms:
            t = C / 'runs' / d / arm / 'calls_genotype.tsv'
            rec = load(t.parent / 'COMPLETE.json')
            assert rec['status'] == 'complete' and rec['native_table_sha256'] == hashlib.sha256(t.read_bytes()).hexdigest()
            tab = {}
            for line in t.read_text().splitlines():
                if line and not line.startswith('#'):
                    g = line.split('\t')[0].removeprefix('HLA-')
                    if g in GENES:
                        tab[g] = (line, decode_line(line))
                        for s in tab[g][1]:
                            need[str(refpath(arm, d))].update(s)
            calls[(d, arm)] = tab
    refs = {}
    for p, names in need.items():
        refs[p] = R.read_fasta_needed(p, names)
        print('ref', Path(p).parent.name, Path(p).name, 'needed', len(names), 'found', len(refs[p]), flush=True)

    def lookup(p, name):
        d = refs[str(p)]
        for n in (name, 'HLA-' + name, name.removeprefix('HLA-')):
            if n in d:
                return R.Seq(n, *d[n])
        return None

    # ---- truth, rows
    rows, frozen = [], []
    gtf_check = collections.Counter()
    genomic_ref = {}
    for d in donors:
        T = load(C / 'truth' / d / 'TRUTH.json')
        genes, feats = gtf_features(C / 'truth' / d / 'immuannot.gtf.gz')
        for gname in GENES:
            if not T['genes'][gname]['eligible']:
                continue
            recs = [r for r in T['records'] if r['gene'] == gname]
            assert len(recs) == 2
            copies = []
            for r in recs:
                ex = truth_exons(r, genes, feats, gname)
                copies.append(R.Seq(r['core_id'], r['sequence'].upper(), ex))
                # validation against IPD header exons for copies identical to an IPD genomic allele
                if r['genomic_labels']:
                    lab = 'HLA-' + r['genomic_labels'][0].removeprefix('HLA-')
                    if lab not in genomic_ref:
                        genomic_ref[lab] = None
                    gtf_check['present_copies'] += 1
            for r, c in zip(recs, copies):
                if not c.exons:
                    gtf_check['truth_without_exons'] += 1
            copies_meta = [{k: r[k] for k in ['haplotype', 'genomic_labels', 'core_id', 'core_sha256']} for r in recs]
            for arm in arms:
                line_dec = calls[(d, arm)].get(gname)
                rp = refpath(arm, d)
                done = True
                if line_dec is None:
                    slots_labels = [[], []]
                else:
                    slots_labels = line_dec[1]
                slot_alleles, unres = [], []
                for labs in slots_labels:
                    al = [lookup(rp, n) for n in labs]
                    unres.append(sum(a is None for a in al))
                    slot_alleles.append([a for a in al if a is not None])
                res = R.assign(slot_alleles, copies)
                for ci, (rr, c, (p, rs)) in enumerate(zip(recs, copies, res)):
                    base = dict(set='confirm', readset='real', donor=d, gene=gname, copy=ci + 1,
                                haplotype=rr['haplotype'], ipd_absent=int(not rr['genomic_labels']),
                                inP=int(d not in excluded))
                    names = slots_labels[p] if p is not None else []
                    row = R.make_row(base, arm, c, rs, names, slot_alleles[p] if p is not None else [],
                                     unres[p] if p is not None else 0)
                    if arm == 'ipd':
                        # the ipd reference holds exon-only / truncated-intron sequences: gene-body values are not available
                        for k in ('call_gb_len', 'gb_dist', 'nc_dist', 'nc_negative', 'gb_dist_best', 'gb_exact', 'gb_exact_best', 'full_exact', 'full_exact_best'):
                            row[k] = ''
                        row['gb_defined'] = 0
                    rows.append(row)
                # frozen endpoints on this donor-gene-arm
                if line_dec is None:
                    pred = [dict(labels=[], label_unknown=True, cores=[], core_unknown=True)] * 2
                else:
                    ak = (arm, fold[d])
                    if ak not in ALIAS:
                        ALIAS[ak] = {} if arm in ('ipd', 'genomic') else load(C / plan['references'][arm][str(fold[d])]['aliases'])
                    aliases = ALIAS[ak]
                    pred = decode_slots(line_dec[0], aliases)
                npair = named_pair(pred, copies_meta)
                un = unmatched_cores(pred, copies_meta)
                frozen.append(dict(donor=d, gene=gname, arm=arm, eligible=int(npair['eligible']),
                                   called=int(npair['called']), correct=int(npair['correct']),
                                   unm_copies=un['eligible_copies'], unm_recovered=un['exact_copies_recovered']))
    cols = R.ROW_COLS + ['inP']
    R.write_rows(out / 'rows_confirm.tsv', rows, cols)
    with open(out / 'frozen_confirm.tsv', 'w') as f:
        f.write('donor\tgene\tarm\teligible\tcalled\tcorrect\tunm_copies\tunm_recovered\n')
        for r in frozen:
            f.write('\t'.join(str(r[k]) for k in ['donor', 'gene', 'arm', 'eligible', 'called', 'correct', 'unm_copies', 'unm_recovered']) + '\n')
    # ---- validation of truth exons against IPD header exons (identical sequences only)
    gen = R.read_fasta_needed(C / 'refs/genomic.fa', None)
    chk = collections.Counter()
    bad = []
    for d in donors:
        T = load(C / 'truth' / d / 'TRUTH.json')
        genes, feats = gtf_features(C / 'truth' / d / 'immuannot.gtf.gz')
        for r in T['records']:
            if r['gene'] not in GENES or not T['genes'][r['gene']]['eligible'] or not r['genomic_labels']:
                continue
            lab = 'HLA-' + r['genomic_labels'][0].removeprefix('HLA-')
            if lab not in gen:
                chk['label_not_in_genomic_fa'] += 1
                continue
            seq, ex = gen[lab]
            if seq != r['sequence'].upper():
                chk['core_differs_from_ipd_allele_sequence'] += 1
                continue
            tex = truth_exons(r, genes, feats, r['gene'])
            if tex == ex:
                chk['exons_equal'] += 1
            elif tex and tex[0][0] == ex[0][0] and tex[-1][1] == ex[-1][1]:
                chk['body_equal_exons_differ'] += 1
                bad.append((d, r['gene'], 'body_equal', len(tex), len(ex)))
            else:
                chk['body_differs'] += 1
                bad.append((d, r['gene'], 'body_differs', len(tex), len(ex), tex[:1] + tex[-1:], ex[:1] + ex[-1:]))
    print('TRUTH_EXON_CHECK', dict(chk))
    for b in bad[:20]:
        print('  mismatch', b)
    print('rows', len(rows), 'frozen', len(frozen), 'truth_without_exons', gtf_check['truth_without_exons'])


if __name__ == '__main__':
    main(sys.argv[1])
