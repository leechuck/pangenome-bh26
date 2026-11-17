"""Score the thin-data runs (stage thin-data-20261114).

  score_thin.py rows     -> results/rows.tsv (long: one row per copy x arm x depth x rule), depth_rows.tsv
Scoring reuses the frozen scorer (mosaic-typer-20261113/src/score.py: Truth, NovelIndex, slot_eval, copy_eval)
unchanged; only the slot-to-call rule differs (prior rules, post-processing of T1K ambiguity lists).

Rules
  P0s  frozen scorer: the whole ambiguity list; a field level is correct only if ALL listed alleles agree on it
  P0f  first listed allele (the brief's P0)
  PM   listed allele most frequent in the donor's own population group (truth labels of the OTHER-family donors)
  PO   same, frequencies from the two other population groups together
"""
import collections, json, re, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib as L
C = L.C
import score as S0          # frozen scorer (imported read-only)

W = L.TS / 'work'
RES = L.TS / 'results'
DEPTHS = [('full', 1), ('1/2', 2), ('1/4', 4), ('1/8', 8), ('1/16', 16)]
RULES = ('P0s', 'P0f', 'PM', 'PO')


def group_of_pop(pop, resource):
    if resource in ('Arab', 'Saudi'):
        return 'ME'
    code = re.search(r'\(([A-Z]{3})\)', pop)
    code = code.group(1) if code else pop
    if code in ('CDX', 'CHB', 'CHS', 'KHV', 'JPT'):
        return 'EAS'
    if code in ('BEB', 'GIH', 'ITU', 'PJL', 'STU'):
        return 'SAS'
    raise ValueError(pop)


def table_path(rs, donor, arm, n):
    """Calls table of an arm at depth 1/n (n=1 -> existing full-depth run, verified against its receipt)."""
    if n == 1:
        run = L.arm_run(rs, donor, arm)
        if run is None:
            return None
        rec = C.load_json(run / 'COMPLETE.json')
        t = run / 'calls_genotype.tsv'
        if rec['status'] != 'complete' or C.sha256_file(t) != rec['native_table_sha256']:
            return None
        return t
    d = W / 'runs' / rs / donor / arm / f'd{n}'
    if not (d / 'COMPLETE.json').exists():
        return None
    rec = C.load_json(d / 'COMPLETE.json')
    t = d / 'calls_genotype.tsv'
    if C.sha256_file(t) != rec['table_sha256']:
        return None
    return t


class Freq:
    """Population-group allele counts from truth labels, own family component excluded."""
    def __init__(self, T, pop_group, excl):
        self.pop_group, self.excl = pop_group, excl
        self.per_donor = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
        for r in T.records:
            d, g = r['donor'], r['gene']
            pd = self.per_donor[d][g]
            for k in (2, 3, 4):
                for nm in {C.normalized(l, k) for l in r['genomic_labels']} - {None}:
                    pd[(k, nm)] += 1
            pd[('sha', C.sha256_str(r['core_sequence']))] += 1
        self.cache = {}

    def table(self, donor, mode):
        key = (donor, mode)
        if key in self.cache:
            return self.cache[key]
        own = self.excl[donor]
        g0 = self.pop_group[donor]
        tot = collections.defaultdict(collections.Counter)
        for d, per in self.per_donor.items():
            if d in own or d not in self.pop_group:
                continue
            same = self.pop_group[d] == g0
            if (mode == 'PM') == same:
                for g, cnt in per.items():
                    tot[g].update(cnt)
        self.cache[key] = tot
        return tot


def label_sha(T, g, label, nidx):
    if S0.is_novel(label):
        n = nidx.get(label)
        return n[0] if n else None
    row = T.msas[g].row_of_allele.get(label)
    return C.sha256_str(T.msas[g].seqs[row]) if row is not None else None


def freq_of(label, k, g, tab, nidx):
    c = tab.get(g)
    if c is None:
        return 0
    if S0.is_novel(label):
        n = nidx.get(label)
        if not n:
            return 0
        if k == 4:
            return c.get(('sha', n[0]), 0)
        names = n[1] if k == 3 else n[2]
        return c.get((k, next(iter(names))), 0) if len(names) == 1 else 0
    nm = C.normalized(label, k)
    return c.get((k, nm), 0) if nm else 0


def choose(labels, g, tab, nidx):
    if not labels:
        return []
    for k in (4, 3, 2):
        f = [freq_of(l, k, g, tab, nidx) for l in labels]
        if max(f) > 0:
            return [labels[f.index(max(f))]]
    return [labels[0]]


def score_gene2(T, g, copies, slots, nidx):
    ev = [S0.slot_eval(s, nidx) for s in slots]
    best = None
    for perm in itertools_perms(len(ev), len(copies)):
        ms = [S0.copy_eval(c, ev[p]) for c, p in zip(copies, perm)]
        sc = sum(1000 * m['exact'] + 100 * m['named4'] + 10 * m['named3_decoded'] + m['named2_decoded'] for m in ms)
        if best is None or sc > best[0]:
            best = (sc, perm, ms)
    gene_shas = set().union(*(e['all_shas'] for e in ev)) if ev else set()
    out = []
    for c, p, m in zip(copies, best[1], best[2]):
        e = ev[p]
        shas = [label_sha(T, g, l, nidx) for l in e['typed']]
        seqok = bool(e['typed']) and all(s == c['sha'] for s in shas)
        out.append(dict(slot=p, n_listed=len(slots[p]), no_call=not slots[p],
                        chosen='|'.join(e['typed']) if e['typed'] else '.',
                        ok2=m['named2_decoded'], ok3=m['named3_decoded'], ok4=m['named4'],
                        seqok=seqok, recovered_any=c['sha'] in gene_shas))
    return out


def itertools_perms(n, k):
    import itertools
    return itertools.permutations(range(n), k)


def main():
    t0 = time.time()
    C.verify_inputs()
    rsets = L.readsets()
    T = S0.Truth(['D-real', 'T-real', 'T-sim'])
    folds = C.load_json(C.S / 'inputs/FOLDS.json')
    pop_group = {d: group_of_pop(v['population'], v['resource']) for d, v in folds['donors'].items()}
    excl = C.exclusion_map()
    fam = {d: ','.join(sorted(excl[d])) for d in excl}
    F = Freq(T, pop_group, excl)
    cols = ['readset', 'donor', 'group', 'family', 'depth', 'arm', 'rule', 'gene', 'haplotype', 'slot', 'n_listed',
            'no_call', 'chosen', 'ok2', 'ok3', 'ok4', 'has_truth2', 'has_truth3', 'present4', 'ipd_absent',
            'seqok', 'recovered_any', 'truth4']
    nidx_cache = {}
    missing = []
    fo = open(RES / 'rows.tsv', 'w')
    fo.write('\t'.join(cols) + '\n')
    fd = open(RES / 'depth_rows.tsv', 'w')
    fd.write('\t'.join(['readset', 'donor', 'group', 'depth', 'arm', 'pairs', 'candidate_fragments',
                        'assigned_fragments_analyzer', 'abundance_depth_median']) + '\n')
    nrows = 0
    for rs, donors in rsets.items():
        for d in donors:
            sub = W / 'reads' / rs / d / 'SUBSAMPLE.json'
            sj = C.load_json(sub) if sub.exists() else None
            for arm in L.ARMS_ALL:
                if L.arm_run(rs, d, arm) is None:
                    continue
                ref = L.arm_info(rs, d, arm)['ref']
                if str(ref) not in nidx_cache:
                    ni = S0.NovelIndex(T.maps)
                    if arm != 'A':
                        ni.add_fasta(ref)
                    nidx_cache[str(ref)] = ni
                nidx = nidx_cache[str(ref)]
                for lab, n in DEPTHS:
                    if n != 1 and n not in L.levels_for_arm(arm):
                        continue
                    t = table_path(rs, d, arm, n)
                    if t is None:
                        missing.append((rs, d, arm, lab))
                        continue
                    calls = C.read_t1k(t)
                    raw = C.read_t1k_raw(t)
                    ab = []
                    for g in C.GENES:
                        r = raw.get(g)
                        if r and r['copies']:
                            ab.append(sum(max(s['abundance'], 0) for s in r['slots'][:r['copies']]) / 2.0)
                    ab = sorted(ab)
                    med = (ab[len(ab) // 2] if len(ab) % 2 else (ab[len(ab) // 2 - 1] + ab[len(ab) // 2]) / 2) if ab else ''
                    if n == 1:
                        lg = L.log_stats(L.arm_run(rs, d, arm) / 'run.log')
                        pairs = sj['full']['pairs'] if sj else ''
                    else:
                        rec = C.load_json(W / 'runs' / rs / d / arm / f'd{n}' / 'COMPLETE.json')
                        lg = {k: rec.get(k, '') for k in ('candidate_fragments', 'assigned_fragments_analyzer')}
                        pairs = sj[str(n)]['pairs'] if sj else ''
                    fd.write('\t'.join(map(str, [rs, d, pop_group[d], lab, arm, pairs,
                                                  lg.get('candidate_fragments', ''),
                                                  lg.get('assigned_fragments_analyzer', ''), med])) + '\n')
                    tabs = {'PM': F.table(d, 'PM'), 'PO': F.table(d, 'PO')}
                    for g in C.GENES:
                        copies = sorted(T.copies.get((rs, d, g), []), key=lambda c: c['haplotype'])
                        if not copies:
                            continue
                        slots0 = calls.get(g, [[], []])
                        ruleslots = {'P0s': slots0,
                                     'P0f': [s[:1] for s in slots0],
                                     'PM': [choose(s, g, tabs['PM'], nidx) for s in slots0],
                                     'PO': [choose(s, g, tabs['PO'], nidx) for s in slots0]}
                        for rule in RULES:
                            res = score_gene2(T, g, copies, ruleslots[rule], nidx)
                            for c, r in zip(copies, res):
                                # ambiguity and no-call refer to the T1K output, not the chosen allele
                                r['n_listed'] = len(slots0[r['slot']])
                                r['no_call'] = not slots0[r['slot']]
                                v = dict(readset=rs, donor=d, group=pop_group[d], family=fam[d], depth=lab,
                                         arm=arm, rule=rule, gene=g, haplotype=c['haplotype'], slot=r['slot'],
                                         n_listed=r['n_listed'], no_call=int(r['no_call']), chosen=r['chosen'],
                                         ok2=int(r['ok2']), ok3=int(r['ok3']), ok4=int(r['ok4']),
                                         has_truth2=int(bool(c['truth2'])), has_truth3=int(bool(c['truth3'])),
                                         present4=int((not c['ipd_absent']) and bool(c['truth4'])),
                                         ipd_absent=int(c['ipd_absent']), seqok=int(r['seqok']),
                                         recovered_any=int(r['recovered_any']), truth4=c['truth4'] or '')
                                fo.write('\t'.join(str(v[k]) for k in cols) + '\n')
                                nrows += 1
    fo.close(); fd.close()
    C.save_json(RES / 'rows_meta.json', dict(rows=nrows, missing_tables=missing, seconds=time.time() - t0,
                                             code_sha256={p.name: C.sha256_file(p) for p in Path(__file__).parent.glob('*.py')}))
    print('rows', nrows, 'missing', len(missing))


if __name__ == '__main__':
    main()
