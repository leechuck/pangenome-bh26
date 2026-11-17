"""Tables from results/rows.tsv: cells with family-clustered bootstrap intervals, contrasts, per-gene, depth."""
import collections, csv, itertools, json, statistics, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib as L
import score as S0
C = L.C
RES = L.TS / 'results'
DEPTHS = ['full', '1/2', '1/4', '1/8', '1/16']
RULES = ['P0s', 'P0f', 'PM', 'PO']
GROUPS = {'ME': {'ME'}, 'EAS': {'EAS'}, 'SAS': {'SAS'}, 'REAL': {'EAS', 'SAS'}}
METRICS = ['acc2', 'acc3', 'acc4', 'nocall', 'ambig', 'listlen', 'seqok_all', 'recov_absent']


def mvals(r):
    n = int(r['n_listed'])
    return {'acc2': (int(r['ok2']), int(r['has_truth2'])), 'acc3': (int(r['ok3']), int(r['has_truth3'])),
            'acc4': (int(r['ok4']), int(r['present4'])), 'nocall': (int(r['no_call']), 1),
            'ambig': (int(n >= 2), 1), 'listlen': (n, int(n >= 2)),
            'seqok_all': (int(r['seqok']), 1), 'recov_absent': (int(r['seqok']) * int(r['ipd_absent']), int(r['ipd_absent']))}


def boot(clusters):
    """clusters: list of (num, den). Percentile interval as in the frozen scorer (4000, seed 20261020)."""
    num = np.array([c[0] for c in clusters], float); den = np.array([c[1] for c in clusters], float)
    if den.sum() == 0:
        return None, None, None
    rng = np.random.default_rng(20261020)
    idx = rng.integers(0, len(num), (4000, len(num)))
    d = den[idx].sum(1)
    ok = d > 0
    v = np.sort(num[idx].sum(1)[ok] / d[ok])
    return 100 * num.sum() / den.sum(), 100 * v[int(0.025 * len(v))], 100 * v[int(0.975 * len(v)) - 1]


def load():
    rows = []
    with open(RES / 'rows.tsv') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            r['rt'] = 'sim' if r['readset'] == 'T-sim' else 'real'
            r['m'] = mvals(r)
            rows.append(r)
    return rows


def cells(rows):
    by = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    gene = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    for r in rows:
        for gname, gs in GROUPS.items():
            if r['group'] not in gs or (gname == 'REAL' and r['rt'] != 'real') or (gname != 'REAL' and False):
                continue
            if (gname == 'ME') != (r['rt'] == 'sim'):
                continue
            for mt in METRICS:
                if mt == 'recov_absent' and r['arm'] == 'A':
                    continue
                num, den = r['m'][mt]
                c = by[(r['rt'], gname, r['depth'], r['arm'], r['rule'], mt)][r['family']]
                c[0] += num; c[1] += den
            if r['rule'] in ('P0s', 'P0f', 'PM'):
                num, den = r['m']['acc2']
                c = gene[(r['rt'], gname, r['depth'], r['arm'], r['rule'], r['gene'])][r['family']]
                c[0] += num; c[1] += den
    return by, gene


def write_cells(by, gene):
    with open(RES / 'figure_table.tsv', 'w') as f:
        f.write('\t'.join(['read_type', 'group', 'depth', 'arm', 'rule', 'metric', 'n', 'value_pct', 'lo95', 'hi95',
                           'family_clusters']) + '\n')
        for k in sorted(by):
            cl = [tuple(v) for v in by[k].values() if v[1] > 0]
            if not cl:
                continue
            if k[5] == 'listlen':   # mean list length among ambiguous slots (not a percentage)
                pass
            v, lo, hi = boot(cl)
            scale = 0.01 if k[5] == 'listlen' else 1.0
            f.write('\t'.join(map(str, list(k) + [sum(c[1] for c in cl), round(v * scale, 3), round(lo * scale, 3),
                                                   round(hi * scale, 3), len(cl)])) + '\n')
    with open(RES / 'pergene_acc2.tsv', 'w') as f:
        f.write('\t'.join(['read_type', 'group', 'depth', 'arm', 'rule', 'gene', 'n', 'acc2_pct', 'lo95', 'hi95']) + '\n')
        for k in sorted(gene):
            cl = [tuple(v) for v in gene[k].values() if v[1] > 0]
            if not cl:
                continue
            v, lo, hi = boot(cl)
            f.write('\t'.join(map(str, list(k) + [sum(c[1] for c in cl), round(v, 2), round(lo, 2), round(hi, 2)])) + '\n')


def depth_summary():
    rows = list(csv.DictReader(open(RES / 'depth_rows.tsv'), delimiter='\t'))
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        if r['arm'] != 'A':
            continue
        k = (('sim' if r['readset'] == 'T-sim' else 'real'), r['group'], r['depth'])
        for fld in ('pairs', 'candidate_fragments', 'assigned_fragments_analyzer', 'abundance_depth_median'):
            if r[fld] != '':
                out[k][fld].append(float(r[fld]))
    with open(RES / 'depth_summary.tsv', 'w') as f:
        f.write('\t'.join(['read_type', 'group', 'depth', 'donors', 'median_pairs', 'median_candidate_fragments',
                           'median_assigned_fragments', 'median_T1K_abundance_depth']) + '\n')
        for k in sorted(out, key=lambda k: (k[0], k[1], DEPTHS.index(k[2]))):
            d = out[k]
            med = lambda x: round(statistics.median(x), 2) if x else ''
            f.write('\t'.join(map(str, list(k) + [len(d['pairs']), med(d['pairs']), med(d['candidate_fragments']),
                                                  med(d['assigned_fragments_analyzer']),
                                                  med(d['abundance_depth_median'])])) + '\n')


# ----------------------------------------------------------------------------- contrasts
G = {}


def init(idx, group_of):
    G['idx'] = idx; G['group_of'] = group_of


def one_contrast(job):
    kind, rt, gname, depth, spec, mt = job
    idx = G['idx']
    def vec(arm, rule, dep):
        return idx.get((rt, arm, dep, rule), {})
    def ok(r):
        return r['group'] in GROUPS[gname]
    def val(r):
        n, d = r['m'][mt]
        return (n, d)
    if kind == 'diff':
        (a1, r1), (a0, r0) = spec
        A = vec(a1, r1, depth); B = vec(a0, r0, depth)
        pairs = []
        for k, x in A.items():
            y = B.get(k)
            if y is None or not ok(x):
                continue
            (n1, d1), (n0, d0) = val(x), val(y)
            if d1 == 0 or d0 == 0:
                continue
            pairs.append(dict(donor=x['donor'], difference=n1 - n0))
    else:   # interaction: (spec at `depth`) minus (spec at full)
        (a1, r1), (a0, r0) = spec
        out = {}
        for dep in (depth, 'full'):
            A = vec(a1, r1, dep); B = vec(a0, r0, dep)
            for k, x in A.items():
                y = B.get(k)
                if y is None or not ok(x):
                    continue
                (n1, d1), (n0, d0) = val(x), val(y)
                if d1 == 0 or d0 == 0:
                    continue
                out.setdefault(k, {})[dep] = (x['donor'], n1 - n0)
        pairs = [dict(donor=v[depth][0], difference=v[depth][1] - v['full'][1]) for v in out.values()
                 if depth in v and 'full' in v]
    u = S0.uncertainty(pairs, G['group_of'])
    ci = u['bootstrap95_pp']
    return [kind, rt, gname, depth, '%s/%s - %s/%s' % (spec[0] + spec[1]), mt, len(pairs), u['family_clusters'],
            None if u['difference_pp'] is None else round(u['difference_pp'], 3),
            None if ci is None else round(ci[0], 3), None if ci is None else round(ci[1], 3),
            u['sign_flip_p'], sum(p['difference'] > 0 for p in pairs), sum(p['difference'] < 0 for p in pairs)]


def main():
    rows = load()
    by, gene = cells(rows)
    write_cells(by, gene)
    depth_summary()
    idx = collections.defaultdict(dict)
    for r in rows:
        idx[(r['rt'], r['arm'], r['depth'], r['rule'])][(r['donor'], r['gene'], r['haplotype'])] = r
    arms_present = collections.defaultdict(set)
    for (rt, arm, dep, rule) in idx:
        arms_present[rt].add((arm, dep))
    jobs = []
    cm = ['acc2', 'acc3', 'acc4', 'seqok_all', 'recov_absent', 'nocall', 'ambig']
    rm = ['acc2', 'acc3', 'acc4', 'seqok_all']
    for rt, gnames in (('real', ['EAS', 'SAS', 'REAL']), ('sim', ['ME'])):
        arms = sorted({a for a, _ in arms_present[rt]})
        for gname in gnames:
            for dep in DEPTHS:
                avail = {a for a, d in arms_present[rt] if d == dep}
                for rule in ('P0f', 'PM', 'P0s'):
                    for new in ('B', 'N', 'M'):
                        if new in avail and 'A' in avail:
                            for mt in cm:
                                jobs.append(('diff', rt, gname, dep, ((new, rule), ('A', rule)), mt))
                for arm in sorted(avail):
                    for hi, lo in (('PM', 'P0f'), ('PM', 'PO'), ('PM', 'P0s'), ('P0f', 'P0s')):
                        for mt in rm:
                            jobs.append(('diff', rt, gname, dep, ((arm, hi), (arm, lo)), mt))
                if dep != 'full':
                    for rule in ('P0f', 'PM'):
                        for new in ('B', 'N', 'M'):
                            if new in avail and 'A' in avail and (new, 'full') in arms_present[rt]:
                                for mt in rm + ['recov_absent']:
                                    jobs.append(('inter', rt, gname, dep, ((new, rule), ('A', rule)), mt))
                    for arm in sorted(avail):
                        for hi, lo in (('PM', 'P0f'), ('PM', 'PO'), ('PM', 'P0s')):
                            for mt in rm:
                                jobs.append(('inter', rt, gname, dep, ((arm, hi), (arm, lo)), mt))
    group_of = {r['donor']: tuple(sorted(r['family'].split(','))) for r in rows}
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4, initializer=init, initargs=(idx, group_of)) as p:
        res = p.map(one_contrast, jobs, chunksize=4)
    with open(RES / 'contrasts.tsv', 'w') as f:
        f.write('\t'.join(['kind', 'read_type', 'group', 'depth', 'contrast', 'metric', 'n_pairs', 'family_clusters',
                           'diff_pp', 'lo95_pp', 'hi95_pp', 'sign_flip_p', 'gains', 'losses']) + '\n')
        for r in res:
            f.write('\t'.join('' if x is None else str(x) for x in r) + '\n')
    print('contrasts', len(res))


if __name__ == '__main__':
    main()
