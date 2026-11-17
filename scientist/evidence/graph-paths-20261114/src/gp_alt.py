"""DEVELOPMENT ONLY (truth read for scoring): alternative path picking on the same clean graphs.
Giraffe maps the gene's reads (input ii) to the graph; each alignment adds 1 to every graph edge it walks.
Two routes P1, P2 through the DAG are chosen to make the expected edge coverage h*([e in P1]+[e in P2]) match
the observed coverage on branch edges (edges not walked by every stored path). Alternating longest-path DP
over the DAG (any route allowed, so new combinations of stored pieces are reachable), started from the arm-A
pair; edges with no read on either end node are uninformative and ties keep the current route.
  python gp_alt.py READSET DONOR        -> work/<rs>/<donor>/alt_rows.tsv (distances and flags only)"""
import collections
import json
import statistics
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gp_lib as G      # noqa: E402
import gp_donor as GD   # noqa: E402
C = G.C
ROUNDS = 6
EPS = 1e-3


def load_graph(gdir):
    seq, walks, out, inn = {}, {}, collections.defaultdict(set), collections.defaultdict(set)
    for line in open(gdir / 'graph.gfa'):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'S':
            seq[int(f[1])] = f[2]
        elif f[0] == 'L':
            a, b = int(f[1]), int(f[3])
            out[a].add(b)
            inn[b].add(a)
        elif f[0] == 'W' and f[1] != 'REF':
            walks[f[1]] = [n for n, _ in G.parse_walk(f[6])]
    return seq, walks, out, inn


def best_route(seq, out, starts, ends, gain, bonus):
    """Longest path in the DAG (node ids are topological) from a start node to an end node."""
    NEG = float('-inf')
    score, prev = {}, {}
    for n in sorted(seq):
        best = 0.0 if n in starts else NEG
        bp = None
        for p in inn_of[n]:
            if score.get(p, NEG) == NEG:
                continue
            v = score[p] + gain.get((p, n), 0.0) + bonus.get((p, n), 0.0)
            if v > best:
                best, bp = v, p
        score[n], prev[n] = best, bp
    e = max((n for n in ends if score.get(n, NEG) > NEG), key=lambda n: score[n])
    path = [e]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    return path[::-1], score[e]


inn_of = {}


def main():
    global inn_of
    rs, donor = sys.argv[1], sys.argv[2]
    out = G.S / 'work' / rs / donor
    assert (out / 'COMPLETE.json').exists()
    info = C.load_json(out / 'coverage.json')
    records = C.load_records()
    classes = C.load_classes()
    msas = C.load_msas(records, genes=G.GENES3)
    fold = info['fold']
    frag = (['--fragment-mean', 350, '--fragment-stdev', 50] if rs == 'D-sim' else
            ['--fragment-mean', info['fragment_estimate']['median'], '--fragment-stdev', max(10, info['fragment_estimate']['robust_sd'])])
    rows = []
    for g in G.GENES3:
        msa = msas[g]
        truth = []
        for r in sorted((r for r in records if r['donor'] == donor and r['gene'] == g), key=lambda r: str(r['haplotype'])):
            core, ex = r['core_sequence'].upper(), r['core_exons']
            bc = msa.base_cols(msa.row_of_source[r['source_id']])
            truth.append(dict(hap=r['haplotype'], body=core[ex[0][0]:ex[-1][1]], c0=bc[ex[0][0]], c1=bc[ex[-1][1] - 1],
                              absent=not r['genomic_labels'], cls=classes[r['source_id']]['ipd_label_L150']))
        called = [msa.row_of_allele[n] for n in info['called'][g] if n in msa.row_of_allele]
        if len(called) == 1:
            called = called * 2
        for variant in G.VARIANTS:
            gdir = G.graph_dir(g, variant, fold)
            wd = out / g / variant / 'alt'
            wd.mkdir(parents=True, exist_ok=True)
            t0 = time.time()
            base = dict(readset=rs, donor=donor, gene=g, variant=variant)
            try:
                seq, walks, outE, inn = load_graph(gdir)
                inn_of = inn
                node_seq, node_cols = G.load_nodes(gdir)
                G.sh([G.VG, 'giraffe', '-Z', gdir / 'graph.gbz', '-d', gdir / 'graph.dist', '-m', gdir / 'graph.min',
                      '-z', gdir / 'graph.zipcodes', '-f', out / g / 'ii_1.fq', '-f', out / g / 'ii_2.fq', '-t', 4] + frag,
                     wd, wd / 'mapping.gam', timeout=GD.GIRAFFE_TIMEOUT)
                G.sh([G.VG, 'view', '-aj', 'mapping.gam'], wd, wd / 'aln.jsonl', timeout=GD.GIRAFFE_TIMEOUT)
                ecov, ncov = collections.Counter(), collections.Counter()
                nal = 0
                for line in open(wd / 'aln.jsonl'):
                    a = json.loads(line)
                    ns = [(int(m['position']['node_id']), bool(m['position'].get('is_reverse'))) for m in a.get('path', {}).get('mapping', [])]
                    if not ns:
                        continue
                    nal += 1
                    if ns[0][1]:
                        ns = [(n, not r) for n, r in ns[::-1]]
                    for n, _ in ns:
                        ncov[n] += 1
                    for (x, _), (y, _) in zip(ns, ns[1:]):
                        if x != y:
                            ecov[(min(x, y), max(x, y))] += 1
                for f in ('mapping.gam', 'aln.jsonl'):
                    (wd / f).unlink()
                use = collections.Counter()
                for w in walks.values():
                    for e in set(zip(w, w[1:])):
                        use[e] += 1
                npaths = len(walks)
                allE = {(a, b) for a in outE for b in outE[a]}
                branch = {e for e in allE if use[e] < npaths}
                # haploid edge depth, truth-free: half the median coverage of edges on both arm-A paths
                wa, wb = walks[G.path_name(msa, called[0])], walks[G.path_name(msa, called[1])]
                both = set(zip(wa, wa[1:])) & set(zip(wb, wb[1:]))
                h = max(1.0, statistics.median([ecov[e] for e in both]) / 2) if both else 1.0
                starts = {w[0] for w in walks.values()}
                ends = {w[-1] for w in walks.values()}
                names = {C.load_json(gdir / 'STATS.json')['ref_row']: None}
                P = [walks[G.path_name(msa, called[0])], walks[G.path_name(msa, called[1])]]
                for it in range(ROUNDS):
                    changed = False
                    for k in (0, 1):
                        other = set(zip(P[1 - k], P[1 - k][1:]))
                        gain = {}
                        for e in branch:
                            if ncov[e[0]] == 0 and ncov[e[1]] == 0:
                                continue
                            pw = h * (e in other)
                            gain[e] = abs(ecov[e] - pw) - abs(ecov[e] - pw - h)
                        bonus = {e: EPS for e in zip(P[k], P[k][1:])}
                        newp, _ = best_route(seq, outE, starts, ends, gain, bonus)
                        if newp != P[k]:
                            P[k], changed = newp, True
                    if not changed:
                        break
                haps = [GD.hap_bases([(n, '+') for n in p], node_seq, node_cols) for p in P]
                stored = {tuple(w) for w in walks.values()}
                B = [[GD.body_of(hh, t['c0'], t['c1']) for hh in haps] for t in truth]
                D = [[GD.ed(truth[i]['body'], b) for b in B[i]] for i in range(2)]
                pick = (0, 1) if D[0][0] + D[1][1] <= D[0][1] + D[1][0] else (1, 0)
                ipd_rows = [r for r in range(len(msa.uids)) if msa.ipd_names[r]]
                for i, t in enumerate(truth):
                    pb = B[i][pick[i]]
                    ipdb = {msa.aligned[r][t['c0']:t['c1'] + 1].replace('-', '') for r in ipd_rows}
                    gr = G.panel_rows(msa, variant, fold)
                    stb = {msa.aligned[r][t['c0']:t['c1'] + 1].replace('-', '') for r in gr}
                    rows.append(dict(base, copy=t['hap'], ipd_absent=t['absent'], cls=t['cls'], body_len=len(t['body']),
                                     dist=D[i][pick[i]], exact=D[i][pick[i]] == 0, picked_named=pb in ipdb,
                                     picked_stored=pb in stb, route_is_stored_path=tuple(P[pick[i]]) in stored,
                                     h=h, alignments=nal, rounds=it + 1, seconds=round(time.time() - t0, 1), error=''))
            except Exception as e:
                for t in truth:
                    rows.append(dict(base, copy=t['hap'], ipd_absent=t['absent'], cls=t['cls'], body_len=len(t['body']),
                                     dist=len(t['body']), exact=False, picked_named=False, picked_stored=False,
                                     route_is_stored_path=False, h='', alignments='', rounds='',
                                     seconds=round(time.time() - t0, 1), error=repr(e)[:200].replace('\t', ' ')))
    cols = list(rows[0])
    with open(out / 'alt_rows.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    C.save_json(out / 'ALT_COMPLETE.json', dict(status='complete', rows=len(rows), errors=sum(bool(r['error']) for r in rows)))
    print('done', rs, donor, len(rows))


if __name__ == '__main__':
    main()
