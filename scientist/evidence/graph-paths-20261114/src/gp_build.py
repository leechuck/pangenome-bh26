"""Build one acyclic graph straight from the per-gene alignment, verify it, and index it.
  python gp_build.py GENE VARIANT FOLD        (VARIANT ipd: FOLD ignored, one graph for all folds)
Writes work/graphs/<gene>/<variant[_foldK]>/ : graph.gfa, node_cols.json, graph.gbz, graph.dist, graph.ri,
graph.min, graph.zipcodes, hapl_<L>.hapl, STATS.json, COMPLETE.json | FAILURE.json.
No truth of development or test donors is read (cohort rows are chosen by FOLDS.json training donors)."""
import collections
import json
import shutil
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gp_lib as G   # noqa: E402
C = G.C


def encoding(msa, node_seq, walks, rows):
    """DRB1 named loci: how the paths store them (drb1_graphtest.enc logic on this graph)."""
    res = {}
    pseq = {r: ''.join(node_seq[n] for n in walks[r]) for r in rows}
    for name, (c0, c1) in G.GT.NAMED.items():
        la, ra = G.GT.anchors(msa, c0, c1)
        strings, walk_c, lens = collections.Counter(), collections.Counter(), set()
        nodes_used, one_base = set(), set()
        for r in rows:
            iv = G.GT.locate(pseq[r], la, ra)
            if iv is None:
                continue
            x = pseq[r][iv[0]:iv[1]]
            strings[x] += 1
            lens.add(len(x))
            w, pos = [], 0
            for n in walks[r]:
                e = pos + len(node_seq[n])
                if e > iv[0] and pos < iv[1]:
                    w.append(n)
                pos = e
            walk_c[tuple(w)] += 1
            nodes_used |= set(w)
        one_base = {n for n in nodes_used if len(node_seq[n]) == 1}
        # column-local branching: distinct (column-span) routes; strings vs walks equal means clean branches
        res[name] = dict(paths_located=sum(strings.values()), distinct_strings=len(strings),
                         distinct_lengths=len(lens), min_len=min(lens) if lens else None,
                         max_len=max(lens) if lens else None, distinct_walks=len(walk_c),
                         walks_repeating_a_node=sum(len(set(w)) != len(w) for w in walk_c),
                         nodes_used=len(nodes_used), one_base_nodes=len(one_base))
    return res


def acyclic(node_seq, edges):
    indeg = collections.Counter(b for _, b in edges)
    out = collections.defaultdict(list)
    for a, b in edges:
        out[a].append(b)
    q = [n for n in node_seq if indeg[n] == 0]
    seen = 0
    while q:
        n = q.pop()
        seen += 1
        for b in out[n]:
            indeg[b] -= 1
            if indeg[b] == 0:
                q.append(b)
    return seen == len(node_seq)


def main():
    gene, variant, fold = sys.argv[1], sys.argv[2], int(sys.argv[3])
    out = G.graph_dir(gene, variant, fold)
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    t0 = time.time()
    st = dict(gene=gene, variant=variant, fold=None if variant == 'ipd' else fold, vg=G.VG, times={})
    try:
        records = C.load_records()
        msa = C.load_msas(records, genes=(gene,))[gene]
        rows = G.panel_rows(msa, variant, fold)
        if variant == 'ipd_cohort':
            # leakage guard: no test donor of this fold (or its family component) carries a cohort row alone
            excl_map = C.exclusion_map()
            test = set(G.folds()[fold]['test'])
            train = set(G.folds()[fold]['train'])
            for d in test:
                if d in excl_map:
                    assert not (excl_map[d] & train), f'component of {d} in training set'
            row_is_ipd = [bool(n) for n in msa.ipd_names]
            for d in test:
                if d not in excl_map:
                    continue
                excl = excl_map[d]
                comp_absent = {msa.row_of_source[r['source_id']] for r in records
                               if r['gene'] == gene and r['donor'] in excl and not r['genomic_labels']}
                PS = G.PS
                PS.check_no_leakage(rows, row_is_ipd, msa.cohort_donors, excl, d, comp_absent, 'D')
        ref_row = rows[0]
        st.update(rows=len(rows), ipd_rows=sum(bool(msa.ipd_names[r]) for r in rows),
                  cohort_rows=sum(not msa.ipd_names[r] for r in rows),
                  ipd_names=sum(len(msa.ipd_names[r]) for r in rows),
                  cohort_copies=sum(len(msa.cohort_sources[r]) for r in rows if not msa.ipd_names[r]),
                  msa_columns=msa.ncol, ref_row=ref_row)
        t = time.time()
        node_seq, node_cols, walks, edges = G.build_gfa(msa, rows, out / 'graph.gfa', ref_row)
        st['times']['build_gfa'] = round(time.time() - t, 2)
        C.save_json(out / 'node_cols.json', {str(k): v for k, v in node_cols.items()})
        bad = [r for r in rows if ''.join(node_seq[n] for n in walks[r]) != msa.seqs[r]]
        assert not bad, f'{len(bad)} paths do not spell their sequence'
        assert acyclic(node_seq, edges), 'cycle'
        lens = collections.Counter(len(s) for s in node_seq.values())
        st.update(nodes=len(node_seq), edges=len(edges), paths=len(rows), graph_bases=sum(map(len, node_seq.values())),
                  one_base_nodes=lens[1], acyclic=True, paths_verified_python=len(rows))
        if gene == 'DRB1':
            st['encoding'] = encoding(msa, node_seq, walks, rows)
        st['body_cols'] = None
        # indexes
        T = st['times']
        T['gbwt_G'] = G.sh([G.VG, 'gbwt', '-G', 'graph.gfa', '--translation', 'translation.tsv', '-g', 'raw.gbz'], out)
        T['set_ref'] = G.sh([G.VG, 'gbwt', '-Z', 'raw.gbz', '--set-reference', 'REF', '-g', 'graph.gbz'], out)
        (out / 'raw.gbz').unlink()
        trf = out / 'translation.tsv'
        tr = [l.split('\t') for l in trf.read_text().splitlines() if l.startswith('T\t')] if trf.exists() else []
        st['translation_identity'] = all(x[1] == x[2].strip() for x in tr)
        assert st['translation_identity'], 'node ids changed in GBZ'
        G.sh([G.VG, 'paths', '-x', 'graph.gbz', '-F'], out, out / 'paths.fa')
        seqs, name = {}, None
        for line in open(out / 'paths.fa'):
            if line.startswith('>'):
                name = line[1:].strip()
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
        by_sample = {n.split('#')[0]: ''.join(v) for n, v in seqs.items()}
        want = {G.path_name(msa, r): msa.seqs[r] for r in rows}
        want['REF'] = msa.seqs[ref_row]
        miss = [k for k, v in want.items() if by_sample.get(k) != v]
        assert not miss, f'GBZ paths differ: {len(miss)}'
        st['paths_verified_gbz'] = len(want)
        (out / 'paths.fa').unlink()
        G.sh([G.VG, 'stats', '-z', '-A', 'graph.gbz'], out, out / 'vg_stats.txt')
        st['vg_stats'] = (out / 'vg_stats.txt').read_text()
        T['dist'] = G.sh([G.VG, 'index', '-j', 'graph.dist', '-P', 'REF', '-t', 4, 'graph.gbz'], out)
        T['ri'] = G.sh([G.VG, 'gbwt', '-p', '--num-threads', 4, '-r', 'graph.ri', '-Z', 'graph.gbz'], out)
        T['minimizer'] = G.sh([G.VG, 'minimizer', '-t', 4, '-d', 'graph.dist', '-z', 'graph.zipcodes',
                               '-o', 'graph.min', 'graph.gbz'], out)
        st['hapl'] = {}
        for Lc in G.SUBCHAINS:
            T[f'hapl_{Lc}'] = G.sh([G.VG, 'haplotypes', '-v', 1, '-t', 4, '--subchain-length', Lc,
                                    '-d', 'graph.dist', '-r', 'graph.ri', '-H', f'hapl_{Lc}.hapl', 'graph.gbz'], out,
                                   log=f'hapl_{Lc}.log')
            st['hapl'][Lc] = dict(bytes=(out / f'hapl_{Lc}.hapl').stat().st_size,
                                  log_tail=(out / f'hapl_{Lc}.log').read_text()[-1500:])
        st['seconds'] = round(time.time() - t0, 1)
        C.save_json(out / 'STATS.json', st)
        C.save_json(out / 'COMPLETE.json', dict(status='complete', **{k: st[k] for k in ('gene', 'variant', 'fold', 'nodes', 'edges', 'paths', 'seconds')}))
        print(json.dumps({k: v for k, v in st.items() if k not in ('vg_stats',)}, indent=1, default=str))
    except Exception:
        C.save_json(out / 'FAILURE.json', dict(status='failed', error=traceback.format_exc(), **st))
        raise


if __name__ == '__main__':
    main()
