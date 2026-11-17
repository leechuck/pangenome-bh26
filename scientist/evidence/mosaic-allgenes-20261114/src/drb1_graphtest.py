"""DEVELOPMENT ONLY (reads truth): how does the cohort pangenome graph store the DRB1 difference sites,
and where does vg giraffe put the reads that carry the true local version?
  python drb1_graphtest.py enc PANEL                 graph encoding, all five folds
  python drb1_graphtest.py map READSET DONOR PANEL   one development donor, graph of the fold that holds it out
Outputs hold positions, lengths, counts and flags only."""
import collections
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L
import drb1_patch2 as PT
C, PR = L.C, L.PR
import panel_select as PS
import run_donor as RD

S2 = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113')
CV = C.ZEN / 'pangenome-reference-cv-20261014'
SIF = C.BENCH / 'containers/ubuntu24.sif'
VG = '/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin/vg'
G = 'DRB1'
K = 15
NAMED = {'homopolymer_5785': (5775, 5795), 'intron2_repeat': (11858, 12081)}
COMP = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(COMP)[::-1]


FLIP = {'+': '-', '-': '+'}


class Graph:
    """Bidirected graph; an oriented node is (id, '+'|'-')."""

    def __init__(self, gfa: Path):
        raw, self.out, self.paths = {}, collections.defaultdict(set), {}
        self.minus = 0
        for line in open(gfa):
            f = line.rstrip('\n').split('\t')
            if f[0] == 'S':
                raw[f[1]] = f[2].upper()
            elif f[0] == 'L':
                self.out[(f[1], f[2])].add((f[3], f[4]))
                self.out[(f[3], FLIP[f[4]])].add((f[1], FLIP[f[2]]))
            elif f[0] == 'P':
                steps = [(x[:-1], x[-1]) for x in f[2].split(',')]
                self.paths[f[1]] = steps
            elif f[0] == 'W':
                steps = [(n, '+' if o == '>' else '-') for o, n in re.findall(r'([<>])([^<>]+)', f[6])]
                self.paths['#'.join(f[1:4])] = steps
        self.raw = raw
        self.seq = {}
        for n, x in raw.items():
            self.seq[(n, '+')] = x
            self.seq[(n, '-')] = rc(x)
        self.minus = sum(o == '-' for st in self.paths.values() for _, o in st)
        self.pseq = {p: ''.join(self.seq[n] for n in st) for p, st in self.paths.items()}
        self.bypos = collections.defaultdict(set)
        for n, x in self.seq.items():
            for o, ch in enumerate(x):
                self.bypos[ch].add((n, o))

    def cyclic(self):
        indeg = collections.Counter()
        for a, bs in self.out.items():
            for b in bs:
                indeg[b] += 1
        q = [n for n in self.seq if indeg[n] == 0]
        seen = 0
        while q:
            n = q.pop()
            seen += 1
            for b in self.out[n]:
                indeg[b] -= 1
                if indeg[b] == 0:
                    q.append(b)
        return seen != len(self.seq)

    def thread(self, s):
        """Fewest restarts needed to spell s as walks over existing nodes and edges (0 = one walk)."""
        breaks, cur = 0, None
        for ch in s:
            if cur is not None:
                nxt = set()
                for n, o in cur:
                    sn = self.seq[n]
                    if o + 1 < len(sn):
                        if sn[o + 1] == ch:
                            nxt.add((n, o + 1))
                    else:
                        for m in self.out[n]:
                            if self.seq[m][0] == ch:
                                nxt.add((m, 0))
                if nxt:
                    cur = nxt
                    continue
                breaks += 1
            cur = set(self.bypos.get(ch, ()))
        return breaks


def period(x):
    for p in range(1, 7):
        if len(x) >= 2 * p and all(x[i] == x[i - p] for i in range(p, len(x))):
            return p
    return 0


def fold_of(donor):
    for f in C.load_json(CV / 'FOLDS.json')['folds']:
        if donor in f['test']:
            return f['fold']
    raise KeyError(donor)


def anchors(msa, c0, c1, k=20):
    la, ra = set(), set()
    for al in msa.aligned:
        l = al[:c0].replace('-', '')[-k:]
        r = al[c1 + 1:].replace('-', '')[:k]
        if len(l) == k and len(r) == k and not set(al[:c0]) <= {'-'} and not set(al[c1 + 1:]) <= {'-'}:
            la.add(l)
            ra.add(r)
    return la, ra


def locate(ps, la, ra, k=20, maxlen=700):
    """(start, end) of the stretch between a left and a right anchor in a path sequence, or None."""
    best = None
    for l in la:
        i = ps.find(l)
        while i >= 0:
            s0 = i + k
            for r in ra:
                j = ps.find(r, s0, s0 + maxlen + k)
                if j >= 0 and (best is None or j - s0 < best[1] - best[0]):
                    best = (s0, j)
            i = ps.find(l, i + 1)
    return best


def enc(panel):
    records = C.load_records()
    msa = C.load_msas(records, genes=(G,))[G]
    anc = {name: anchors(msa, c0, c1) for name, (c0, c1) in NAMED.items()}
    out = ['fold\tpanel\tnodes\tedges\tpaths\treverse_steps\tdirected_cycle\tpaths_repeating_a_node_same_direction\tpaths_repeating_a_node_other_direction_only']
    loc = ['fold\tlocus\tpaths_located\tdistinct_strings\tdistinct_lengths\tmin_len\tmax_len\tdistinct_walks\twalks_repeating_a_node\tnodes_used']
    loops = ['fold\tnode_len\trepeat_unit_len\tpaths_passing_more_than_once\tmax_passes\tmedian_path_offset\tinside_locus']
    for k in range(5):
        g = Graph(CV / f'fold{k}/graphs/{panel}/{G}/graph.gfa')
        same = sum(len(set(st)) != len(st) for st in g.paths.values())
        other = sum(len(set(st)) == len(st) and len({n for n, _ in st}) != len(st) for st in g.paths.values())
        out.append('\t'.join(map(str, (k, panel, len(g.raw), sum(map(len, g.out.values())) // 2, len(g.paths), g.minus,
                                       g.cyclic(), same, other))))
        span = {}
        for name in NAMED:
            strings, walks = collections.Counter(), collections.Counter()
            for p, ps in g.pseq.items():
                iv = locate(ps, *anc[name])
                if iv is None:
                    continue
                span[(name, p)] = iv
                strings[ps[iv[0]:iv[1]]] += 1
                w, pos = [], 0
                for n in g.paths[p]:
                    e = pos + len(g.seq[n])
                    if e > iv[0] and pos < iv[1]:
                        w.append(n)
                    pos = e
                walks[tuple(w)] += 1
            lens = sorted({len(x) for x in strings})
            loc.append('\t'.join(map(str, (k, name, sum(strings.values()), len(strings), len(lens), lens[0] if lens else '',
                                           lens[-1] if lens else '', len(walks),
                                           sum(len(set(w)) != len(w) for w in walks),
                                           len({n for w in walks for n in w})))))
        multi = collections.defaultdict(list)
        for p, st in g.paths.items():
            cnt = collections.Counter(st)
            pos, first = 0, {}
            for n in st:
                first.setdefault(n, pos)
                pos += len(g.seq[n])
            for n, c in cnt.items():
                if c > 1:
                    where = next((name for name in NAMED if (name, p) in span
                                  and span[(name, p)][0] - 30 <= first[n] <= span[(name, p)][1] + 30), '')
                    multi[n].append((c, first[n], where))
        for n, v in sorted(multi.items(), key=lambda kv: np.median([x[1] for x in kv[1]])):
            w = collections.Counter(x[2] for x in v).most_common(1)[0][0]
            loops.append('\t'.join(map(str, (k, len(g.seq[n]), period(g.seq[n]), len(v), max(x[0] for x in v),
                                             int(np.median([x[1] for x in v])), w))))
    txt = '\n'.join(out) + '\n\n' + '\n'.join(loc) + '\n\n' + '\n'.join(loops) + '\n'
    (S2 / f'results/graphtest_encoding_{panel}.tsv').write_text(txt)
    print(txt)


def vg(args, cwd, stdout=None):
    binds = f'{C.BENCH}/hla:/home/leechuck/hla,/encrypted/e3008:/encrypted/e3008,/ibex/scratch/projects/c2014:/ibex/scratch/projects/c2014'
    cmd = ['singularity', 'exec', '--cleanenv', '--no-home', '--pwd', str(cwd), '--bind', binds, str(SIF), VG] + list(map(str, args))
    with open(cwd / 'vg.log', 'ab') as log:
        if stdout:
            with open(stdout, 'wb') as o:
                subprocess.run(cmd, stdout=o, stderr=log, check=True)
        else:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)


def sites(P_t, P_a, gap=6):
    """Column clusters where the true copy differs from the called allele, with K true bases either side."""
    inside = (P_t != L.OUT) & (P_a != L.OUT)
    d = np.nonzero(inside & (P_t != P_a))[0]
    cl = []
    for c in d:
        if cl and c - cl[-1][1] <= gap:
            cl[-1][1] = int(c)
        else:
            cl.append([int(c), int(c)])
    isb = P_t < 4
    res = []
    for a, b in cl:
        l, n = a, 0
        while l > 0 and n < K:
            l -= 1
            if not inside[l]:
                break
            n += isb[l]
        r, m = b, 0
        while r < len(P_t) - 1 and m < K:
            r += 1
            if not inside[r]:
                break
            m += isb[r]
        if n < K or m < K:
            continue
        res.append((l, r, a, b))
    return res


def mapd(readset, donor, panel):
    assert readset in ('D-sim', 'D-real')
    out = S2 / 'work/graphtest' / panel / readset / donor
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    records = C.load_records()
    msa = C.load_msas(records, genes=(G,))[G]
    shas = [C.sha256_str(s) for s in msa.seqs]
    rows = PS.select_rows('C', [bool(n) for n in msa.ipd_names], shas, None)
    names = [msa.ipd_names[r][0] for r in rows]
    ipd = PR.Panel(names, [msa.aligned[r] for r in rows])
    rown = {}
    for i, r in enumerate(rows):
        for n in msa.ipd_names[r]:
            rown.setdefault(n, i)
    k = fold_of(donor)
    gdir = CV / f'fold{k}/graphs/{panel}/{G}'
    g = Graph(gdir / 'graph.gfa')
    pathset = set(g.pseq.values())
    tr = [r for r in records if r['donor'] == donor and r['gene'] == G]
    tp = PR.Panel([r['source_id'] for r in tr], [msa.aligned[msa.row_of_source[r['source_id']]] for r in tr])
    cls = C.load_classes()
    prop = C.load_json(S2 / 'work' / readset / donor / 'C3/proposals.json')
    p1 = [h['pass1_allele'] for h in prop['haplotypes'] if 'pass1_allele' in h]
    # pair each true copy with the nearest pass-1 allele
    todo = []
    for i, r in enumerate(tr):
        a = min(p1, key=lambda n: L.mism(ipd.P[rown[n]], tp.P[i]))
        todo.append((i, r, a))
    # giraffe, same commands as the earlier graph mapping stage (unsampled graph)
    spec = RD.pass1_spec(readset, donor, 'C')
    shutil.copy2(gdir / 'graph.gbz', out / 'graph.gbz')
    backbone = C.load_json(gdir / 'manifest.json')['indexing_backbone']
    vg(['index', '-j', 'mapping.dist', '-t', 4, '-P', backbone, 'graph.gbz'], out)
    vg(['giraffe', '-Z', 'graph.gbz', '-d', 'mapping.dist', '-f', spec['run'] / 'calls_candidate_1.fq',
        '-f', spec['run'] / 'calls_candidate_2.fq', '-t', 4], out, out / 'mapping.gam')
    vg(['view', '-aj', 'mapping.gam'], out, out / 'aln.jsonl')
    reads = []
    for line in open(out / 'aln.jsonl'):
        a = json.loads(line)
        s = a.get('sequence', '').upper()
        w = []
        for m in a.get('path', {}).get('mapping', []):
            p = m['position']
            ns = g.seq[(str(p['node_id']), '-' if p.get('is_reverse') else '+')]
            o = int(p.get('offset', 0))
            fl = sum(int(e.get('from_length', 0)) for e in m.get('edit', []))
            w.append(ns[o:o + fl])
        reads.append((s, ''.join(w)))
    n_mapped = sum(bool(w) for _, w in reads)
    copy_rows, site_rows = [], []
    for i, r, a in todo:
        c = cls[r['source_id']]['ipd_label_L150']
        absent = not r['genomic_labels']
        seq = r['core_sequence'].upper()
        copy_rows.append('\t'.join(map(str, (readset, donor, panel, k, i, c, absent, a, len(seq), seq in pathset,
                                             g.thread(seq), len(reads), n_mapped))))
        for l, rr, a0, b0 in sites(tp.P[i], ipd.P[rown[a]]):
            t = PT.row_string(tp.P[i], l, rr)
            cu = PT.row_string(ipd.P[rown[a]], l, rr)
            if t is None or cu is None:
                continue
            kind = 'length' if len(t) != len(cu) else 'substitution'
            named = next((n for n, (c0, c1) in NAMED.items() if a0 <= c1 and b0 >= c0), '')
            in_path = any(t in s for s in pathset)
            walk = g.thread(t) == 0
            ix = seq.find(t)
            ext = seq[max(0, ix - 60):ix + len(t) + 60] if ix >= 0 else ''
            win = bool(ext) and any(ext in s for s in pathset)
            win_walk = bool(ext) and g.thread(ext) == 0
            cnt = collections.Counter()
            if len(t) <= 120:
                trc, crc = rc(t), rc(cu)
                for s, w in reads:
                    if t in s or trc in s:
                        x = t if t in s else trc
                        cnt['t'] += 1
                        cnt['t_exact' if x in w else 't_edit' if w else 't_unmapped'] += 1
                    elif cu in s or crc in s:
                        cnt['c'] += 1
            site_rows.append('\t'.join(map(str, (readset, donor, panel, i, c, absent, a0, b0, named, kind, len(t), len(cu),
                                                 len(t) <= 120, in_path, walk, cnt['t'], cnt['t_exact'], cnt['t_edit'],
                                                 cnt['t_unmapped'], cnt['c'], win, win_walk))))
    (out / 'copies.tsv').write_text('\n'.join(copy_rows) + '\n')
    (out / 'sites.tsv').write_text('\n'.join(site_rows) + ('\n' if site_rows else ''))
    for f in ('mapping.gam', 'aln.jsonl', 'graph.gbz', 'mapping.dist'):
        (out / f).unlink(missing_ok=True)
    for f in out.glob('graph.*'):
        f.unlink()
    C.save_json(out / 'COMPLETE.json', dict(donor=donor, readset=readset, panel=panel, fold=k, reads=len(reads), mapped=n_mapped))
    print(donor, readset, 'copies', len(copy_rows), 'sites', len(site_rows), 'mapped', n_mapped, 'of', len(reads))


if __name__ == '__main__':
    if sys.argv[1] == 'enc':
        enc(sys.argv[2])
    else:
        mapd(sys.argv[2], sys.argv[3], sys.argv[4])
