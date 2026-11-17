"""Shared helpers for stage graph-paths-20261114 (step 5 of route-20261114).
Imports the earlier stages' modules read-only (common, proposer, drb1_lib, drb1_graphtest, drb1_patch2,
run_donor, panel_select). No truth is read here."""
from __future__ import annotations

import collections
import json
import os
import subprocess
import sys
import time
from pathlib import Path

E = Path('/encrypted/e3008/hohndor-dogohla-20261020')
S = E / 'graph-paths-20261114'
sys.path.insert(0, str(E / 'mosaic-drb1-20261113/src'))
import drb1_lib as L            # noqa: E402  (puts mosaic-typer src on sys.path)
import drb1_graphtest as GT     # noqa: E402
import drb1_patch2 as PT        # noqa: E402
import run_donor as RD          # noqa: E402
import panel_select as PS       # noqa: E402
C, PR = L.C, L.PR

CV = C.ZEN / 'pangenome-reference-cv-20261014'
SIF = C.BENCH / 'containers/ubuntu24.sif'
VG = '/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin/vg'
KMC = '/home/leechuck/hla/t1k-pangenome/tools/kmc-3.2.4/kmc'
GENES3 = ('DRB1', 'DQA1', 'B')
VARIANTS = ('ipd', 'ipd_cohort')
SUBCHAINS = (10000, 3000, 1000, 500, 300)
MAX_NODE = 32
GDIR = S / 'work/graphs'


def graph_dir(gene, variant, fold):
    return GDIR / gene / (variant if variant == 'ipd' else f'{variant}_fold{fold}')


def folds():
    return {f['fold']: f for f in C.load_json(CV / 'FOLDS.json')['folds']}


def sh(cmd, cwd: Path, stdout=None, log='run.log', timeout=None):
    """Run a command; container commands go through singularity. Logs stderr to cwd/log."""
    cmd = [str(x) for x in cmd]
    if cmd[0] in (VG, KMC):
        binds = (f'{C.BENCH}/hla:/home/leechuck/hla,/encrypted/e3008:/encrypted/e3008,'
                 f'/ibex/scratch/projects/c2014:/ibex/scratch/projects/c2014')
        cmd = ['singularity', 'exec', '--cleanenv', '--no-home', '--pwd', str(cwd), '--bind', binds, str(SIF)] + cmd
    if timeout:
        cmd = ['timeout', '-s', 'KILL', str(int(timeout))] + cmd
    t0 = time.time()
    with open(Path(cwd) / log, 'ab') as lg:
        lg.write(('$ ' + ' '.join(cmd) + '\n').encode())
        lg.flush()
        if stdout:
            with open(stdout, 'wb') as o:
                subprocess.run(cmd, stdout=o, stderr=lg, check=True, cwd=cwd)
        else:
            subprocess.run(cmd, stdout=lg, stderr=lg, check=True, cwd=cwd)
    return round(time.time() - t0, 2)


def panel_rows(msa, variant, fold):
    """IPD rows (one per distinct IPD genomic sequence) plus, for ipd_cohort, every cohort row carried by
    a training donor of the fold (identical cohort copies share one row)."""
    rows = [r for r in range(len(msa.uids)) if msa.ipd_names[r]]
    if variant == 'ipd_cohort':
        train = set(folds()[fold]['train'])
        rows += [r for r in range(len(msa.uids)) if not msa.ipd_names[r] and msa.cohort_donors[r] & train]
    return rows


def path_name(msa, r):
    return f'IPD{r}' if msa.ipd_names[r] else f'COH{r}'


def build_gfa(msa, rows, gfa: Path, ref_row: int):
    """Acyclic graph straight from the alignment: one element per (column, base); runs of elements with a
    single successor/predecessor are merged into nodes of <= MAX_NODE bases. Each row becomes a W line
    (sample = path_name, haplotype 1, contig = gene); a REF sample (haplotype 0) repeats ref_row's walk.
    Returns (node_seq, node_cols, walks)."""
    out_e, in_e = collections.defaultdict(set), collections.defaultdict(set)
    starts, ends, elems = set(), set(), {}
    for r in rows:
        u = [(c, ch) for c, ch in enumerate(msa.aligned[r]) if ch != '-']
        elems[r] = u
        starts.add(u[0])
        ends.add(u[-1])
        for a, b in zip(u, u[1:]):
            out_e[a].add(b)
            in_e[b].add(a)
        for x in u:
            out_e.setdefault(x, set())
            in_e.setdefault(x, set())
    nid_of, node_seq, node_cols = {}, {}, {}
    nxt = 1
    for u in sorted(out_e):
        if u in nid_of:
            continue
        n = nxt
        nxt += 1
        chars, cols, cur = [u[1]], [u[0]], u
        nid_of[u] = n
        while len(chars) < MAX_NODE and cur not in ends and len(out_e[cur]) == 1:
            v = next(iter(out_e[cur]))
            if v in starts or len(in_e[v]) != 1 or v in nid_of:
                break
            chars.append(v[1])
            cols.append(v[0])
            nid_of[v] = n
            cur = v
        node_seq[n] = ''.join(chars)
        node_cols[n] = cols
    walks, edges = {}, set()
    for r in rows:
        w = []
        for x in elems[r]:
            n = nid_of[x]
            if not w or w[-1] != n:
                if w:
                    edges.add((w[-1], n))
                w.append(n)
        walks[r] = w
    with open(gfa, 'w') as f:
        f.write('H\tVN:Z:1.1\n')
        for n in range(1, nxt):
            f.write(f'S\t{n}\t{node_seq[n]}\n')
        for a, b in sorted(edges):
            f.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
        g = msa.gene
        for r in rows:
            w = walks[r]
            ln = sum(len(node_seq[n]) for n in w)
            f.write(f'W\t{path_name(msa, r)}\t1\t{g}\t0\t{ln}\t' + ''.join(f'>{n}' for n in w) + '\n')
        w = walks[ref_row]
        ln = sum(len(node_seq[n]) for n in w)
        f.write(f'W\tREF\t0\t{g}\t0\t{ln}\t' + ''.join(f'>{n}' for n in w) + '\n')
    return node_seq, node_cols, walks, edges


def load_nodes(gdir: Path):
    seq, cols = {}, {}
    with open(gdir / 'graph.gfa') as f:
        for line in f:
            if line[0] == 'S':
                _, n, s = line.rstrip('\n').split('\t')[:3]
                seq[int(n)] = s
    for k, v in C.load_json(gdir / 'node_cols.json').items():
        cols[int(k)] = v
    return seq, cols


def parse_walk(s):
    """GAF/W walk '>1>2<3' -> [(1,'+'),(2,'+'),(3,'-')]"""
    out, i = [], 0
    import re
    for o, n in re.findall(r'([<>])(\d+)', s):
        out.append((int(n), '+' if o == '>' else '-'))
    return out


COMP = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(COMP)[::-1]


def kmers(seq, k=29):
    """canonical k-mer strings"""
    seq = seq.upper()
    out = set()
    for i in range(len(seq) - k + 1):
        x = seq[i:i + k]
        if 'N' in x:
            continue
        y = rc(x)
        out.add(x if x < y else y)
    return out


def read_fastq_seqs(path):
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().rstrip('\n')
            f.readline()
            f.readline()
            yield h, s


def qname(header):
    q = header[1:].split()[0]
    if q.endswith('/1') or q.endswith('/2'):
        q = q[:-2]
    return q
