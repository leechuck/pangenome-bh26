"""Partial-credit re-scoring helpers (stage rescore-partial-20261114). Scoring only, no typing.

Definitions (fixed by the step brief):
  gene body   = sequence from the first coding-exon start to the last coding-exon end (T1K header
                exon coordinates; UTR excluded).
  coding      = concatenated coding exons.
  distance    = edlib global edit distance (mode NW).
  non-coding  = max(0, gene-body distance - coding distance)   (negatives are counted).
A slot with several tied alleles is represented by its WORST allele (largest gene-body distance, the
frozen scorers' "all tied labels must be right" rule) and by its BEST allele (sensitivity).
"""
from __future__ import annotations

import collections
import itertools
import random
import sys

import edlib

sys.dont_write_bytecode = True
_ED = {}


def ed(a: str, b: str) -> int:
    if a == b:
        return 0
    if not a or not b:
        return max(len(a), len(b))
    k = (a, b)
    v = _ED.get(k)
    if v is None:
        v = edlib.align(a, b, mode='NW', task='distance')['editDistance']
        _ED[k] = v
    return v


class Seq:
    """A sequence with coding-exon intervals (half-open, relative to seq)."""
    __slots__ = ('name', 'seq', 'exons', '_body', '_cod')

    def __init__(self, name, seq, exons):
        self.name = name
        self.seq = seq
        self.exons = [tuple(e) for e in exons] if exons else []
        self._body = None
        self._cod = None

    @property
    def body(self):
        if self._body is None:
            self._body = self.seq[self.exons[0][0]:self.exons[-1][1]] if self.exons else None
        return self._body

    @property
    def cod(self):
        if self._cod is None:
            self._cod = ''.join(self.seq[s:e] for s, e in self.exons) if self.exons else None
        return self._cod

    @property
    def has_introns(self):
        return bool(self.exons) and len(self.body) > len(self.cod)


def dist_pair(a: Seq, t: Seq):
    """distances of called sequence a to true copy t: dict(gb, cod, full)"""
    cod = ed(a.cod, t.cod) if (a.cod is not None and t.cod is not None) else None
    gb = None
    if (a.exons and t.exons and a.has_introns and t.has_introns and len(a.exons) == len(t.exons)):
        gb = ed(a.body, t.body)
    return dict(gb=gb, cod=cod, full=(a.seq == t.seq))


def slot_vs_copy(alleles, t: Seq):
    """worst/best representative over the tied alleles of one slot vs true copy t."""
    ds = [dist_pair(a, t) for a in alleles]

    def key(d):
        return ((d['gb'] if d['gb'] is not None else 0), (d['cod'] if d['cod'] is not None else 0))
    gb_def = all(d['gb'] is not None for d in ds)
    cod_def = all(d['cod'] is not None for d in ds)
    worst = max(ds, key=key)
    best = min(ds, key=key)
    return dict(gb_w=worst['gb'] if gb_def else None, cod_w=worst['cod'] if cod_def else None,
                gb_b=best['gb'] if gb_def else None, cod_b=best['cod'] if cod_def else None,
                full_w=all(d['full'] for d in ds), full_b=any(d['full'] for d in ds),
                gb_all0=gb_def and all(d['gb'] == 0 for d in ds),
                gb_any0=gb_def and any(d['gb'] == 0 for d in ds),
                cod_all0=cod_def and all(d['cod'] == 0 for d in ds),
                cod_any0=cod_def and any(d['cod'] == 0 for d in ds))


def assign(slots, copies):
    """slots: list (len 2) of list-of-Seq or None/[] (missing). copies: list of true Seq.
    Return list (per copy) of (slot index or None, dict from slot_vs_copy or None)."""
    n = len(copies)
    best = None
    idxs = list(range(len(slots)))
    perms = list(itertools.permutations(idxs, n)) if len(slots) >= n else [tuple(idxs) + (None,) * (n - len(idxs))]
    cache = {}
    for perm in perms:
        tot = 0
        totc = 0
        res = []
        for c, p in zip(copies, perm):
            if p is None or not slots[p]:
                res.append((p, None))
                continue
            if (p, id(c)) not in cache:
                cache[(p, id(c))] = slot_vs_copy(slots[p], c)
            r = cache[(p, id(c))]
            res.append((p, r))
            tot += (r['gb_w'] if r['gb_w'] is not None else r['cod_w'] or 0)
            totc += (r['cod_w'] or 0)
        # primary: smallest total gene-body distance (coding distance for exon-only calls); tie-break: more
        # gene-body-identical copies, then smaller total coding distance, then permutation order
        sc = (tot, -sum(1 for _, r in res if r and (r['gb_all0'] or (r['gb_w'] is None and r['cod_all0']))), totc)
        if best is None or sc < best[0]:
            best = (sc, res)
    return best[1]


# ----------------------------------------------------------------------------- statistics

def uncertainty(pairs, group_of, scale=100.0, nboot=4000, seed=20261020):
    """Family-clustered paired contrast: copy of score.py uncertainty() with a scale argument
    (scale=100 gives percentage points as in the frozen scorers; scale=1 gives mean difference).
    pairs: dicts with 'donor' and numeric 'difference'."""
    clusters = collections.defaultdict(list)
    for r in pairs:
        clusters[group_of.get(r['donor'], (r['donor'],))].append(r['difference'])
    values = [(sum(v), len(v)) for v in clusters.values()]
    if not values:
        return dict(n=0, families=0, diff=None, lo=None, hi=None, p=None)
    rng = random.Random(seed)
    boots = []
    for _ in range(nboot):
        draw = [rng.choice(values) for _ in values]
        boots.append(scale * sum(x[0] for x in draw) / sum(x[1] for x in draw))
    boots.sort()
    differences = [v[0] for v in values if v[0]]
    observed = abs(sum(differences))
    if len(differences) <= 18:
        outcomes = [sum(s * d for s, d in zip(signs, differences))
                    for signs in itertools.product((-1, 1), repeat=len(differences))]
        p = sum(abs(x) >= observed - 1e-9 for x in outcomes) / len(outcomes)
    else:
        hits = sum(abs(sum(rng.choice((-1, 1)) * d for d in differences)) >= observed - 1e-9
                   for _ in range(20000))
        p = (hits + 1) / 20001
    return dict(n=sum(v[1] for v in values), families=len(values),
                diff=scale * sum(x[0] for x in values) / sum(x[1] for x in values),
                lo=boots[100], hi=boots[3899], p=p)


ROW_COLS = ['set', 'readset', 'donor', 'gene', 'copy', 'haplotype', 'ipd_absent', 'arm', 'call',
            'n_tied', 'missing', 'unresolved', 'true_gb_len', 'true_cod_len', 'call_gb_len',
            'gb_dist', 'cod_dist', 'nc_dist', 'nc_negative', 'gb_dist_best', 'cod_dist_best',
            'full_exact', 'full_exact_best', 'gb_exact', 'cod_exact', 'gb_exact_best',
            'gb_defined', 'novel_final']


def make_row(base, arm, copy, res, call_names, alleles, unresolved, novel_final=False):
    """base: dict with set/readset/donor/gene/copy/haplotype/ipd_absent."""
    r = dict(base)
    r['arm'] = arm
    r['call'] = '|'.join(call_names) if call_names else '.'
    r['n_tied'] = len(alleles) if alleles else 0
    r['unresolved'] = unresolved
    r['true_gb_len'] = len(copy.body) if copy.body is not None else ''
    r['true_cod_len'] = len(copy.cod) if copy.cod is not None else ''
    r['novel_final'] = int(novel_final)
    if res is None:
        r['missing'] = 1
        for k in ('call_gb_len', 'gb_dist', 'cod_dist', 'nc_dist', 'nc_negative', 'gb_dist_best',
                  'cod_dist_best', 'full_exact', 'full_exact_best', 'gb_exact', 'cod_exact',
                  'gb_exact_best', 'gb_defined'):
            r[k] = ''
        return r
    r['missing'] = 0
    r['call_gb_len'] = max((len(a.body) for a in alleles if a.body is not None), default='')
    gbd = res['gb_w']
    r['gb_defined'] = int(gbd is not None)
    r['gb_dist'] = '' if gbd is None else gbd
    r['cod_dist'] = '' if res['cod_w'] is None else res['cod_w']
    if gbd is not None and res['cod_w'] is not None:
        diff = gbd - res['cod_w']
        r['nc_negative'] = int(diff < 0)
        r['nc_dist'] = max(0, diff)
    else:
        r['nc_negative'] = ''
        r['nc_dist'] = ''
    r['gb_dist_best'] = '' if res['gb_b'] is None else res['gb_b']
    r['cod_dist_best'] = '' if res['cod_b'] is None else res['cod_b']
    r['full_exact'] = int(res['full_w'])
    r['full_exact_best'] = int(res['full_b'])
    r['gb_exact'] = '' if gbd is None else int(res['gb_all0'])
    r['gb_exact_best'] = '' if gbd is None else int(res['gb_any0'])
    r['cod_exact'] = '' if res['cod_w'] is None else int(res['cod_all0'])
    return r


def write_rows(path, rows, cols=ROW_COLS):
    with open(path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')


def read_fasta_needed(path, names=None):
    """{name: (seq, exons)} for headers whose first token is in names (or all if None)."""
    out = {}
    name, ints, parts, keep = None, [], [], False

    def flush():
        if name is not None and keep:
            n = ints[0] if ints else 0
            ex = [(ints[1 + 2 * i], ints[2 + 2 * i] + 1) for i in range(n)]
            out[name] = (''.join(parts).upper(), ex)
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                flush()
                tok = line[1:].split()
                name = tok[0]
                ints = [int(x) for x in tok[1:]]
                parts = []
                keep = names is None or name in names
            elif keep:
                parts.append(line.strip())
    flush()
    return out
