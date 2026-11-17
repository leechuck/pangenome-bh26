"""Local-version voting and patching (candidate builder v3, stage mosaic-drb1-20261113).

Diagnosis: most DRB1 recombinant ("mosaic") copies differ from the nearest IPD allele by the length
of a short repeat or homopolymer, or by a few sparse sites. Column-wise evidence is least reliable
exactly there. This builder never projects reads onto columns. For each short variable stretch
("block") among the allele-family rows it counts the reads that span the block and contain one of
the known local versions exactly (with invariant anchors on both sides), and patches the pass-1
allele where its own local version has no reads and another known version has clear support.
No truth is read here.
"""
from __future__ import annotations

import collections
from typing import Dict, List

import numpy as np

import drb1_lib as L

PR = L.PR
OUT = L.OUT
BASES = 'ACGT'


def row_string(P_row: np.ndarray, c0: int, c1: int):
    """Bases of a row over columns [c0, c1]; None if the row has no sequence there (OUTSIDE)."""
    seg = P_row[c0:c1 + 1]
    if (seg == OUT).any() or (seg == PR.OTHER).any():
        return None
    return ''.join(BASES[x] for x in seg[seg < 4])


def build_blocks(fp: PR.Panel, k: int = 15, max_len: int = 110, max_k: int = 45):
    """Variable stretches among the family rows, with anchors of >= k invariant bases each side.
    Anchors are lengthened until no local version is a substring of another."""
    P = fp.P
    inside = P != OUT
    mx = np.where(inside, P, 0).max(0)
    mn = np.where(inside, P, 9).min(0)
    var = (mx != mn) & (inside.sum(0) >= 2)
    maj = ((P < 4).sum(0) * 2 >= np.maximum(inside.sum(0), 1)) & (inside.sum(0) > 0)
    cb = np.cumsum(maj)
    cols = np.nonzero(var)[0]
    blocks = []
    if len(cols) == 0:
        return blocks
    groups = [[int(cols[0]), int(cols[0])]]
    for c in cols[1:].tolist():
        if cb[c] - cb[groups[-1][1]] <= k:
            groups[-1][1] = c
        else:
            groups.append([c, c])
    C = P.shape[1]
    for f, l in groups:
        kk = k
        done = None
        while kk <= max_k:
            # anchor columns: kk majority-base columns before f and after l
            cl = int(np.searchsorted(cb, cb[f] - kk - (1 if maj[f] else 0), side='right'))
            cr = int(np.searchsorted(cb, cb[l] + kk, side='left'))
            cl, cr = max(0, cl), min(C - 1, cr)
            alle = collections.defaultdict(list)
            for r in range(fp.R):
                s = row_string(P[r], cl, cr)
                if s is not None:
                    alle[s].append(r)
            keys = sorted(alle, key=len)
            sub = any(a != b and a in b for i, a in enumerate(keys) for b in keys[i + 1:])
            longest = max((len(s) for s in keys), default=0)
            done = dict(first=f, last=l, cl=cl, cr=cr, alleles=dict(alle), anchor=kk,
                        substring_conflict=sub, longest=longest)
            if not sub or longest > max_len:
                break
            kk += 10
        done['resolvable'] = (not done['substring_conflict']) and done['longest'] <= max_len and len(done['alleles']) >= 2
        blocks.append(done)
    return blocks


def count_blocks(blocks, reads, pad_cols: int = 0):
    """reads: dicts with seq (reference orientation), reverse, c0, c1 (MSA window of the placement).
    Adds per block counts[s] = [fwd, rev] of reads containing the local version s."""
    order = sorted(range(len(blocks)), key=lambda i: blocks[i]['cl'])
    cls = np.array([blocks[i]['cl'] for i in order])
    for b in blocks:
        b['counts'] = {s: [0, 0] for s in b['alleles']}
        b['spanning'] = 0
    for rd in reads:
        lo = int(np.searchsorted(cls, rd['c0'] - pad_cols, side='left'))
        hi = int(np.searchsorted(cls, rd['c1'], side='right'))
        for j in range(lo, hi):
            b = blocks[order[j]]
            if not b['resolvable'] or b['cr'] > rd['c1'] + pad_cols:
                continue
            b['spanning'] += 1
            for s in b['alleles']:
                if s in rd['seq']:
                    b['counts'][s][1 if rd['reverse'] else 0] += 1
    return blocks


def decide(blocks, fp: PR.Panel, row_a: int, partner_codes, D: float, read_len: int = 150,
           same_family: bool = False, row_b: int = None, sup_share: float = 0.3, min_sup: int = 3):
    """Patches for the haplotype whose pass-1 allele is row_a (index in fp). Returns (patches, notes).
    A patch = (block index, version string, donor row)."""
    patches, amb = [], []
    for i, b in enumerate(blocks):
        if not b['resolvable']:
            continue
        cur = row_string(fp.P[row_a], b['cl'], b['cr'])
        if cur is None:
            continue
        n = {s: sum(v) for s, v in b['counts'].items()}
        exp_h = (D / 2.0) * max(0.0, (read_len - b['longest']) / read_len)
        thr = max(min_sup, sup_share * exp_h)
        low = max(1.0, 0.1 * exp_h)
        sup = {s for s, v in b['counts'].items() if sum(v) >= thr and v[0] >= 1 and v[1] >= 1}
        pstr = None
        if partner_codes is not None:
            seg = partner_codes[b['cl']:b['cr'] + 1]
            if not ((seg == OUT).any() or (seg == PR.OTHER).any()):
                pstr = ''.join(BASES[x] for x in seg[seg < 4])
        shared = pstr is not None and pstr in b['alleles']
        if n.get(cur, 0) > low:
            # own version has reads. Heterozygous extra version when both haplotypes show the same string
            if shared and pstr == cur:
                extra = [s for s in sup if s != cur and n[s] >= 0.3 * (n[s] + n[cur])]
                if len(extra) == 1:
                    amb.append(dict(block=i, kind='het_extra_version_phase_unknown', version_len=len(extra[0]),
                                    n_cur=n[cur], n_new=n[extra[0]], new=extra[0]))
            continue
        cands = [s for s in sup if s != cur]
        if shared and pstr != cur and pstr in cands:
            # reads of the partner explain its own version; keep it only if it looks present twice
            others = [s for s in cands if s != pstr]
            if others or n[pstr] < 1.6 * max(exp_h, 1.0):
                cands = others
        if not cands:
            continue
        cands.sort(key=lambda s: -n[s])
        if len(cands) > 1 and n[cands[0]] < 2 * n[cands[1]]:
            amb.append(dict(block=i, kind='two_supported_versions', n=[n[cands[0]], n[cands[1]]]))
            continue
        x = cands[0]
        patches.append(dict(block=i, new=x, donor_row=b['alleles'][x][0], n_cur=n.get(cur, 0), n_new=n[x],
                            fwd_rev=b['counts'][x], exp_h=round(exp_h, 1), cl=b['cl'], cr=b['cr'],
                            len_cur=len(cur), len_new=len(x), anchor=b['anchor'], n_versions=len(b['alleles'])))
    return patches, amb


def patched_path(fp: PR.Panel, row_a: int, patches, blocks) -> PR.Path:
    segs, pos = [], 0
    for p in sorted(patches, key=lambda p: p['cl']):
        b = blocks[p['block']]
        if b['cl'] < pos:
            continue
        if b['cl'] > pos:
            segs.append((row_a, pos, b['cl']))
        segs.append((p['donor_row'], b['cl'], b['cr'] + 1))
        pos = b['cr'] + 1
    if pos < fp.C:
        segs.append((row_a, pos, fp.C))
    return PR.merge_path(segs)


def spell(fp: PR.Panel, codes: np.ndarray, prm: PR.Params):
    """Fewest-switch spelling of a sequence (given as MSA codes) from family rows with >= min_seg_bp
    segments: (mismatching columns, switches)."""
    T = np.where((np.arange(PR.NCODE)[None, :] != codes[:, None]), 1000.0, 0.0).astype(np.float32)
    T[:, OUT] = 0
    T[codes == OUT, :] = 0
    c, path = PR.viterbi(fp, T, prm, set(), True)
    return int(round(c // 1000)), path.switches
