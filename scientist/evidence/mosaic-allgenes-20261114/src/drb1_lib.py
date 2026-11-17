"""DRB1-focused recombinant proposer, version 2 (stage mosaic-drb1-20261113).

Reuses the frozen stage's modules read-only (common, evidence, proposer, run_donor) and changes
four things, each aimed at one diagnosed cause of the DRB1 failure:
  1. reads      projection against the allele-family rows of the pass-1 alleles only, small edit
                distance, and no mate placed on another gene (removes look-alike-gene reads);
  2. evidence   a character is supported only at a share of the expected depth, pile-up columns
                carry no evidence (removes low-level noise);
  3. search     each haplotype may only copy from its own family rows, no read-count term inside
                the search, higher switch cost;
  4. gate       relative gain over the best no-switch pair instead of an absolute residual cap.
No truth is read here.
"""
from __future__ import annotations

import collections
import sys
from pathlib import Path

import numpy as np

OLD = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src')
sys.path.insert(0, str(OLD))
import common as C          # noqa: E402
import evidence as EV       # noqa: E402
import proposer as PR       # noqa: E402

OUT = PR.OUTSIDE


# ----------------------------------------------------------------------------- reads

def map_reads2(target: dict, fq1, fq2, work: Path, threads: int, env: Path):
    """bwa mem to the pass-1 alleles of all genes; keep per read its gene and its mate's gene."""
    import subprocess
    import pysam
    work.mkdir(parents=True, exist_ok=True)
    fa = work / 'target.fa'
    with open(fa, 'w') as f:
        for n, s in target.items():
            f.write(f'>{n}\n{s}\n')
    subprocess.run([str(env / 'bwa'), 'index', str(fa)], check=True, capture_output=True)
    sam = work / 'aln.sam'
    with open(sam, 'w') as f:
        subprocess.run([str(env / 'bwa'), 'mem', '-t', str(threads), '-v', '1', str(fa), str(fq1),
                        str(fq2)], stdout=f, stderr=subprocess.DEVNULL, check=True)
    reads = collections.defaultdict(list)
    gene_of_q = collections.defaultdict(set)
    stats = collections.Counter()
    with pysam.AlignmentFile(str(sam), 'r') as af:
        lens = dict(zip(af.references, af.lengths))
        for rec in af.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary:
                continue
            stats['reads'] += 1
            if rec.is_unmapped:
                stats['unmapped'] += 1
                continue
            g = C.gene_of_name(rec.reference_name)
            stats['gene:' + g] += 1
            gene_of_q[rec.query_name].add(g)
            if g not in C.GENES:
                continue
            cig = rec.cigartuples or []
            lclip = cig[0][1] if cig and cig[0][0] == 4 else 0
            rclip = cig[-1][1] if cig and cig[-1][0] == 4 else 0
            reads[g].append(dict(qname=rec.query_name, seq=rec.query_sequence,
                                 reverse=rec.is_reverse, contig=rec.reference_name,
                                 ref_start=rec.reference_start, ref_end=rec.reference_end,
                                 lclip=lclip, rclip=rclip, contig_len=lens[rec.reference_name],
                                 nm=rec.get_tag('NM') if rec.has_tag('NM') else 99))
    sam.unlink()
    return reads, dict(stats), gene_of_q


def windows(panel: PR.Panel, reads, row_of_contig):
    """Read dicts -> projection inputs (same rules as run_donor.proposer_arm)."""
    rd, unplaced = [], 0
    for x in reads:
        r = row_of_contig(x['contig'])
        if r is None:
            unplaced += 1
            continue
        bcols = np.nonzero(panel.isbase[r])[0]
        c0, c1, tl, tr = EV.window_from_alignment(bcols, x['ref_start'], x['ref_end'], x['lclip'],
                                                 x['rclip'], x['contig_len'])
        seq = x['seq']
        if tr:
            seq = seq[:len(seq) - tr]
        if tl:
            seq = seq[tl:]
        if len(seq) < 30:
            unplaced += 1
            continue
        rd.append(dict(qname=x['qname'], seq=seq, reverse=x['reverse'], c0=c0, c1=c1, prefer_row=r,
                       nm=x.get('nm', 99), clip=x['lclip'] + x['rclip'],
                       cs=int(bcols[x['ref_start']]), ce=int(bcols[x['ref_end'] - 1])))
    return rd, unplaced


def project_all(panel: PR.Panel, rd, max_keep: int = 6):
    """Per read: (qname, reverse, edit distance to the best panel row, cols, codes, overflow)."""
    widx = EV.WindowIndex(panel)
    out = []
    for x in rd:
        w0, w1 = EV.quantize(x['c0'], x['c1'], panel.C)
        r, d = EV.best_row(x['seq'], widx, w0, w1, x.get('prefer_row'))
        if r is None:
            out.append((x['qname'], x['reverse'], None, None, None, None))
            continue
        if d > max_keep:
            out.append((x['qname'], x['reverse'], d, None, None, None))
            continue
        pr = EV.project(x['seq'], panel, r, w0, w1)
        if pr is None:
            out.append((x['qname'], x['reverse'], None, None, None, None))
            continue
        cols, codes, ovf, ed = pr
        out.append((x['qname'], x['reverse'], ed, cols, codes, ovf))
    return out


def make_evidence(ncol: int, proj, max_edit: int, bad_q=frozenset(), whole_fragment=True) -> PR.Evidence:
    """Evidence from projected reads with edit distance <= max_edit. A fragment is dropped when
    its qname is in bad_q (mate on another gene) or, if whole_fragment, when any of its DRB1-placed
    reads fails the edit filter."""
    fail = set()
    if whole_fragment:
        for q, rev, d, cols, codes, ovf in proj:
            if d is None or d > max_edit or cols is None:
                fail.add(q)
    counts = np.zeros((ncol, PR.NCODE, 2), dtype=np.int32)
    overflow = np.zeros(ncol, dtype=np.int32)
    frag = collections.defaultdict(list)
    st = collections.Counter()
    for q, rev, d, cols, codes, ovf in proj:
        if d is None or cols is None or d > max_edit:
            st['excluded_edit'] += 1
            continue
        if q in bad_q:
            st['excluded_mate_other_gene'] += 1
            continue
        if q in fail:
            st['excluded_mate_failed'] += 1
            continue
        st['used'] += 1
        np.add.at(counts, (cols, codes, np.full(len(cols), 1 if rev else 0)), 1)
        for c in ovf:
            overflow[c] += 1
        frag[q].append((cols, codes))
    fc, fk = [], []
    for q, parts in frag.items():
        if len(parts) == 1:
            cols, codes = parts[0]
        else:
            m, bad = {}, set()
            for cols, codes in parts:
                for c, k in zip(cols.tolist(), codes.tolist()):
                    if c in m and m[c] != k:
                        bad.add(c)
                    m[c] = k
            keys = sorted(c for c in m if c not in bad)
            cols = np.array(keys, dtype=np.int64)
            codes = np.array([m[c] for c in keys], dtype=np.uint8)
        fc.append(cols)
        fk.append(codes)
    st['fragments'] = len(fc)
    return PR.Evidence(ncol=ncol, counts=counts, frag_cols=fc, frag_codes=fk, overflow=overflow,
                       stats=dict(st))


# ----------------------------------------------------------------------------- families

def subst_rate(P: np.ndarray, a: int) -> np.ndarray:
    """Per row: share of columns where both that row and row a hold a base and the bases differ."""
    both = (P < 4) & (P[a] < 4)[None, :]
    diff = both & (P != P[a][None, :])
    return diff.sum(1) / np.maximum(both.sum(1), 1)


def first_field(name: str):
    """First field of an IPD allele name; None for cohort / candidate rows."""
    if '*' not in name:
        return None
    f = name.split('*')[1].split(':')[0]
    return f if f.isdigit() else None


def family_rows(panel: PR.Panel, a: int, theta: float):
    rate = subst_rate(panel.P, a)
    ff = first_field(panel.names[a])
    return sorted({r for r in range(panel.R)
                   if rate[r] <= theta or (ff is not None and first_field(panel.names[r]) == ff)})


def subpanel(panel: PR.Panel, rows):
    sp = PR.Panel.__new__(PR.Panel)
    sp.names = [panel.names[r] for r in rows]
    sp.P = panel.P[rows]
    sp.R, sp.C = sp.P.shape
    sp.isbase = panel.isbase[rows]
    sp.Pbase = panel.Pbase[rows]
    sp.seq = [panel.seq[r] for r in rows]
    sp._tcache = {}
    sp.rows_in_parent = list(rows)
    return sp


# ----------------------------------------------------------------------------- evidence rules

def tables_frozen(ev: PR.Evidence, prm: PR.Params):
    sup, contra, depth = PR.column_tables(ev, prm)
    return sup, contra, depth, None


def tables_quant(ev: PR.Evidence, a_sup=0.2, c_dep=0.5, over=1.6, min_reads=3):
    n = ev.counts.sum(axis=2)
    depth = n.sum(axis=1)
    pos = depth[depth > 0]
    D = float(np.median(pos)) if len(pos) else 0.0
    fwd, rev = ev.counts[:, :, 0], ev.counts[:, :, 1]
    sup = (n >= max(min_reads, a_sup * D)) & (fwd >= 1) & (rev >= 1)
    sup[:, PR.OTHER] = False
    sup[:, OUT] = False
    contra = (depth[:, None] >= c_dep * D) & (n <= 1) & ~sup
    contra[:, OUT] = False
    mask = depth > over * D
    sup[mask] = False
    contra[mask] = False
    return sup, contra, depth, dict(D=D, masked_overdeep=int(mask.sum()))


def tables_nogap(ev: PR.Evidence, fam_P: np.ndarray, widen: int = 0, **kw):
    """tables_quant, with no evidence at columns where any family row has a gap (insertion/deletion
    positions, optionally widened), and gap observations never counted as supported."""
    sup, contra, depth, info = tables_quant(ev, **kw)
    gapcol = (fam_P == 4).any(axis=0)
    if widen:
        idx = np.nonzero(gapcol)[0]
        for k in range(1, widen + 1):
            gapcol[np.clip(idx - k, 0, len(gapcol) - 1)] = True
            gapcol[np.clip(idx + k, 0, len(gapcol) - 1)] = True
    sup[:, 4] = False
    contra[:, 4] = False
    sup[gapcol] = False
    contra[gapcol] = False
    info = dict(info, masked_gap_columns=int(gapcol.sum()))
    return sup, contra, depth, info


# ----------------------------------------------------------------------------- search

def solve2(pa: PR.Panel, pb: PR.Panel, sup, contra, ra: int, rb: int, prm: PR.Params,
           prior_a, prior_b, allow_switch: bool, first='b'):
    """Alternating exact Viterbi; haplotype A copies from panel pa only, B from pb only."""
    ha, hb = PR.Path([(ra, 0, pa.C)]), PR.Path([(rb, 0, pb.C)])
    seen = set()
    for _ in range(prm.max_iter):
        order = ('b', 'a') if first == 'b' else ('a', 'b')
        for which in order:
            if which == 'b':
                T = PR.pair_cost_table(sup, contra, ha.codes(pa.P))
                _, hb = PR.viterbi(pb, T, prm, prior_b, allow_switch)
            else:
                T = PR.pair_cost_table(sup, contra, hb.codes(pb.P))
                _, ha = PR.viterbi(pa, T, prm, prior_a, allow_switch)
        key = (ha.key(), hb.key())
        if key in seen:
            break
        seen.add(key)
    cost = PR.pair_evidence_cost(sup, contra, ha.codes(pa.P), hb.codes(pb.P))
    return dict(ha=ha, hb=hb, cost=cost, switches=ha.switches + hb.switches)


def total(s, lam):
    return s['cost'] + lam * s['switches']


def search(pa, pb, sup, contra, starts_a, starts_b, prm: PR.Params):
    """Best no-switch pair and best switching pair. starts_*: rows (in pa / pb) of pass-1 labels."""
    pri_a, pri_b = set(starts_a), set(starts_b)
    zero = []
    for ra in starts_a[:2]:
        for rb in starts_b[:2]:
            c = PR.pair_evidence_cost(sup, contra, pa.P[ra], pb.P[rb])
            zero.append(dict(ha=PR.Path([(ra, 0, pa.C)]), hb=PR.Path([(rb, 0, pb.C)]), cost=c, switches=0))
            for first in 'ab':
                zero.append(solve2(pa, pb, sup, contra, ra, rb, prm, pri_a, pri_b, False, first))
    best0 = min(zero, key=lambda z: z['cost'])
    sols = []
    for ra, rb in {(starts_a[0], starts_b[0]), (best0['ha'].segs[0][0], best0['hb'].segs[0][0])}:
        for first in 'ab':
            sols.append(solve2(pa, pb, sup, contra, ra, rb, prm, pri_a, pri_b, True, first))
    # no-switch refinement from the rows the switching solutions use
    for so in list(sols):
        for ra in {r for r, _, _ in so['ha'].segs}:
            for rb in {r for r, _, _ in so['hb'].segs}:
                c = PR.pair_evidence_cost(sup, contra, pa.P[ra], pb.P[rb])
                if c < best0['cost']:
                    best0 = dict(ha=PR.Path([(ra, 0, pa.C)]), hb=PR.Path([(rb, 0, pb.C)]), cost=c, switches=0)
    best = min(sols + [best0], key=lambda z: (total(z, prm.lam), z['switches']))
    return best0, best


def blind_intervals(panel: PR.Panel, row: int, spans, min_run: int = 100, edge: int = 150):
    """Stretches of >= min_run bases of the allele `row` covered by at most one well-fitting read.
    spans: (first column, last column) of reads placed on this haplotype with few edits.
    Returns [(first column, last column)] and the per-column depth."""
    C_ = panel.C
    d = np.zeros(C_ + 1, dtype=np.int32)
    for cs, ce in spans:
        d[cs] += 1
        d[ce + 1] -= 1
    depth = np.cumsum(d)[:C_]
    bcols = np.nonzero(panel.isbase[row])[0]
    flags = depth[bcols] <= 1
    out, i, n = [], 0, len(bcols)
    while i < n:
        if not flags[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and flags[j + 1]:
            j += 1
        if j - i + 1 >= min_run and i >= edge and j < n - edge:
            out.append((int(bcols[i]), int(bcols[j])))
        i = j + 1
    return out, depth


def near_blind(panel: PR.Panel, row: int, cl: int, cr: int, blind, near: int = 150) -> bool:
    pb = panel.Pbase[row]
    a, b = int(pb[cl]), int(pb[cr + 1])
    for bs, be in blind:
        x, y = int(pb[bs]), int(pb[be + 1])
        if a - near <= y and b + near >= x:
            return True
    return False


def mism(a: np.ndarray, b: np.ndarray) -> int:
    v = (a != OUT) & (b != OUT)
    return int((v & (a != b)).sum())
