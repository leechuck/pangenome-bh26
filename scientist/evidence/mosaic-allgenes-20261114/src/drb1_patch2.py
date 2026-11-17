"""Local-version voting and patching, version 3.1 (stage mosaic-drb1-20261113).

Changes from drb1_patch.py, each from a measured failure on development donors:
  * events are small clusters of variable columns (gap <= 6 bases), so dense families no longer
    chain into unreadable blocks; anchors may contain variation;
  * a read supports a local version only if the WHOLE read fits the haplotype carrying that version
    (edit distance <= max_edit) and does not fit the partner haplotype at least as well. This removes
    reads from repeat copies elsewhere that share a 31-base string;
  * support thresholds are relative to the expected number of spanning reads, with an upper bound
    (pile-ups) and a ratio rule instead of "own version has no reads" (stutter in real reads).
No truth is read here.
"""
from __future__ import annotations

import collections

import numpy as np

import drb1_lib as L

PR = L.PR
OUT = L.OUT
BASES = 'ACGT'

# Deterministic work bounds (stage mosaic-allgenes, item 1). Exceeding one raises BoundExceeded; the
# driver then records a FAILURE for that gene and keeps the pass-1 call (fallback). Counters are
# reset per gene by drb1_run.run() and written to proposals.json.
WORK = collections.Counter()
WORK_CAP = dict(family_rows=20_000, events=20_000, align_calls=20_000_000)   # dev max: 5.4k rows, 0.4k events, 14k aligns


class BoundExceeded(RuntimeError):
    pass


def _work(key, n=1):
    WORK[key] += n
    if WORK[key] > WORK_CAP[key]:
        raise BoundExceeded(f'{key} > {WORK_CAP[key]}')


def row_string(P_row, c0, c1):
    seg = P_row[c0:c1 + 1]
    if (seg == OUT).any() or (seg == PR.OTHER).any():
        return None
    return ''.join(BASES[x] for x in seg[seg < 4])


def build_events(fp: PR.Panel, row_a: int, k: int = 15, gap: int = 6, max_len: int = 120):
    P = fp.P
    inside = P != OUT
    mx = np.where(inside, P, 0).max(0)
    mn = np.where(inside, P, 9).min(0)
    var = (mx != mn) & (inside.sum(0) >= 2)
    maj = ((P < 4).sum(0) * 2 >= np.maximum(inside.sum(0), 1)) & (inside.sum(0) > 0)
    cb = np.cumsum(maj)
    cols = np.nonzero(var)[0]
    C = P.shape[1]
    out = []
    if len(cols) == 0:
        return out
    groups = [[int(cols[0]), int(cols[0])]]
    for c in cols[1:].tolist():
        if cb[c] - cb[groups[-1][1]] <= gap:
            groups[-1][1] = c
        else:
            groups.append([c, c])
    steps = [k, k + 10, k + 20, k + 30]
    combos = sorted(((a, b) for a in steps for b in steps), key=lambda t: (t[0] + t[1], abs(t[0] - t[1])))
    _work('events', len(groups))
    for f, l in groups:
        ev = None
        for kl, kr in combos:
            cl = int(np.searchsorted(cb, cb[f] - kl - (1 if maj[f] else 0), side='right'))
            cr = int(np.searchsorted(cb, cb[l] + kr, side='left'))
            cl, cr = max(0, cl), min(C - 1, cr)
            cur = row_string(P[row_a], cl, cr)
            alle = collections.defaultdict(list)
            for r in range(fp.R):
                s = row_string(P[r], cl, cr)
                if s is not None and len(s) >= 12:
                    alle[s].append(r)
            keys = sorted(alle, key=len)
            sub = any(a != b and a in b for i, a in enumerate(keys) for b in keys[i + 1:])
            longest = max((len(s) for s in keys), default=0)
            ev = dict(first=f, last=l, cl=cl, cr=cr, alleles=dict(alle), anchor=[kl, kr], cur=cur,
                      substring_conflict=sub, longest=longest)
            if not sub or longest > max_len + 40:
                break
        ev['resolvable'] = (ev['cur'] is not None and not ev['substring_conflict'] and
                            ev['longest'] <= max_len and len(ev['alleles']) >= 2)
        out.append(ev)
    return out


def local_ref(fp: PR.Panel, row_a: int, cl: int, cr: int, donor: int = None, flank: int = 170):
    """Sequence of row_a around columns [cl, cr], with that window taken from the donor row."""
    pb = fp.Pbase
    sa = fp.seq[row_a]
    a_cl, a_cr = int(pb[row_a, cl]), int(pb[row_a, cr + 1])
    left = sa[max(0, a_cl - flank):a_cl]
    right = sa[a_cr:a_cr + flank]
    if donor is None or donor == row_a:
        mid = sa[a_cl:a_cr]
    else:
        mid = fp.seq[donor][int(pb[donor, cl]):int(pb[donor, cr + 1])]
    return left + mid + right


def count_events(events, fp: PR.Panel, row_a: int, partner_seq_of, reads, max_edit: int = 2):
    """reads: dicts with seq, reverse, c0, c1. partner_seq_of(cl, cr) -> partner haplotype's local
    sequence (or None). Fills per event: sup[s] = [fwd, rev] reads that contain s, fit haplotype A
    carrying s, and are not explained at least as well by the partner; cur_total = reads containing
    the pass-1 version and fitting A (partner ignored)."""
    import edlib
    order = sorted(range(len(events)), key=lambda i: events[i]['cl'])
    cls = np.array([events[i]['cl'] for i in order])
    for e in events:
        e['sup'] = {}
        e['cur_total'] = 0
        e['raw'] = collections.Counter()
        e['_ref'] = {}
        e['_pref'] = False
    for rd in reads:
        lo = int(np.searchsorted(cls, rd['c0'], side='left'))
        hi = int(np.searchsorted(cls, rd['c1'], side='right'))
        seq = rd['seq']
        for j in range(lo, hi):
            e = events[order[j]]
            if not e['resolvable'] or e['cr'] > rd['c1']:
                continue
            for s, rows in e['alleles'].items():
                if s not in seq:
                    continue
                e['raw'][len(s) if s != e['cur'] else -1] += 1
                if s not in e['_ref']:
                    e['_ref'][s] = local_ref(fp, row_a, e['cl'], e['cr'], None if s == e['cur'] else rows[0])
                _work('align_calls')
                da = edlib.align(seq, e['_ref'][s], mode='HW', task='distance', k=max_edit)['editDistance']
                if da < 0:
                    continue
                if s == e['cur']:
                    e['cur_total'] += 1
                if e['_pref'] is False:
                    e['_pref'] = partner_seq_of(e['cl'], e['cr'])
                if e['_pref']:
                    _work('align_calls')
                    db = edlib.align(seq, e['_pref'], mode='HW', task='distance', k=da)['editDistance']
                    if db >= 0:
                        continue
                v = e['sup'].setdefault(s, [0, 0])
                v[1 if rd['reverse'] else 0] += 1
    for e in events:
        e.pop('_ref', None)
        e.pop('_pref', None)
    return events


def decide(events, fp: PR.Panel, row_a: int, partner_string_of, D: float, read_len: int = 150,
           sup_share: float = 0.25, min_sup: int = 3, cur_ratio: float = 0.25, het_share: float = 0.3,
           pile: float = 2.5):
    patches, notes = [], []
    for i, e in enumerate(events):
        if not e['resolvable']:
            continue
        cur = e['cur']
        exp_h = (D / 2.0) * max(0.0, (read_len - e['longest']) / read_len)
        thr = max(min_sup, sup_share * exp_h)
        cands = [(sum(v), s) for s, v in e['sup'].items()
                 if s != cur and sum(v) >= thr and v[0] >= 1 and v[1] >= 1]
        if not cands:
            continue
        cands.sort(reverse=True)
        n_new, x = cands[0]
        rec = dict(event=i, cl=e['cl'], cr=e['cr'], len_cur=len(cur), len_new=len(x), n_new=n_new,
                   fwd_rev=e['sup'][x], cur_total=e['cur_total'], exp_h=round(exp_h, 1),
                   n_versions=len(e['alleles']), anchor=e['anchor'],
                   new_sha=L.C.sha256_str(x)[:16], cur_sha=L.C.sha256_str(cur)[:16],
                   n_rows_new=len(e['alleles'][x]), n_rows_cur=len(e['alleles'].get(cur, [])),
                   thr=round(thr, 1))
        if len(cands) > 1 and cands[1][0] * 2 > n_new:
            notes.append(dict(rec, kind='two_supported_versions', second=cands[1][0]))
            continue
        if n_new > pile * exp_h + 3:
            notes.append(dict(rec, kind='pile_up'))
            continue
        shared = partner_string_of(e['cl'], e['cr']) == cur
        rec['shared'] = shared
        if shared:
            ok = n_new >= het_share * (n_new + e['cur_total'])
            kind = 'shared_with_partner'
        else:
            ok = e['cur_total'] <= max(1.0, cur_ratio * n_new)
            kind = 'own_version_unsupported'
        if not ok:
            notes.append(dict(rec, kind='own_version_still_supported', shared=shared))
            continue
        patches.append(dict(rec, kind=kind, new=x, donor_row=e['alleles'][x][0]))
    return patches, notes


def patched_path(fp: PR.Panel, row_a: int, patches, events) -> PR.Path:
    segs, pos = [], 0
    for p in sorted(patches, key=lambda p: p['cl']):
        e = events[p['event']]
        if e['cl'] < pos:
            continue
        if e['cl'] > pos:
            segs.append((row_a, pos, e['cl']))
        segs.append((p['donor_row'], e['cl'], e['cr'] + 1))
        pos = e['cr'] + 1
    if pos < fp.C:
        segs.append((row_a, pos, fp.C))
    return PR.merge_path(segs)


def spell(fp: PR.Panel, codes: np.ndarray, prm: PR.Params):
    T = np.where((np.arange(PR.NCODE)[None, :] != codes[:, None]), 1000.0, 0.0).astype(np.float32)
    T[:, OUT] = 0
    T[codes == OUT, :] = 0
    c, path = PR.viterbi(fp, T, prm, set(), True)
    return int(round(c // 1000)), path.switches


def structural_filter(panel: PR.Panel, fam, row_p1: int, max_indel_cols: int = 150):
    """Family rows whose insertion/deletion layout is close to the pass-1 allele (rows with a large
    structural difference would make every column of a long stretch 'variable')."""
    P = panel.P
    a = P[row_p1]
    keep = []
    for r in fam:
        b = P[r]
        inside = (a != OUT) & (b != OUT)
        d = int((inside & ((a < 4) != (b < 4))).sum())
        if d <= max_indel_cols or r == row_p1:
            keep.append(r)
    return keep


def build_haplotype(panel: PR.Panel, fam, row_p1: int, partner_row: int, reads, D: float, **kw):
    """Patch the pass-1 allele (panel row row_p1) using its family rows. Returns dict with the
    sub-panel, path, patches and notes."""
    fam = structural_filter(panel, list(fam), row_p1, kw.pop('max_indel_cols', 150))
    _work('family_rows', len(fam))
    fp = L.subpanel(panel, fam)
    ra = fam.index(row_p1)
    events = build_events(fp, ra)
    pb = panel.Pbase
    pseq = panel.seq[partner_row]

    def partner_seq_of(cl, cr, flank=170):
        seg = panel.P[partner_row, cl:cr + 1]
        if (seg == OUT).all():
            return None
        a, b = int(pb[partner_row, cl]), int(pb[partner_row, cr + 1])
        return pseq[max(0, a - flank):b + flank] or None

    def partner_string_of(cl, cr):
        return row_string(panel.P[partner_row], cl, cr)
    max_edit = kw.pop('max_edit', 2)
    count_events(events, fp, ra, partner_seq_of, reads, max_edit)
    patches, notes = decide(events, fp, ra, partner_string_of, D, **kw)
    path = patched_path(fp, ra, patches, events)
    return dict(fp=fp, fam=fam, ra=ra, events=events, patches=patches, notes=notes, path=path)
