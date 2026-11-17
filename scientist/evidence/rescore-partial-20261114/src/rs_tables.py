"""Tables, paired contrasts and validation from the row-level files. Usage: rs_tables.py RESULTS_DIR"""
import collections
import json
import statistics
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rs_common as R

INTS = ['copy', 'ipd_absent', 'n_tied', 'missing', 'unresolved', 'true_gb_len', 'true_cod_len', 'call_gb_len',
        'gb_dist', 'cod_dist', 'nc_dist', 'nc_negative', 'gb_dist_best', 'cod_dist_best', 'full_exact',
        'full_exact_best', 'gb_exact', 'cod_exact', 'gb_exact_best', 'gb_defined', 'novel_final', 'inP',
        'arm_complete', 'body_final_any', 'fz_recovered', 'fz_recovered_final', 'fz_named4', 'fz_body_recovered', 'fz_final_mos']
GENES = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']


def load_rows(p):
    with open(p) as f:
        head = f.readline().rstrip('\n').split('\t')
        out = []
        for line in f:
            v = line.rstrip('\n').split('\t')
            r = dict(zip(head, v))
            for k in INTS:
                if k in r:
                    r[k] = int(r[k]) if r[k] != '' else None
            out.append(r)
    return out


def pct(n, N):
    return f'{n} ({100 * n / N:.1f}%)' if N else f'{n} (-)'


def mm(vals):
    return f'{statistics.mean(vals):.1f} / {statistics.median(vals):g}' if vals else 'NA'


def md(header, rows):
    s = '| ' + ' | '.join(header) + ' |\n|' + '|'.join('---' for _ in header) + '|\n'
    for r in rows:
        s += '| ' + ' | '.join(str(x) for x in r) + ' |\n'
    return s


def summ(rows, gb_na=False):
    n = len(rows)
    d = dict(n=n, missing=sum(1 for r in rows if r['missing']))
    d['gb_exact'] = None if gb_na else sum(1 for r in rows if r['gb_exact'] == 1)
    d['cod_exact'] = sum(1 for r in rows if r['cod_exact'] == 1)
    for k in (1, 5, 20):
        d[f'gb_le{k}'] = None if gb_na else sum(1 for r in rows if r['gb_dist'] is not None and r['gb_dist'] <= k)
        d[f'cod_le{k}'] = sum(1 for r in rows if r['cod_dist'] is not None and r['cod_dist'] <= k)
    for f in ('gb_dist', 'cod_dist', 'nc_dist'):
        d[f] = None if (gb_na and f != 'cod_dist') else [r[f] for r in rows if r['missing'] == 0 and r[f] is not None]
    d['nc_neg'] = sum(1 for r in rows if r['nc_negative'] == 1)
    return d


def share_tables(rows, arms, gb_na_arms=()):
    subsets = [('all copies', lambda r: True), ('IPD-absent copies', lambda r: r['ipd_absent'] == 1),
               ('IPD-present copies', lambda r: r['ipd_absent'] == 0)]
    out = ''
    for title, f in subsets:
        t1, t2 = [], []
        for a in arms:
            sel = [r for r in rows if r['arm'] == a and f(r)]
            na = a in gb_na_arms
            s = summ(sel, na)
            N = s['n']
            NA = lambda x: 'NA' if x is None else pct(x, N)
            t1.append([a, N, s['missing'], NA(s['gb_exact']), pct(s['cod_exact'], N), NA(s['gb_le1']), NA(s['gb_le5']),
                       NA(s['gb_le20']), pct(s['cod_le1'], N), pct(s['cod_le5'], N), pct(s['cod_le20'], N)])
            t2.append([a, N - s['missing'], mm(s['gb_dist']) if s['gb_dist'] is not None else 'NA',
                       mm(s['cod_dist']), mm(s['nc_dist']) if s['nc_dist'] is not None else 'NA', s['nc_neg']])
        out += f'\n**{title}: shares (denominator = all copies; a missing call counts as not correct)**\n\n'
        out += md(['arm', 'copies', 'missing calls', 'gene body exact', 'coding exact', 'gb <=1', 'gb <=5', 'gb <=20',
                   'cod <=1', 'cod <=5', 'cod <=20'], t1)
        out += f'\n**{title}: wrong bases per called copy, mean / median (missing calls excluded)**\n\n'
        out += md(['arm', 'called copies', 'gene body', 'coding', 'non-coding', 'nc negative (clipped to 0)'], t2)
    return out


def gene_table(rows, arms, gb_na_arms=()):
    hdr = ['gene'] + [f'{a}: copies; gb exact; mean gb / mean cod' for a in arms]
    body = []
    for g in GENES:
        line = [g]
        for a in arms:
            sel = [r for r in rows if r['arm'] == a and r['gene'] == g]
            s = summ(sel, a in gb_na_arms)
            gbm = 'NA' if s['gb_dist'] is None or not s['gb_dist'] else f"{statistics.mean(s['gb_dist']):.1f}"
            cm = f"{statistics.mean(s['cod_dist']):.1f}" if s['cod_dist'] else 'NA'
            ge = 'NA' if s['gb_exact'] is None else s['gb_exact']
            line.append(f"{s['n']}; {ge}; {gbm} / {cm}")
        body.append(line)
    return md(hdr, body)


def sensitivity_table(rows, arms, gb_na_arms=()):
    t = []
    for a in arms:
        sel = [r for r in rows if r['arm'] == a]
        N = len(sel)
        if a in gb_na_arms:
            continue
        w = sum(1 for r in sel if r['gb_exact'] == 1)
        b = sum(1 for r in sel if r['gb_exact_best'] == 1)
        wd = [r['gb_dist'] for r in sel if r['missing'] == 0 and r['gb_dist'] is not None]
        bd = [r['gb_dist_best'] for r in sel if r['missing'] == 0 and r['gb_dist_best'] is not None]
        tied = sum(1 for r in sel if (r['n_tied'] or 0) > 1)
        t.append([a, tied, pct(w, N), pct(b, N), f'{statistics.mean(wd):.2f}' if wd else 'NA',
                  f'{statistics.mean(bd):.2f}' if bd else 'NA'])
    return md(['arm', 'copies with tied alleles', 'gb exact (worst tied allele; primary)', 'gb exact (best tied allele)',
               'mean gb wrong bases (worst)', 'mean gb wrong bases (best)'], t)


def contrast(rows, a, b, subset, metric, group_of):
    ka = {(r['readset'], r['donor'], r['gene'], r['copy']): r for r in rows if r['arm'] == a and subset(r)}
    kb = {(r['readset'], r['donor'], r['gene'], r['copy']): r for r in rows if r['arm'] == b and subset(r)}
    keys = sorted(set(ka) & set(kb))
    pairs = []
    ma, mb = [], []
    if metric in ('gb_exact', 'cod_exact'):
        for k in keys:
            va = 1 if ka[k][metric] == 1 else 0
            vb = 1 if kb[k][metric] == 1 else 0
            pairs.append(dict(donor=k[1], difference=va - vb))
            ma.append(va)
            mb.append(vb)
        u = R.uncertainty(pairs, group_of, scale=100.0)
        gains = sum(1 for p in pairs if p['difference'] > 0)
        losses = sum(1 for p in pairs if p['difference'] < 0)
        return dict(n=len(pairs), a=sum(ma), b=sum(mb), gains=gains, losses=losses, **{k: u[k] for k in ('diff', 'lo', 'hi', 'p', 'families')})
    for k in keys:
        x, y = ka[k], kb[k]
        if x['missing'] or y['missing'] or x[metric] is None or y[metric] is None:
            continue
        pairs.append(dict(donor=k[1], difference=x[metric] - y[metric]))
        ma.append(x[metric])
        mb.append(y[metric])
    u = R.uncertainty(pairs, group_of, scale=1.0)
    return dict(n=len(pairs), a=statistics.mean(ma) if ma else None, b=statistics.mean(mb) if mb else None,
                gains=sum(1 for p in pairs if p['difference'] < 0), losses=sum(1 for p in pairs if p['difference'] > 0),
                **{k: u[k] for k in ('diff', 'lo', 'hi', 'p', 'families')})


def fmt_contrast(c, share):
    if c['n'] == 0 or c['diff'] is None:
        return ['0', 'NA', 'NA', 'NA', 'NA']
    if share:
        return [c['n'], f"{c['a']} vs {c['b']} (+{c['gains']}/-{c['losses']})", f"{c['diff']:+.2f} pp",
                f"[{c['lo']:+.2f}, {c['hi']:+.2f}]", f"{c['p']:.3g}"]
    return [c['n'], f"{c['a']:.2f} vs {c['b']:.2f} (closer {c['gains']}/farther {c['losses']})", f"{c['diff']:+.2f}",
            f"[{c['lo']:+.2f}, {c['hi']:+.2f}]", f"{c['p']:.3g}"]


def contrast_tables(rows, pairs_list, group_of, tag, metrics=('gb_exact', 'cod_exact', 'gb_dist', 'cod_dist')):
    subsets = [('all copies', lambda r: True), ('IPD-absent', lambda r: r['ipd_absent'] == 1),
               ('IPD-present', lambda r: r['ipd_absent'] == 0)]
    out = ''
    res = {}
    lines = []
    for a, b, ms in pairs_list:
        for sname, f in subsets:
            for m in ms:
                c = contrast(rows, a, b, f, m, group_of)
                res[f'{a}-{b}|{sname}|{m}'] = c
                lines.append([f'{a} - {b}', sname, {'gb_exact': 'share gene body exact (pp)', 'cod_exact': 'share coding exact (pp)',
                                                    'gb_dist': 'mean wrong bases, gene body', 'cod_dist': 'mean wrong bases, coding'}[m]]
                             + fmt_contrast(c, m.endswith('exact')))
    out += md(['contrast', 'subset', 'metric', 'n paired', 'values (a vs b)', 'difference', '95% CI (family bootstrap)', 'sign-flip p'], lines)
    return out, res


def drb1_special(rows, readsets):
    out = ''
    for novel_arm, base_arm in (('C4', 'A'), ('D4', 'B')):
        lines = []
        dist_lines = []
        for rs in readsets:
            R_ = [r for r in rows if r['readset'] == rs]
            ka = {(r['donor'], r['gene'], r['copy']): r for r in R_ if r['arm'] == base_arm}
            kn = {(r['donor'], r['gene'], r['copy']): r for r in R_ if r['arm'] == novel_arm}
            if not kn:
                continue
            for lab, f in (('IPD-absent', 1), ('IPD-present', 0)):
                deltas = []
                n_nov = n_inc = 0
                for k, r in kn.items():
                    if r['ipd_absent'] != f or r['novel_final'] != 1 or r['missing']:
                        continue
                    n_nov += 1
                    b = ka.get(k)
                    if b is None or b['missing'] or b['gb_dist'] is None or r['gb_dist'] is None:
                        n_inc += 1
                        continue
                    deltas.append(r['gb_dist'] - b['gb_dist'])
                sm = sum(1 for x in deltas if x < 0)
                eq = sum(1 for x in deltas if x == 0)
                lg = sum(1 for x in deltas if x > 0)
                lines.append([rs, lab, n_nov, len(deltas), n_inc, sm, eq, lg])
                if deltas:
                    ds = sorted(deltas)
                    dist_lines.append([rs, lab, len(ds), ds[0], statistics.median(ds), ds[-1],
                                       sum(1 for x in ds if x <= -20), sum(1 for x in ds if -20 < x <= -5),
                                       sum(1 for x in ds if -5 < x < 0), sum(1 for x in ds if x == 0),
                                       sum(1 for x in ds if 0 < x < 5), sum(1 for x in ds if 5 <= x < 20),
                                       sum(1 for x in ds if x >= 20)])
        out += f'\n**{novel_arm} final novel (patched) calls vs the {base_arm} call for the same copy (gene-body distance to truth)**\n\n'
        out += md(['read set', 'copy type', 'copies with novel final call', 'comparable', 'not comparable (arm {} call missing)'.format(base_arm),
                   f'closer than {base_arm}', 'equal', 'farther'], lines)
        out += f'\nChange in gene-body wrong bases ({novel_arm} minus {base_arm}); negative = patched call is closer:\n\n'
        out += md(['read set', 'copy type', 'n', 'min', 'median', 'max', '<=-20', '-19..-5', '-4..-1', '0', '+1..+4', '+5..+19', '>=+20'], dist_lines)
    return out


def validate_confirm(rows, frozen, sets):
    out = ''
    fz = {(r['donor'], r['gene'], r['arm']): r for r in frozen}
    by_dg = collections.defaultdict(list)
    for r in rows:
        by_dg[(r['donor'], r['gene'], r['arm'])].append(r)
    for sname, filt in sets:
        t = []
        mism = []
        ipd_called = [0]
        for arm in sorted({r['arm'] for r in rows}):
            # IPD-absent exact recovery
            fr_unm = mine_full = mine_gb = tot = 0
            for (d, g, a), rs in by_dg.items():
                if a != arm or not filt(rs[0]):
                    continue
                f = fz[(d, g, a)]
                fr_unm += int(f['unm_recovered'])
                tot += int(f['unm_copies'])
                mf = sum(1 for r in rs if r['ipd_absent'] == 1 and r['full_exact'] == 1)
                mg = sum(1 for r in rs if r['ipd_absent'] == 1 and r['gb_exact'] == 1)
                mine_full += mf
                mine_gb += mg
                if mf != int(f['unm_recovered']):
                    mism.append((arm, d, g, 'absent', f['unm_recovered'], mf))
            # four-field named correct on eligible donor-genes
            fe = fc = m_full = m_gb = m_cod = both = 0
            for (d, g, a), rs in by_dg.items():
                if a != arm or not filt(rs[0]):
                    continue
                f = fz[(d, g, a)]
                if not int(f['eligible']):
                    continue
                fe += 1
                fc += int(f['correct'])
                full = all(r['full_exact'] == 1 for r in rs)
                gb = all(r['gb_exact'] == 1 for r in rs)
                cod = all(r['cod_exact'] == 1 for r in rs)
                m_full += full
                m_gb += gb
                m_cod += cod
                if arm != 'ipd':
                    ntie = max((r['n_tied'] or 0) for r in rs)
                    if gb != bool(int(f['correct'])):
                        kind = ('frozen_right_mine_not' if int(f['correct']) else
                                ('frozen_not_called(tie/non-4-field)' if not int(f['called']) else 'frozen_wrong_name_but_identical_gene_body'))
                        mism.append((arm, d, g, kind, f['correct'], int(gb), ntie))
                if arm == 'ipd':
                    ipd_called[0] += int(f['called'])
            t.append([arm, f'{fr_unm}/{tot}', f'{mine_full}', f'{mine_gb}', f'{fc}/{fe}', m_full,
                      m_gb if arm != 'ipd' else 'NA', m_cod])
        out += f'\n**Set {sname}: frozen scorer vs this scorer (distance = 0)**\n\n'
        out += md(['arm', 'frozen IPD-absent recovered / copies', 'mine: whole-sequence identical', 'mine: gene body identical',
                   'frozen four-field correct / eligible donor-genes', 'mine: both copies whole-sequence identical (core differs from IPD allele ends in many copies)',
                   'mine: both copies gene body identical', 'mine: both copies coding identical'], t)
        cnt = collections.Counter((m[0], m[3]) for m in mism)
        out += f"\nipd arm: frozen four-field endpoint counts a donor-gene as called only if both slots carry one four-field name: {ipd_called[0]} of the eligible donor-genes were called that way but far fewer match the truth name; the ipd reference has no or little intron/UTR content, so four-field (non-coding) differences cannot be resolved there. Coding-identical is a different, weaker test and is not expected to match the frozen count.\n"
        out += '\nGene-body-identical vs frozen four-field correct, mismatches at donor-gene level, counts by arm and kind: ' + \
            (', '.join(f'{k[0]}/{k[1]}: {v}' for k, v in sorted(cnt.items())) or 'none') + '\n'
        out += '\nFirst mismatches (arm, donor, gene, kind, frozen, mine): ' + '; '.join(str(m) for m in mism[:25]) + '\n'
    return out


def validate_typer(rows, rs, arms):
    t = []
    mism = []
    for a in arms:
        sel = [r for r in rows if r['arm'] == a and r['readset'] == rs]
        if not sel:
            continue
        ab = [r for r in sel if r['ipd_absent'] == 1]
        pr = [r for r in sel if r['ipd_absent'] == 0 and r['truth4']]
        pr_all = [r for r in sel if r['ipd_absent'] == 0]
        t.append([a, len(ab), sum(r['fz_recovered_final'] for r in ab), sum(r['fz_recovered'] for r in ab),
                  sum(1 for r in ab if r['full_exact'] == 1), sum(1 for r in ab if r['full_exact_best'] == 1),
                  sum(1 for r in ab if r['gb_exact'] == 1), len(pr), sum(r['fz_named4'] for r in pr),
                  sum(1 for r in pr if r['full_exact'] == 1), sum(1 for r in pr if r['gb_exact'] == 1),
                  len(pr_all), sum(1 for r in pr_all if r['gb_exact'] == 1)])
        for r in ab:
            if (r['fz_recovered_final'] == 1) != (r['full_exact_best'] == 1):
                mism.append((a, r['donor'], r['gene'], r['copy'], 'absent', r['fz_recovered_final'], r['full_exact_best'], r['n_tied']))
        for r in pr:
            if (r['fz_named4'] == 1) != (r['gb_exact'] == 1):
                kind = 'frozen_right_gb_not' if r['fz_named4'] == 1 else (
                    'gb_identical_frozen_wrong(tied)' if (r['n_tied'] or 0) > 1 else 'gb_identical_other_4field_name_or_ends')
                mism.append((a, r['donor'], r['gene'], r['copy'], kind, r['fz_named4'], r['gb_exact'], r['n_tied']))
    out = md(['arm', 'IPD-absent copies', 'frozen recovered_final', 'frozen recovered (any slot)', 'mine whole-seq exact (worst tied)',
              'mine whole-seq exact (any tied)', 'mine gene-body exact', 'IPD-present copies with one 4-field truth',
              'frozen named4 correct', 'mine whole-seq exact', 'mine gene-body exact', 'all IPD-present copies', 'mine gene-body exact (all present)'], t)
    cnt = collections.Counter((m[0], m[4]) for m in mism)
    out += '\nPer-copy disagreements (frozen four-field correct vs gene-body identical on IPD-present copies with one 4-field truth; frozen recovered_final vs whole-sequence identical on IPD-absent copies): count by arm and kind\n\n' + \
        ('\n'.join(f'- {k}: {v}' for k, v in sorted(cnt.items())) or '- none') + '\n'
    out += '\nFirst disagreements (arm, donor, gene, copy, kind, frozen, mine, n_tied): ' + '; '.join(str(m) for m in mism[:20]) + '\n'
    return out


def main(resdir):
    res = Path(resdir)
    sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src')
    import common as C
    gmap = {d: tuple(sorted(v)) for d, v in C.exclusion_map().items()}
    summary = {}
    # ---------------- confirmatory
    if (res / 'rows_confirm.tsv').exists():
        rows = load_rows(res / 'rows_confirm.tsv')
        g = json.loads((res / 'groups_confirm.json').read_text())
        gconf = {d: tuple(v) for d, v in g['groups'].items()}
        arms = ['ipd', 'genomic', 'genomic_core_hprc_asian', 'genomic_core_nonasian']
        frozen = [dict(zip(['donor', 'gene', 'arm', 'eligible', 'called', 'correct', 'unm_copies', 'unm_recovered'], l.rstrip('\n').split('\t')))
                  for l in open(res / 'frozen_confirm.tsv').read().splitlines()[1:]]
        sets = [('P', lambda r: r['inP'] == 1), ('A', lambda r: True)]
        txt = '# Confirmatory test, real reads (partial-credit re-scoring)\n'
        for sname, f in sets:
            sub = [r for r in rows if f(r)]
            nd = len({r['donor'] for r in sub})
            txt += f'\n## Set {sname}: {nd} donors, {len(sub) // len(arms)} true copies per arm\n'
            gbna = ('ipd',)
            txt += share_tables(sub, arms, gbna)
            txt += '\n**By gene**\n\n' + gene_table(sub, arms, gbna)
            txt += '\n**Sensitivity to the tie rule (gene body, arms with a full gene body)**\n\n' + sensitivity_table(sub, arms, gbna)
            pl = [('genomic', 'ipd', ('cod_exact', 'cod_dist')),
                  ('genomic_core_hprc_asian', 'genomic', ('gb_exact', 'cod_exact', 'gb_dist', 'cod_dist')),
                  ('genomic_core_hprc_asian', 'genomic_core_nonasian', ('gb_exact', 'cod_exact', 'gb_dist', 'cod_dist'))]
            ct, cres = contrast_tables(sub, pl, gconf, sname)
            txt += f'\n**Paired contrasts (set {sname}); for exactness a - b in percentage points, for wrong bases a - b (negative = a is closer)**\n\n' + ct
            summary[f'confirm_{sname}'] = cres
        txt += '\n## Validation against the frozen scorer\n' + validate_confirm(rows, frozen, sets)
        # truth-without-exons and unresolved
        txt += f"\nTruth copies without Immuannot coding exons (coding/gene-body undefined): {sum(1 for r in rows if r['true_cod_len'] is None) // len(arms)}\n"
        txt += f"Rows with unresolved call labels: {sum(1 for r in rows if (r['unresolved'] or 0) > 0)}\n"
        txt += f"Gene-body defined share per arm: " + ', '.join(
            f"{a} {sum(1 for r in rows if r['arm'] == a and r['gb_defined'] == 1)}/{sum(1 for r in rows if r['arm'] == a and r['missing'] == 0)}" for a in arms) + '\n'
        (res / 'TABLES_confirm.md').write_text(txt)
    # ---------------- typer
    for rs in ('T-sim', 'T-real'):
        p = res / f'rows_typer_{rs}.tsv'
        if not p.exists():
            continue
        rows = load_rows(p)
        arms = [a for a in 'ABCDE' if any(r['arm'] == a for r in rows)]
        txt = f'# Mosaic-typer test, read set {rs} (all 8 genes), partial-credit re-scoring\n\n{len(rows) // len(arms)} true copies per arm.\n'
        txt += share_tables(rows, arms)
        txt += '\n**By gene**\n\n' + gene_table(rows, arms)
        txt += '\n**Sensitivity to the tie rule**\n\n' + sensitivity_table(rows, arms)
        pl = [('B', 'A', ('gb_exact', 'cod_exact', 'gb_dist', 'cod_dist')), ('D', 'B', ('gb_exact', 'gb_dist')), ('E', 'A', ('gb_exact', 'gb_dist'))]
        ct, cres = contrast_tables(rows, pl, gmap, rs)
        txt += '\n**Paired contrasts**\n\n' + ct
        summary[f'typer_{rs}'] = cres
        txt += '\n## Validation against the frozen scorer\n\n' + validate_typer(rows, rs, arms)
        txt += f"\nRows with unresolved call labels: {sum(1 for r in rows if (r['unresolved'] or 0) > 0)}\n"
        (res / f'TABLES_typer_{rs}.md').write_text(txt)
    # ---------------- drb1
    allrows = []
    for rs in ('T-sim', 'T-real', 'D-sim', 'D-real'):
        p = res / f'rows_drb1_{rs}.tsv'
        if not p.exists():
            continue
        rows = load_rows(p)
        allrows += rows
        arms = [a for a in ('A', 'B', 'C4', 'D4') if any(r['arm'] == a for r in rows)]
        lab = 'test' if rs.startswith('T') else 'DEVELOPMENT'
        txt = f'# DRB1 patch-builder, read set {rs} ({lab}), partial-credit re-scoring\n\n{len(rows) // len(arms)} true DRB1 copies per arm.\n'
        txt += share_tables(rows, arms)
        txt += '\n**Sensitivity to the tie rule**\n\n' + sensitivity_table(rows, arms)
        pl = [('C4', 'A', ('gb_exact', 'cod_exact', 'gb_dist', 'cod_dist')), ('D4', 'B', ('gb_exact', 'cod_exact', 'gb_dist', 'cod_dist')),
              ('B', 'A', ('gb_exact', 'gb_dist'))]
        pl = [x for x in pl if x[0] in arms and x[1] in arms]
        ct, cres = contrast_tables(rows, pl, gmap, rs)
        txt += '\n**Paired contrasts**\n\n' + ct
        summary[f'drb1_{rs}'] = cres
        txt += drb1_special(rows, [rs])
        txt += '\n## Validation against the frozen scorer (DRB1 only)\n\n' + validate_typer(rows, rs, arms)
        # frozen gene-body column
        t = []
        for a in arms:
            ab = [r for r in rows if r['arm'] == a and r['ipd_absent'] == 1]
            t.append([a, len(ab), sum(r['fz_body_recovered'] for r in ab), sum(r['body_final_any'] for r in ab),
                      sum(1 for r in ab if r['gb_exact'] == 1), sum(1 for r in ab if r['gb_exact_best'] == 1)])
        txt += '\nGene-body recovery on IPD-absent DRB1 copies (frozen "gene body exact" column vs mine):\n\n' + md(
            ['arm', 'absent copies', 'frozen gene body exact (candidate in ANY slot)', 'same, restricted to FINAL calls (any slot)', 'mine gene-body exact (assigned slot, worst tied)', 'mine (any tied)'], t)
        (res / f'TABLES_drb1_{rs}.md').write_text(txt)
    (res / 'contrasts.json').write_text(json.dumps(summary, indent=1))
    print('done')


if __name__ == '__main__':
    main(sys.argv[1])
