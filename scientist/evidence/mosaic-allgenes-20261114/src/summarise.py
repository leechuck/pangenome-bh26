"""Local summary of the row-level outputs (counts only; runs on the workstation on copied TSVs).
  python3 summarise.py RESULTS_DIR TAG ARMS > tables.md
Reads score_rows_TAG.tsv (per copy x arm) and haplotypes_TAG.tsv (per builder haplotype)."""
import collections
import csv
import sys
from pathlib import Path

res, tag, arms = Path(sys.argv[1]), sys.argv[2], sys.argv[3].split(',')
GENES = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')
T = lambda x: x == 'True'
rows = list(csv.DictReader(open(res / f'score_rows_{tag}.tsv'), delimiter='\t'))
haps = list(csv.DictReader(open(res / f'haplotypes_{tag}.tsv'), delimiter='\t')) \
    if (res / f'haplotypes_{tag}.tsv').exists() else []
by = collections.defaultdict(dict)
for r in rows:
    by[(r['readset'], r['arm'])][(r['donor'], r['gene'], r['haplotype'])] = r


def ed(r):
    return None if r['ed_final'] in ('None', '') else int(r['ed_final'])


def base_arm(arm):
    return 'B' if arm[0] in ('B', 'D') else 'A'


def rule_of(arm):
    return arm[2:4] if len(arm) > 2 else None


def unresolved(rs, gene, arm):
    """(haplotypes flagged under the arm's rule, of which with a patched call, all builder haplotypes)"""
    r = rule_of(arm)
    hs = [h for h in haps if h['readset'] == rs and h['arm'] == arm[0] and (gene is None or h['gene'] == gene)
          and h['hap'] != '']
    if not r:
        return '', '', len(hs)
    fl = [h for h in hs if r in h['rules'].split(',')]
    pc = [h for h in fl if h['decision'] in ('proposed', 'rename_to_known_allele')]
    return len(fl), len(pc), len(hs)


def line(rs, arm, gene):
    R = [r for k, r in by[(rs, arm)].items() if gene is None or k[1] == gene]
    if not R:
        return None
    ab = [r for r in R if T(r['ipd_absent'])]
    pr = [r for r in R if not T(r['ipd_absent']) and r['truth4']]
    right = lambda r: T(r['recovered_final']) or T(r['final_body_ok'])
    eds = [ed(r) for r in R if ed(r) is not None]
    eda = [ed(r) for r in ab if ed(r) is not None]
    base = by.get((rs, base_arm(arm)), {})
    closer = farther = same = 0
    if arm not in ('A', 'B'):
        for k, r in by[(rs, arm)].items():
            if (gene is None or k[1] == gene) and T(r['final_mos']) and k in base:
                a, b = ed(r), ed(base[k])
                if a is None or b is None:
                    continue
                closer += a < b
                farther += a > b
                same += a == b
    fl, pc, nh = unresolved(rs, gene, arm)
    na = sum(ed(r) is None for r in R)
    return (f"| {gene or 'all'} | {arm} | {len(ab)} | {sum(T(r['body_recovered']) or T(r['recovered']) for r in ab)} | "
            f"{sum(T(r['recovered']) for r in ab)} / {sum(T(r['recovered_final']) for r in ab)} | "
            f"{sum(T(r['named4']) for r in pr)}/{len(pr)} | {sum(T(r['final_mos']) for r in R if not T(r['ipd_absent']))} | "
            f"{sum(T(r['final_mos']) and right(r) for r in ab)} / {sum(T(r['final_mos']) and not right(r) for r in ab)} | "
            f"{fl}{'' if fl == '' else f' ({pc})'} | "
            f"{sum(eds) / max(len(eds), 1):.1f} | {sum(eda) / max(len(eda), 1):.1f} | "
            f"{'' if arm in ('A', 'B') else f'{closer} / {farther} / {same}'} | {na} |")


HDR = ('| gene | arm | IPD-absent copies | gene body recovered | byte-exact any slot / final | 4-field named correct (IPD-present) '
       '| novel final on IPD-present | novel final on IPD-absent: right body / not right | unresolved haplotypes (with a patched call) '
       '| mean wrong bases per copy (all) | mean wrong bases (IPD-absent) | novel calls closer / farther / same vs base | copies without edit distance |')
SEP = '|' + '---|' * HDR.count('|')[:-0] if False else '|' + '|'.join(['---'] * (HDR.count('|') - 1)) + '|'
for rs in ('D-sim', 'D-real'):
    if not any(k[0] == rs for k in by):
        continue
    print(f'\n### {rs}\n')
    print(HDR)
    print(SEP)
    for gene in list(GENES) + [None]:
        for arm in arms:
            if (rs, arm) not in by:
                continue
            ln = line(rs, arm, gene)
            if ln:
                print(ln)
