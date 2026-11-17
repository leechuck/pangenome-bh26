"""Summarise COV.json files (counts only): per donor and true DRB1 copy, own-copy coverage in the
intron-1 hole columns (3300-5600) and elsewhere, recruited reads vs T1K candidates, sim vs real.
  python cov_table.py COV_DIR READSET_REAL READSET_SIM > table.md"""
import json
import statistics
import sys
from pathlib import Path

d = Path(sys.argv[1])
real, sim = sys.argv[2], sys.argv[3]
H0, H1 = 3300, 5600


def stats(c, key='own_unique'):
    hole = [w[key] for w in c['windows'] if H0 <= w['col0'] <= H1]
    rest = [w[key] for w in c['windows'] if not (H0 <= w['col0'] <= H1)]
    low = sum(1 for w in c['windows'] if H0 <= w['col0'] <= H1 and w[key] <= 1)
    return (statistics.median(hole) if hole else float('nan'), statistics.median(rest) if rest else float('nan'), low, len(hole))


def runs(cr, cs, key='own_unique'):
    """column runs where real <= 1 and sim >= 5"""
    out, cur = [], None
    for wr, ws in zip(cr['windows'], cs['windows']):
        if wr[key] <= 1 and ws[key] >= 5:
            if cur and wr['col0'] <= cur[1] + 250:
                cur[1] = wr['col1']
            else:
                cur = [wr['col0'], wr['col1']]
                out.append(cur)
    return [f'{a}-{b}' for a, b in out if b - a >= 250]


print('| donor | copy | IPD-absent | pass-1 (real) | real recruited: hole median / rest median / hole windows <=1 | real T1K candidates: hole median / windows <=1 | sim recruited: hole median / rest median / windows <=1 | columns with real <=1 and sim >=5 (runs >= 250 cols) |')
print('|---|---|---|---|---|---|---|---|')
for p in sorted(d.glob(f'{real}_*/COV.json')):
    donor = p.parent.name.split('_', 1)[1]
    R = json.loads(p.read_text())
    sp = d / f'{sim}_{donor}/COV.json'
    S = json.loads(sp.read_text()) if sp.exists() else None
    for i, c in enumerate(R['sources']['recruited']['copies']):
        a = stats(c)
        b = stats(R['sources']['t1k_candidates']['copies'][i])
        s = stats(S['sources']['recruited']['copies'][i]) if S else None
        rr = runs(c, S['sources']['recruited']['copies'][i]) if S else []
        print(f"| {donor} | {c['hap']} | {R['copies'][i]['ipd_absent']} | {','.join(x.split('*')[1] for x in R['pass1'])} | "
              f"{a[0]:.0f} / {a[1]:.0f} / {a[2]} of {a[3]} | {b[0]:.0f} / {b[2]} | "
              + (f"{s[0]:.0f} / {s[1]:.0f} / {s[2]}" if s else 'n/a') + f" | {', '.join(rr) or '-'} |")
