"""Compare two score_rows TSV files restricted to one read set (row-level identity check). Counts only.
  python cmp_rows.py NEW.tsv OLD.tsv READSET OUT.json"""
import csv
import json
import sys
new, old, rs, out = sys.argv[1:5]
key = lambda r: (r['readset'], r['donor'], r['gene'], r['haplotype'], r['arm'])
N = {key(r): r for r in csv.DictReader(open(new), delimiter='\t') if r['readset'] == rs}
O = {key(r): r for r in csv.DictReader(open(old), delimiter='\t') if r['readset'] == rs}
diff = [k for k in N if k in O and N[k] != O[k]]
res = dict(rows_new=len(N), rows_old=len(O), only_new=len(set(N) - set(O)), only_old=len(set(O) - set(N)),
           differing=len(diff), differing_examples=[list(k) for k in diff[:10]],
           identical=not diff and set(N) == set(O))
json.dump(res, open(out, 'w'), indent=1)
print(res)
