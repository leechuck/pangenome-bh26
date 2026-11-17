"""Names, first/last MSA column and length of every IPD DRB1 row (no sequence)."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L
C, PR = L.C, L.PR
recs = C.load_records()
msa = C.load_msas(recs, genes=('DRB1',))['DRB1']
out = ['name\tall_names\tstart_col\tend_col\tlength\tfirst_exon_start_col\tlast_exon_end_col']
for r, names in enumerate(msa.ipd_names):
    if not names:
        continue
    a = msa.aligned[r]
    bc = [i for i, ch in enumerate(a) if ch != '-']
    ex = msa.exons(r)
    out.append(f"{names[0]}\t{','.join(names)}\t{bc[0]}\t{bc[-1]}\t{len(bc)}\t{bc[ex[0][0]]}\t{bc[ex[-1][1]-1]}")
Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-drb1-20261113/results/drb1_row_extents.tsv').write_text('\n'.join(out) + '\n')
