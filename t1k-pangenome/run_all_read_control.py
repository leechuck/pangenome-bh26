"""Development cohort control: bypass T1K candidate extraction without a graph."""
import argparse
import csv
import hashlib
from pathlib import Path
from graph_recruitment import run

ROOT=Path('/home/leechuck/hla/t1k-pangenome')
COHORT_SHA='b378e30b28ea4b934dc8eed1c6d5c6c13336504bb6e80dd2f90870d73cad8c9d'


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--index',type=int,required=True)
    p.add_argument('--cohort',type=Path,required=True)
    a=p.parse_args()
    if hashlib.sha256(a.cohort.read_bytes()).hexdigest()!=COHORT_SHA:
        raise ValueError('Development cohort changed')
    rows=list(csv.DictReader(a.cohort.open(),delimiter='\t'))
    if not 0<=a.index<len(rows):raise ValueError('Invalid donor index')
    donor=rows[a.index]['donor']
    run(ROOT/'development/linear-v1/ipd_genome'/donor,
        Path('/home/leechuck/hla/hla-typer/reads')/donor,
        ROOT/'t1k-references/genome-v2',ROOT/'tools/t1k-evidence-v1',
        ROOT/'development/graph-recruitment-v1'/donor/'all','all',None,4)
