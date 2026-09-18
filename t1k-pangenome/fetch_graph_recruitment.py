"""Transfer small graph recruitment results and score development outcomes."""
from pathlib import Path
import subprocess
from score_graph_recruitment import run


if __name__=='__main__':
    here=Path(__file__).resolve().parent
    snapshot=here/'work/graph-recruitment-snapshot';snapshot.mkdir(parents=True,exist_ok=True)
    source=('hohndor@ilogin.ibex.kaust.edu.sa:'
            '/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/'
            't1k-pangenome/development/graph-recruitment-v1/')
    subprocess.run(['rsync','-az','--include=*/','--include=manifest.json',
                    '--include=COMPLETE.json','--include=calls.json',
                    '--include=*_genotype.tsv','--exclude=*',source,str(snapshot)+'/'],check=True)
    run(snapshot,here/'development/graph-recruitment-pilot','HG00658')
