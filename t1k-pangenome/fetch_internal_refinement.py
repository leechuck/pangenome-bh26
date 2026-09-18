"""Collect the gene-internal refinement comparison on exposed donors only."""
import subprocess
import argparse
from score_graph_pair_refinement import HERE,run


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pair-specific',action='store_true');a=p.parse_args()
    name='graph-pairwise-refinement' if a.pair_specific else 'graph-internal-refinement'
    version='v3-pairwise' if a.pair_specific else 'v2-internal'
    snapshot=HERE/'work'/(name+'-snapshot');snapshot.mkdir(parents=True,exist_ok=True)
    source=('hohndor@ilogin.ibex.kaust.edu.sa:'
            '/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/t1k-pangenome/development/graph-pair-refinement-'+version+'/')
    subprocess.run(['rsync','-az','--include=*/','--include=manifest.json','--include=COMPLETE.json',
                    '--include=calls.json','--include=decisions.json','--exclude=*',source,str(snapshot)+'/'],check=True)
    run(snapshot,HERE/'development'/name)
